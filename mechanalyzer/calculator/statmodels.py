""" Compute product energy distribution according to
    different statistical models
"""

import sys
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from phydat import phycon
from autoparse import find


MW_dct_elements = {
    'C': 12e-3,
    'N': 14e-3,
    'O': 16e-3,
    'H': 1e-3,
    'S': 32e-3,
    'P': 31e-3,
    'F': 19e-3,
    'Cl': 35.45e-3
}  # kg/mol


def get_dof_info(block, ask_for_ts=False):
    """ Gets the N of degrees of freedom and MW of each species
        :param block: bimol species of which you want the dofs
        :type block: list(str1, str2)
        :param ask_for_ts: build the dof info also for the ts
        :type ask_for_ts: bool
        :return dof_info: dataframe with vibrat/rot degrees of freedom
            and molecular weight
        :rtype: dataframe(index=species, columns=['vib dof', 'rot dof', 'mw'])
    """
    info_array = np.zeros((2+int(ask_for_ts), 3))
    keys = []
    atoms_ts = 0
    # extract N of dofs and MW
    for i, block_i in enumerate(block):
        info = block_i.splitlines()
        where_name = find.where_in('Species', info)[0]
        where_hind = find.where_in('Hindered', info)
        where_geom = find.where_in('Geometry', info)[0]
        num_atoms = int(info[where_geom].strip().split()[1])
        atoms_ts += num_atoms

        key = info[where_name].strip().split()[1]
        keys.append(key)
        try:
            where_freq = find.where_in('Frequencies', info)[0]
            num_dof = (
                int(info[where_freq].strip().split()[1]) + len(where_hind)
            )
            if 3*num_atoms - num_dof == 6:
                rot_dof = 3
            else:
                rot_dof = 2
        except IndexError:
            # if 1 atom only: no 'Frequencies', set to 0
            num_dof = 0
            rot_dof = 0
        # this allows to get 3N-5 or 3N-6 without analyzing the geometry
        info_array[i, 0] = num_dof
        info_array[i, 1] = rot_dof

        # MW from type of atoms:
        geom_in = where_geom+1
        geom_fin = geom_in+num_atoms
        atoms_array = np.array([geomline.strip().split()[0]
                                for geomline in info[geom_in:geom_fin]])

        info_array[i, 2] = np.sum(np.array([MW_dct_elements[at]
                                  for at in atoms_array], dtype=float))

    # if ask for ts: assume first 2 blocks are 2 reactants of bimol reaction
    # and derive the DOFs of the TS
    if ask_for_ts:
        keys.append('TS')

        # assume there are no linear TSs
        info_array[2, :] = [3*atoms_ts - 7, 3,
                            info_array[0, 2]+info_array[1, 2]]

    dof_info = pd.DataFrame(info_array, index=keys, columns=[
                            'vib dof', 'rot dof', 'mw'])

    return dof_info


def dos_trasl(mass1, mass2, ene_grid, pressure, temp):
    """ Compute the translational density of states per unit volume
        mass1, mass2: MW in kg/mol
        ene_grid: energy grid in kcal/mol (array)
        P, T: pressure and temperature in SI units [Pa, K]
        :return dos_tr_series: dos in mol/kcal
        :rtype: array
    """

    # conversions
    ene_grid_j = ene_grid*4184/phycon.NAVO  # kcal/mol*J/kcal/NAVO=J
    mass1 /= phycon.NAVO  # kg/mol/(molecule/mol)
    mass2 /= phycon.NAVO
    red_mass = (mass1 * mass2) / (mass1 + mass2)    # kg

    # 1 molecule
    # rhotr/V = pi/4*(8m/h2)^3/2*e^1/2
    rho = (np.pi/4*np.power(8*red_mass/phycon.H**2, 3/2) *
           np.power(ene_grid_j, 1/2))  # unit: 1/m3/J

    # consider the molar volume RT/P [m3/1]
    v_mol = 1.38e-23*temp/pressure  # m3/1 for ideal gases
    # 1/m3*(m3/1)/J*(J/kcal)*(mol) = mol/kcal
    rho_kcal_mol = rho*v_mol*4184/phycon.NAVO

    return rho_kcal_mol


class PEDModels:
    """ Statistical models for Product Energy Distributions
    """

    def __init__(self, ped_df, hotfrg, otherfrg,
                 dos_df=None, dof_info=None, ene_bw=None):
        """ initialize variables
            :param ped_df: dataframe(columns:P, rows:T)
                with the Series of energy distrib
            :type ped_df: dataframe(series(float))
            :param hotfrg: selected hot fragment between frag1 and frag2
            :type hotfrg: str
            :param otherfrg: the other fragment
            :type otherfrg: str
            :param dos_df: rovibr dos for each fragment
            :type dos_df: dataframe()
            :param dof_dct: dct with dofs {prodi: Ni}
            :type dof_dct: dct
            :param prod1: fragment of the energy distribution we want
            :type prod1: str
        """
        self.ped_df = ped_df
        self.dos_df = dos_df

        self.ene_bw = ene_bw
        self.prod1 = hotfrg
        self.prod2 = otherfrg
        self.dof_info = dof_info

        self.phi = None
        self.ene1_step = None
        self.ene1_vect = None
        self.rho_rovib_prod1 = None
        self.rho_non1 = None
        self.ene_dos0 = None
        self.f_rho_rovib_prod1 = None
        self.f_rho_rovib_prod2 = None

        try:
            self.mw_dct = dof_info['mw']
            self.vibdof_dct = dof_info['vib dof']
            self.rotdof_dct = dof_info['rot dof']
        except TypeError:
            pass

        self.models_dct = {
            'equip_simple': self.equip_simple,
            'equip_phi': self.equip_phi,
            'rovib_dos': self.rovib_dos,
            'beta_phi1a': self.beta_phi1a,
            'beta_phi2a': self.beta_phi2a,
            'beta_phi3a': self.beta_phi3a
        }

    def compute_ped(self, modeltype):
        """ compute ped according to the desired model
        """

        try:
            ped_df_prod = self.models_dct[modeltype]()
        except KeyError:
            ped_df_prod = None
            _mods = '\n'.join(self.models_dct.keys())
            print('*Error: model not available. '
                  f'Please select among \n {_mods} \n')
            sys.exit()

        return ped_df_prod

    def get_dofs(self):
        """ get dofs from dof_info
        """
        try:
            self.dof_info.empty
        except AttributeError:
            print('Error: DOFs not defined, now exiting\n')
            sys.exit()

        # derive the energy fraction from the equipartition theorem
        try:
            vibdof_prod1, vibdof_prod2 = self.vibdof_dct[[
                self.prod1, self.prod2]]
            rotdof_prod1, rotdof_prod2 = self.rotdof_dct[[
                self.prod1, self.prod2]]
            vibdof_ts = self.vibdof_dct['TS']
            rotdof_ts = self.rotdof_dct['TS']
        except KeyError:
            _str1 = ' '.join([self.prod1, self.prod2])
            _str2 = ' '.join(self.vibdof_dct.keys())
            print(f'incomplete degrees of freedom info: species are \n {_str1}'
                  f'\n while dof info is available for \n {_str2} \n '
                  'Exiting ...')
        return (vibdof_prod1, rotdof_prod1,
                vibdof_prod2, rotdof_prod2,
                vibdof_ts, rotdof_ts)

    def equip_simple(self):
        """ Derive the energy distribution of 1 product from the
            energy equipartition theorem

            :return ped_df_prod: energy distribution of the product prod
            :rtype: dataframe(series(float))
        """
        dofs = self.get_dofs()
        vibdof_prod1, rotdof_prod1, vibdof_prod2, rotdof_prod2, _, _ = dofs

        beta_prod = (vibdof_prod1+rotdof_prod1/2) / \
            (vibdof_prod1+vibdof_prod2+(3+rotdof_prod1+rotdof_prod2)/2)
        # 3/2: 1/2kbT for each rotation, no trasl (ts trasl energy preserved)
        # 9/2: 1/2kbT*6 rotational dofs for products, +3 for relative trasl
        print(f'Fraction of energy transferred to products: {beta_prod:.2f}')

        # rescale all energies with beta: allocate values in new dataframe
        ped_df_prod = pd.DataFrame(index=self.ped_df.index,
                                   columns=self.ped_df.columns, dtype=object)
        for pressure in self.ped_df.columns:
            for temp in self.ped_df.sort_index().index:
                idx_new = self.ped_df[pressure][temp].index * beta_prod
                norm_factor = np.trapz(
                    self.ped_df[pressure][temp].values, x=idx_new)
                vals = self.ped_df[pressure][temp].values/norm_factor
                ped_df_prod[pressure][temp] = pd.Series(vals, index=idx_new)

        return ped_df_prod

    def prob_ene1_fct(self, distr_type):
        """ Derive the energy distribution of 1 product from one
            of the statistical models by Danilack Goldsmith PROCI 2020
            phi is the average fraction of energy transferred to the products

            calculate P(E') = probability of fragment 1 to have energy E'
            a. alfa(E) = PED(E)
            b. select E' and calculate P(E',E) according to normal distribution
            c. integrate over dE: int(P(E',E)*PED(E)dE) = P(E')
        """

        def norm_distr(ene1, ene, phi, ene_bw):
            """ P(ene1; ene) = exp(
                    -(ene1-phi*ene)^2/(2^0.5*sigma(ene_bw))) /
                    ((2*pi)^0.5*sigma(ene_bw)
                )
                mi = phi*ene
                sigma = f(ene_bw, phi*ene)
                ene1 is a number
                ene is a vector
            """

            mtermi = np.array(phi*ene, dtype=float)
            # correlation from Danilack-Goldsmith
            # I add additional fraction of energy transferred to the products
            sigma = np.array(0.87+0.04*(ene_bw+phi*(ene-ene_bw)), dtype=float)
            num = np.exp(-((ene1-mtermi)/(2**0.5)/sigma)**2)
            den = np.power(2*np.pi, 0.5)*sigma

            prob_e1e = num/den

            return prob_e1e

        def init_dos(pressure, temp):
            """ initialize variables for DOS calculation
            """
            rho_rovib_prod1 = self.f_rho_rovib_prod1(self.ene1_vect)
            ene1_vect_w0 = np.concatenate((np.array([0]), self.ene1_vect))
            rho_rovib_prod2 = self.f_rho_rovib_prod2(ene1_vect_w0)
            rho_trasl = dos_trasl(
                self.mw_dct[self.prod1],
                self.mw_dct[self.prod2],
                ene1_vect_w0,
                pressure*101325,
                temp)

            # calculate rho_non1(ene1_vect)
            rho_non1 = []
            for idx in self.ene1_vect:
                # first iter should be 0
                # the sum of the energies in rhovib_prod2 and
                # rho_trasl is always ene1
                idx_ene_int = np.arange(0, idx+1, dtype=int)
                idx_ene_minus_ene_int = idx_ene_int[::-1]
                rho_non1_integrand = (
                    rho_rovib_prod2[idx_ene_int] *
                    rho_trasl[idx_ene_minus_ene_int]
                )
                rho_non1.append(np.trapz(rho_non1_integrand,
                                         x=self.ene1_vect[idx_ene_int]))

            rho_non1 = np.array(rho_non1)

            return rho_rovib_prod1, rho_non1

        def dos(idx_ene1, idx_ene_new):
            """ Calculate density of states """

            prob_ene1ene = []
            # loops over total energy:
            # loop over the idxs where you find enetot in ene1vect
            for idx_ene in idx_ene_new:

                if idx_ene1 == idx_ene:
                    # rho1(ene1)*rhonon1(ene-ene1) = 0
                    prob_ene1ene.append(0)
                    continue
                # index of ene1<ene (fixed ene)
                idx_ene1_array = np.arange(0, idx_ene)
                # index of ene-ene1 (fixed ene)
                idx_ene_minus_ene1_array = idx_ene1_array[::-1]
                rho1_ene1 = self.rho_rovib_prod1[idx_ene1]
                rho_non1 = self.rho_non1[idx_ene_minus_ene1_array[idx_ene1]]
                # rho1(ene1) with ene1<ene (fixed ene)
                rho1_ene1_array = self.rho_rovib_prod1[idx_ene1_array]
                # rhonon1(ene-ene1) with ene1<ene (fixed ene)
                rho_non1_array = self.rho_non1[idx_ene_minus_ene1_array]
                num = rho1_ene1 * rho_non1
                den = np.trapz(rho1_ene1_array*rho_non1_array,
                               x=self.ene1_vect[idx_ene1_array])
                prob_ene1ene.append(num/den)

            prob_ene1ene = np.array(prob_ene1ene)
            # print(prob_ene1ene)
            return prob_ene1ene

        # preallocations
        ped_df_prod = pd.DataFrame(index=self.ped_df.index,
                                   columns=self.ped_df.columns, dtype=object)

        for pressure in self.ped_df.columns:
            for temp in self.ped_df.sort_index().index:
                ped_series = self.ped_df[pressure][temp].sort_index()
                ene = ped_series.index
                ped_step = ene[1]-ene[0]
                self.ene1_step = ped_step/3
                # ene1_vect_low = np.arange(
                # ene[0], self.ene1_step, -self.ene1_step)[1:]
                # ene1_vect_high = np.arange(
                # ene[0], max(ene), self.ene1_step)
                steps_to_zero = round((ene[0]-0)/self.ene1_step)
                ene1_vect_low = np.linspace(
                    ene[0],
                    ene[0]-steps_to_zero*self.ene1_step,
                    steps_to_zero+1)[1:-1]
                ene1_vect_high = np.linspace(
                    ene[0],
                    ene[-1],
                    num=round((ene[-1]-ene[0])/self.ene1_step)+1)
                # ^ includes enemax
                self.ene1_vect = np.sort(
                    np.concatenate((ene1_vect_low, ene1_vect_high)))

                # idx_ene_vect: indices in ene1_vect corresponding to values
                # of ene: ene[0]=self.ene1_vect[idx_ene_vect[0]]
                idx_ene_vect = np.arange(
                    len(ene1_vect_low), len(self.ene1_vect), 3)
                if distr_type == 'dos':
                    self.rho_rovib_prod1, self.rho_non1 = init_dos(
                        pressure, temp)

                prob_ene1_vect = []
                for idx_ene1, ene1 in enumerate(self.ene1_vect):
                    # indexes from self.ene1_vect: almost identical
                    # to energies in ene, but more consistent
                    idx_ene_new = idx_ene_vect[idx_ene_vect >= idx_ene1]
                    ene_new = self.ene1_vect[idx_ene_new]

                    if distr_type == 'phi':
                        prob_ene1ene = norm_distr(
                            ene1, ene_new, self.phi, self.ene_bw)
                    elif distr_type == 'dos':
                        prob_ene1ene = dos(idx_ene1, idx_ene_new)
                    prob_ene1ene_tot_pressure_ped = (
                        prob_ene1ene *
                        ped_series.values[idx_ene_vect >= idx_ene1]
                    )
                    prob_ene1 = np.trapz(
                        prob_ene1ene_tot_pressure_ped, ene_new)
                    prob_ene1_vect.append(prob_ene1)

                norm_factor_prob_ene1 = np.trapz(
                    prob_ene1_vect, x=self.ene1_vect)
                prob_ene1_norm = prob_ene1_vect/norm_factor_prob_ene1
                ped_df_prod[pressure][temp] = pd.Series(
                    prob_ene1_norm, index=self.ene1_vect)

                # remove comments to print P(E1)|T,P
                # prob_ene1_df = ped_df_prod[P][T].reset_index()
                # header_label = np.array(prob_ene1_df.columns, dtype=str)
                # header_label[0] = 'E [kcal/mol]'
                # labels = '\t\t'.join(header_label)
                # np.savetxt('PE1_{}_{}.txt'.format(P,T), prob_ene1_df.values,
                #         delimiter='\t', header=labels, fmt='%1.3e')
                # print(P, T, ped_df_prod[P][T].idxmax(), '\n')

        return ped_df_prod

    def equip_phi(self):
        """ Derive the energy distribution of 1 product from the
            energy equipartition theorem

            :return ped_df_prod: energy distribution of the product prod
            :rtype: dataframe(series(float))
        """

        dofs = self.get_dofs()
        vibdof_prod1, rotdof_prod1, vibdof_prod2, rotdof_prod2, _, _ = dofs

        phi_prod = (vibdof_prod1+rotdof_prod1/2) / \
            (vibdof_prod1+vibdof_prod2+(3+rotdof_prod1+rotdof_prod2)/2)
        print('Fraction of energy transferred to '
              f'products phi: {phi_prod:.2f}')
        self.phi = phi_prod
        ped_df_prod = self.prob_ene1_fct('phi')

        return ped_df_prod

    def beta_phi1a(self):
        """ Derive the energy distribution of 1 product from
            statistical model phi1a Danilack Goldsmith PROCI 2020

            :return ped_df_prod: energy distribution of the product prod
            :rtype: dataframe(series(float))
        """
        dofs = self.get_dofs()
        vibdof_prod1, _, _, _, vibdof_ts, _ = dofs

        # derive the energy fraction phi
        phi1a = vibdof_prod1/vibdof_ts
        print('Fraction of energy transferred to '
              f'products phi1a: {phi1a:.2f}')

        # rescale all energies with beta: allocate values in new dataframe
        self.phi = phi1a
        ped_df_prod = self.prob_ene1_fct('phi')

        return ped_df_prod

    def beta_phi2a(self):
        """ Derive the energy distribution of 1 product from
            statistical model phi2a Danilack Goldsmith PROCI 2020

            :return ped_df_prod: energy distribution of the product prod
            :rtype: dataframe(series(float))
        """

        dofs = self.get_dofs()
        vibdof_prod1, rotdof_prod1, _, _, vibdof_ts, rotdof_ts = dofs

        # derive the energy fraction phi
        phi2a = (vibdof_prod1+(3+rotdof_prod1)/2)/(vibdof_ts+(3+rotdof_ts)/2)
        print('Fraction of energy transferred to '
              f'products phi2a: {phi2a:.2f}')

        self.phi = phi2a
        ped_df_prod = self.prob_ene1_fct('phi')

        return ped_df_prod

    def beta_phi3a(self):
        """ Derive the energy distribution of 1 product from
            statistical model phi2a Danilack Goldsmith PROCI 2020

            :return ped_df_prod: energy distribution of the product prod
            :rtype: dataframe(series(float))
        """
        # derive the energy fraction phi
        dofs = self.get_dofs()
        vibdof_prod1, rotdof_prod1, vibdof_prod2, rotdof_prod2, _, _ = dofs

        # derive the energy fraction phi
        phi3a = (vibdof_prod1+(3+rotdof_prod1)/2) / \
            (vibdof_prod1+vibdof_prod2+(3+3+rotdof_prod1+rotdof_prod2)/2+3)
        print('Fraction of energy transferred to '
              f'products phi3a: {phi3a:.2f}')

        self.phi = phi3a
        ped_df_prod = self.prob_ene1_fct('phi')

        return ped_df_prod

    def rovib_dos(self):
        """ Derive the energy distribution of 1 product from the
            convolution of the density of states

            :return ped_df_prod: probability distribution of energy of prod1
            :rtype: dataframe(series(float))
        """
        # checks on input
        try:
            self.dos_df.empty
        except AttributeError:
            print('*Error: dos not defined, exiting now \n')
            sys.exit()

        if (
            self.prod1 not in self.dos_df.columns or
            self.prod2 not in self.dos_df.columns
        ):
            print('*Error: rovibrational density of states unavailable '
                  'for prod1/prod2. add "Fragment" next to frag name')
            sys.exit()

        # preallocations
        ped_df_prod = pd.DataFrame(index=self.ped_df.index,
                                   columns=self.ped_df.columns, dtype=object)

        ene_dos0 = self.dos_df.index
        self.ene_dos0 = ene_dos0

        # dos functions for prod1, prod2: more convenient because
        # many values are duplicate
        self.f_rho_rovib_prod1 = interp1d(
            ene_dos0, self.dos_df[self.prod1][ene_dos0].values, kind='cubic',
            fill_value='extrapolate')
        self.f_rho_rovib_prod2 = interp1d(
            ene_dos0, self.dos_df[self.prod2][ene_dos0].values, kind='cubic',
            fill_value='extrapolate')

        ped_df_prod = self.prob_ene1_fct('dos')

        return ped_df_prod