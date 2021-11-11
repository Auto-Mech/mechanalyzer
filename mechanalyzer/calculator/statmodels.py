"""
    Compute product energy distribution according to different statistical models
"""
import sys
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
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
        :return dof_info: dataframe with vibrat/rot degrees of freedom and molecular weight
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
        N_atoms = int(info[where_geom].strip().split()[1])
        atoms_ts += N_atoms

        key = info[where_name].strip().split()[1]
        keys.append(key)
        try:
            where_freq = find.where_in('Frequencies', info)[0]
            N_dof = int(info[where_freq].strip().split()[1]) + len(where_hind)
            if 3*N_atoms - N_dof == 6:
                rot_dof = 3
            else:
                rot_dof = 2
        except IndexError:
            # if 1 atom only: no 'Frequencies', set to 0
            N_dof = 0
            rot_dof = 0
        # this allows to get 3N-5 or 3N-6 without analyzing the geometry
        info_array[i, 0] = N_dof
        info_array[i, 1] = rot_dof

        # MW from type of atoms:
        geom_in = where_geom+1
        geom_fin = geom_in+N_atoms
        atoms_array = np.array([geomline.strip().split()[0]
                                for geomline in info[geom_in:geom_fin]])

        mw = np.sum(np.array([MW_dct_elements[at]
                              for at in atoms_array], dtype=float))
        info_array[i, 2] = mw

    # if ask for ts: assume the first 2 blocks are 2 reactants of a bimol reaction
    # and derive the DOFs of the TS
    if ask_for_ts == True:
        keys.append('TS')

        # assume there are no linear TSs
        info_array[2, :] = [3*atoms_ts - 7, 3,
                            info_array[0, 2]+info_array[1, 2]]

    dof_info = pd.DataFrame(info_array, index=keys, columns=[
                            'vib dof', 'rot dof', 'mw'])

    return dof_info


def dos_trasl(m1, m2, E_grid, P, T):
    """ Compute the translational density of states per unit volume
        m1, m2: MW in kg/mol
        E_grid: energy grid in kcal/mol (array)
        P, T: pressure and temperature in SI units [Pa, K]
        :return dos_tr_series: dos in mol/kcal
        :rtype: array
    """
    # conversions
    NAVO = 6.022e+23  # molec/mol
    E_grid_J = E_grid*4184/NAVO  # kcal/mol*J/kcal/NAVO=J
    m1 /= NAVO  # kg/mol/(molecule/mol)
    m2 /= NAVO
    red_mass = (m1 * m2) / (m1 + m2)    # kg
    h = 6.626e-34   # Planck's constant, J*s
    # 1 molecule
    # rhotr/V = pi/4*(8m/h2)^3/2*e^1/2
    rho = (np.pi/4*np.power(8*red_mass/h**2, 3/2) *
           np.power(E_grid_J, 1/2))  # unit: 1/m3/J

    # consider the molar volume RT/P [m3/1]
    V_mol = 1.38e-23*T/P  # m3/1 for ideal gases
    # 1/m3*(m3/1)/J*(J/kcal)*(mol) = mol/kcal
    rho_kcal_mol = rho*V_mol*4184/NAVO

    return rho_kcal_mol


class ped_models:

    def __init__(self, ped_df, hotfrg, otherfrg, dos_df=None, dof_info=None, E_BW=None):
        """ initialize variables
            :param ped_df: dataframe(columns:P, rows:T) with the Series of energy distrib
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

        self.E_BW = E_BW
        self.prod1 = hotfrg
        self.prod2 = otherfrg
        self.dof_info = dof_info

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
            return ped_df_prod

        except KeyError:
            print('*Error: model not available. Please select among \n {} \n'.format(
                '\n'.join(self.models_dct.keys())))
            exit()

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
            print('incomplete degrees of freedom info: species are \n {} \
              \n while dof info is available for \n  {} \n Exiting ...'.format(
                ' '.join([self.prod1, self.prod2]), ' '.join(self.vibdof_dct.keys())))
        return vibdof_prod1, rotdof_prod1, vibdof_prod2, rotdof_prod2, vibdof_ts, rotdof_ts

    def equip_simple(self):
        """ Derive the energy distribution of 1 product from the energy equipartition theorem

            :return ped_df_prod: energy distribution of the product prod
            :rtype: dataframe(series(float))
        """
        vibdof_prod1, rotdof_prod1, vibdof_prod2, rotdof_prod2, _, _ = self.get_dofs()

        beta_prod = (vibdof_prod1+rotdof_prod1/2) / \
            (vibdof_prod1+vibdof_prod2+(3+rotdof_prod1+rotdof_prod2)/2)
        # 3/2: 1/2kbT for each rotation, no trasl (ts trasl energy preserved)
        # 9/2: 1/2kbT*6 rotational dofs for products, +3 for relative trasl
        print(
            'fraction of energy transferred to products: {:.2f}'.format(beta_prod))
        # rescale all energies with beta: allocate values in new dataframe
        ped_df_prod = pd.DataFrame(index=self.ped_df.index,
                                   columns=self.ped_df.columns, dtype=object)
        for P in self.ped_df.columns:
            for T in self.ped_df.sort_index().index:
                idx_new = self.ped_df[P][T].index * beta_prod
                norm_factor = np.trapz(self.ped_df[P][T].values, x=idx_new)
                vals = self.ped_df[P][T].values/norm_factor
                ped_df_prod[P][T] = pd.Series(vals, index=idx_new)

        return ped_df_prod

    def P_E1_fct(self, distr_type):
        """ Derive the energy distribution of 1 product from one
            of the statistical models by Danilack Goldsmith PROCI 2020
            phi is the average fraction of energy transferred to the products

            calculate P(E') = probability of fragment 1 to have energy E'
            a. alfa(E) = PED(E)
            b. select E' and calculate P(E',E) according to normal distribution
            c. integrate over dE: int(P(E',E)*PED(E)dE) = P(E')
        """

        def norm_distr(E1, E, phi, E_BW):
            """ P(E1; E) = exp(-(E1-phi*E)^2/(2^0.5*sigma(E_BW)))/((2*pi)^0.5*sigma(E_BW))
                mi = phi*E
                sigma = f(E_BW, phi*E)
                E1 is a number
                E is a vector
            """

            mi = np.array(phi*E, dtype=float)
            # correlation from Danilack Goldsmith - I add the additional fraction of energy transferred to the products
            sigma = np.array(0.87+0.04*(E_BW+phi*(E-E_BW)), dtype=float)
            num = np.exp(-((E1-mi)/(2**0.5)/sigma)**2)
            den = np.power(2*np.pi, 0.5)*sigma

            P_E1E = num/den

            return P_E1E

        def init_dos(P, T):
            """ initialize variables for DOS calculation
            """
            rho_rovib_prod1 = self.f_rho_rovib_prod1(self.E1_vect)
            E1_vect_w0 = np.concatenate((np.array([0]), self.E1_vect))
            rho_rovib_prod2 = self.f_rho_rovib_prod2(E1_vect_w0)
            rho_trasl = dos_trasl(self.mw_dct[self.prod1],
                self.mw_dct[self.prod2], E1_vect_w0, P*101325, T)     

            # calculate rho_non1(E1_vect)
            rho_non1 = []
            for idx in self.E1_vect:
                # first iter should be 0
                # the sum of the energies in rhovib_prod2 and rho_trasl is always E1
                idx_Eint = np.arange(0, idx+1, dtype=int)
                idx_EminusEint = idx_Eint[::-1]
                rho_non1_integrand = rho_rovib_prod2[idx_Eint] * rho_trasl[idx_EminusEint]
                rho_non1.append(np.trapz(rho_non1_integrand, x = self.E1_vect[idx_Eint]))

            rho_non1 = np.array(rho_non1)

            return rho_rovib_prod1, rho_non1

        def dos(idx_E1, idx_E_new):
            
            P_E1E = []
            # loops over total energy: loop over the indexes where you find Etot in E1vect
            for idx_E in idx_E_new: 
                
                if idx_E1 == idx_E:
                    # rho1(E1)*rhonon1(E-E1) = 0
                    P_E1E.append(0)
                    continue
                idx_E1_array = np.arange(0, idx_E) # index of E1<E (fixed E)
                idx_E_minus_E1_array = idx_E1_array[::-1] # index of E-E1 (fixed E)
                rho1_E1 = self.rho_rovib_prod1[idx_E1] # rho1(E1)
                rho_non1 = self.rho_non1[idx_E_minus_E1_array[idx_E1]] # rhonon1(E-E1)
                rho1_E1_array = self.rho_rovib_prod1[idx_E1_array] # rho1(E1) with E1<E (fixed E)
                rho_non1_array = self.rho_non1[idx_E_minus_E1_array] # rhonon1(E-E1) with E1<E (fixed E)
                num = rho1_E1 * rho_non1
                den = np.trapz(rho1_E1_array * rho_non1_array, x = self.E1_vect[idx_E1_array])
                P_E1E.append(num/den)

            P_E1E = np.array(P_E1E)
            # print(P_E1E)
            return P_E1E

        # preallocations
        ped_df_prod = pd.DataFrame(index=self.ped_df.index,
                                   columns=self.ped_df.columns, dtype=object)

        for P in self.ped_df.columns:
            for T in self.ped_df.sort_index().index:
                ped_series = self.ped_df[P][T].sort_index()
                E = ped_series.index
                ped_step = E[1]-E[0]
                self.E1_step = ped_step/3
                #E1_vect_low = np.arange(E[0], self.E1_step, -self.E1_step)[1:]
                #E1_vect_high = np.arange(E[0], max(E), self.E1_step)
                steps_to_zero = round((E[0]-0)/self.E1_step)
                E1_vect_low = np.linspace(E[0], E[0]-steps_to_zero*self.E1_step, steps_to_zero+1)[1:-1]
                E1_vect_high = np.linspace(E[0], E[-1], num=round((E[-1]-E[0])/self.E1_step)+1) # includes Emax
                self.E1_vect = np.sort(np.concatenate((E1_vect_low, E1_vect_high)))
                
                # idx_E_vect: indices in E1_vect corresponding to values of E: E[0]=self.E1_vect[idx_E_vect[0]]
                idx_E_vect = np.arange(len(E1_vect_low), len(self.E1_vect), 3) # 
                if distr_type == 'dos':
                    self.rho_rovib_prod1, self.rho_non1 = init_dos(P, T)

                P_E1_vect = []
                for idx_E1, E1 in enumerate(self.E1_vect):
                    #indexes from self.E1_vect: almost identical to energies in E, but more consistent
                    idx_E_new = idx_E_vect[idx_E_vect >= idx_E1]
                    E_new = self.E1_vect[idx_E_new]
                    
                    if distr_type == 'phi':
                        P_E1E = norm_distr(E1, E_new, self.phi, self.E_BW)
                    elif distr_type == 'dos':
                        P_E1E = dos(idx_E1, idx_E_new)
                    P_E1Etot_P_ped = P_E1E*ped_series.values[idx_E_vect >= idx_E1]
                    P_E1 = np.trapz(P_E1Etot_P_ped, E_new)
                    P_E1_vect.append(P_E1)

                norm_factor_P_E1 = np.trapz(
                    P_E1_vect, x=self.E1_vect)
                P_E1_norm = P_E1_vect/norm_factor_P_E1
                ped_df_prod[P][T] = pd.Series(P_E1_norm, index=self.E1_vect)

                # remove comments to print P(E1)|T,P
                #P_E1_df = ped_df_prod[P][T].reset_index()
                #header_label = np.array(P_E1_df.columns, dtype=str)
                #header_label[0] = 'E [kcal/mol]'
                #labels = '\t\t'.join(header_label)
                #np.savetxt('PE1_{}_{}.txt'.format(P,T), P_E1_df.values,
                #        delimiter='\t', header=labels, fmt='%1.3e')
                #print(P, T, ped_df_prod[P][T].idxmax(), '\n')

        return ped_df_prod

    def equip_phi(self):
        """ Derive the energy distribution of 1 product from the energy equipartition theorem

            :return ped_df_prod: energy distribution of the product prod
            :rtype: dataframe(series(float))
        """

        vibdof_prod1, rotdof_prod1, vibdof_prod2, rotdof_prod2, _, _ = self.get_dofs()

        phi_prod = (vibdof_prod1+rotdof_prod1/2) / \
            (vibdof_prod1+vibdof_prod2+(3+rotdof_prod1+rotdof_prod2)/2)
        print(
            'fraction of energy transferred to products phi: {:.2f}'.format(phi_prod))
        self.phi = phi_prod
        ped_df_prod = self.P_E1_fct('phi')

        return ped_df_prod

    def beta_phi1a(self):
        """ Derive the energy distribution of 1 product from 
            statistical model phi1a Danilack Goldsmith PROCI 2020

            :return ped_df_prod: energy distribution of the product prod
            :rtype: dataframe(series(float))
        """
        vibdof_prod1, _, _, _, vibdof_ts, _ = self.get_dofs()
        # derive the energy fraction phi
        phi1a = vibdof_prod1/vibdof_ts
        print(
            'fraction of energy transferred to products phi1a: {:.2f}'.format(phi1a))
        # rescale all energies with beta: allocate values in new dataframe
        self.phi = phi1a
        ped_df_prod = self.P_E1_fct('phi')


        return ped_df_prod

    def beta_phi2a(self):
        """ Derive the energy distribution of 1 product from 
            statistical model phi2a Danilack Goldsmith PROCI 2020

            :return ped_df_prod: energy distribution of the product prod
            :rtype: dataframe(series(float))
        """
        vibdof_prod1, rotdof_prod1, _, _, vibdof_ts, rotdof_ts = self.get_dofs()
        # derive the energy fraction phi
        phi2a = (vibdof_prod1+(3+rotdof_prod1)/2)/(vibdof_ts+(3+rotdof_ts)/2)
        print(
            'fraction of energy transferred to products phi2a: {:.2f}'.format(phi2a))
        self.phi = phi2a
        ped_df_prod = self.P_E1_fct('phi')

        return ped_df_prod

    def beta_phi3a(self):
        """ Derive the energy distribution of 1 product from 
            statistical model phi2a Danilack Goldsmith PROCI 2020

            :return ped_df_prod: energy distribution of the product prod
            :rtype: dataframe(series(float))
        """
        # derive the energy fraction phi
        vibdof_prod1, rotdof_prod1, vibdof_prod2, rotdof_prod2, _, _ = self.get_dofs()
        phi3a = (vibdof_prod1+(3+rotdof_prod1)/2) / \
            (vibdof_prod1+vibdof_prod2+(3+3+rotdof_prod1+rotdof_prod2)/2+3)
        print(
            'fraction of energy transferred to products phi3a: {:.2f}'.format(phi3a))
        self.phi = phi3a
        ped_df_prod = self.P_E1_fct('phi')

        return ped_df_prod

    def rovib_dos(self):
        """ Derive the energy distribution of 1 product from the 
            convolution of the density of states

            :return ped_df_prod: probability distribution of the energy of prod1
            :rtype: dataframe(series(float))
        """
        # checks on input
        try:
            self.dos_df.empty
        except AttributeError:
            print('*Error: dos not defined, exiting now \n')
            sys.exit()

        if self.prod1 not in self.dos_df.columns or self.prod2 not in self.dos_df.columns:
            print('*Error: rovibrational density of states unavailable for prod1/prod2. add "Fragment" next to frag name')
            sys.exit()

        # preallocations
        ped_df_prod = pd.DataFrame(index=self.ped_df.index,
                                   columns=self.ped_df.columns, dtype=object)

        E_dos0 = self.dos_df.index
        self.E_dos0 = E_dos0
        # dos functions for prod1, prod2: more convenient because many values are duplicate
        self.f_rho_rovib_prod1 = interp1d(E_dos0, self.dos_df[self.prod1][E_dos0].values, kind='cubic',
                                     fill_value='extrapolate')
        self.f_rho_rovib_prod2 = interp1d(E_dos0, self.dos_df[self.prod2][E_dos0].values, kind='cubic',
                                     fill_value='extrapolate')

        ped_df_prod = self.P_E1_fct('dos')

        return ped_df_prod

