""" Extract branching fractions of the products from
    hotenergies and pedspecies
    Fit and write the rates
"""

import sys
from tokenize import Single
from tracemalloc import start
from scipy.interpolate import interp1d
import numpy as np
import pandas as pd


def ped_df_rescale_test(starthot_df, energy_scale, save=False):
    for temp in starthot_df.index:
        for pressure in starthot_df.columns:
            vals = starthot_df[pressure][temp].values
            dfnew_index = starthot_df[pressure][temp].index + energy_scale
            starthot_df[pressure][temp] = pd.Series(vals, index=dfnew_index)
            starthot_df[pressure][temp] = starthot_df[pressure][temp][starthot_df[pressure][temp].index > 0]
            if pressure == 1 and temp in [500, 1000, 2000] and save == True:
                tosave = starthot_df[pressure][temp].reset_index()
                np.savetxt('hotdf_shift_{}_{}.txt'.format(
                    pressure, temp), tosave, fmt='%1.3e')

    return starthot_df


def ped_df_rescale(starthot_df, ped_df_fromhot):
    """ obtain a new energy distribution for ped_df_fromhot
        based on the energy distribution of hot_df

        :param starthot_df: dataframe[P][T] with the
            Series of energy distrib [en: prob(en)]
        :type starthot_df: dataframe(series(float))
        :param ped_df_fromhot: (dataframe(columns:P, rows:T))
                            with series(hoten: series([en: prob(en)]))
        :type ped_df_fromhot: df[P][T]:series[energies: df[allspecies][energies]: prob]}
        :return ped_df: ped_df weighted on hot_df distribution
        :rtype: df[P][T]:series(en: prob(en))
        sum(ped_df_fromhot(E';E)*starthot_df(E)), E is the starting hoten, E' is the prod en
    """
    # sort indexes
    starthot_df = starthot_df.sort_index()
    ped_df_fromhot = ped_df_fromhot.sort_index()
    starthot_df, ped_df_fromhot = checks_temp_pressure_and_extend(
        starthot_df, ped_df_fromhot)
    temps, pressures = [ped_df_fromhot.index, ped_df_fromhot.columns]

    ped_df = pd.DataFrame(index=temps, columns=pressures, dtype=object)
    T_del = []
    for temp in temps:
        for pressure in pressures:

            # initial distribution: sort and fit
            starthot = starthot_df[pressure][temp].sort_index()
            # rescale values based on probability - too low probability excluded
            starthot = starthot[starthot > max(starthot)*1e-4]  # 99.99%
            # refit starthot to derive the weight factors later

            f_starthot = interp1d(
                starthot.index, starthot.values, bounds_error=False,
                kind='cubic', fill_value=(starthot.values[0], starthot.values[-1]))

            # reduce the energy range of ped_fromhot

            ped_fromhot = ped_df_fromhot[pressure][temp].sort_index()
            # print('before: ', ped_fromhot, '\n')
            ped_fromhot = ped_fromhot.iloc[(
                starthot.index[0] <= ped_fromhot.index)*(ped_fromhot.index <= starthot.index[-1])]
            # print('after: ', ped_fromhot, '\n')
            # set new energy vector from min and max energies in peden_fromhot
            min_en_fromhot = min([min(ped_fromhot.iloc[i].index)
                                  for i in np.arange(0, len(ped_fromhot))])
            max_en_fromhot = min([max(ped_fromhot.iloc[i].index)
                                  for i in np.arange(0, len(ped_fromhot))])
            ene_vect = np.arange(min_en_fromhot, max_en_fromhot, 0.5)
            prob_vect = np.zeros(ene_vect.shape)
            # if temp == 1300 and pressure == 0.01:
            #    print(temp, pressure, starthot, '\n')
            #    print(ped_fromhot_0, '\n', ped_fromhot)

            # weight factor from fitted starthot
            for starten in ped_fromhot.index:
                weightfactor = f_starthot(starten)
                hoten = ped_fromhot[starten].index
                pedhot = ped_fromhot[starten].values
                # interpolate values
                if len(hoten) > 3:
                    f_ped_fromhot = interp1d(
                        hoten, pedhot, bounds_error=False,
                        kind='cubic', fill_value=(0., 0.))
                    prob_vect += f_ped_fromhot(ene_vect)*weightfactor

                elif 1 >= len(hoten) >= 3:
                    # find max val of hoten and set that one
                    idx_max = np.argmax(pedhot)
                    # find where hoten is closest and add 1
                    idx = np.argmin(ene_vect - hoten[idx_max])
                    single_1 = np.zeros(ene_vect.shape)
                    single_1[idx] = pedhot[idx_max]
                    prob_vect += single_1*weightfactor

            # renormalize and put in dataframe
            prob_vect /= np.trapz(prob_vect, x=ene_vect)
            ped_df[pressure][temp] = pd.Series(prob_vect, index=ene_vect)
            if ped_df[pressure][temp].empty:
                T_del.append(temp)
            # if temp == 1300 and pressure == 0.01:
            #    print(ped_df[pressure][temp], '\n')

    ped_df = ped_df.drop(index=list(set(T_del)))

    return ped_df


def bf_tp_df_full(ped_df, hotbf_df):
    """ Build a branching fractions dataframe as a
        function of temprature and pressure containing the BFs of
        each product of the PES

        :param ped_df: dataframe[P][T] with the
            Series of energy distrib [en: prob(en)]
        :type ped_df: dataframe(series(float))
        :param hotbf_df: hot branching fractions for hotspecies
        :type hotbf_df: df[P][T]:df[allspecies][energies]}
        :return bf_tp_df: branching fractions at T,P
            for each product for the selected hotspecies
        :rtype: dataframe of series df[P][T]:series[species],
                dataframe(series(float))
    """

    # sort indexes
    ped_df = ped_df.sort_index()
    hotbf_df = hotbf_df.sort_index()
    ped_df, hotbf_df = checks_temp_pressure_and_extend(ped_df, hotbf_df)
    # compute branching fractions
    # derive keys for dicts and complete set of T,P (should be identical)
    # temp_ped contains the desired T range; temp_hot may have a larger range
    temps, pressures = [ped_df.index, ped_df.columns]
    allspecies = hotbf_df[pressures[0]][temps[0]].columns  # extract species
    # for each T, P: compute BF
    bf_tp_df = pd.DataFrame(index=temps, columns=pressures, dtype=object)
    for temp in temps:
        for pressure in pressures:
            # extract ped and hoten by increasing index

            ped = ped_df[pressure][temp].sort_index()
            hoten = hotbf_df[pressure][temp].sort_index().index

            # reduce the energy range of hoten and ped
            hoten = hoten[(ped.index[0] <= hoten)*(hoten <= ped.index[-1])]
            # ped = ped[(min_en <= ped.index)*(ped.index <= max_en)]
            if len(hoten) > 3:
                ene_vect = ped.index
                ped_vect = ped.values
                # Series to allocate
                bf_series = pd.Series(0, index=allspecies, dtype=float)
                for spc in allspecies:
                    hoten_spc = hotbf_df[pressure][temp][spc][hoten]
                    f_hoten = interp1d(
                        hoten_spc.index, hoten_spc.values, bounds_error=False,
                        kind='cubic', fill_value=(hoten_spc.values[0], hoten_spc.values[-1]))
                    hoten_vect = f_hoten(ene_vect)
                    # recompute in an appropriate range
                    bf_series[spc] = np.trapz(ped_vect*hoten_vect, x=ene_vect)
                # renormalize for all species and put in dataframe
                bf_tp_df[pressure][temp] = bf_series/np.sum(bf_series.values)

        # if any nan: delete column
        if any(bf_tp_df.loc[temp].isnull()):
            bf_tp_df = bf_tp_df.drop(index=[temp])

    return bf_tp_df


def bf_tp_df_todct(bf_tp_df, bf_threshold, savefile=False, rxn='', model=''):
    """ Converts the dataframe of hot branching fractions to dictionary and
        excludes invalid BFs

        :param bf_tp_df: branching fractions at T,P for each product
            for the selected hotspecies
        :type bf_tp_df: dataframe of series df[P][T]:series[species],
            dataframe(series(float))
        :param bf_threshold: threshold to filter out all species with lower BFs
        :type bf_threshold: float
        :param rxn: name of the reaction producing hot species - only for file saving
        :type rxn: str
        :param model: type of model - only required for filename
        :type model: str
        :param savefile: save file with branching fractions?
        :type savefile: bool
        :return bf_tp_dct: branching fractions at T,P for each product
        :rtype: dct{species: {pressure: (array(T), array(BF))}} -
            same type as ktp dct
    """
    # get species
    temps, pressures = bf_tp_df.index, bf_tp_df.columns
    allspecies = bf_tp_df.iloc[0, 0].index
    bf_tp_dct_out = {}
    # fill the dictionary:
    for spc in allspecies:
        num_data_highenough = 0
        bf_tp_dct_i = {}
        bf_df_sp_i = pd.DataFrame(
            np.zeros((len(temps), len(pressures))),
            index=temps, columns=pressures)

        for pressure in pressures:
            bf_temp = []
            temp_new = []
            for temp in temps:
                bfrac = bf_tp_df[pressure][temp][spc]
                if 1 > bfrac >= 1e-30:
                    # avoid too small values and 1 to avoid discontinuities
                    temp_new.append(temp)
                    bf_temp.append(bfrac)
                    bf_df_sp_i[pressure][temp] = bfrac
                # check if values is high enough
                if bfrac >= bf_threshold:
                    num_data_highenough += 1
                # condition on BF and temperature
            if bf_temp:
                bf_tp_dct_i[pressure] = (np.array(temp_new), np.array(bf_temp))

        if num_data_highenough > 0:
            bf_tp_dct_out[spc] = bf_tp_dct_i

        # write file with the BFs
        if num_data_highenough > 0 and savefile:
            bf_df_sp_i = bf_df_sp_i.reset_index()
            header_label = np.array(bf_df_sp_i.columns, dtype=str)
            header_label[0] = 'T [K]'
            labels = '\t\t'.join(header_label)
            np.savetxt(f'bf_{model}_{rxn}_{spc}.txt', bf_df_sp_i.values,
                       delimiter='\t', header=labels, fmt='%1.2e')

    return bf_tp_dct_out


def merge_bf_rates(bf_tp_dct, ktp_dct):
    """ Read the branching fractions of the products and
        derive final rate constants

        :param bf_tp_dct: branching fractions at T,P for each product
            for the selected hotspecies
        :type: dct{species: {pressure: (array(T), array(BF))}}}
        :param ktp_dct: set of rates for the "total" rate constant
        :type ktp_dct: dictionary {P: (T, k)}
        :return new_ktp_dct: rates for each species multiplied by bf: bf_i*k
        :rtype: dct {species: {P: (T, k)}}
    """
    # drop "high" from ktp_dct
    if 'high' in ktp_dct.keys():
        ktp_dct.pop('high')

    # get all species and preallocate new dct
    all_species = bf_tp_dct.keys()
    new_ktp_dct = dict.fromkeys(all_species)

    # loop over all species
    for spc, bf_tp_dct_sp in bf_tp_dct.items():

        pressure_all = list(bf_tp_dct_sp.keys())
        pressure_ktp = list(ktp_dct.keys())
        # extend ktp_dct to the full range of pressures
        if not all(pressure in pressure_ktp for pressure in pressure_all):
            print('Warning: P range of ktp dictionary extended to match BF:\n')
            print(' Values at other pressures approximated from available \n')
            ktp_dct = extend_dct_with_pressure(
                pressure_all, pressure_ktp, ktp_dct)

        if not all(pressure in pressure_all for pressure in pressure_ktp):
            print('Warning: P range of BF dictionary extended to match ktp:\n')
            print(' Values at other pressures approximated from available \n')
            bf_tp_dct_sp = extend_dct_with_pressure(
                pressure_ktp, pressure_all, bf_tp_dct_sp)

        pressure_all = list(set(pressure_all + pressure_ktp))
        bf_ktp_dct = dict.fromkeys(pressure_all)

        for pressure in pressure_all:
            temp_bf = bf_tp_dct_sp[pressure][0]
            temp_ktp = np.array(ktp_dct[pressure][0])
            # choose the smallest vector
            temp_common = temp_bf[
                [temp_bf_i in temp_ktp for temp_bf_i in temp_bf]]
            # reconstruct appropriate datasets
            bf_series = pd.Series(bf_tp_dct_sp[pressure][1],
                                  index=temp_bf)
            ktp_series = pd.Series(np.array(ktp_dct[pressure][1]),
                                   index=temp_ktp)
            # reconstruct k and allocate
            kvals = (bf_series[temp_common].values *
                     ktp_series[temp_common].values)
            bf_ktp_dct[pressure] = (temp_common, kvals)
        # allocate the new dictionary
        new_ktp_dct[spc] = bf_ktp_dct

    return new_ktp_dct

################# useful functions  ###########################################


def checks_temp_pressure_and_extend(ped_df, hotbf_df):
    """ compare T,P """
    temp_ped, pressure_ped = [ped_df.index, ped_df.columns]
    temp_hot, pressure_hot = [hotbf_df.index, hotbf_df.columns]
    # check that T of temp_hot are at least as many as those of temp_ped
    # if they're not, it doesn't make sense to continue
    if not all(temp_ped_i in temp_hot for temp_ped_i in temp_ped):
        print('*Error: temperature range in HOTenergies '
              'does not cover the full range')
        sys.exit()

    # if in pressure_ped not all pressures of hoten are available:
    # extend the range of pressures
    # ex.: for H abstractions, they will be pressure independent
    # but probably the successive decomposition is not
    if not all(pressure_hot_i in pressure_ped
               for pressure_hot_i in pressure_hot):
        print('*Warning: P range of PedOutput smaller than HOTenergies:\n')
        print('Energy distribution at other pressure '
              'approximated from available values \n')
        ped_df = extend_dct_with_pressure(
            pressure_hot, pressure_ped, ped_df)

    # check that pressures of pressure_ped are contained in hot energies
    # if they are not: extend pressure range assuming behavior is the same
    # ELIF no good here - might be that pressure_ped and pressure_hot
    #   only intersect in a limited range
    if not all(pressure_ped_i in pressure_hot
               for pressure_ped_i in pressure_ped):
        print('Warning: P range of HOTenergies smaller than PEDoutput:')
        print('Energy distribution at other pressures '
              'approximated from available values \n')
        hotbf_df = extend_dct_with_pressure(
            pressure_ped, pressure_hot, hotbf_df)

    return ped_df, hotbf_df


def extend_dct_with_pressure(pressure_all, pressure_red, dct_toextend):
    """ Takes a dictionary(dataframe) with pressure keys(columns) and extends it:
        approximate the values at missing pressures with available values

        :param pressure_all: all pressures
        :type pressure_all: list
        :param dct_toextend: dictionary (dataframe) w/ pressure_reduced as keys
        :type dct_toextend: dictionary (dataframe)
        :return dct_extended: dictionary with pressure_all as keys (columns)
        :rtype: dictionary (dataframe)
    """

    for pressure in pressure_all:
        if pressure not in pressure_red:
            # approximate pressure:
            # provides the minimum difference with pressure_ped_i
            pressure_approx = pressure_red[
                np.argmin(
                    [np.log(abs(pressure_red_i-pressure))
                     for pressure_red_i in pressure_red])]
            # extend original dataframe
            dct_toextend[pressure] = dct_toextend[pressure_approx]

    return dct_toextend
