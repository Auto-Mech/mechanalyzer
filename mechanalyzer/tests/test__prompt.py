""" test mechanalyzer.calculator.statmodels
    test mechanalyzer.calculator.bf
    called by
    mechanalyzer.builder.bf
    mechanalyzer.builder.ped
    similar structure to script in mechanalyzer_bin
"""

import os
import numpy as np
from autofile.io_ import read_file
import mechanalyzer
import mess_io
import autoparse.pattern as app
from ioformat import remove_comment_lines

PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH = os.path.join(PATH, 'test_nB', 'C3H8_OH')
PED_INP = read_file(os.path.join(INP_PATH, 'me_ktp_ped.inp'))
PED_INP = remove_comment_lines(PED_INP, delim_pattern=app.escape('!'))
PED_INP = remove_comment_lines(PED_INP, delim_pattern=app.escape('#'))
PED_OUT = read_file(os.path.join(INP_PATH, 'ped.out'))
KE_PED_OUT = read_file(os.path.join(INP_PATH, 'ke_ped.out'))
HOT_OUT = read_file(os.path.join(INP_PATH, 'me_ktp_hoten.log'))

T = [400.0, 600.0, 800.0, 1200.0, 1800.0, 2000.0]
P = [0.1, 1.0, 100.0]
PEDSPECIES = [['RH', 'NC3H7'], ['RH', 'IC3H7']]

ENERGY_DCT = {'W0': -6.0, 'RH': -2.2, 'NC3H7': -20.14, 'IC3H7': -23.18,
              'B0': -2.2, 'B1': 0.0, 'B2': -0.89}
{'B0': ('W0', 'RH'), 'B1': ('W0', 'NC3H7'), 'B2': ('W0', 'IC3H7')}

ene_bw_dct = {'RH->NC3H7': 17.94, 'RH->IC3H7': 20.98}

KTP_DCT = {'RH->NC3H7':
           {
               1.0: ((400.0, 600.0, 800.0, 1200.0, 1800.0, 2000.0),
                     (5.14854e-13, 1.49542e-12, 2.95517e-12, 7.33206e-12, 1.72203e-11, 2.12587e-11))
           },
           'RH->IC3H7':
           {
               1.0: ((400.0, 600.0, 800.0, 1200.0, 1800.0, 2000.0),
                     (1.8988e-13, 4.35041e-13, 7.7321e-13, 1.73144e-12, 3.78948e-12, 4.61016e-12))
           }}

LABELS = list(KTP_DCT.keys())

FRAGMENTS_DCT = {'RH->NC3H7': ['CH3CH2CH2',
                               'H2O'], 'RH->IC3H7': ['CH3CHCH3', 'H2O']}

FRAG_REACS = ['C3H8', 'OH']

HOT_FRAG_DCT = {'CH3CH2CH2': ['CH3CH2CH2'], 'CH3CHCH3': [
    'CH3CHCH3'], 'P1': ['C2H4', 'CH3'], 'P2': ['CH3CHCH2', 'H']}

HOTSPECIES = ['CH3CH2CH2', 'CH3CHCH3']


class test_prompt:

    def __init__(self):
        # HERE OPERATIONS TESTED IN OTHER FUNCTIONS - REQUIRED TO GET DATA
        # get dof info
        species_blocks_ped = mess_io.reader.get_species(PED_INP)
        self.DOF_DCT = {}
        for label in LABELS:
            prods = label.split('->')[1]
            # NB FCT TESTED IN TEST__CALC_STATMODELS
            self.DOF_DCT[label] = mechanalyzer.calculator.statmodels.get_dof_info(
                species_blocks_ped[prods], ask_for_ts=True)
        # GET PED
        self.PED_DCT = mess_io.reader.ped.get_ped(
            PED_OUT, PEDSPECIES, ENERGY_DCT)
        # GET DOS
        self.DOS_ROVIB = mess_io.reader.rates.dos_rovib(KE_PED_OUT)
        # HOTEN
        self.HOTEN_DCT = mess_io.reader.hoten.extract_hot_branching(
            HOT_OUT, HOTSPECIES, list(HOT_FRAG_DCT.keys()), T, P)

    # test different models
    def test_equip_simple(self):
        """ test statmodels.pedmodels.equip_simple
        """
        ped_df_frag1_dct = mechanalyzer.builder.ped.ped_frag1(
            self.PED_DCT['RH->NC3H7'], 'CH3CH2CH2', 'H2O', ['equip_simple'], dof_info=self.DOF_DCT['RH->NC3H7'], ene_bw=ene_bw_dct['RH->NC3H7'])
        ped_600 = ped_df_frag1_dct['equip_simple'][1.0][600]
        ped_1200 = ped_df_frag1_dct['equip_simple'][1.0][1200]
        assert np.isclose((ped_600.iloc[100]), 0.10436991291189263)
        assert np.isclose((ped_1200.iloc[100]), 0.020849788247301364)
        assert np.isclose(np.trapz(ped_600.values, x=ped_600.index), 1)
        assert np.isclose(np.trapz(ped_1200.values, x=ped_1200.index), 1)

    def test_equip_phi(self):
        """ test statmodels.pedmodels.equip_phi
        """
        ped_df_frag1_dct = mechanalyzer.builder.ped.ped_frag1(
            self.PED_DCT['RH->NC3H7'], 'CH3CH2CH2', 'H2O', ['equip_phi'], dof_info=self.DOF_DCT['RH->NC3H7'], ene_bw=ene_bw_dct['RH->NC3H7'])
        ped_600 = ped_df_frag1_dct['equip_phi'][1.0][600]
        ped_1200 = ped_df_frag1_dct['equip_phi'][1.0][1200]
        assert np.isclose((ped_600.iloc[500]), 0.046399200133238984)
        assert np.isclose((ped_1200.iloc[500]), 0.037261506158156384)
        assert np.isclose(np.trapz(ped_600.values, x=ped_600.index), 1)
        assert np.isclose(np.trapz(ped_1200.values, x=ped_1200.index), 1)

    def test_beta_phi1a(self):
        """ test statmodels.pedmodels.beta_phi1a
        """
        ped_df_frag1_dct = mechanalyzer.builder.ped.ped_frag1(
            self.PED_DCT['RH->IC3H7'], 'CH3CHCH3', 'H2O', ['beta_phi1a'], dof_info=self.DOF_DCT['RH->IC3H7'], ene_bw=ene_bw_dct['RH->IC3H7'])
        ped_400 = ped_df_frag1_dct['beta_phi1a'][1.0][400]
        ped_800 = ped_df_frag1_dct['beta_phi1a'][1.0][800]
        assert np.isclose((ped_400.iloc[500]), 0.0009891036253060362)
        assert np.isclose((ped_800.iloc[500]), 0.059832564512289556)
        assert np.isclose(np.trapz(ped_400.values, x=ped_400.index), 1)
        assert np.isclose(np.trapz(ped_800.values, x=ped_800.index), 1)

    def test_beta_phi2a(self):
        """ test statmodels.pedmodels.beta_phi2a
        """
        ped_df_frag1_dct = mechanalyzer.builder.ped.ped_frag1(
            self.PED_DCT['RH->IC3H7'], 'CH3CHCH3', 'H2O', ['beta_phi2a'], dof_info=self.DOF_DCT['RH->IC3H7'], ene_bw=ene_bw_dct['RH->IC3H7'])
        ped_400 = ped_df_frag1_dct['beta_phi2a'][1.0][400]
        ped_800 = ped_df_frag1_dct['beta_phi2a'][1.0][800]
        assert np.isclose((ped_400.iloc[500]), 0.0004357549913514037)
        assert np.isclose((ped_800.iloc[500]), 0.05299289492753652)
        assert np.isclose(np.trapz(ped_400.values, x=ped_400.index), 1)
        assert np.isclose(np.trapz(ped_800.values, x=ped_800.index), 1)

    def test_beta_phi3a(self):
        """ test statmodels.pedmodels.beta_phi3a
        """
        ped_df_frag1_dct = mechanalyzer.builder.ped.ped_frag1(
            self.PED_DCT['RH->IC3H7'], 'CH3CHCH3', 'H2O', ['beta_phi3a'], dof_info=self.DOF_DCT['RH->IC3H7'], ene_bw=ene_bw_dct['RH->IC3H7'])
        ped_400 = ped_df_frag1_dct['beta_phi3a'][1.0][400]
        ped_800 = ped_df_frag1_dct['beta_phi3a'][1.0][800]
        assert np.isclose((ped_400.iloc[500]), 0.0009891036253060362)
        assert np.isclose((ped_800.iloc[500]), 0.059832564512289556)
        assert np.isclose(np.trapz(ped_400.values, x=ped_400.index), 1)
        assert np.isclose(np.trapz(ped_800.values, x=ped_800.index), 1)

    def test_rovib_dos(self):
        """ test statmodels.pedmodels.beta_rovib_dos
        """
        ped_df_frag1_dct = mechanalyzer.builder.ped.ped_frag1(
            self.PED_DCT['RH->IC3H7'], 'CH3CHCH3', 'H2O', ['rovib_dos'], dos_df=self.DOS_ROVIB, dof_info=self.DOF_DCT['RH->IC3H7'], ene_bw=ene_bw_dct['RH->IC3H7'])
        ped_1800 = ped_df_frag1_dct['rovib_dos'][1.0][1800]
        ped_2000 = ped_df_frag1_dct['rovib_dos'][1.0][2000]
        assert np.isclose((ped_1800.iloc[500]), 0.016860110020388747)
        assert np.isclose((ped_2000.iloc[500]), 0.013911651621907656)
        assert np.isclose(np.trapz(ped_1800.values, x=ped_1800.index), 1)
        assert np.isclose(np.trapz(ped_2000.values, x=ped_2000.index), 1)

    def test_bf_from_phi1a(self):
        """ test builder.bf.bf_tp_dct
            calls calculator.bf.bf_tp_df_full, bf_tp_df_todct
        """
        ped_df_frag1_dct = mechanalyzer.builder.ped.ped_frag1(
            self.PED_DCT['RH->IC3H7'], 'CH3CHCH3', 'H2O', ['beta_phi1a'], dof_info=self.DOF_DCT['RH->IC3H7'], ene_bw=ene_bw_dct['RH->IC3H7'])
        bf_tp_dct = mechanalyzer.builder.bf.bf_tp_dct(['beta_phi1a'], ped_df_frag1_dct, self.HOTEN_DCT['CH3CHCH3'], 0.1)
        assert np.allclose(bf_tp_dct['beta_phi1a']['CH3CHCH3'][1.0][1], 
            np.array([1., 0.99999461, 0.99904393, 0.88909341, 0.25846626, 0.01535025]))
        assert np.allclose(bf_tp_dct['beta_phi1a']['CH3CHCH3'][100.0][1], 
            np.array([1., 0.99999997, 0.99998977, 0.99576998, 0.82944956, 0.5149821]))
        assert np.allclose(bf_tp_dct['beta_phi1a']['P2'][1.0][1], 
            np.array([3.57703898e-10, 5.38693230e-06, 9.53732832e-04, 1.10030973e-01, 7.29985636e-01, 9.67187058e-01]))
        assert np.allclose(bf_tp_dct['beta_phi1a']['P2'][100.0][1], 
            np.array([3.03637884e-08, 1.02118304e-05, 4.17895042e-03, 1.66495546e-01, 4.73210178e-01]))
        self.bf_tp_dct = bf_tp_dct

    def test_new_ktp_dct(self):
        """ test builder.bf.merge_bf_ktp
            calls calculator.bf.merge_bf_rates
            calls builder.bf.rename_ktp_dct
        """
        rxn_ktp_dct = mechanalyzer.builder.bf.merge_bf_ktp(self.bf_tp_dct, KTP_DCT['RH->IC3H7'], FRAG_REACS, 'CH3CHCH3', 'H2O', HOT_FRAG_DCT)
        assert np.allclose(rxn_ktp_dct['beta_phi1a']['C3H8+OH=CH3CHCH3+H2O'][1.0][1], np.array([1.89880000e-13, 4.35038655e-13, 7.72470759e-13, 1.53941189e-12, 9.79452717e-13, 7.07670899e-14]))
        assert np.allclose(rxn_ktp_dct['beta_phi1a']['C3H8+OH=>CH3CHCH2+H+H2O'][1.0][1], np.array([6.79208161e-23, 2.34353641e-18, 7.37435763e-16, 1.90512027e-13, 2.76626597e-12, 4.45888709e-12]))


if __name__ == '__main__':
    tests = test_prompt()
    tests.test_equip_simple()
    tests.test_equip_phi()
    tests.test_beta_phi1a()
    tests.test_beta_phi2a()
    tests.test_beta_phi3a()
    tests.test_bf_from_phi1a()
    tests.test_new_ktp_dct()
    tests.test_rovib_dos()
