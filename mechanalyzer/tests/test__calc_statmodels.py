""" test mechanalyzer.calculator.statmodels
    mess_io required for input generation - too large input strings otherwise
    statmodels.ped_models is a class, tested in test__prompt.py
"""

import os
import numpy as np
import ioformat.pathtools
import mess_io
import mechanalyzer


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH_SGL = os.path.join(PATH, 'data', 'prompt', 'CH2O_H')
INP_PATH_DBL = os.path.join(PATH, 'data', 'prompt', 'C3H8_OH')
PED_INP_SGL = ioformat.pathtools.read_file(INP_PATH_SGL, 'me_ktp_ped.inp')
PED_INP_DBL = ioformat.pathtools.read_file(INP_PATH_DBL, 'me_ktp_ped.inp')


def test_get_dof_info():
    """ test mechanalyzer.calculator.statmodels
    """

    # PES1
    species_blocks_ped1 = mess_io.reader.get_species(PED_INP_SGL)
    dof1 = mechanalyzer.calculator.statmodels.get_dof_info(
        species_blocks_ped1['R'], ask_for_ts=True)

    assert all(np.isclose(list(dof1.loc['HCO'].values), [
               3., 3., 0.029], atol=1e-4))
    assert all(np.isclose(list(dof1.loc['H2'].values), [
               1.0, 2.0, 0.002], atol=1e-4))
    assert all(np.isclose(list(dof1.loc['TS'].values), [
               8.0, 3.0, 0.031], atol=1e-4))

    # PES2
    species_blocks_ped2 = mess_io.reader.get_species(PED_INP_DBL)
    dof2 = mechanalyzer.calculator.statmodels.get_dof_info(
        species_blocks_ped2['NC3H7'])

    assert all(np.isclose(list(dof2.loc['CH3CH2CH2'].values), [
               24., 3., 0.043], atol=1e-4, rtol=1e-3))
    assert all(np.isclose(list(dof2.loc['H2O'].values), [
               3., 3., 0.018], atol=1e-4))


if __name__ == '__main__':
    test_get_dof_info()
