""" test mechanalyzer.calculator.statmodels
    mess_io required for input generation - too large input strings otherwise
    statmodels.ped_models is a class, tested in test__prompt.py
"""

import os
import numpy as np
from ioformat import pathtools
from ioformat import remove_comment_lines
import autoparse.pattern as app
import mess_io
from mechanalyzer.calculator.spinfo_frommess import get_info


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH_SGL = os.path.join(PATH, 'data', 'prompt', 'CH2O_H')
INP_PATH_DBL = os.path.join(PATH, 'data', 'prompt', 'C3H8_OH')
PED_INP_SGL = pathtools.read_file(INP_PATH_SGL, 'me_ktp_ped.inp')
PED_INP_DBL = pathtools.read_file(INP_PATH_DBL, 'me_ktp_ped.inp')

PED_INP_SGL = remove_comment_lines(
    PED_INP_SGL, delim_pattern=app.escape('#'))

def test_get_dof_info():
    """ test mechanalyzer.calculator.spinfo_frommess.get_info
    """

    # PES1
    species_blocks_ped1 = mess_io.reader.get_species(PED_INP_SGL)
    dof1 = get_info(
        species_blocks_ped1['R'])

    assert all(np.isclose(list(dof1.loc['HCO',['n_atoms', 'vib dof','rot dof', 'mw']].values), [
               3., 3., 3., 0.029], atol=1e-4))
    assert all(np.isclose(list(dof1.loc['H2',['n_atoms', 'vib dof','rot dof', 'mw']].values), [
               2.0, 1.0, 2.0, 0.002], atol=1e-4))
    assert all(np.isclose(list(dof1.loc['TS',['n_atoms', 'vib dof','rot dof', 'mw']].values), [
               5.0, 8.0, 3.0, 0.031], atol=1e-4))

    # PES2
    species_blocks_ped2 = mess_io.reader.get_species(PED_INP_DBL)
    dof2 = get_info(
        species_blocks_ped2['CH3CH2CH2+H2O'])

    assert all(np.isclose(list(dof2.loc['CH3CH2CH2',['vib dof','rot dof', 'mw']].values), [
               24., 3., 0.043], atol=1e-4, rtol=1e-3))
    assert all(np.isclose(list(dof2.loc['H2O',['vib dof','rot dof', 'mw']].values), [
               3., 3., 0.018], atol=1e-4))


if __name__ == '__main__':
    test_get_dof_info()
