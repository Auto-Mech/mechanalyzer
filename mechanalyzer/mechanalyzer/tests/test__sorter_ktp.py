"""
Read the mechanism file
"""

import os
import numpy as np
import mechanalyzer

BIG_ARRAY = np.array([1e15, 1e15, 1e15])
MIDDLE_ARRAY = np.array([1e14, 1e14, 1e14]) 
LITTLE_ARRAY = np.array([1e13, 1e13, 1e13]) 

al_ktp_dct = {
(('H2', 'O'), ('OH', 'H'), (None,)): [
    {'high': (np.array([500, 1000, 1500]), np.array([0.157572885e+134, 2.79926202e+143, 1.72670689e+149])),
     1: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
     10: (np.array([500, 1000, 1500]), np.array([6.57572885e+134, 8.79926202e+143, 4.72670689e+149]))},
    {'high': (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
     1: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
     10: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))}],
 (('H', 'O2'), ('OH', 'O'), (None,)): [
     {'high': (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
      1: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
      10: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))},
     {'high': (np.array([500, 1000, 1500]), [4.6548277231154764e+45, 8.556998184634325e+52, 4.662500917095324e+56]),
      1: (np.array([500, 1000, 1500]), [4.6548277231154764e+45, 8.556998184634325e+52, 4.662500917095324e+56]),
      10: (np.array([500, 1000, 1500]), [4.6548277231154764e+45, 8.556998184634325e+52, 4.662500917095324e+56])}],
 (('H2', 'O'), ('OH', 'OH'), (None,)): [
     {'high': (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
      1: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
      10: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))}, None],
 (('H', 'O'), ('OH',), (None,)): [
     {'high': (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
      1: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
      10: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))},
     {'high': (np.array([500, 1000, 1500]), [1.420849319576619e+96, 2.8405686431169553e+77, 3.4922934313599517e+72]),
      1: (np.array([500, 1000, 1500]), [5.8295576381190475e+100, 2.3308958102634265e+82, 4.2985260083885116e+77]),
      10: (np.array([500, 1000, 1500]), [2.3917907260059445e+105, 1.912671707993609e+87, 5.2908858341829314e+82])}],
 (('H', 'O'), ('OH',), ('(+M)',)): [
     {'high': (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
      1: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
      10: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))},
     {'high': (np.array([500, 1000, 1500]), [9.813202359645695e+109, 1.569488025355258e+92, 6.512342336821681e+87]),
      1: (np.array([500, 1000, 1500]), [4.0262276922599165e+114, 1.2878805345625882e+97, 8.015784887656628e+92]),
      10: (np.array([500, 1000, 1500]), [1.6519081983453455e+119, 1.0568008449314422e+102, 9.866312924289953e+97])}],
# (('H', 'O'), ('OH',), ('+O(S)',)): [
#     {'high': (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
#      1: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
#      10: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))},
#     {'high': (np.array([500, 1000, 1500]), [6.777561788188122e+123, 8.671829380720677e+106, 1.2144054772466455e+103]),
#      1: (np.array([500, 1000, 1500]), [2.7807443439483792e+128, 7.115874780854176e+111, 1.4947637222572564e+108]),
#      10: (np.array([500, 1000, 1500]), [1.1409027830446496e+133, 5.839099418788192e+116, 1.8398456094270225e+113])}],
 (('H2', 'O(S)'), ('OH', 'H'), (None,)): [
     {'high': (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+139])),
      1: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+139])),
      10: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+139]))},
     {'high': (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+139])),
      1: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+139])),
      10: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+139]))}],
 (('H2', 'O2'), ('HO2V', 'H'), (None,)): [None, {
     'high': (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
     1: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
     10: (np.array([500, 1000, 1500]), np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))}]}

# start here
CWD = os.getcwd()

# MODIFY THIS SECTION WITH INPUT NAMES AND SORTING OPTIONS

SPC_NAMES = ['data/spc2.csv', 'data/spc1B.csv']
MECH_NAMES = ['data/mech2.txt', 'data/mech1.txt']
SORTMECH_NAME = 'sorted.txt'
ISOLATE_SPECIES = []
SORT_STR = ['molecularity','rxn_max_vals','rxn_max_ratio','rxn_class_broad',0] 

############ input reading ####################

# READ FILE# READ FILE AND BUILD DICTIONARIES
for i, SPC_NAME in enumerate(SPC_NAMES):
        spc_dct_full,rxn_param_dct,elem_tuple = mechanalyzer.parser.mech.readfiles(os.path.join(CWD,SPC_NAME),os.path.join(CWD,MECH_NAMES[i]))

# BUILD  MECH INFORMATION
mech_info = mechanalyzer.parser.mech.build_dct(spc_dct_full,al_ktp_dct)

# SORTING: sort the mech and build the sorted rxn param dct
sorted_idx,cmts_dct,spc_dct = mechanalyzer.parser.mech.sort_mechanism(mech_info,spc_dct_full,SORT_STR,ISOLATE_SPECIES)
al_ktp_dct_sorted = mechanalyzer.parser.mech.reordered_mech(al_ktp_dct,sorted_idx)
rxn_param_dct_sorted =  mechanalyzer.parser.mech.reordered_mech(rxn_param_dct,sorted_idx)
print(al_ktp_dct_sorted)
# WRITE THE NEW MECHANISM
spc_dct = mechanalyzer.parser.spc.order_species_by_atomcount(spc_dct)
chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=elem_tuple, spc_dct=spc_dct, rxn_param_dct=rxn_param_dct_sorted,
    comments=None)


