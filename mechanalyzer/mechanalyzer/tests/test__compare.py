from mechanalyzer.calculator import compare
from mechanalyzer.calculator.sample_dcts import SPC_IDENT_DCT1
from mechanalyzer.calculator.sample_dcts import SPC_IDENT_DCT2
from mechanalyzer.calculator.sample_dcts import SPC_THERMO_DCT1
from mechanalyzer.calculator.sample_dcts import SPC_THERMO_DCT2
from mechanalyzer.calculator.sample_dcts import RXN_KTP_DCT1
from mechanalyzer.calculator.sample_dcts import RXN_KTP_DCT2
from mechanalyzer.calculator.sample_dcts import TEMPS




CORRECT_SPC_KEYS = ('H', 'OH', 'O', 'O2', 'H2', 'HO2V', 'O(S)')
CORRECT_RENAMED_RXN_KEYS = (
    (('H2', 'O'), ('OH', 'H'), (None,)),
    (('O', 'OH'), ('O2', 'H'), (None,)),
    (('H2', 'O2'), ('HO2V', 'H'), (None,)),
    (('OH',), ('H', 'O'), (None,)),
    (('OH',), ('H', 'O'), ('(+M)',)),
    (('OH',), ('H', 'O'), ('+O(S)',)),
    (('H2', 'O(S)'), ('O', 'OH'), (None,))
)

# For the test case of when rev_rates=True
# (These are seemingly out of order due to the fact that they are popped out and then readded)
CORRECT_REVERSED_RXN_KEYS = (
    (('H2', 'O'), ('OH', 'H'), (None,)), 
    (('H2', 'O2'), ('HO2V', 'H'), (None,)), 
    (('H', 'O2'), ('OH', 'O'), (None,)), 
    (('H', 'O'), ('OH',), (None,)), 
    (('H', 'O'), ('OH',), ('(+M)',)), 
    (('H', 'O'), ('OH',), ('+O(S)',)),
    (('H2', 'O(S)'), ('OH', 'O'), (None,))
)

# For the test case of when rev_rates=False 
# (Only the last reaction has had the products reordered)
CORRECT_PARTIALLY_REVERSED_RXN_KEYS = (
    (('H2', 'O'), ('OH', 'H'), (None,)), 
    (('O', 'OH'), ('O2', 'H'), (None,)), 
    (('H2', 'O2'), ('HO2V', 'H'), (None,)), 
    (('OH',), ('H', 'O'), (None,)), 
    (('OH',), ('H', 'O'), ('(+M)',)), 
    (('OH',), ('H', 'O'), ('+O(S)',)), 
    (('H2', 'O(S)'), ('OH', 'O'), (None,))
)

# Without removing loners
CORRECT_ALIGNED_RXN_KTP_KEYS = (
    (('H2', 'O'), ('OH', 'H'), (None,)), 
    (('H', 'O2'), ('OH', 'O'), (None,)), 
    (('H2', 'O'), ('OH', 'OH'), (None,)), 
    (('H', 'O'), ('OH',), (None,)), 
    (('H', 'O'), ('OH',), ('(+M)',)), 
    (('H', 'O'), ('OH',), ('+O(S)',)), 
    (('H2', 'O(S)'), ('OH', 'O'), (None,)), 
    (('H2', 'O2'), ('HO2V', 'H'), (None,))
)

# With removing loners
CORRECT_ALIGNED_RXN_KTP_KEYS_NO_LONERS = (
    (('H2', 'O'), ('OH', 'H'), (None,)), 
    (('H', 'O2'), ('OH', 'O'), (None,)), 
    (('H', 'O'), ('OH',), (None,)), 
    (('H', 'O'), ('OH',), ('(+M)',)), 
    (('H', 'O'), ('OH',), ('+O(S)',)),
    (('H2', 'O(S)'), ('OH', 'O'), (None,))
)


def test_rename_spc_ident_dct():
    rename_instructions, _ = compare.get_rename_instructions(SPC_IDENT_DCT1, SPC_IDENT_DCT2)
    renamed_dct = compare.rename_species(SPC_IDENT_DCT2, rename_instructions, target_type='spc')
    assert tuple(renamed_dct.keys()) == CORRECT_SPC_KEYS
    assert tuple(renamed_dct.values()) == tuple(SPC_IDENT_DCT2.values())


def test_rename_spc_thermo_dct():
    rename_instructions, _ = compare.get_rename_instructions(SPC_IDENT_DCT1, SPC_IDENT_DCT2)
    renamed_dct = compare.rename_species(SPC_THERMO_DCT2, rename_instructions, target_type='spc')
    assert tuple(renamed_dct.keys()) == CORRECT_SPC_KEYS
         

def test_rename_rxn_ktp_dct():
    rename_instructions, _ = compare.get_rename_instructions(SPC_IDENT_DCT1, SPC_IDENT_DCT2)
    renamed_dct = compare.rename_species(RXN_KTP_DCT2, rename_instructions, target_type='rxn')
    assert tuple(renamed_dct.keys()) == CORRECT_RENAMED_RXN_KEYS


def test_rename_dcts():
    rxn_ktp_dcts = [RXN_KTP_DCT1, RXN_KTP_DCT2] 
    spc_thermo_dcts = [SPC_THERMO_DCT1, SPC_THERMO_DCT2]
    spc_ident_dcts = [SPC_IDENT_DCT1, SPC_IDENT_DCT2]

    # Rename the rxn_ktp and spc_thermo dcts
    renamed_rxn_ktp_dcts, _ = compare.rename_dcts(rxn_ktp_dcts, spc_ident_dcts, target_type='rxn')
    renamed_spc_thermo_dcts, _ = compare.rename_dcts(spc_thermo_dcts, spc_ident_dcts, 
                                                     target_type='spc')
    assert tuple(renamed_rxn_ktp_dcts[1].keys()) == CORRECT_RENAMED_RXN_KEYS
    assert tuple(renamed_spc_thermo_dcts[1].keys()) == CORRECT_SPC_KEYS
    

def test_reverse_rxn_ktp_dcts():
    rxn_ktp_dcts = [RXN_KTP_DCT1, RXN_KTP_DCT2]
    spc_thermo_dcts = [SPC_THERMO_DCT1, SPC_THERMO_DCT2]
    spc_ident_dcts = [SPC_IDENT_DCT1, SPC_IDENT_DCT2]

    # Rename the rxn_ktp and spc_thermo dcts
    renamed_rxn_ktp_dcts, _ = compare.rename_dcts(rxn_ktp_dcts.copy(), spc_ident_dcts.copy(), 
                                                  target_type='rxn')
    renamed_spc_thermo_dcts, _ = compare.rename_dcts(spc_thermo_dcts.copy(), spc_ident_dcts.copy(),
                                                     target_type='spc')
 #   print('before\n', rxn_ktp_dcts) 
    # Check functionality of the full reversal
    reversed_rxn_ktp_dcts = compare.reverse_rxn_ktp_dcts(
        renamed_rxn_ktp_dcts.copy(), renamed_spc_thermo_dcts.copy(), TEMPS, rev_rates=True
    )
#    print('after\n', rxn_ktp_dcts) 
    
    # Check functionality of flipping the ordering of reactants and/or products but NOT reversing
    partially_reversed_rxn_ktp_dcts = compare.reverse_rxn_ktp_dcts(
        renamed_rxn_ktp_dcts, renamed_spc_thermo_dcts, TEMPS, rev_rates=False
    )

    # Check that keys are correct after being reversed
    assert tuple(reversed_rxn_ktp_dcts[0].keys()) == tuple(RXN_KTP_DCT1.keys())  
    assert tuple(reversed_rxn_ktp_dcts[1].keys()) == CORRECT_REVERSED_RXN_KEYS
    assert tuple(partially_reversed_rxn_ktp_dcts[1].keys()) == CORRECT_PARTIALLY_REVERSED_RXN_KEYS

    # Check a numerical case
    


def test_align_rxn_ktp_dcts():
    rxn_ktp_dcts = [RXN_KTP_DCT1, RXN_KTP_DCT2]
    spc_thermo_dcts = [SPC_THERMO_DCT1, SPC_THERMO_DCT2]
    spc_ident_dcts = [SPC_IDENT_DCT1, SPC_IDENT_DCT2]

    # Try without removing loners
    aligned_rxn_ktp_dct = compare.get_aligned_rxn_ktp_dct(
        rxn_ktp_dcts, spc_thermo_dcts, spc_ident_dcts, TEMPS, rev_rates=True, remove_loners=False,
        write_file=False
    )
    print('aligned:\n', aligned_rxn_ktp_dct)
    # Try with removing loners
    no_loners_aligned_rxn_ktp_dct = compare.get_aligned_rxn_ktp_dct(
        rxn_ktp_dcts, spc_thermo_dcts, spc_ident_dcts, TEMPS, rev_rates=True, remove_loners=True,
        write_file=False
    )
    assert tuple(aligned_rxn_ktp_dct.keys()) == CORRECT_ALIGNED_RXN_KTP_KEYS
    assert tuple(no_loners_aligned_rxn_ktp_dct.keys()) == CORRECT_ALIGNED_RXN_KTP_KEYS_NO_LONERS


if __name__ == '__main__':
    test_rename_spc_ident_dct()
    test_rename_spc_thermo_dct()
    test_rename_rxn_ktp_dct()
    test_rename_dcts()
    test_reverse_rxn_ktp_dcts()
    test_align_rxn_ktp_dcts()

