""" Test connections code
"""

import mechanalyzer

PES_DCT = {
    'C3H7O2': (
        (('R', 'O2'), ('RO2',)),
        (('RO2',), ('QOOH1',)),
        (('RO2',), ('QOOH2',)),
        (('RO2',), ('ALK', 'HO2')),
        (('QOOH1',), ('ALK', 'HO2')),
        (('QOOH1',), ('EPO', 'OH')),
        (('QOOH2',), ('EPO2', 'OH'))
    ),
    'C3H7O4(1)': (
        (('QOOH1', 'O2'), ('OOQOOH1',)),
        (('OOQOOH1',), ('HOOQhOOH1',)),
        (('HOOQhOOH1',), ('OQhOOH1', 'OH'))
    ),
    'C3H7O4(2)': (
        (('QOOH2', 'O2'), ('OOQOOH2',)),
        (('OOQOOH2',), ('HOOQhOOH2',)),
        (('HOOQhOOH2',), ('OQhOOH2', 'OH'))
    ),
    'C3H7O6(1)': (
        (('HOOQhOOH1', 'O2'), ('HOOQh(O2)OOH1',)),
    )
}


def test__():
    """ conn
    """

    conn_lst = mechanalyzer.parser.conn_pes(PES_DCT)
    print('CONN LST')
    print(conn_lst)


if __name__ == '__main__':
    test__()
