""" test mechanalyzer.builder.checker
"""

from mechanalyzer.builder import checker

RXN_PARAM_DCT1 = {
    (('H2', 'O'), ('OH', 'H'), (None,)): None,
    (('H', 'O2'), ('OH', 'O'), (None,)): None,
    (('H2', 'O'), ('OH', 'OH'), (None,)): None,
    (('H', 'O'), ('OH',), (None,)): None,
    (('H', 'O'), ('OH',), ('(+M)',)): None,
    (('H', 'O'), ('OH',), ('+O(S)',)): None,
    (('H2', 'O(S)'), ('OH', 'O'), (None,)): None,
    (('H2', 'O2'), ('HO2', 'H'), (None,)): None
}

RXN_PARAM_DCT2 = {
    (('H', 'O'), ('OH',), (None,)): None,
    (('OH',), ('H', 'O'), (None,)): None
}


def test__source_and_sink_species():
    """ test mechanalyzer.builder.checker.source_and_sink_species
    """

    sources1, sinks1 = checker.source_and_sink_species(RXN_PARAM_DCT1)
    sources2, sinks2 = checker.source_and_sink_species(RXN_PARAM_DCT2)
    assert list(set(sources1) - set(['O2', 'O(S)', 'H2'])) == []
    assert list(set(sinks1) - set(['OH', 'HO2'])) == []
    assert sources2 == sinks2 == []


if __name__ == '__main__':
    test__source_and_sink_species()
