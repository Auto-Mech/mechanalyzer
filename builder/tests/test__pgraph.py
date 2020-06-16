"""
test the code
"""

import pgraph


# List of reactions in the mechanism file
RCT_NAMES = [
    'CH2CH2+OH',
    'CH2CH2+OH',
    'C2H4OH',
    'C2H5O',
    'C2H4OH',
]
PRD_NAMES = [
    'CH2CH+H2O',
    'C2H4OH',
    'C2H5O',
    'CH3CHO+H',
    'CH2CH+H2O',
]

# Build a graph
PES_GRAPH = pgraph.build_pes_graph(RCT_NAMES, PRD_NAMES)


def test__graph_structure():
    """ test pgraph.build_PES_GRAPH
        test pgraph.species
        test pgraph.channels
    """
    assert PES_GRAPH == (
        ('CH2CH2+OH', 'CH2CH+H2O', 'C2H4OH', 'C2H5O', 'CH3CHO+H'),
        (frozenset({0, 1}), frozenset({0, 2}), frozenset({2, 3}),
         frozenset({3, 4}), frozenset({1, 2})))
    assert pgraph.species(PES_GRAPH) == (
        ('CH2CH2+OH', 'CH2CH+H2O', 'C2H4OH', 'C2H5O', 'CH3CHO+H'))
    assert pgraph.channels(PES_GRAPH) == (
        (frozenset({0, 1}), frozenset({0, 2}), frozenset({2, 3}),
         frozenset({3, 4}), frozenset({1, 2})))
    print('\npes graph')
    print(PES_GRAPH)


def test__paths():
    """ test pgraph.find_pathways
    """
    paths_idxs = pgraph.find_pathways(
        PES_GRAPH, 'CH2CH2+OH', 'CH2CH+H2O', pathval='idx')
    paths_names = pgraph.find_pathways(
        PES_GRAPH, 'CH2CH2+OH', 'CH2CH+H2O', pathval='name')
    assert paths_idxs == [[0, 1], [0, 2, 1]]
    assert paths_names == [
        ['CH2CH2+OH', 'CH2CH+H2O'],
        ['CH2CH2+OH', 'C2H4OH', 'CH2CH+H2O']]
    print('\npaths w/ idxs')
    print(paths_idxs)
    print('\npaths w/ names')
    print(paths_names)


def test__channel_manipulations():
    """ test pgraph.replace_species
    """
    new_graph = pgraph.replace_species(PES_GRAPH, 'C2H4OH', 'AAAA')
    print('\nreplaced graph')
    print(new_graph)


def test__check_connectivity():
    """ test pgraph.isolated_species
    """
    iso_spc_idxs = pgraph.isolated_species(PES_GRAPH, pathval='idx')
    iso_spc_names = pgraph.isolated_species(PES_GRAPH, pathval='name')
    assert not iso_spc_idxs
    assert not iso_spc_names
    print('\nisolated species w/ idxs')
    print(iso_spc_idxs)
    print('\nisolated species w/ names')
    print(iso_spc_names)


if __name__ == '__main__':
    test__graph_structure()
    test__paths()
    test__channel_manipulations()
    test__check_connectivity()
