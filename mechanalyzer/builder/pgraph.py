"""
Looks at handling PES graphs
"""

import networkx


def ngx_from_graph(xgr):
    """ networkx graph object from a molecular graph
    """
    spcs, chnls = xgr
    nxg = networkx.Graph()
    for spc in spcs:
        nxg.add_node(spcs.index(spc))
    for chnl in chnls:
        nxg.add_edge(*chnl)

    return nxg


def build_pes_graph(rct_names, prd_names):
    """ Build a PES graph using the reactant and product names
    """

    # Get a list of the species
    spc = []
    for rct, prd in zip(rct_names, prd_names):
        if rct not in spc:
            spc.append(rct)
        if prd not in spc:
            spc.append(prd)

    # Get a list of the reactions
    chns = []
    for rct, prd in zip(rct_names, prd_names):
        rct_idx = spc.index(rct)
        prd_idx = spc.index(prd)
        chns.append(frozenset({rct_idx, prd_idx}))
        # Add directonality if wanted

    return (tuple(spc), tuple(chns))


def species(pes_graph):
    """ Get the species from a PES graph
    """
    return pes_graph[0]


def channels(pes_graph):
    """ Get the channels from a PES graph
    """
    return pes_graph[1]


def find_pathways(pes_graph, rct, prd, pathval='idx'):
    """ For a given reactant and product find any and all
        pathways that connect the two species
    """
    assert pathval in ('idx', 'name')

    # Get the indices for the rct and prd
    rct_idx = get_species_index(pes_graph, rct)
    prd_idx = get_species_index(pes_graph, prd)

    # Convert to ngx graph
    ngx = ngx_from_graph(pes_graph)

    # Get the paths and return them
    paths = networkx.all_simple_paths(ngx, source=rct_idx, target=prd_idx)

    # Convert paths using idxs to names if requested
    if pathval == 'name':
        paths = [[species(pes_graph)[idx] for idx in path] for path in paths]

    return list(paths)


def replace_species(pes_graph, in_spc, new_spc):
    """ Find a species and replace it with an alternative species
        Useful for changing CH2OOH to CH2O+OH
    """
    spcs, chnls = list(species(pes_graph)), channels(pes_graph)
    new_spcs = [spc if spc != in_spc else new_spc
                for spc in spcs]
    return (tuple(new_spcs), chnls)


def get_species_index(pes_graph, spc):
    """ Get the index in the PES species for a
    """
    return pes_graph[0].index(spc)


def isolated_species(pes_graph, pathval='idx'):
    """ Check connectivity of graph to ensure all species connected by rxn
    """
    # Get the isolated species
    iso_spc = list(networkx.isolates(ngx_from_graph(pes_graph)))

    # Convert paths using idxs to names if requested
    if pathval == 'name':
        iso_spc = [species(pes_graph)[idx] for idx in iso_spc]

    return iso_spc
# def remove_unphysical_channels_widx():
#     """ Take out the channels that of have been deemed to be unphysical
#         Uses the reactant and product index
#     """
#     pass
#
#
# def remove_unphysical_channels_wspc():
#     """ Take out the channels that of have been deemed to be unphysical
#         Remove any channel that contains a spc as a rct, prd, or either
#     """
#     pass
#
#
# def alter_species(channels, spc_idx, idx1, idx2):
#     """ Alter the channel indices
#     """
#     spc_idx = get_species_index(pes_graph, spc)
#     return None
#
#
# def alter_channel_index(channels, spc, idx1, idx2):
#     """ Alter the channel indices
#     """
#     spc_idx = get_species_index(pes_graph, spc)
#     new_channels = tuple(
#         [channel
#          if spc_idx not in channel
#          else frozenset({idx1, idx2})
#          for channel in channels]
#     )
#     return new_channels
