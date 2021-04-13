def source_and_sink_species(rxn_param_dct):
    """ Produce two separate lists: (1) species that only appear as reactants (sources) and (2) 
        species that only appear as products (sinks). Either or both list can be empty if no 
        species meet these criteria.

        :param rxn_param_dct:
        :type rxn_param_dct: dct {rxn1: (param_tuple1, param_tuple2, ...), rxn2: ...}
        :return source_species: list of species that only appear as reactants
        :rtype: lst [spc1, spc2, ...]
        :return sink_species: list of species that only appear as products
        :rtype: lst [spc1, spc2, ...]
    """
    all_rcts = []
    all_prds = []
    for rxn in rxn_param_dct.keys():
        rcts, prds, _ = rxn
        for rct in rcts:
            all_rcts.append(rct)
        for prd in prds:
            all_prds.append(prd)
    source_species = list(set(all_rcts) - set(all_prds))
    sink_species = list(set(all_prds) - set(all_rcts))

    return source_species, sink_species


    
