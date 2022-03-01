""" Handles names for mechanisms
"""

import automol.inchi
import automol.geom
import automol.graph
import automol.formula


# Name remaping function
def remap_mechanism_names(mech_spc_dct, rxn_param_dct, map_dct):
    """ Change all of the names in a spc_dct and rxn_param_dct
        according to the provided dictionary
    """

    # Fill the map dict if any names are missing
    for name in mech_spc_dct:
        if name not in map_dct:
            map_dct[name] = name

    # Alter the names of the dictionary
    re_mech_spc_dct = {}
    for name, dct in mech_spc_dct.items():
        re_mech_spc_dct[map_dct[name]] = dct

    # Alter the names of the reactions
    re_rxn_param_dct = {}
    for rxn, params in rxn_param_dct.items():
        rcts, prds, thrd = rxn[0], rxn[1], rxn[2]
        re_rxn = (
            tuple(map_dct[rct] for rct in rcts),
            tuple(map_dct[prd] for prd in prds),
            thrd
        )
        re_rxn_param_dct[re_rxn] = params

    return re_mech_spc_dct, re_rxn_param_dct


# Build various dictionaries
def functional_group_name_dct(mech_spc_dct, rename_rule_dct):
    """ Build a dictionary to map the names of a mechanism according
        to its functional groups and number of carbon atoms.

        The function detects which functional groups are in a species
        and then uses the list of these groups to assign a mechanism
        name according to a rule set bu the user in a supplied dictionary.

        Example:
            OOQOOH is name assigned to species with -OOH and -OO groups
            RO2 is name assigned to species with JUST -OO groups

        The hierarchy of names is determined by the order of the supplied
        dictionary.

        So the final name should correspond to
            <NCarbon>-<FUNC-GRP-RULE-NAME>

        In cases where there are less that 1 carbons or ID'd functional groups,
        no naming rule established

        :param rename_rule_dct: assigns string name to set of func. groups
        :type rename_rule_dct: dct[frozenset(grp1, grp2): str]
    """

    def _fgrp_name_string(fgrps, rule_dct, name):
        """ Determines what the new name for a species should
            be according the functional groups it has and the
            rule set by the user.

            grps and dct must be sets
        """
        if fgrps:
            _name = None
            for fgrp_lst, fgrp_name in rule_dct.items():
                if fgrp_lst <= fgrps:
                    _name = fgrp_name
                    break
            if _name is None:
                _name = name
        else:
            _name = name

        return _name

    fgrp_map_dct = {}
    count_dct = {}  # For adding numbers
    for name, dct in mech_spc_dct.items():
        # Initialize re_name to None in case we don't wish to establish a rule
        re_name = None
        print(name)
        if 'cbh' not in name:
            # Grab SPC inchi to use to get other info
            _ich = dct['inchi']

            # Get the number of atoms
            geo = automol.inchi.geometry(_ich)
            c_cnt = automol.geom.atom_count(geo, 'C', match=True)

            # Set the name based on the counts and functional groups
            # Do the formula
            if c_cnt > 1:

                # Get the formula
                fml_str = automol.formula.string(automol.inchi.formula(_ich))

                # Get func groups
                gra = automol.inchi.graph(_ich)
                fgrp_dct = automol.graph.functional_group_dct(gra)
                fgrps = frozenset(fgrp for fgrp, grp_idxs in fgrp_dct.items()
                                  if grp_idxs)
                print(fgrps)
                _fgrp_name = _fgrp_name_string(fgrps, rename_rule_dct, name)

                re_name = f'{fml_str}-{_fgrp_name}'

        # Fill the dictionary, increasing the count as needed
        if re_name is not None:
            if re_name in count_dct:
                count_dct[re_name] += 1
                fgrp_map_dct[name] = re_name + f'({count_dct[re_name]})'
            else:
                count_dct[re_name] = 1
                fgrp_map_dct[name] = re_name + '(1)'
            print(fgrp_map_dct[name])
        else:
            print(name)
        print()

    return fgrp_map_dct
