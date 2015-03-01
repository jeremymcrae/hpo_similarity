""" simple functions to load and parse data for HPO terms analyses
"""

import json

def load_participants_hpo_terms(path, alt_ids, obsolete):
    """ loads patient HPO terms
    
    Args:
        path: path to JSON-encoded patient phenotype file.
        alt_ids: dict to map HPO terms from their alt_id, to their current ID
        obsolete: set of obsolete HPO IDs
    
    Returns:
        dictionary of HPO term lists indexed by proband ID e.g. {DDD01:
        [HP:01, HP:02], "DDD02": [HP:03, HP:03]}
    """
    
    # load the phenotype data for each participant
    with open(path, "r") as handle:
        hpo = json.load(handle)
    
    for proband in hpo:
        terms = hpo[proband]
        
        # strip out the obsolete terms, currently there are two probands (out of
        # >4000) who each have an obsolete term, so it's not worth converting the
        # obsolete terms to a more appropriate term
        terms = [term for term in terms if term not in obsolete]
        
        # convert each term to it's standard HPO ID if the term is in the HPO IDs,
        # otherwise just assume it is a standard HPO ID already.
        terms = [alt_ids[term] if term in alt_ids else term for term in terms]
        
        hpo[proband] = terms
    
    return hpo

def load_variants(path):
    """ loads dictionary of probands indexed by HGNC symbol
    
    Args:
        path: path to JSON-encoded dictionary of probands per gene
    
    Returns:
        dictionary of probands indexed by HGNC symbol, e.g. {ADNP: [DDD01,
        DDD02], ANKRD11: [DDD03, DDD04]}.
    """
    
    with open(path, "r") as handle:
        genes = json.load(handle)
    
    return genes
            
