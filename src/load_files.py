""" simple functions to load and parse data for HPO terms analyses
"""

from src.family_hpo_terms import FamilyHPO

def load_participants_hpo_terms(pheno_path, alt_id_path, alt_node_ids, obsolete_ids):
    """ loads patient HPO terms
    
    Args:
        pheno_path: path to patient pheotype file, containing one line per
            proband, with HPO codes as a field in the line.
        alt_id_path: path to set of alternate IDs for each individual
        alt_ids: dict to map HPO terms from their alt_id, to their current ID
        obsolete_ids: set of obsolete HPO IDs
    
    Returns:
        dictionary of FamilyHPO objects (containing HPO codes for trio members),
        indexed by proband ID.
    """
    
    # set the alt HPO terms and obsolete node IDs in the FamilyHPO class, to be
    # shared among all instances
    FamilyHPO.alt_ids = alt_node_ids
    FamilyHPO.obsolete_ids = obsolete_ids
    
    alt_ids = load_alt_id_map(alt_id_path)
    
    # load the phenotype data for each participant
    handle = open(pheno_path, "r")
    
    # get the positions of the columns in the list of header labels
    header = handle.readline().strip().split("\t")
    proband_column = header.index("patient_id")
    child_hpo_column = header.index("child_hpo")
    maternal_hpo_column = header.index("maternal_hpo")
    paternal_hpo_column = header.index("paternal_hpo")
    
    participant_hpo = {}
    for line in handle:
        line = line.split("\t")
        proband_id = line[proband_column]
        child = line[child_hpo_column]
        mother = line[maternal_hpo_column]
        father = line[paternal_hpo_column]
        
        # don't use probands who lack HPO terms
        if child == "NA":
            continue
        
        # swap the proband across to the DDD ID if it exists
        if proband_id in alt_ids:
            proband_id = alt_ids[proband_id]
        
        participant_hpo[proband_id] = FamilyHPO(child, mother, father)
    
    return participant_hpo

def load_alt_id_map(alt_id_path):
    """ loads the decipher to DDD ID mapping file.
    
    Args:
        alt_id_path: path to file containing alternate IDs (ie DECIPHER IDs) for
            each proband.
    
    Returns:
        dictionary of DDD IDs, indexed by their DECIPHER ID.
    """
    
    alt_ids = {}
    
    with open(alt_id_path) as f:
        for line in f:
            line = line.split("\t")
            ref_id = line[0]
            alt_id = line[1]
            
            alt_ids[alt_id] = ref_id
    
    return alt_ids

def load_variants(path):
    """ load variants found in probands, organised per gene.
    
    Args:
        path: path to variant containing file.
    
    Returns:
        dictionary of proband ID lists, indexed by HGNC symbol
    """
    
    functional = set(["stop_gained", "splice_acceptor_variant",
        "splice_donor_variant", "frameshift_variant", "missense_variant",
        "initiator_codon_variant", "stop_lost", "inframe_deletion",
        "inframe_insertion", "splice_region_variant"])
    
    variants = {}
    with open(path) as handle:
        header = handle.readline().strip().split("\t")
        
        for line in handle:
            line = line.strip().split("\t")
            
            proband_id = line[0]
            hgnc_symbol = line[8]
            consequence = line[10]
            
            if consequence not in functional:
                continue
            
            if hgnc_symbol not in variants:
                variants[hgnc_symbol] = set()
            
            variants[hgnc_symbol].add(proband_id)
    
    return variants
            
