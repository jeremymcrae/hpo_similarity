""" simple functions to load and parse data for HPO terms analyses
"""

import sys

from src.family_hpo_terms import FamilyHPO

def load_ddg2p(ddg2p_path):
    """ load a DDG2P gene file, so we can extract HPO terms for each gene
    
    Args:
        ddg2p_path: path to DDG2P file
    
    Returns:
        dictionary of hpo terms indexed by gene name, then inheritance type
    """
    
    f = open(ddg2p_path)
    
    # get the positions of the columns in the header
    header = f.readline().strip().split("\t")
    gene_column = header.index("gencode_gene_name")
    inheritance_column = header.index("Allelic_requirement")
    hpo_column = header.index("HPO_ids")
    
    genes = {}
    for line in f:
        line = line.split("\t")
        
        gene = line[gene_column]
        inheritance = line[inheritance_column]
        hpo_terms = line[hpo_column].strip()
        
        if gene not in genes:
            genes[gene] = {}
        
        if inheritance not in genes[gene]:
            genes[gene][inheritance] = set()
        
        # split the hpo terms, then update the gene/inheritance set
        hpo_terms = hpo_terms.split(";")
        hpo_terms = [x.strip() for x in hpo_terms]
        genes[gene][inheritance].update(hpo_terms)
        
    return genes

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
    
    # use the alt HPO terms and obsolete node Ids in the FamilyHPO class, to be
    # shared among all instances
    FamilyHPO.alt_ids = alt_node_ids
    FamilyHPO.obsolete_ids = obsolete_ids
    
    alt_ids = load_alt_id_map(alt_id_path)
    
    # load the phenotype data for each participant
    f = open(pheno_path)
    
    # get the positions of the columns in the list of header labels
    header = f.readline().strip().split("\t")
    proband_column = header.index("patient_id")
    child_hpo_column = header.index("child_hpo")
    maternal_hpo_column = header.index("maternal_hpo")
    paternal_hpo_column = header.index("paternal_hpo")
    
    participant_hpo = {}
    for line in f:
        line = line.split("\t")
        proband_id = line[proband_column]
        child_hpo = line[child_hpo_column]
        maternal_hpo = line[maternal_hpo_column]
        paternal_hpo = line[paternal_hpo_column]
        
        # don't bother to use probands who lack HPO terms
        if child_hpo == "NA":
            continue
        
        # swap the proband across to the DDD ID if it exists
        if proband_id in alt_ids:
            proband_id = alt_ids[proband_id]
        
        participant_hpo[proband_id] = FamilyHPO(child_hpo, maternal_hpo, paternal_hpo)
    
    return participant_hpo

def load_alt_id_map(alt_id_path):
    """ loads the decipher to DDD ID mapping file
    
    Args:
        alt_id_path: path to file containing alternate IDs (ie DECIPHER IDs) for
            each proband
    
    Returns:
        dictionary of DDD Ids, indexed by their DECIPHER ID.
    """
    
    alt_ids = {}
    
    with open(alt_id_path) as f:
        for line in f:
            line = line.split("\t")
            ref_id = line[0]
            alt_id = line[1]
            
            # if ":" in alt_id:
            #     alt_id = alt_id.split(":")[0]
            
            alt_ids[alt_id] = ref_id
    
    return alt_ids

def load_clinical_filter_variants(path):
    """ loads the genes for variants passing clinical filtering for each proband
    
    Args:
        path: path to clinical filtering output
    
    Returns:
        dictionary of [(gene, inheritance)] tuple lists, indexed by proband
    """
    
    f = open(path)
    header = f.readline().strip().split("\t")
    
    proband_column = header.index("proband")
    gene_column = header.index("gene")
    inheritance_column = header.index("inheritance")
    
    probands = {}
    
    for line in f:
        if line == "\n": # ignore blank lines
            continue
        
        line = line.strip().split("\t")
        proband_ID = line[proband_column]
        gene = line[gene_column]
        inheritance = line[inheritance_column]
        
        if proband_ID not in probands:
            probands[proband_ID] = set()
        
        probands[proband_ID].add((gene, inheritance))
    
    return probands

def load_de_novos(path):
    """ load de novos found in probands
    
    Args:
        path: path to de novo containing file
    
    Returns:
        dictionary of proband ID lists, indexed by HGNC symbol
    """
    
    functional = set(["stop_gained", "splice_acceptor_variant",
        "splice_donor_variant", "frameshift_variant", "missense_variant",
        "initiator_codon_variant", "stop_lost", "inframe_deletion",
        "inframe_insertion", "splice_region_variant"])
    
    de_novos = {}
    with open(path) as handle:
        header = handle.readline().strip().split("\t")
        
        for line in handle:
            line = line.strip().split("\t")
            
            proband_id = line[0]
            hgnc_symbol = line[8]
            consequence = line[10]
            
            if consequence not in functional:
                continue
            
            if hgnc_symbol not in de_novos:
                de_novos[hgnc_symbol] = set()
            
            de_novos[hgnc_symbol].add(proband_id)
    
    return de_novos
            
