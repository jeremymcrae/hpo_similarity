""" convenience script to prepare input files for hpo_similarity.py, using DDD
data sources.
"""

import os
import json

DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/"
PHENOTYPES_PATH = os.path.join(DATAFREEZE_DIR, "phenotypes_and_patient_info.txt")
ALTERNATE_IDS_PATH = os.path.join(DATAFREEZE_DIR, "person_sanger_decipher.txt")

USER_DIR = "/nfs/users/nfs_j/jm33/"
VARIANTS_PATH = os.path.join(USER_DIR, "apps/mupit/data-raw/de_novo_datasets/de_novos.ddd_4k.ddd_only.txt")
PHENOTYPES_OUT = os.path.join(USER_DIR, "apps/hpo_similarity/data/phenotypes_by_proband.json")
GENES_OUT = os.path.join(USER_DIR, "apps/hpo_similarity/data/probands_by_gene.json")


def load_participants_hpo_terms(pheno_path, alt_id_path, output_path):
    """ loads patient HPO terms
    
    Args:
        pheno_path: path to patient pheotype file, containing one line per
            proband, with HPO codes as a field in the line.
        alt_id_path: path to set of alternate IDs for each individual
        output_path: path to save HPO terms per proband as JSON-encoded file.
    """
    
    alt_ids = load_alt_id_map(alt_id_path)
    
    # load the phenotype data for each participant
    handle = open(pheno_path, "r")
    
    # get the positions of the columns in the list of header labels
    header = handle.readline().strip().split("\t")
    proband_column = header.index("patient_id")
    child_hpo_column = header.index("child_hpo")
    
    hpo_by_proband = {}
    for line in handle:
        line = line.split("\t")
        proband_id = line[proband_column]
        child_terms = line[child_hpo_column]
        
        # don't use probands who lack HPO terms
        if child_terms == "NA":
            continue
        
        # swap the proband across to the DDD ID if it exists
        if proband_id in alt_ids:
            proband_id = alt_ids[proband_id]
        
        if "|" in child_terms:
            child_terms = child_terms.split("|")
        else:
            child_terms = [child_terms]
        
        hpo_by_proband[proband_id] = child_terms
    
    with open(output_path, "w") as output:
        json.dump(hpo_by_proband, output, indent=4, sort_keys=True)

def load_alt_id_map(alt_id_path):
    """ loads the decipher to DDD ID mapping file.
    
    Args:
        alt_id_path: path to file containing alternate IDs (ie DECIPHER IDs) for
            each proband.
    
    Returns:
        dictionary of DDD IDs, indexed by their DECIPHER ID.
    """
    
    alt_ids = {}
    
    with open(alt_id_path) as handle:
        for line in handle:
            line = line.split("\t")
            ref_id = line[0]
            alt_id = line[1]
            
            alt_ids[alt_id] = ref_id
    
    return alt_ids

def load_variants(path, output_path):
    """ load variants found in probands, organised per gene.
    
    Args:
        path: path to variant containing file.
        output_path: path to save probands per gene as JSON-encoded file.
    """
    
    functional = set(["stop_gained", "splice_acceptor_variant",
        "splice_donor_variant", "frameshift_variant", "missense_variant",
        "initiator_codon_variant", "stop_lost", "inframe_deletion",
        "inframe_insertion", "splice_region_variant"])
    
    variants = {}
    with open(path) as handle:
        header = handle.readline()
        
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
    
    for hgnc in variants:
        variants[hgnc] = list(variants[hgnc])
    
    with open(output_path, "w") as output:
        json.dump(variants, output, indent=4, sort_keys=True)

def main():
    
    load_participants_hpo_terms(PHENOTYPES_PATH, ALTERNATE_IDS_PATH, PHENOTYPES_OUT)
    load_variants(VARIANTS_PATH, GENES_OUT)

if __name__ == "__main__":
    main()
