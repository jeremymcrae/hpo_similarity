""" some graph analyses of HPO terms and their usage in patients and gene hits
"""

# dependencies: 
#    go-parser (http://bazaar.launchpad.net/~ntamas/+junk/go-parser/files)
#    networkx (http://networkx.github.io/, pip install --user networkx)
#    matplotlib (http://matplotlib.org/, pip install --user matplotlib)
#    HPO (http://compbio.charite.de/hudson/job/hpo/ or http://www.human-phenotype-ontology.org/)
#    DDG2P

# TODO: quantify success against reported variants
# TODO: set up information content calculations
# TODO: investigate graph similarity analyses
# TODO: syne1  - add in specific term matching, eg ataxia
# TODO: growth parameters - phenotype matching with percentiles - eg if a
#       child has microcephaly, exclude them if their head circumference
#       is above the 50th percentile


import networkx as nx
import os
import sys

import load_files
import create_hpo_graph as hpo
import match_hpo

# load the HPO database as graph database
    # convert terms between database versions, if necessary

# load clinical reporting file, find frequently identified genes, check their HPO terms from DDG2P, 
# load the patients HPO terms, check for matches between top genes and 

USER_PATH = "/nfs/users/nfs_j/jm33/"
HPO_FOLDER = os.path.join(USER_PATH, "apps", "hpo_filter")
HPO_PATH = os.path.join(HPO_FOLDER, "hpo_data", "hp.obo")
OBLIGATE_GENES_PATH = os.path.join(HPO_FOLDER, "obligate_terms", "obligate_hpo_terms.frequent_genes.txt")
ORGAN_TO_HPO_MAPPER_PATH = os.path.join(HPO_FOLDER, "obligate_terms", "ddg2p_organ_code_to_hpo_term.txt")
DDG2P_ORGAN_PATH = os.path.join(HPO_FOLDER, "obligate_terms", "ddg2p_organ_specificity.txt")
CANDIDATE_VARIANTS_PATH = os.path.join(USER_PATH, "clinical_reporting.txt")

ddd_freeze = "/nfs/ddd0/Data/datafreeze/1139trios_20131030/"
DDG2P_PATH = os.path.join(ddd_freeze, "DDG2P_with_genomic_coordinates_20131107.tsv")
PHENOTYPES_PATH = os.path.join(ddd_freeze, "phenotypes.shared.version.20131129.txt")
ALTERNATE_IDS_PATH = os.path.join(ddd_freeze, "person_sanger_decipher.private.txt")

def count_genes(genes_index):
    """ count the number of times each gene was found in the candidate genes
    
    Args:
        genes_index: dictionary of proband IDs and inheritance tuples indexed by gene
    
    Returns:
        list of genes, sorted by number of occurences
    """
    
    genes_count = []
    for gene in genes_index:
        genes_count.append((len(genes_index[gene]), gene))
    
    genes_count.sort()
    genes_count.reverse()
    
    return genes_count

def add_matches_to_report(path, matches, searched_genes, column_label):
    """ annotates a file of candidate variants using dicts of proband matches indexed by gene
    
    Args:
        path: path to file listing variants for each proband
        matches: fict of probands that were matched for each gene
        searched_genes: list of searched genes, so we can add whether a gene 
            was not found because it was not searched for
        column_label: label to use in the file header
    """
    
    f = open(path)
    header = f.readline().strip().split("\t")
    proband_label = "proband"
    gene_label = "gene"
    proband_column = header.index(proband_label)
    gene_column = header.index(gene_label)
    
    header.append(column_label + "_searched")
    header.append(column_label + "_passed")
    header = "\t".join(header) + "\n"
    
    parsed_lines = [header]
    for line in f:
        if line == "\n":
            parsed_lines.append(line)
            continue
        
        line = line.strip().split("\t")
        proband = line[proband_column]
        gene = line[gene_column]
        
        searched = "not_checked"
        if gene in searched_genes:
            searched = "gene_checked"
        
        if searched == "not_checked":
            matched = "NA"
        else:
            matched = "FAIL"
        
        if proband in matches:
            if gene in matches[proband]:
                matched = "PASS"
        
        line.append(searched)
        line.append(matched)
        
        line = "\t".join(line) + "\n"
        parsed_lines.append(line)
    
    f.close()
    
    output = open(path, "w")
    output.writelines(parsed_lines)
    output.close()

def main():
    # build a graph of DDG2P terms, so we can trace paths between terms
    hpo_file = hpo.loadHPONetwork(HPO_PATH)
    graph = hpo_file.get_graph()
    
    # load gene HPO terms, proband HPO terms, probands with candidates in these genes
    ddg2p_genes = load_files.load_ddg2p(DDG2P_PATH)
    family_hpo_terms = load_files.load_participants_hpo_terms(PHENOTYPES_PATH, ALTERNATE_IDS_PATH)
    genes_index = load_files.load_candidate_genes(CANDIDATE_VARIANTS_PATH)
    
    # load hpo terms that are required for specific genes
    obligate_hpo_terms = load_files.load_obligate_terms(OBLIGATE_GENES_PATH)
    ddg2p_hpo_organ_terms = load_files.load_organ_terms(ORGAN_TO_HPO_MAPPER_PATH, DDG2P_ORGAN_PATH)
    
    # # report the number of probands found for each gene
    # for gene in count_genes(genes_index):
    #     print str(gene[1]) + "\t" + str(gene[0])
    
    # find the probands that have hpo terms that match the terms required for certain genes
    matcher = match_hpo.checkHPOMatches(family_hpo_terms, graph, genes_index)
    obligate_matches = matcher.find_matches(obligate_hpo_terms)
    ddg2p_organ_matches = matcher.find_matches(ddg2p_hpo_organ_terms)
    
    # add the gene matches to the reporting file
    add_matches_to_report(CANDIDATE_VARIANTS_PATH, obligate_matches, obligate_hpo_terms, "obligate_terms")
    add_matches_to_report(CANDIDATE_VARIANTS_PATH, ddg2p_organ_matches, ddg2p_hpo_organ_terms, "ddg2p_organ_terms")
    

if __name__ == '__main__':
    main()


