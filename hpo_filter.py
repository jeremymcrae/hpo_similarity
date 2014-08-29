""" some graph analyses of HPO terms and their usage in patients and gene hits
"""

# dependencies: 
#    go-parser (http://bazaar.launchpad.net/~ntamas/+junk/go-parser/files)
#    networkx (http://networkx.github.io/, pip install --user networkx)
#    matplotlib (http://matplotlib.org/, pip install --user matplotlib)
#    HPO (http://compbio.charite.de/hudson/job/hpo/ or http://www.human-phenotype-ontology.org/)
#    DDG2P

# TODO: quantify success against reported variants
# TODO: investigate graph similarity analyses
# TODO: growth parameters - phenotype matching with percentiles - eg if a
#       child has microcephaly, exclude them if their head circumference
#       is above the 50th percentile


import os
import sys

import load_files
import create_hpo_graph as hpo
import match_hpo
import hpo_similarity as similarity
import hpo_filter_reporting as reporting

USER_PATH = "/nfs/users/nfs_j/jm33/"
HPO_FOLDER = os.path.join(USER_PATH, "apps", "hpo_filter")
HPO_PATH = os.path.join(HPO_FOLDER, "hpo_data", "hp.obo")
ALL_PROBAND_HPO_TERMS_PATH = os.path.join(HPO_FOLDER, "hpo_data", "patient_hpo_terms.txt")
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
    
    for gene in genes_count:
        print(str(gene[1]) + "\t" + str(gene[0]))
    
    return genes_count

def main():
    # build a graph of DDG2P terms, so we can trace paths between terms
    hpo_file = hpo.loadHPONetwork(HPO_PATH)
    graph = hpo_file.get_graph()
    alt_node_ids = hpo_file.get_alt_ids()
    
    # load gene HPO terms, proband HPO terms, probands with candidates in these genes
    ddg2p_genes = load_files.load_ddg2p(DDG2P_PATH)
    family_hpo_terms = load_files.load_participants_hpo_terms(PHENOTYPES_PATH, ALTERNATE_IDS_PATH)
    genes_index, probands_index = load_files.load_candidate_genes(CANDIDATE_VARIANTS_PATH)
    all_proband_hpo_terms = load_files.load_full_proband_hpo_list(ALL_PROBAND_HPO_TERMS_PATH)
    
    # load hpo terms that are required for specific genes
    obligate_hpo_terms = load_files.load_obligate_terms(OBLIGATE_GENES_PATH)
    ddg2p_hpo_organ_terms = load_files.load_organ_terms(ORGAN_TO_HPO_MAPPER_PATH, DDG2P_ORGAN_PATH)
    
    # # report the number of probands found for each gene
    # genes_count = count_genes(genes_index)
    
    # find the probands that have hpo terms that match the terms required for certain genes
    matcher = match_hpo.CheckHPOMatches(family_hpo_terms, graph, genes_index)
    obligate_matches = matcher.find_matches(obligate_hpo_terms)
    ddg2p_organ_matches = matcher.find_matches(ddg2p_hpo_organ_terms)
    
    # add the gene matches to the reporting file
    reporting_path = CANDIDATE_VARIANTS_PATH[:-4] + ".with_hpo_matches.txt"
    report = reporting.Reporting(CANDIDATE_VARIANTS_PATH, reporting_path)
    
    report.add_matches_to_report(obligate_matches, obligate_hpo_terms, "obligate_terms")
    report.add_matches_to_report(ddg2p_organ_matches, ddg2p_hpo_organ_terms, "ddg2p_organ_terms")
    
    # now look for similarity scores
    matcher = similarity.ICSimilarity(family_hpo_terms, ddg2p_genes, graph, alt_node_ids, all_proband_hpo_terms)
    scores_IC = matcher.get_similarity_scores(probands_index)
    # print(matcher.find_most_informative_term())
    
    matcher = similarity.ICDistanceSimilarity(family_hpo_terms, ddg2p_genes, graph, alt_node_ids, all_proband_hpo_terms)
    scores_distance = matcher.get_similarity_scores(probands_index)
    
    matcher = similarity.PathLengthSimilarity(family_hpo_terms, ddg2p_genes, graph, alt_node_ids, all_proband_hpo_terms)
    scores_length = matcher.get_similarity_scores(probands_index)
    
    matcher = similarity.JaccardSimilarity(family_hpo_terms, ddg2p_genes, graph, alt_node_ids, all_proband_hpo_terms)
    scores_jaccard = matcher.get_similarity_scores(probands_index)
    
    report.add_scores_to_report(scores_IC, "similarity_max_IC")
    report.add_scores_to_report(scores_distance, "similarity_max_distance")
    report.add_scores_to_report(scores_length, "similarity_path_length")
    report.add_scores_to_report(scores_jaccard, "similarity_jaccard")
    
    # found_terms = []
    # for node in graph:
    #     # print(node, graph.node[node])
    #     name = graph.node[node]["name"]
    #     found = False
    #     if any("tall" in s.lower() for s in name):
    #         found = True
        
    #     definition = [""]
    #     if "def" in graph.node[node]:
    #         definition = graph.node[node]["def"]
        
    #     if found:
    #         found_terms.append([node, name, definition])
    
    # for term in sorted(found_terms):
    #     print(term)

if __name__ == '__main__':
    main()


