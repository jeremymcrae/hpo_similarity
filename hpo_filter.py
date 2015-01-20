""" some graph analyses of HPO terms and their usage in patients and gene hits
"""

# dependencies:
#    go-parser (http://bazaar.launchpad.net/~ntamas/+junk/go-parser/files)
#    networkx (http://networkx.github.io/, pip install --user networkx)
#    matplotlib (http://matplotlib.org/, pip install --user matplotlib)
#    HPO (http://compbio.charite.de/hudson/job/hpo/ or http://www.human-phenotype-ontology.org/)
#    DDG2P

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import sys

from src.load_files import load_ddg2p, load_participants_hpo_terms, \
    load_clinical_filter_variants
from src.create_hpo_graph import loadHPONetwork
from src.similarity import ICSimilarity, ICDistanceSimilarity, \
    PathLengthSimilarity, JaccardSimilarity
from src.reporting import Reporting

USER_PATH = "/nfs/users/nfs_j/jm33/"
HPO_FOLDER = os.path.join(USER_PATH, "apps", "hpo_filter")
HPO_PATH = os.path.join(HPO_FOLDER, "hpo_data", "hp.obo")
CANDIDATE_VARIANTS_PATH = os.path.join(USER_PATH, "clinical_reporting.2015-01-15.txt")

ddd_freeze = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/"
DDG2P_PATH = os.path.join(ddd_freeze, "DDG2P_freeze_with_gencode19_genomic_coordinates_20141118_fixed.txt")
PHENOTYPES_PATH = os.path.join(ddd_freeze, "phenotypes_and_patient_info.txt")
ALTERNATE_IDS_PATH = os.path.join(ddd_freeze, "person_sanger_decipher.txt")

def count_genes(probands):
    """ count the number of times each gene was found in the candidate genes
    
    Args:
        probands: dictionary of genes and inheritance tuples indexed by proband
    
    Returns:
        list of genes, sorted by number of occurences
    """
    
    # find how many times a gene was found within the probands
    genes = {}
    for proband in probands:
        for (gene, inh) in probands[proband]:
            if gene not in genes:
                genes[gene] = 0
            genes[gene] += 1
    
    counts = genes.values()
    gene_names = genes.keys()
    
    # get a gene list sorted by number of times a gene was present, highest first
    genes_count = zip(counts, gene_names)
    genes_count.sort()
    genes_count.reverse()
    
    for gene in genes_count:
        print("{0}\t{1}".format(gene[1], gene[0]))
    
    return genes_count

def main():
    # build a graph of DDG2P terms, so we can trace paths between terms
    hpo_file = loadHPONetwork(HPO_PATH)
    graph = hpo_file.get_graph()
    alt_node_ids = hpo_file.get_alt_ids()
    obsolete_ids = hpo_file.get_obsolete_ids()
    
    # load gene HPO terms, proband HPO terms, probands with candidates in these genes
    ddg2p_genes = load_ddg2p(DDG2P_PATH)
    family_hpo_terms = load_participants_hpo_terms(PHENOTYPES_PATH, ALTERNATE_IDS_PATH, alt_node_ids, obsolete_ids)
    probands = load_clinical_filter_variants(CANDIDATE_VARIANTS_PATH)
    
    # add the gene matches to the reporting file
    reporting_path = CANDIDATE_VARIANTS_PATH[:-4] + ".with_hpo_matches.txt"
    report = Reporting(CANDIDATE_VARIANTS_PATH, reporting_path)
    
    # # now look for similarity scores
    # matcher = ICSimilarity(family_hpo_terms, ddg2p_genes, graph, alt_node_ids)
    # matcher.tally_hpo_terms(ddg2p_genes, source="ddg2p")
    # scores_IC = matcher.get_similarity_scores(ddg2p_genes, probands)
    #
    # matcher = ICDistanceSimilarity(family_hpo_terms, ddg2p_genes, graph, alt_node_ids)
    # matcher.tally_hpo_terms(ddg2p_genes, source="ddg2p")
    # scores_distance = matcher.get_similarity_scores(ddg2p_genes, probands)
    
    matcher = PathLengthSimilarity(family_hpo_terms, graph, alt_node_ids)
    matcher.tally_hpo_terms(ddg2p_genes, source="ddg2p")
    scores_length = matcher.get_similarity_scores(ddg2p_genes, probands)
    
    # matcher = JaccardSimilarity(family_hpo_terms, ddg2p_genes, graph, alt_node_ids)
    # matcher.tally_hpo_terms(ddg2p_genes, source="ddg2p")
    # scores_jaccard = matcher.get_similarity_scores(ddg2p_genes, probands)
    
    # report.add_scores_to_report(scores_IC, "similarity_max_IC")
    # report.add_scores_to_report(scores_distance, "similarity_max_distance")
    report.add_scores_to_report(scores_length, "similarity_path_length")
    # report.add_scores_to_report(scores_jaccard, "similarity_jaccard")
    

if __name__ == '__main__':
    main()
