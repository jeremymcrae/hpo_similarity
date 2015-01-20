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
import argparse

from src.load_files import load_participants_hpo_terms, load_de_novos
from src.create_hpo_graph import loadHPONetwork
from src.similarity import PathLengthSimilarity
from src.weighted_choice import WeightedChoice

USER_PATH = "/nfs/users/nfs_j/jm33/"
HPO_FOLDER = os.path.join(USER_PATH, "apps", "hpo_filter")
HPO_PATH = os.path.join(HPO_FOLDER, "hpo_data", "hp.obo")

ddd_freeze = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/"
PHENOTYPES_PATH = os.path.join(ddd_freeze, "phenotypes_and_patient_info.txt")
ALTERNATE_IDS_PATH = os.path.join(ddd_freeze, "person_sanger_decipher.txt")

def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description="Examines the likelihood of \
        obtaining similar HPO terms in probands with de novos in the same gene.")
    parser.add_argument(dest="input", \
        default="/nfs/users/nfs_j/jm33/apps/mupit/data-raw/de_novo_datasets/de_novos.ddd_4k.ddd_only.txt", \
        help="Path to file listing known mutations in genes. See example file \
            in data folder for format.")
    
    args = parser.parse_args()
    
    return args.input

def product(values):
    
    total = 1
    for x in values:
        total *= x
    
    return total

def find_closest_ancestor_prob(matcher, proband_term, nonproband_terms, family_hpo_terms):
    
    min_distance = 1e9
    best_terms = []
    
    for term in nonproband_terms:
        
        # some HPO terms are obsolete, and cannot be placed on the graph eg
        # HP:0001113. We shall simply ignore these HPO terms
        try:
            path = matcher.get_shortest_path(proband_term, term)
        except KeyError:
            continue
        
        if len(path) < min_distance:
            min_distance = len(path)
            best_terms = [term]
        elif len(path) == min_distance:
            best_terms.append(term)
    
    ic = []
    for term in set(best_terms):
        common_terms = matcher.find_common_ancestors(proband_term, term)
        ancestor = matcher.find_closest_ancestor(proband_term, common_terms)
        ic.append(matcher.calculate_information_content(ancestor))
    
    return max(ic)

def analyse_gene(matcher, sampler, family_hpo_terms, probands):
    
    if len(probands) == 1:
        return None
    
    hpo_terms = [family_hpo_terms[proband].get_child_hpo() for proband in probands if proband in family_hpo_terms]
    
    simulated_list = []
    for proband_list in hpo_terms:
        a = [sampler.choice() for x in range(0, len(proband_list))]
        simulated_list.append(a)
    
    log_probs = []
    for pos in range(len(hpo_terms)):
        proband = hpo_terms[pos]
        nonprobands = [x for i,x in enumerate(hpo_terms) if i > pos]
        for term in proband:
            for nonproband in nonprobands:
                ic = find_closest_ancestor_prob(matcher, term, nonproband, family_hpo_terms)
                log_probs.append(ic)
    
    log_prob = sum(log_probs)
    # if prob == 0:
    #     print(prob, probs)
    
    return log_prob

def analyse_proband_similarity(matcher, alt_node_ids, family_hpo_terms, de_novos):
    
    # get the prior probability that each HPO term will be seleceted
    max_rate = 0
    rates = []
    for term in matcher.hpo_counts:
        rate = matcher.get_term_count(term)/len(family_hpo_terms)
        rates.append((term, rate))
    
    # get weighted sampler
    sampler = WeightedChoice(rates)
    
    max_prob = 0
    for gene in sorted(de_novos):
        probands = de_novos[gene]
        log_prob = analyse_gene(matcher, sampler, family_hpo_terms, probands)
        
        if log_prob > 500:
            max_prob = log_prob
            print(gene)

def main():
    
    de_novo_path = get_options()
    
    # build a graph of DDG2P terms, so we can trace paths between terms
    hpo_file = loadHPONetwork(HPO_PATH)
    graph = hpo_file.get_graph()
    alt_node_ids = hpo_file.get_alt_ids()
    obsolete_ids = hpo_file.get_obsolete_ids()
    
    print(obsolete_ids)
    
    # load HPO terms and de novos per gene
    family_hpo_terms = load_participants_hpo_terms(PHENOTYPES_PATH, ALTERNATE_IDS_PATH, alt_node_ids, obsolete_ids)
    de_novos = load_de_novos(de_novo_path)
    
    # print(family_hpo_terms)
    #
    # for family in family_hpo_terms:
    #     print(family_hpo_terms[family].get_child_hpo())
    
    matcher = PathLengthSimilarity(family_hpo_terms, graph, alt_node_ids)
    matcher.tally_hpo_terms(family_hpo_terms, source="child_hpo")
    analyse_proband_similarity(matcher, alt_node_ids, family_hpo_terms, de_novos)
    
    # scores_length = matcher.get_similarity_scores(probands)

if __name__ == '__main__':
    main()
