""" assess similarity of HPO terms in sets of probands.

If we have a set of probands who we know share some genetic feature in a gene,
we would like to know what is the probability of them sharing sharing their
Human Phenotype Ontology (HPO; standardised phenotypes) terms.

The HPO terms form a graph, so in order to estimate similarity, we use common
ancestors of sets of terms.
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import argparse
import bisect
import math
import random

from src.load_files import load_participants_hpo_terms, load_variants
from src.create_hpo_graph import loadHPONetwork
from src.shared_term_plots import plot_shared_terms
from src.similarity import PathLengthSimilarity

HPO_PATH = os.path.join(os.path.dirname(__file__), "data", "hp.obo")
DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/"
PHENOTYPES_PATH = os.path.join(DATAFREEZE_DIR, "phenotypes_and_patient_info.txt")
ALTERNATE_IDS_PATH = os.path.join(DATAFREEZE_DIR, "person_sanger_decipher.txt")

def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description="Examines the likelihood of \
        obtaining similar HPO terms in probands with variants in the same gene.")
    parser.add_argument("--variants", dest="variants_path", required=True, \
        default="/nfs/users/nfs_j/jm33/apps/mupit/data-raw/de_novo_datasets/ \
            de_novos.ddd_4k.ddd_only.txt", \
        help="Path to file listing known variants in genes. See example file \
            in data folder for format.")
    parse.add_argument("--output", required=True, help="path to output file")
    
    args = parser.parse_args()
    
    return args

def get_ic(matcher, proband_term, other_terms):
    """ find the information content of the ancestor of the closest term
    
    For a given proband HPO term, we search the terms from another proband for
    the terms that are closest to the probands term. We use the path distance
    between terms to assess closeness. Sometimes we get multiple terms from the
    other proband that are equally close to the first probands term. We want
    something like the probability of obtaining this matching term.
    
    Rather than using the matching term itself, instead we use the closest
    common ancestor to the probands term and the matching HPO terms from the
    other proband.
    
    We check the ancestors to the matching terms for the one that has the
    highest information content (ie most rare). We return the information
    content.
    
    Args:
        matcher: PathLengthSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        proband_term: HPO term from a proband (e.g. "HP:000118").
        other_terms: list of HPO terms for a single proband e.g. ["HP:000110",
            "HP:000220", "HP000330"].
    
    Returns:
        The information content for the ancestor of the HPO term that is closest
        to the probands HPO term. Where multiple HPO terms are equally close to
        the probands term, use the one with the highest information content.
    """
    
    min_distance = 1e9
    best_terms = []
    
    for term in other_terms:
        path = matcher.get_shortest_path(proband_term, term)
        
        if len(path) < min_distance:
            min_distance = len(path)
            best_terms = [term]
        elif len(path) == min_distance:
            best_terms.append(term)
    
    # if we have multiple HPO terms that are equally close to the probands term,
    # then we find the highest information content for the ancestor closest to
    # the probands term
    ic = []
    for term in set(best_terms):
        common_terms = matcher.find_common_ancestors(proband_term, term)
        ancestor = matcher.find_closest_ancestor(proband_term, common_terms)
        ic.append(matcher.calculate_information_content(ancestor))
    
    return max(ic)

def geomean(values):
    """ calculate the geometric mean of a list of floats.
    
    Args:
        values: list of values, none of which are zero or less. The function
            will raise an error if the list is empty.
    
    Returns:
        The geometric mean as float.
    """
    
    values = [math.log10(x) for x in values]
    mean = sum(values)/float(len(values))
    
    return 10**mean
    
def get_proband_similarity(matcher, hpo_terms):
    """ calculate the similarity of HPO terms across different individuals.
    
    We start with a list of lists for each individual, where the list for each
    individual is a list of HPO terms. For each term in each individual, we find
    the closest matching terms from the other probands. Given the matching term,
    we get the information content of the term. We return the sum of the
    information content values from the different matches.
    
    Args:
        matcher: PathLengthSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        hpo_terms: list of HPO terms found for each proband with variants for
            the current gene.
    
    Returns:
        the summed negative log probability of observing the HPO terms in the
        population.
    """
    
    ic_scores = []
    for pos in range(len(hpo_terms)):
        proband = hpo_terms[pos]
        
        # remove the proband, so we don't match to itself
        others = hpo_terms[:]
        others.pop(pos)
        
        proband_ic_scores = []
        for other in others:
            
            # for each term in the proband, measure how well it matches the
            # terms in another proband
            temp_ic = []
            for term in proband:
                ic = get_ic(matcher, term, other)
                temp_ic.append(ic)
            
            # # TODO:  check out the impact of using max IC for a proband pair
            # proband_ic_scores.append(max(temp_ic))
            proband_ic_scores.append(geomean(temp_ic))
        
        ic_scores.append(geomean(proband_ic_scores))
    
    return sum(ic_scores)

def test_similarity(matcher, family_hpo, probands):
    """ find if groups of probands per gene share HPO terms more than by chance.
    
    Args:
        matcher: PathLengthSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        family_hpo: list of FamilyHPO objects for all probands.
        probands: list of proband IDs.
    
    Returns:
        The probability that the HPO terms used in the probands match as well as
        they do.
    """
    
    hpo_terms = [family_hpo[x].get_child_hpo() for x in probands if x in family_hpo]
    other_probands = [x for x in family_hpo if x not in probands]
    
    # We can't test similarity from a single proband. We don't call this
    # function for genes with a single proband, however, sometimes only one of
    # the probands has HPO terms recorded. We cannot estimate the phenotypic
    # similarity between probands in this case, so return None instead.
    if len(hpo_terms) == 1:
        return None
    
    observed = get_proband_similarity(matcher, hpo_terms)
    
    # get a distribution of scores for randomly sampled HPO terms
    distribution = []
    for x in range(1000):
        sampled = random.sample(other_probands, len(probands))
        simulated = [family_hpo[n].get_child_hpo() for n in sampled]
        predicted = get_proband_similarity(matcher, simulated)
        distribution.append(predicted)
    distribution = sorted(distribution)
    
    # figure out where in the distribution the observed value occurs
    pos = bisect.bisect_left(distribution, observed)
    sim_prob = (abs(pos - len(distribution)))/(1 + len(distribution))
    
    return sim_prob

def analyse_genes(matcher, family_hpo, probands_by_gene, output_path):
    """ tests genes to see if their probands share HPO terms more than by chance.
    
    Args:
        matcher: PathLengthSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        family_hpo: list of FamilyHPO objects for all probands
        probands_by_gene: dictionary of genes, to the probands who have variants
            in those genes.
        output_path: path to file to write the results to.
    """
    
    output = open(output_path, "w")
    output.write("hgnc\thpo_similarity_p_value\n")
    
    for gene in sorted(probands_by_gene):
        probands = probands_by_gene[gene]
        
        p_value = "NA"
        if len(probands) > 1:
            p_value = analyse_probands(matcher, family_hpo, probands)
        
        if p_value is None:
            p_value = "NA"
            
        output.write("{0}\t{1}\n".format(gene, p_value))
    
    output.close()

def main():
    
    options = get_options()
    
    # build a graph of DDG2P terms, so we can trace paths between terms
    hpo_file = loadHPONetwork(HPO_PATH)
    graph = hpo_file.get_graph()
    alt_node_ids = hpo_file.get_alt_ids()
    obsolete_ids = hpo_file.get_obsolete_ids()
    
    # load HPO terms and probands for each gene
    family_hpo_terms = load_participants_hpo_terms(PHENOTYPES_PATH, ALTERNATE_IDS_PATH, alt_node_ids, obsolete_ids)
    probands_by_gene = load_variants(options.variants_path)
    
    matcher = PathLengthSimilarity(family_hpo_terms, graph, alt_node_ids)
    matcher.tally_hpo_terms(family_hpo_terms, source="child_hpo")
    
    analyse_genes(matcher, family_hpo_terms, probands_by_gene, options.output)

if __name__ == '__main__':
    main()
