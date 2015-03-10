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
import sys
import argparse
import bisect
import math
import random

from src.load_files import load_participants_hpo_terms, load_genes
from src.ontology import Ontology
from src.similarity import ICSimilarity
from src.permute_probands import permute_probands


def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description="Examines the likelihood of \
        obtaining similar HPO terms in probands with variants in the same gene.")
    parser.add_argument("--genes", dest="genes_path", required=True, \
        help="Path to JSON file listing probands per gene. See \
            data/example_genes.json for format.")
    parser.add_argument("--phenotypes", dest="phenotypes_path", required=True, \
        help="Path to JSON file listing phenotypes per proband. See \
            data/example_phenotypes.json for format.")
    parser.add_argument("--ontology", \
        default=os.path.join(os.path.dirname(__file__), "data", "hp.obo"), \
        help="path to HPO ontology obo file, see http://human-phenotype-ontology.org")
    parser.add_argument("--output", default=sys.stdout, \
        help="path to output file, defaults to standard out.")
    parser.add_argument("--permute", action="store_true", default=False,
        help="whether to permute the probands across genes, in order to assess \
            method robustness.")
    parser.add_argument("--iterations", type=int, default=100000,
        help="whether to permute the probands across genes, in order to assess \
            method robustness.")
    
    args = parser.parse_args()
    
    return args

def get_score_for_pair(hpo_graph, proband_1, proband_2):
    """ Calculate the similarity in HPO terms between terms for two probands.
    
    This runs through the pairs of HPO terms from the two probands and finds
    the information content for the most informative common ancestor for each
    pair. We return the largest of these IC scores, known as the maxIC.
    
    Args:
        hpo_graph: ICSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        proband_1: list of HPO terms for one proband
        proband_2: list of HPO terms for the other proband
    
    Returns:
        A score for how similar the terms are between the two probands.
    """
    
    ic = []
    for term_1 in proband_1:
        for term_2 in proband_2:
            ic.append(hpo_graph.get_most_informative_ic(term_1, term_2))
    
    return max(ic)
    
def get_proband_similarity(hpo_graph, probands):
    """ calculate the similarity of HPO terms across different individuals.
    
    We start with a list of HPO lists e.g. [[HP:01, HP:02], [HP:02, HP:03]],
    and calculate a matrix of similarity scores for each pair of probands in the
    HPO lists. We collapse that to a single score that estimates the similarity
    across all the probands.
    
    Args:
        hpo_graph: ICSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        probands: List of HPO terms found for each proband with variants for
            the current gene e.g. [[HP:01, HP:02], [HP:02, HP:03]].
    
    Returns:
        The summed similarity score across the HPO terms for each proband.
    """
    
    ic_scores = []
    for x in range(len(probands)):
        for y in range(x, len(probands)):
            # don't match a proband to itself
            if x == y:
                continue
            
            # for each term in the proband, measure how well it matches the
            # terms in another proband
            score = get_score_for_pair(hpo_graph, probands[x], probands[y])
            ic_scores.append(score)
    
    return sum(ic_scores)

def test_similarity(hpo_graph, hpo_by_proband, probands, n_sims):
    """ find if groups of probands per gene share HPO terms more than by chance.
    
    We simulate a distribution of similarity scores by randomly sampling groups
    of probands. I tried matching the number of sampled HPO terms to the numbers
    in the probands for the gene. For that, I gave each term the chance of being
    sampled as the rate at which it was observed in all the probands. However,
    these sampled terms gave abberant QQ plots, with excessive numbers of
    extremely signficant P values. I suspect this is due to underlying
    relationships between HPO terms.
    
    Args:
        hpo_graph: ICSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        hpo_by_proband: dictionary of HPO terms per proband
        probands: list of proband IDs.
        n_sims: number of simulations to run.
    
    Returns:
        The probability that the HPO terms used in the probands match as well as
        they do.
    """
    
    probands = [hpo_by_proband[x] for x in probands if x in hpo_by_proband]
    other_probands = [x for x in hpo_by_proband if x not in probands]
    
    # We can't test similarity from a single proband. We don't call this
    # function for genes with a single proband, however, sometimes only one of
    # the probands has HPO terms recorded. We cannot estimate the phenotypic
    # similarity between probands in this case, so return None instead.
    if len(probands) < 2:
        return None
    
    observed = get_proband_similarity(hpo_graph, probands)
    
    # get a distribution of scores for randomly sampled HPO terms
    distribution = []
    for x in range(n_sims):
        sampled = random.sample(other_probands, len(probands))
        simulated = [hpo_by_proband[n] for n in sampled]
        predicted = get_proband_similarity(hpo_graph, simulated)
        distribution.append(predicted)
    
    distribution = sorted(distribution)
    
    # figure out where in the distribution the observed value occurs
    pos = bisect.bisect_left(distribution, observed)
    sim_prob = (abs(pos - len(distribution)))/(1 + len(distribution))
    
    if sim_prob == 0:
        sim_prob = 1 / (1 + len(distribution))
    
    return sim_prob

def analyse_genes(hpo_graph, hpo_by_proband, probands_by_gene, output_path, iterations):
    """ tests genes to see if their probands share HPO terms more than by chance.
    
    Args:
        hpo_graph: ICSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        hpo_by_proband: dictionary of HPO terms per proband
        probands_by_gene: dictionary of genes, to the probands who have variants
            in those genes.
        output_path: path to file to write the results to, or sys.stdout object.
        iterations: number of iterations to run.
    """
    
    # Sometimes output_path is actually sys.stdout, other times it is a path.
    try:
        output = open(output_path, "w")
    except TypeError:
        output = output_path
    
    output.write("hgnc\thpo_similarity_p_value\n")
    
    for gene in sorted(probands_by_gene):
        probands = probands_by_gene[gene]
        
        p_value = None
        if len(probands) > 1:
            p_value = test_similarity(hpo_graph, hpo_by_proband, probands, iterations)
        
        if p_value is None:
            continue
        
        output.write("{0}\t{1}\n".format(gene, p_value))
    
    output.close()

def main():
    
    options = get_options()
    
    # build a graph of HPO terms, so we can trace paths between terms
    hpo_ontology = Ontology(options.ontology)
    hpo_graph = hpo_ontology.get_graph()
    alt_node_ids = hpo_ontology.get_alt_ids()
    obsolete_ids = hpo_ontology.get_obsolete_ids()
    
    # load HPO terms and probands for each gene
    print("loading HPO terms and probands by gene")
    hpo_by_proband = load_participants_hpo_terms(options.phenotypes_path, \
        alt_node_ids, obsolete_ids)
    probands_by_gene = load_genes(options.genes_path)
    
    if options.permute:
        probands_by_gene = permute_probands(probands_by_gene)
    
    hpo_graph = ICSimilarity(hpo_by_proband, hpo_graph, alt_node_ids)
    
    print("analysing similarity")
    try:
        analyse_genes(hpo_graph, hpo_by_proband, probands_by_gene, \
            options.output, options.iterations)
    except KeyboardInterrupt:
        sys.exit("HPO similarity exited.")

if __name__ == '__main__':
    main()
