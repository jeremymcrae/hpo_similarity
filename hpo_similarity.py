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
import itertools

from src.load_files import load_participants_hpo_terms, load_variants
from src.ontology import Ontology
from src.similarity import ICSimilarity
from src.permute_probands import permute_probands


def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description="Examines the likelihood of \
        obtaining similar HPO terms in probands with variants in the same gene.")
    parser.add_argument("--variants", dest="variants_path", required=True, \
        help="Path to file listing probands per gene. See \
            data/example_variants.json for format.")
    parser.add_argument("--phenotypes", dest="phenotypes_path", required=True, \
        help="Path to file listing phenotypes per proband. See \
            data/example_phenotypes.json for format.")
    parser.add_argument("--ontology", \
        default=os.path.join(os.path.dirname(__file__), "data", "hp.obo"), \
        help="path to HPO ontology obo file, see http://human-phenotype-ontology.org/.")
    parser.add_argument("--output", required=True, help="path to output file")
    parser.add_argument("--permute", action="store_true", default=False,
        help="whether to permute the probands across genes, in order to assess \
            method robustness.")
    
    args = parser.parse_args()
    
    return args

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

def get_score_for_pair(matcher, proband_1, proband_2):
    """ Calculate the similarity in HPO terms between terms for two probands.
    
    Currently we use the geometric mean of the max IC values for all the pairs
    of terms between the two probands. I trialled the maximum rather than the
    geometric mean, but the QQ plots were distorted, and the P values were
    excessively correlated with P values from gene enrichment testing, that is,
    the P values were not independent of the enrichment tests.
    
    Args:
        matcher: ICSimilarity object for the HPO term graph, with
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
            ic.append(matcher.get_max_ic(term_1, term_2))
    
    return geomean(ic)
    
def get_proband_similarity(matcher, probands):
    """ calculate the similarity of HPO terms across different individuals.
    
    We start with a list of HPO lists e.g. [[HP:01, HP:02], [HP:02, HP:03]],
    and calculate a matrix of similarity scores for each pair of probands in the
    HPO lists. We collapse that to a single score that estimates the similarity
    across all the probands.
    
    Args:
        matcher: ICSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        probands: List of HPO terms found for each proband with variants for
            the current gene e.g. [[HP:01, HP:02], [HP:02, HP:03]].
    
    Returns:
        The summed similarity score across the HPO terms for each proband.
    """
    
    ic_scores = []
    for pos in range(len(probands)):
        proband = probands[pos]
        
        # remove the proband, so we don't match to itself
        others = probands[:]
        others.pop(pos)
        
        for other in others:
            # for each term in the proband, measure how well it matches the
            # terms in another proband
            score = get_score_for_pair(matcher, proband, other)
            ic_scores.append(score)
    
    return sum(ic_scores)

def test_similarity(matcher, hpo_by_proband, probands, n_sims=1000):
    """ find if groups of probands per gene share HPO terms more than by chance.
    
    We simulate a distribution of similarity scores by randomly sampling groups
    of probands. I tried matching the number of sampled HPO terms to the numbers
    in the probands for the gene. For that, I gave each term the chance of being
    sampled as the rate at which it was observed in all the probands. However,
    these sampled terms gave abberant QQ plots, with excessive numbers of
    extremely signficant P values. I suspect this is due to underlying
    relationships between HPO terms.
    
    Args:
        matcher: ICSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        hpo_by_proband: dictionary of HPO terms per proband
        probands: list of proband IDs.
    
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
    
    observed = get_proband_similarity(matcher, probands)
    
    # get a distribution of scores for randomly sampled HPO terms
    distribution = []
    for x in range(n_sims):
        sampled = random.sample(other_probands, len(probands))
        simulated = [hpo_by_proband[n] for n in sampled]
        predicted = get_proband_similarity(matcher, simulated)
        distribution.append(predicted)
    
    distribution = sorted(distribution)
    
    # figure out where in the distribution the observed value occurs
    pos = bisect.bisect_left(distribution, observed)
    sim_prob = (abs(pos - len(distribution)))/(1 + len(distribution))
    
    return sim_prob

def analyse_genes(matcher, hpo_by_proband, probands_by_gene, output_path):
    """ tests genes to see if their probands share HPO terms more than by chance.
    
    Args:
        matcher: ICSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        hpo_by_proband: dictionary of HPO terms per proband
        probands_by_gene: dictionary of genes, to the probands who have variants
            in those genes.
        output_path: path to file to write the results to.
    """
    
    output = open(output_path, "w")
    output.write("hgnc\thpo_similarity_p_value\n")
    
    for gene in sorted(probands_by_gene):
        probands = probands_by_gene[gene]
        
        p_value = None
        if len(probands) > 1:
            p_value = test_similarity(matcher, hpo_by_proband, probands)
        
        if p_value is None:
            continue
        
        print(gene)
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
    hpo_by_proband = load_participants_hpo_terms(options.phenotypes_path, alt_node_ids, obsolete_ids)
    probands_by_gene = load_variants(options.variants_path)
    
    if options.permute:
        probands_by_gene = permute_probands(probands_by_gene)
    
    matcher = ICSimilarity(hpo_by_proband, hpo_graph, alt_node_ids)
    
    print("analysing similarity")
    analyse_genes(matcher, hpo_by_proband, probands_by_gene, options.output)

if __name__ == '__main__':
    main()
