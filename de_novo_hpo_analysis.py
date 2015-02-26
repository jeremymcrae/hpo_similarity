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
import bisect
import itertools
import math
import random

from matplotlib import use
use("Agg")
from matplotlib import pyplot
from matplotlib import lines
import networkx

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

def get_ic_from_closest_ancestor(matcher, proband_term, other_terms):
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
        values: list of values, none of which are zero or less. The list must
            contain at least some values.
    
    Returns:
        The geometric mean as float.
    """
    
    values = [math.log10(x) for x in values]
    mean = sum(values)/float(len(values))
    
    return 10**mean
    
def calculate_hpo_prob(matcher, hpo_terms):
    """ calculate the total information content across a list of HPO lists.
    
    We start with a list of lists for each proband, where the list for each
    proband is a list of HPO terms. For each term in each proband, we find the
    closest matching terms from the other probands. Given the matching term, we
    get the information content of the term. We return the sum of the
    information content values from the different matches.
    
    Args:
        matcher: PathLengthSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        hpo_terms: HPO terms found for each proband with de novos for the
            current gene.
    
    Returns:
        the summed negative log probability of observing the HPO terms in the
        population.
    """
    
    log_probs = []
    for pos in range(len(hpo_terms)):
        proband = hpo_terms[pos]
        
        # remove the proband, and and preceeding probands, so we
        #   a) don't try to check for matches from an identical list
        #   b) don't check probands for whom we already have results
        others = [x for i,x in enumerate(hpo_terms) if i != pos]
        proband_probs = []
        for other in others:
            other_probs = []
            for term in proband:
                ic = get_ic_from_closest_ancestor(matcher, term, other)
                other_probs.append(ic)
            proband_probs.append(geomean(other_probs))
        
        log_probs.append(geomean(proband_probs))
    
    return sum(log_probs)

def analyse_probands(matcher, family_hpo, probands):
    """ get the probability of HPO terms matching in a set of probands
    
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
    
    # Sometimes only one of the probands has HPO terms recorded. We cannot
    # estimate the phenotypic similarity between probands when we only have one
    # proband, so return None instead.
    if len(hpo_terms) == 1:
        return None
    else:
        return 1
    
    observed = calculate_hpo_prob(matcher, hpo_terms)
    
    # get a distribution of scores for randomly sampled HPO terms
    distribution = []
    for x in range(100):
        sampled = random.sample(other_probands, len(probands))
        simulated = [family_hpo[n].get_child_hpo() for n in sampled]
        predicted = calculate_hpo_prob(matcher, simulated)
        distribution.append(predicted)
    distribution = sorted(distribution)
    
    # figure out where in the distribution the observed value occurs
    pos = bisect.bisect_left(distribution, observed)
    sim_prob = (abs(pos - len(distribution)))/(1 + len(distribution))
    
    return sim_prob

def analyse_de_novos(matcher, family_hpo, de_novos):
    """ find if probands with de novos in each gene share HPO terms more than by chance.
    
    In order to assess a gene
    
    Args:
        matcher: PathLengthSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        family_hpo: list of FamilyHPO objects for all probands
        de_novos: dictionary of genes, to the probands who have de novos in
            those genes
    """
    
    output = open("results/de_novo_analysis.test.txt", "w")
    for gene in sorted(de_novos):
        probands = de_novos[gene]
        
        p_value = "NA"
        if len(probands) > 1:
            p_value = analyse_probands(matcher, family_hpo, probands)
        
        if p_value is None:
            p_value = "NA"
            
        output.write("{0}\t{1}\n".format(gene, p_value))
    
    output.close()

def plot_graph(gene, matcher, hpo):
    """ plots the HPO terms used in a set of probands as a tree
    
    Args:
        gene: gene name as string
        matcher: PathLengthSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        hpo: list of hpo lists for probands
    """
    g = matcher.graph.copy()
    
    all_hpo = [item for sublist in hpo for item in sublist]
    all_paths = get_all_paths_between_hpo(matcher, hpo)
    all_terms = [item for sublist in all_paths for item in sublist]
    
    # get a graph that has unneeded terms removed
    nodes = matcher.graph.nodes()
    [g.remove_node(node) for node in nodes if node not in set(all_terms)]
    
    # plot the network as a hierarchy (ie tree) of nodes,
    pos = networkx.graphviz_layout(g, prog='dot')
    networkx.draw_networkx_nodes(g, pos, with_labels=False, node_color="white", \
        node_size=50, alpha=0.5)
    networkx.draw_networkx_edges(g, pos, width=0.4, alpha=0.5, arrows=False)
    
    # scale the size of the plotted point by how many times the term was used
    terms = sorted(set(all_hpo))
    sizes = [50 + 50 * math.log(all_hpo.count(x), 2) for x in terms]
    
    # shade the plotted points so that rarer terms are more intensely shaded
    ic = [matcher.calculate_information_content(x) for x in terms]
    colors = ["#{0:02x}{0:02x}{0:02x}".format(int(255 - (23 * x))) for x in ic]
    
    # now draw the HPO terms that were used in the probands
    networkx.draw_networkx_nodes(g, pos, nodelist=terms, node_size=sizes, \
        node_color=colors)
    
    # label the HPO terms that were used in the probands, adjust the y position
    # so that the label will sit above the plotted node
    labels = dict(zip(terms, range(1, len(terms) + 1)))
    label_pos = {x: (pos[x][0], pos[x][1] + 20) for x in pos}
    networkx.draw_networkx_labels(g, label_pos, nodelist=terms, labels=labels, \
        font_size=7, font_color="red")
    
    path = os.path.join("results", "{0}_test.pdf".format(gene))
    
    # get the label parameters
    definitions = ["{0} - {1} ({2})".format(terms.index(x) + 1, g.node[x]["name"][0], x) for x in terms]
    artists = [lines.Line2D(range(1), range(1), marker="", color="white") for x in terms]
    
    pyplot.legend(artists, definitions, loc=2, fontsize=4, frameon=False)
    pyplot.axis("off")
    pyplot.title(gene, loc="center")
    pyplot.savefig(path, bbox_inches="tight", pad_inches=0)
    pyplot.close()

def get_all_paths_between_hpo(matcher, hpo):
    """ gets lists of HPO terms for all paths between probands HPO terms
    
    Args:
        matcher: PathLengthSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        hpo: list of hpo lists for probands
    """
    
    terms = [item for sublist in hpo for item in sublist]
    
    pairs = itertools.combinations(terms, 2)
    paths = []
    for (term1, term2) in pairs:
        path = matcher.get_shortest_path(term1, term2)
        paths.append(path)
    
    return paths

def main():
    
    de_novo_path = get_options()
    
    # build a graph of DDG2P terms, so we can trace paths between terms
    hpo_file = loadHPONetwork(HPO_PATH)
    graph = hpo_file.get_graph()
    alt_node_ids = hpo_file.get_alt_ids()
    obsolete_ids = hpo_file.get_obsolete_ids()
    
    # load HPO terms and de novos per gene
    family_hpo_terms = load_participants_hpo_terms(PHENOTYPES_PATH, ALTERNATE_IDS_PATH, alt_node_ids, obsolete_ids)
    de_novos = load_de_novos(de_novo_path)
    
    matcher = PathLengthSimilarity(family_hpo_terms, graph, alt_node_ids)
    matcher.tally_hpo_terms(family_hpo_terms, source="child_hpo")
    analyse_de_novos(matcher, family_hpo_terms, de_novos)

if __name__ == '__main__':
    main()
