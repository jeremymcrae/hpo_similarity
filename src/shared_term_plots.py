""" plots a graph of HPO terms used within a group of probands.
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals


import os
import math
import itertools

from matplotlib import use
use("Agg")
from matplotlib import pyplot
from matplotlib import lines
import networkx


def plot_shared_terms(gene, matcher, hpo):
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
    
    terms = sorted(set(all_hpo))
    shade_used_nodes(g, terms, all_hpo, pos)
    
    # label the HPO terms that were used in the probands, adjust the y position
    # so that the label will sit above the plotted node
    labels = dict(zip(terms, range(1, len(terms) + 1)))
    label_pos = {x: (pos[x][0], pos[x][1] + 20) for x in pos}
    networkx.draw_networkx_labels(g, label_pos, nodelist=terms, labels=labels, \
        font_size=7, font_color="red")
    
    save_figure(g, terms, gene)

def shade_used_nodes(g, terms, all_hpo, pos):
    """ replot the nodes for the terms used in the probands
    
    Plot the used nodes, so the size represents how many times the term was used
    among the group of probands, and the intensity of shading represents the
    rarity of the term in the population.
    
    Args:
        g: networkx graph object, with unsed terms removed.
        terms: list of terms used in the probands.
        all_hpo: list of all HPO terms used in the population.
    """
    
    # scale the size of each plotted node by how many times the term was used
    sizes = [50 + 50 * math.log(all_hpo.count(x), 2) for x in terms]
    
    # shade the plotted points so that rarer terms are more intensely shaded
    ic = [matcher.calculate_information_content(x) for x in terms]
    colors = ["#{0:02x}{0:02x}{0:02x}".format(int(255 - (23 * x))) for x in ic]
    
    # now draw the HPO terms that were used in the probands
    networkx.draw_networkx_nodes(g, pos, nodelist=terms, node_size=sizes, \
        node_color=colors)

def save_figure(g, terms, gene):
    """ exports a file, and makes a legend of HPO terms and their definitions.
    
    Args:
        g: networkx graph object, with unsed terms removed.
        terms: list of terms used in the probands.
        gene: HGNC symbol for the current gene.
    """
    
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
