""" determines similarity scores between proband HPO terms and the DDG2P
genes predicted to be relevant to them.
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import math
import networkx

class CalculateSimilarity(object):
    """ calculate graph similarity scores
    """
    
    def __init__(self, family_hpos, hpo_graph, alt_node_ids):
        """
        
        Args:
            family_hpos: hpo term object, to get terms for probands, or parents
            graph: graph of hpo terms, as networkx object
        """
        
        self.graph = hpo_graph
        self.family_hpos = family_hpos
        self.alt_node_ids = alt_node_ids
        
        self.graph_depth = networkx.eccentricity(self.graph, v="HP:0000001")
        
        self.descendant_cache = {}
        self.ancestor_cache = {}
        
        self.hpo_counts = {}
        self.total_freq = 0
    
    def fix_alternate_id(self, term):
        """ converts HPO terms using alternate IDs to the standard term
        
        some of the HPO terms recorded for the probands, or in the DDG2P
        database are alternate IDs for the HPO terms. If so, swap over to the
        standard HPO term, as these are the node names in the HPO graph.
        """
        
        if self.graph.has_node(term):
            pass
        elif term in self.alt_node_ids:
            term = self.alt_node_ids[term]
        else:
            raise KeyError
        
        return term
    
    def tally_hpo_terms(self, hpo_terms, source="ddg2p"):
        """ tallies each HPO term across the DDG2P genes
        
        Args:
            hpo_terms: hpo term dictionary
        """
        
        for item in hpo_terms:
            child_terms = hpo_terms[item].get_child_hpo()
            for term in child_terms:
                self.add_hpo(term)
    
    def add_hpo(self, term):
        """ increments the count for an HPO term
        
        This increments a) the count for the specific term, and b) the total
        count of all terms.
        
        Args:
            term: string for HPO term
        """
        
        if term not in self.hpo_counts:
            # don't use terms which cannot be placed on the graph
            if not self.graph.has_node(term):
                return
            
            if term not in self.hpo_counts:
                self.hpo_counts[term] = 0
            
            self.hpo_counts[term] += 1
        
        self.total_freq += 1
    
    def get_subterms(self, top_term):
        """ finds the set of subterms that descend from a top level HPO term
        
        Args:
            top_term: hpo term to find descendants of
        
        Returns:
            set of descendant HPO terms
        """
        
        if top_term not in self.descendant_cache:
            self.descendant_cache[top_term] = networkx.descendants(self.graph, top_term)
        
        return self.descendant_cache[top_term]
    
    def get_ancestors(self, bottom_term):
        """ finds the set of subterms that are ancestors of a HPO term
        
        Args:
            bottom_term: hpo term to find ancestors of
        
        Returns:
            set of ancestor HPO terms
        """
        
        if bottom_term not in self.ancestor_cache:
            subterms = networkx.ancestors(self.graph, bottom_term)
            subterms.add(bottom_term)
            self.ancestor_cache[bottom_term] = subterms
        
        return self.ancestor_cache[bottom_term]
    
    def find_common_ancestors(self, term_1, term_2):
        """ finds the common ancestors of two hpo terms
        
        Args:
            term_1: hpo term, eg HP:0000002
            term_2: hpo term, eg HP:0000003
        
        Returns:
            a list of all the common ancestors for the two terms
        """
        
        # ignore terms that are obsolete (ie are not in the graph)
        if term_1 not in self.graph or term_2 not in self.graph:
            return set()
        
        return set(self.get_ancestors(term_1)) & set(self.get_ancestors(term_2))


class ICSimilarity(CalculateSimilarity):
    """ calculate similarity by IC score
    """
    
    counts_cache = {}
    ic_cache = {}
    
    def calculate_information_content(self, term):
        """ calculates the information content for an hpo term
        
        For discussion of information content and similarity scores, see:
        Van der Aalst et al., (2007) Data & Knowledge Engineering 61:137-152
        
        Args:
            term: hpo term, eg "HP:0000001"
        
        Returns:
            the information content value for a single hpo term
        """
        
        if term not in self.ic_cache:
            term_count = self.get_term_count(term)
            
            if term not in self.graph:
                return 0
            
            # cache the IC, so we don't have to recalculate for the term
            self.ic_cache[term] = -math.log(term_count/self.total_freq)
        
        return self.ic_cache[term]
    
    def get_term_count(self, term):
        """ Count how many times a term (or its subterms) was used.
        
        Args:
            term: hpo term, eg "HP:0000001"
        
        Returns:
            the number of times a term (or its subterms) was used.
        """
        
        if term not in self.counts_cache:
            if term not in self.graph:
                return 0
            
            descendants = self.get_subterms(term)
            
            count = 0
            if term in self.hpo_counts:
                count += self.hpo_counts[term]
            for subterm in descendants:
                if subterm in self.hpo_counts:
                    count += self.hpo_counts[subterm]
            
            self.counts_cache[term] = count
        
        return self.counts_cache[term]

class PathLengthSimilarity(ICSimilarity):
    """ calculate similarity score by path length
    """
    
    path_cache = {}
    
    def get_shortest_path(self, term_1, term_2):
        """ finds the shortest path between two terms (and caches the result)
        
        Args:
            term_1: HPO ID for graph node
            term_2: HPO ID for graph node
        
        Returns:
            list of nodes for path
        """
        
        term_1 = self.fix_alternate_id(term_1)
        term_2 = self.fix_alternate_id(term_2)
        
        if (term_1, term_2) not in self.path_cache:
            try:
                path = networkx.shortest_path(self.graph, term_1, term_2)
            except networkx.exception.NetworkXNoPath:
                path = self.get_path_between_nondescendants(term_1, term_2)
                
            self.path_cache[(term_1, term_2)] = path
            self.path_cache[(term_2, term_1)] = path
        
        return self.path_cache[(term_1, term_2)]
    
    def get_path_between_nondescendants(self, term_1, term_2):
        """ gets the shortest path between terms not from the same branch
        
        Args:
            term_1: HPO ID for graph node
            term_2: HPO ID for graph node
        
        Returns:
            shortest_path: list of nodes for path
        """
        
        common_terms = self.find_common_ancestors(term_1, term_2)
        ancestor = self.find_closest_ancestor(term_1, common_terms)
        
        # get the paths from the two terms to their closest ancestor
        path_1 = self.get_shortest_path(ancestor, term_1)[::-1]
        path_2 = self.get_shortest_path(ancestor, term_2)
        
        return path_1[:-1] + path_2
    
    def find_closest_ancestor(self, node, ancestors):
        """ finds the closest ancestor of a term from a list of ancestor terms
        
        Args:
            node: node ID to search from
            ancestors: list of ancestor node IDs
        
        Returns:
            ancestor_to_use: closest ancestor node ID
        """
        
         # find the path to the closest common ancestor
        shortest = None
        ancestor_to_use = ""
        for ancestor in ancestors:
            length = len(self.get_shortest_path(ancestor, node))
            if shortest is None or length < shortest:
                shortest = length
                ancestor_to_use = ancestor
        
        return ancestor_to_use
