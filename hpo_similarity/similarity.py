"""
Copyright (c) 2015 Wellcome Trust Sanger Institute

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import math
from networkx import DiGraph

from hpo_similarity.check_proband_terms import check_terms_in_graph

class CalculateSimilarity(DiGraph):
    """ calculate graph similarity scores
    """
    
    def __init__(self):
        self.descendant_cache = {}
        self.ancestor_cache = {}
        
        self.total_freq = 0
        
        super(CalculateSimilarity, self).__init__()
    
    def tally_hpo_terms(self, hpo_terms):
        """ tallies each HPO term across the DDG2P genes
        
        Args:
            hpo_terms: dictionary of HPO terms for each individual
        """
        
        check_terms_in_graph(self, hpo_terms)
        
        for proband in hpo_terms:
            child_terms = hpo_terms[proband]
            for term in child_terms:
                self.add_proband_term(term, proband)
            
            self.total_freq += 1
    
    def get_ids_per_term(self, term):
        ''' pull out the sample IDs annotated with a specific term
        
        NOTE: this doesn't identify all of the individuals with ancestral terms,
        just the sparse set of sample IDs with that specific term.
        
        Args:
            term: HPO term
        
        Returns:
            set of sample IDs for individiuals with a specific term
        '''
        
        if 'sample_ids' not in self.node[term]:
            return set([])
        
        return self.node[term]['sample_ids']
    
    def add_proband_term(self, term, proband):
        """ adds a proband to the list of probands for an HPO term
        
        We need to keep track of which individuals had which terms. This allows
        us to later count probands with that specific term, or an ancestral term
        
        Args:
            term: HPO term (e.g. "HP:0000001")
            proband: sample ID for an individual
        """
        
        if term not in self:
            return
        
        if 'sample_ids' not in self.node[term]:
            self.node[term]['sample_ids'] = set([])
        
        self.node[term]['sample_ids'].add(proband)
    
    def get_descendants(self, term):
        """ finds the set of subterms that descend from a top level HPO term
        
        Args:
            term: hpo term to find descendants of
        
        Returns:
            set of descendant HPO terms
        """
        
        if term not in self.descendant_cache:
            terms = set(self.successors(term))
            
            # recursively extend up the graph until we hit the top node
            extra = [ self.get_descendants(x) for x in terms ]
            terms |= set([item for sublist in extra for item in sublist])
            
            self.descendant_cache[term] = terms
        
        return self.descendant_cache[term]
    
    def get_ancestors(self, bottom_term):
        """ finds the set of subterms that are ancestors of a HPO term
        
        NOTE: this also includes the search node in the list of ancestors. This
        is so that when we look for matches of common ancestors between two
        nodes, and the two node terms are for the same node, we also include the
        common node in the list. That was awkwardly phrased.
        
        Args:
            bottom_term: hpo term to find ancestors of
        
        Returns:
            set of ancestor HPO terms
        """
        
        if bottom_term not in self.ancestor_cache:
            subterms = set(self.predecessors(bottom_term))
            
            # recursively extend down the graph until we hit the bottom nodes
            extra = [ self.get_ancestors(x) for x in subterms ]
            subterms |= set([item for sublist in extra for item in sublist])
            
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
        if term_1 not in self or term_2 not in self:
            return set()
        
        return set(self.get_ancestors(term_1)) & set(self.get_ancestors(term_2))


class ICSimilarity(CalculateSimilarity):
    """ calculate similarity by IC score
    """
    
    most_informative = {}
    
    def get_most_informative_ic(self, term_1, term_2):
        """ calculate the information content between two HPO terms using the most informative common ancestor
        
        Args:
            term_1: hpo term, eg HP:0000003
            term_2: hpo term, eg HP:0000002
        
        Returns:
            the maximum information content value from the common ancestors of
            term_1 and term_2.
        """
        
        terms = (term_1, term_2)
        
        if terms not in self.most_informative:
            
            ancestors = self.find_common_ancestors(term_1, term_2)
            ic_values = [self.calculate_information_content(x) for x in ancestors]
            
            # cache the most informative IC value, so we only compute this once
            # per pair of HPO terms.
            most_informative = max(ic_values)
            self.most_informative[terms] = most_informative
            self.most_informative[(term_2, term_1)] = most_informative
        
        return self.most_informative[terms]
    
    def calculate_information_content(self, term):
        """ calculates the information content for an hpo term
        
        For discussion of information content and similarity scores, see:
        Van der Aalst et al., (2007) Data & Knowledge Engineering 61:137-152
        
        Args:
            term: hpo term, eg "HP:0000001"
        
        Returns:
            the information content value for a single hpo term
        """
        
        if term not in self:
            return 0
        
        if 'info_content' not in self.node[term]:
            term_count = self.get_term_count(term)
            
            # cache the IC, so we don't have to recalculate for the term
            self.node[term]['info_content'] = -math.log(term_count/self.total_freq)
        
        return self.node[term]['info_content']
    
    def get_term_count(self, term):
        """ Count how many times a term (or its subterms) was used.
        
        Args:
            term: hpo term, eg "HP:0000001"
        
        Returns:
            the number of times a term (or its subterms) was used.
        """
        
        if term not in self:
            return 0
        
        if 'count' not in self.node[term]:
            sample_ids = self.get_ids_per_term(term)
            for subterm in self.get_descendants(term):
                sample_ids |= self.get_ids_per_term(subterm)
            
            # cache the count, so we only have to calculate this once
            self.node[term]['count'] = len(sample_ids)
        
        return self.node[term]['count']
