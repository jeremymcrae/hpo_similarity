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

import os
import unittest

from hpo_similarity.ontology import Ontology
from hpo_similarity.similarity import CalculateSimilarity

class TestCalculateSimilarityPy(unittest.TestCase):
    """ class to test CalculateSimilarity
    """
    
    def setUp(self):
        """ construct a CalculateSimilarity object for unit tests
        """
        
        path = os.path.join(os.path.dirname(__file__), "data", "obo.txt")
        ontology = Ontology(path)
        self.graph = ontology.get_graph()
        
        self.hpo_terms = {
            "person_01": ["HP:0000924"],
            "person_02": ["HP:0000118", "HP:0002011"],
            "person_03": ["HP:0000707", "HP:0002011"]
        }
        
        self.graph.tally_hpo_terms(self.hpo_terms)
        
    def test_setup(self):
        """ test that the class initialised correctly.
        
        Mainly I want to check that when the class initialised, it ran
        tally_hpo_terms() correctly. Check that the counts of the HPO terms
        used in the probands match what is expected.
        """
        
        self.assertEqual(self.graph.total_freq, 3)
        self.assertEqual(self.graph.get_ids_per_term("HP:0002011"),
            {'person_02', 'person_03'} )
        
        # check that a redundant term has been added, even though a more specific
        # descendant term was included
        self.assertTrue('sample_ids' in self.graph.node['HP:0000118'])
        
        # Check that we get an error if we look for probands with a term that was
        # not used in the probands.
        with self.assertRaises(KeyError):
            self.graph.node["HP:0000001"]['sample_ids']
        
        # but a similar check using the official method returns an empty set
        self.assertEqual(self.graph.get_ids_per_term("HP:0000001"), set([]))
    
    def test_add_proband_term(self):
        """ check that HPO counting works correctly
        """
        
        # check the baseline count for a term
        self.assertEqual(self.graph.get_ids_per_term("HP:0002011"),
            {'person_02', 'person_03'})
        
        # add a term, and check that the count for the term increases, but
        # the total frequency doesn't change
        self.graph.add_proband_term("HP:0002011", 'person_01')
        self.assertEqual(self.graph.get_ids_per_term("HP:0002011"),
            {'person_01', 'person_02', 'person_03'})
        self.assertEqual(self.graph.total_freq, 3)
        
        # add a term for a proband which has already been included, and check
        # that the count has not changed
        self.graph.add_proband_term("HP:0002011", 'person_01')
        self.assertEqual(self.graph.get_ids_per_term("HP:0002011"),
            {'person_01', 'person_02', 'person_03'})
        
        # check that if we try to add a term that isn't in the HPO ontology, we
        # don't increment any counts
        self.graph.add_proband_term("unknown_term", 'person_01')
        self.assertEqual(self.graph.total_freq, 3)
        
        # Check that if we add a term that currently doesn't have a tallied
        # count then the term gets inserted correctly, and the counts increment
        # appropriately.
        with self.assertRaises(KeyError):
            self.assertEqual(self.graph.node["HP:0000001"]['sample_ids'])
        
        self.graph.add_proband_term("HP:0000001", 'person_01')
        self.assertEqual(self.graph.get_ids_per_term("HP:0000001"), {'person_01'})
        self.assertEqual(self.graph.total_freq, 3)
    
    def test_get_descendants(self):
        """ check that get_descendants works correctly
        """
        
        # check that a high-level node returns the expected set of nodes
        self.assertEqual(self.graph.get_descendants("HP:0000118"), \
            set(['HP:0000707', 'HP:0002011', 'HP:0000924']))
        
        # check that a terminal node doesn't have any descendants
        self.assertEqual(self.graph.get_descendants("HP:0000924"), \
            set([]))
    
    def test_get_ancestors(self):
        """ check that get_ancestors works correctly
        """
        
        # check that we get an appropriate set of ancestor tersm for a termina
        # node
        self.assertEqual(self.graph.get_ancestors("HP:0000924"), \
            set(['HP:0000001', 'HP:0000118', 'HP:0000924']))
        
        # check that even the top node returns itself as a ancestor node
        self.assertEqual(self.graph.get_ancestors("HP:0000001"), \
            set(['HP:0000001']))
    
    def test_find_common_ancestors(self):
        """ check that find_common_ancestors works correctly
        """
        
        # check that two terms on different arms only return their common
        # ancestors
        self.assertEqual(self.graph.find_common_ancestors('HP:0000924', \
            'HP:0000707'), set(["HP:0000001", "HP:0000118"]))
        
        # check that two identical terms return their list of ancestors
        self.assertEqual(self.graph.find_common_ancestors('HP:0000707', \
            'HP:0000707'), set(["HP:0000001", "HP:0000118", "HP:0000707"]))
        
        # check that if one of the two terms is not in the HPO graqph, then we
        # return an empty set
        self.assertEqual(self.graph.find_common_ancestors('HP:9999999', \
            'HP:0000707'), set([]))
