""" class to test ICSimilarity class
"""

import os
import unittest

from hpo_similarity.ontology import Ontology
from hpo_similarity.similarity import ICSimilarity

class TestICSimilarityPy(unittest.TestCase):
    """ class to test ICSimilarity
    """
    
    def setUp(self):
        """ construct a ICSimilarity object for unit tests
        """
        
        path = os.path.join(os.path.dirname(__file__), "data", "obo.txt")
        ontology = Ontology(path)
        graph = ontology.get_graph()
        
        self.hpo_terms = {
            "person_01": ["HP:0000924"],
            "person_02": ["HP:0000118", "HP:0002011"],
            "person_03": ["HP:0000707", "HP:0002011"]
        }
        
        self.hpo_graph = ICSimilarity(self.hpo_terms, graph)
    
    def test_get_term_count(self):
        """ check that get_term_count works correctly
        
        All of the counts here are derived from their usage in self.hpo_terms
        """
        
        # check that we count the term usage (and subterms correctly)
        self.assertEqual(self.hpo_graph.get_term_count("HP:0000118"), 5)
        self.assertEqual(self.hpo_graph.get_term_count("HP:0000707"), 3)
        self.assertEqual(self.hpo_graph.get_term_count("HP:0002011"), 2)
        
        # check that a terminal node, only used once in the probands, has a
        # count of 1
        self.assertEqual(self.hpo_graph.get_term_count("HP:0000924"), 1)
        
        # check the term/subterm count for a term that isn't used within any of
        # he probands, but which all of the used terms descend from.
        self.assertEqual(self.hpo_graph.get_term_count("HP:0000001"), 5)
    
    def test_calculate_information_content(self):
        """ check that calculate_information_content works correctly
        """
        
        # check that the top node has an information content of 0
        self.assertEqual(self.hpo_graph.calculate_information_content("HP:0000001"), \
            0)
        
        # check the information content for a terminal node
        self.assertAlmostEqual(self.hpo_graph.calculate_information_content("HP:0000924"), \
            1.6094379)
        
        # check the information content for a node that is somewhat distant, but
        # which has some descendant nodes that need to be included in the term
        # count
        self.assertAlmostEqual(self.hpo_graph.calculate_information_content("HP:0000707"), \
            0.5108256)
        
    def test_get_most_informative_ic(self):
        """ check that get_most_informative_ic works correctly
        """
        
        # check the most informative information content for two nodes where
        # every common ancestor is the ancestor of all terms used in the probands
        self.assertAlmostEqual(self.hpo_graph.get_most_informative_ic("HP:0000707", \
            "HP:0000924"), 0)
        
        # check the most informative information content for two nodes where
        # both nodes are somewhat down the HPO graph
        self.assertAlmostEqual(self.hpo_graph.get_most_informative_ic("HP:0000707", \
            "HP:0002011"), 0.5108256)
            
        # check the most informative information content for two identical nodes
        self.assertAlmostEqual(self.hpo_graph.get_most_informative_ic("HP:0000924", \
            "HP:0000924"), 1.6094379)
