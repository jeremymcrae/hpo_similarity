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

from hpo_similarity.ontology import open_ontology
from hpo_similarity.check_proband_terms import check_terms_in_graph

class TestCheckProbandTermsPy(unittest.TestCase):
    """ class to test HPO term sanity checking
    """
    
    def setUp(self):
        """ construct a ICSimilarity object and terms per proband for unit tests
        """
        
        path = os.path.join(os.path.dirname(__file__), "data", "obo.txt")
        self.hpo_graph, _, _ = open_ontology(path)
        
        self.hpo_terms = {
            "person_01": ["HP:0000924"],
            "person_02": ["HP:0000118", "HP:0002011"],
            "person_03": ["HP:0000707", "HP:0002011"]
        }
        
        self.hpo_graph.tally_hpo_terms(self.hpo_terms)
    
    def test_check_tally_term_error(self):
        ''' check we get an error if we tally a missing term
        '''
        
        # add a term not in the graph
        self.hpo_terms['person_04'] = ["HP:0012759"]
        
        with self.assertRaises(ValueError):
            self.hpo_graph.tally_hpo_terms(self.hpo_terms)
    
    def test_check_terms_in_graph(self):
        """ check that check_terms_in_graph works correctly
        """
        
        self.hpo_terms['person_04'] = ["HP:0012759"]
        
        # check that we raise an error when an indiviaul has a term that is not
        # present in the ontology grpah
        with self.assertRaises(ValueError):
            check_terms_in_graph(self.hpo_graph, self.hpo_terms)
        
        # check when the errant individual is removed, the function runs without errors
        del self.hpo_terms['person_04']
        self.assertEqual(check_terms_in_graph(self.hpo_graph, self.hpo_terms), None)
    
