""" class to test the Ontology class
"""

import os
import unittest

import networkx

from hpo_similarity.ontology import Ontology
from hpo_similarity.obo import Stanza, Value

# define the header that will be parsed from the test obo dataset
HEADER = {'default-namespace': ['human_phenotype'], 'ontology': ['hp'],\
    'format-version': ['1.2'], 'data-version': ['2013-11-30']}

# define the list of parsed objects from the test obo dataset
HPO_LIST = [
    Stanza('Term', \
        {'comment': [Value('Root of all terms in the Human Phenotype Ontology.', None)],\
        'id': [Value('HP:0000001', None)], \
        'name': [Value('All', None)]}),
    Stanza('Term', \
        {'comment': [Value('This is the root of the phenotypic abnormality subontology of the HPO.', None)], \
        'synonym': [Value('Organ abnormality', ('EXACT []',))], \
        'name': [Value('Phenotypic abnormality', None)], \
        'is_a': [Value('HP:0000001', None)], \
        'id': [Value('HP:0000118', None)], \
        'def': [Value('A phenotypic abnormality.', ('[HPO:probinson]',))]}), \
    Stanza('Term', \
        {'is_obsolete': [Value('true', None)], \
        'id': [Value('HP:0000489', None)], \
        'name': [Value('Abnormality of globe location or size', None)]}), \
    Stanza('Term', \
        {'comment': [Value('The nervous system comprises the neuraxis (brain, spinal cord, and ventricles), the autonomic nervous system, the enteric nervous system, and the peripheral nervous system.', None)], \
        'synonym': [Value('Neurologic abnormalities', ('EXACT []',)), \
                    Value('Neurological abnormality', ('EXACT []',))], \
        'name': [Value('Abnormality of the nervous system', None)], \
        'id': [Value('HP:0000707', None)], \
        'is_a': [Value('HP:0000118', None)], \
        'alt_id': [Value('HP:0001333', None), Value('HP:0006987', None)], \
        'def': [Value('An abnormality of the `nervous system` (FMA:7157).', ('[HPO:probinson]',))]}), \
    Stanza('Term', \
        {'is_a': [Value('HP:0000118', None)], \
        'synonym': [Value('Skeletal abnormalities', ('EXACT []',)), \
                    Value('Skeletal anomalies', ('EXACT []',))], \
        'id': [Value('HP:0000924', None)], \
        'def': [Value('An abnormality of the `skeletal system` (FMA:23881).', ('[HPO:probinson]',))], \
        'name': [Value('Abnormality of the skeletal system', None)]}), \
    Stanza('Term', \
        {'synonym': [Value('Central nervous system disease', ('RELATED []',))], \
        'name': [Value('Abnormality of the central nervous system', None)], \
        'id': [Value('HP:0002011', None)], \
        'is_a': [Value('HP:0000707', None)], \
        'alt_id': [Value('HP:0002405', None), Value('HP:0002413', None), \
                   Value('HP:0002481', None)], \
        'def': [Value('An abnormality of the `central nervous system` (FMA:55675).', ('[HPO:curators]',))]})]


class TestOntologyPy(unittest.TestCase):
    """ test the Ontology class
    """
    
    def setUp(self):
        """ construct an Ontology object for unit tests
        """
        
        path = os.path.join(os.path.dirname(__file__), "data", "obo.txt")
        self.ontology = Ontology(path)
    
    def test_setup(self):
        """ test that the Ontology class has loaded the test obo file correctly
        """
        
        self.assertEqual(self.ontology.hpo_header, HEADER)
        self.assertEqual(self.ontology.hpo_list, HPO_LIST)
    
    def test_load_hpo_database(self):
        """ check that load_hpo_database works (also in the class initialisation)
        """
        
        path = os.path.join(os.path.dirname(__file__), "data", "obo.txt")
        temp_header, temp_hpo = self.ontology.load_hpo_database(path)
        
        self.assertEqual(temp_header, HEADER)
        self.assertEqual(temp_hpo, HPO_LIST)
    
    def test_add_hpo_attributes_to_node(self):
        """ test that add_hpo_attributes_to_node works correctly
        """
        
        # construct a node without any alternate IDs
        node_id = "HP:0000118"
        tags = {'comment': [Value('This is the root of the phenotypic abnormality subontology of the HPO.', None)], \
            'synonym': [Value('Organ abnormality', ('EXACT []',))], \
            'name': [Value('Phenotypic abnormality', None)], \
            'is_a': [Value('HP:0000001', None)], \
            'id': [Value('HP:0000118', None)], \
            'def': [Value('A phenotypic abnormality.', ('[HPO:probinson]',))]}
        
        graph = networkx.DiGraph()
        graph.add_node(node_id)
        
        self.ontology.add_hpo_attributes_to_node(graph, node_id, tags)
        
        self.assertEqual(set(graph.node[node_id].keys()), \
            set(["comment", "synonym", "name", "def", "is_a", "id"]))
        
        self.assertEqual(graph.node[node_id]["comment"], \
            "This is the root of the phenotypic abnormality subontology of the HPO.")
        self.assertEqual(graph.node[node_id]["synonym"], "Organ abnormality")
        self.assertEqual(graph.node[node_id]["name"], "Phenotypic abnormality")
        self.assertEqual(graph.node[node_id]["def"], "A phenotypic abnormality.")
        self.assertEqual(graph.node[node_id]["is_a"], "HP:0000001")
        self.assertEqual(graph.node[node_id]["id"], "HP:0000118")
    
    def test_track_alt_ids(self):
        """ check that track_alt_ids works correctly
        """
        
        # construct a node without any alternate IDs
        node_id = "HP:0000118"
        tags = {'comment': [Value('This is the root of the phenotypic abnormality subontology of the HPO.', None)], \
            'synonym': [Value('Organ abnormality', ('EXACT []',))], \
            'name': [Value('Phenotypic abnormality', None)], \
            'is_a': [Value('HP:0000001', None)], \
            'id': [Value('HP:0000118', None)], \
            'def': [Value('A phenotypic abnormality.', ('[HPO:probinson]',))]}
        
        # if we track the node, it shouldn't add anything to the empty ref_ids,
        # or alt_ids
        self.ontology.track_alt_ids(tags, node_id)
        self.assertEqual(self.ontology.ref_ids, {})
        self.assertEqual(self.ontology.alt_ids, {})
        
        # construct a node with some alternate IDs
        node_id = "HP:0000707"
        tags = {'comment': [Value('The nervous system comprises the neuraxis (brain, spinal cord, and ventricles), the autonomic nervous system, the enteric nervous system, and the peripheral nervous system.', None)], \
            'synonym': [Value('Neurologic abnormalities', ('EXACT []',)), \
                        Value('Neurological abnormality', ('EXACT []',))], \
            'name': [Value('Abnormality of the nervous system', None)], \
            'id': [Value('HP:0000707', None)], \
            'is_a': [Value('HP:0000118', None)], \
            'alt_id': [Value('HP:0001333', None), Value('HP:0006987', None)], \
            'def': [Value('An abnormality of the `nervous system` (FMA:7157).', ('[HPO:probinson]',))]}
        
        # now when we track the node, it should add in the alternate IDs
        self.ontology.track_alt_ids(tags, node_id)
        self.assertEqual(self.ontology.ref_ids, {node_id: ['HP:0001333', 'HP:0006987']})
        self.assertEqual(self.ontology.alt_ids["HP:0001333"], node_id)
        self.assertEqual(self.ontology.alt_ids["HP:0006987"], node_id)
    
    def test_is_obsolete(self):
        """ check that is_obsolete works correctly
        """
        
        # check that some tags with a True is_obsolete entry is True
        tags = {'is_obsolete': [Value('true', None)], \
            'id': [Value('HP:0000489', None)], \
            'name': [Value('Abnormality of globe location or size', None)]}
        self.assertTrue(self.ontology.is_obsolete(tags))
        
        # check that some tags with a False is_obsolete entry is False
        tags = {'is_obsolete': [Value('false', None)], \
            'id': [Value('HP:0000489', None)], \
            'name': [Value('Abnormality of globe location or size', None)]}
        self.assertFalse(self.ontology.is_obsolete(tags))
        
        # check that some tags without an is_obsolete entry is False
        tags = {'id': [Value('HP:0000707', None)], \
            'name': [Value('Abnormality of the nervous system', None)]}
        self.assertFalse(self.ontology.is_obsolete(tags))
    
    def test_get_graph(self):
        """ test that building a graph with get_graph works correctly
        """
        
        graph = self.ontology.get_graph()
        
        # check that all the nodes (aside from the obsolete node) have been
        # included
        self.assertEqual(set(graph.nodes()), \
            set(['HP:0000001', 'HP:0000118', 'HP:0000707', 'HP:0000924', 'HP:0002011',]))
        
        # check that the graph info is set correctly from the HPO header
        self.assertEqual(graph.graph, HEADER)
        
        # check that the edges match what is expected foir each node
        self.assertEqual(graph.edges('HP:0000001'), [("HP:0000001", "HP:0000118")])
        self.assertEqual(set(graph.edges('HP:0000118')), \
            set([("HP:0000118", "HP:0000707"), ("HP:0000118", "HP:0000924")]))
        self.assertEqual(graph.edges('HP:0000924'), [])
        
        # check that, as part of setting the graph up, we have constructed the
        # correct set of obsolete IDs
        self.assertEqual(self.ontology.obsolete_ids, set(["HP:0000489"]))
