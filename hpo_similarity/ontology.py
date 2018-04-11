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

from pkg_resources import resource_filename

from hpo_similarity.similarity import ICSimilarity
from hpo_similarity.obo import Parser


class Ontology(object):
    """ load a HPO obo file into a graph network for subsequent analysis
    """
    
    def __init__(self, hpo_path=None):
        """ load the hpo file, and process it into a network
        """
        
        self.hpo_header, self.hpo_list = self.load_hpo_database(hpo_path)
        
        # track alternate HPO IDs (since we use HPO IDs as node IDs)
        self.alt_ids = {}
        self.ref_ids = {}
        self.obsolete_ids = set()
    
    def load_hpo_database(self, hpo_path):
        """ load the human phenotype ontology (HPO) database in obo format
        
        Args:
            hpo_path: path to HPO obo formatted file
            
        Returns:
            parser.header: obo header for file
            hpo_entries: list of entries in HPO database
        """
        
        if hpo_path is None:
            hpo_path = resource_filename(__name__, "data/hp.obo")
        
        parser = Parser(hpo_path)
        hpo_entries = []
        for entry in parser:
            hpo_entries.append(entry)
        
        return parser.headers, hpo_entries

    def add_hpo_attributes_to_node(self, graph, node_id, obo_tags):
        """ add hpo attributes to a graph node
        
        Args:
            graph: networkx graph object
            node_id: string ID for a graph node
            obo_tags: tags for an obo entry
            
        Returns:
            nothing, updates the graph node within this function
        """
        
        for key in obo_tags:
            graph.node[node_id][key] = str(obo_tags[key][0])

    def track_alt_ids(self, obo_tags, node_id):
        """ track alternate HPO IDs, to map between alternate and canonical
        
        Args:
            obo_tags: tags for an obo entry
            node_id: string ID for a graph node
        
        Returns:
            nothing, updates the dictionaries within this function
        """
        
        # make sure we can convert between the alternate IDs and their HPO ID
        if "alt_id" in obo_tags:
            for alt_id in obo_tags["alt_id"]:
                alt_id = str(alt_id)
                self.alt_ids[alt_id] = node_id
                
                if node_id not in self.ref_ids:
                    self.ref_ids[node_id] = []
                
                self.ref_ids[node_id].append(alt_id)

    def is_obsolete(self, obo_tags):
        """ checks if an "is_obsolete" flag is in the tags for an obo entry
        
        Args:
            obo_tags: tags for an obo entry
        
        Returns:
            True/False for whether the entry is obsolete
        """
        
        return "is_obsolete" in obo_tags and str(obo_tags["is_obsolete"][0]) == "true"

    def get_graph(self):
        """ builds a networkx graph from obo parsed data
        
        Returns:
            networkx graph object
        """
        
        graph = ICSimilarity()
        
        # add the hpo header values as attributes for the graph
        for header_id in self.hpo_header:
            graph.graph[header_id] = self.hpo_header[header_id]
        
        for entry in self.hpo_list:
            tags = entry.tags
            # ignore obsolete HPO entries
            if self.is_obsolete(tags):
                self.obsolete_ids.add(str(tags["id"][0]))
                continue
            
            node_id = str(tags["id"][0])
            graph.add_node(node_id)
            
            # make sure we can convert between HPO ID and their alternate IDs
            self.track_alt_ids(tags, node_id)
            
            # include the attribute data for the node
            self.add_hpo_attributes_to_node(graph, node_id, tags)
            
            # add the predecessors to the node
            if "is_a" in tags:
                for predecessor in tags["is_a"]:
                    predecessor = str(predecessor)
                    graph.add_edge(predecessor, node_id)
        
        return graph
    
    def get_alt_ids(self):
        """ passes the alt ID dictionary for external functions
        """
        
        return self.alt_ids
    
    def get_obsolete_ids(self):
        """ returns the set of obsolete tags
        """
        
        return self.obsolete_ids
