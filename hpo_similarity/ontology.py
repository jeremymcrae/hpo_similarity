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

def load_hpo_database(hpo_path):
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

def add_hpo_attributes_to_node(graph, node_id, obo_tags):
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

def is_obsolete(obo_tags):
    """ checks if an "is_obsolete" flag is in the tags for an obo entry
    
    Args:
        obo_tags: tags for an obo entry
    
    Returns:
        True/False for whether the entry is obsolete
    """
    
    return "is_obsolete" in obo_tags and str(obo_tags["is_obsolete"][0]) == "true"

def track_alt_ids(alt_ids, obo_tags, node_id):
    """ track alternate HPO IDs, to map between alternate and canonical
    
    Args:
        obo_tags: tags for an obo entry
        node_id: string ID for a graph node
    
    Returns:
        nothing, just updates dictionary passed in
    """
    
    # make sure we can convert between the alternate IDs and their HPO ID
    if "alt_id" in obo_tags:
        for alt_id in obo_tags["alt_id"]:
            alt_id = str(alt_id)
            alt_ids[alt_id] = node_id

def add_entry(graph, entry, alt_ids, obsolete_ids):
    """ add a node to the graph
    
    Args:
        graph: a networkx DiGraph
        entry: HPO ontology object, to be added to the graph
        alt_ids: dictionary of alt IDs mapping to correct IDs
        obsolete_ids: set of odsolete IDs
    """
    tags = entry.tags
    # ignore obsolete HPO entries
    if is_obsolete(tags):
        obsolete_ids.add(str(tags["id"][0]))
        return
    
    node_id = str(tags["id"][0])
    graph.add_node(node_id)
    
    # make sure we can convert between HPO ID and their alternate IDs
    track_alt_ids(alt_ids, tags, node_id)
    
    # include the attribute data for the node
    add_hpo_attributes_to_node(graph, node_id, tags)
    
    # add the predecessors to the node
    if "is_a" in tags:
        for predecessor in tags["is_a"]:
            predecessor = str(predecessor)
            graph.add_edge(predecessor, node_id)

def open_ontology(path=None):
    """ builds a networkx graph from obo parsed data
    
    Returns:
        networkx graph object, alt IDs and obsolete IDs
    """
    
    header, entries = load_hpo_database(path)
    graph = ICSimilarity()
    
    # track alternate HPO IDs (since we use HPO IDs as node IDs)
    alt_ids = {}
    obsolete_ids = set()
    
    # add the hpo header values as attributes for the graph
    for header_id in header:
        graph.graph[header_id] = header[header_id]
    
    for entry in entries:
        add_entry(graph, entry, alt_ids, obsolete_ids)
    
    return graph, alt_ids, obsolete_ids
