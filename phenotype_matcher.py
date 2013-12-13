"""
"""

# dependencies: 
#    go-parser (http://bazaar.launchpad.net/~ntamas/+junk/go-parser/files)
#    networkx (http://networkx.github.io/, pip install --user networkx)
#    matplotlib (, pip install --user matplotlib)
#    HPO (http://compbio.charite.de/hudson/job/hpo/ or http://www.human-phenotype-ontology.org/)
#    DDG2P

import obo
import networkx as nx
import matplotlib as plt


# load the HPO database as graph database
    # convert terms between database versions, if necessary

# load the DDG2P gene list with HPO terms

# load clinical reporting file, find frequently identified genes, check their HPO terms from DDG2P, 
# load the patients HPO terms, check for matches between top genes and 


def load_hpo_database(hpo_path):
    """ load the human phenotype ontology (HPO) database in obo format
    
    Args:
        hpo_path: path to HPO obo formatted file
        
    Returns:
        parser.header: obo header for file
        hpo_entries: list of entries in HPO database
    """
    
    parser = obo.Parser(hpo_path)
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
    
    hpo_keys = ['comment', 'subset', 'xref', 'synonym', 'name', 'created_by', 'creation_date', 'def', 'alt_id']
    
    for key in hpo_keys:
        if key in obo_tags:
            graph.node[node_id] = obo_tags[key]

def track_alternate_ids(obo_tags, node_id, alt_mapper, ref_mapper):
    """ keep track of alternate IDs for HPO terms, so we can map between alternate and canonical
    
    Args:
        obo_tags: tags for an obo entry
        node_id: string ID for a graph node
        alt_mapper: dict to map from an alternate HPO ID to a canonical HPO IDs
        ref_mapper: dict to map from a canonical HPO ID to its potential alternate IDs
    
    Returns:
        nothing, updates the dictionaries within this function
    """
    
    # make sure we can convert between the alternate IDs and their HPO ID
    if "alt_id" in obo_tags:
        for alt_id in obo_tags["alt_id"]:
            alt_id  = str(alt_id)
            alt_mapper[alt_id] = node_id
            
            if node_id not in ref_mapper:
                ref_mapper[node_id] = []
            
            ref_mapper[node_id].append(alt_id)

def check_if_obsolete(obo_tags):
    """ checks if an "is_obsolete" flag is in the tags for an obo entry
    
    Args:
        obo_tags: tags for an obo entry
    
    Returns:
        True/False for whether the entry is obsolete
    """
    
    if "is_obsolete" in obo_tags:
        if str(obo_tags["is_obsolete"][0]) == "true":
            return True 
    
    return False

def create_hpo_graph(hpo_header, hpo_list):
    """builds a networkx graph from obo parsed data
    
    Args:
        hpo_header: header for obo file
        hpo_list: list of obo entries
    
    Returns:
        g: netwrokx graph object for hpo file
        alt_id_mapper: dict to map from alternate HPO IDs to canonical HPO IDs
        ref_id_mapper: dict to map from canonical HPO IDs to alternate HPO IDs
    """
    
    g = nx.DiGraph()
    
    # add the hpo header values as attributes for the graph
    for header_id in hpo_header:
        g.graph[header_id] = hpo_header[header_id]
    
    # track alternate HPO IDs (since we use HPO IDs as node IDs)
    alt_id_mapper = {}
    ref_id_mapper = {}
    
    for entry in hpo_list:
        # ignore obsolete HPO entries
        if check_if_obsolete(entry.tags):
            continue
        
        node_id = str(entry.tags["id"][0])
        g.add_node(node_id)
        
        # make sure we can convert between HPO ID and their alternate IDs
        track_alternate_ids(entry.tags, node_id, alt_id_mapper, ref_id_mapper)
        
        # include the attribute data for the node
        add_hpo_attributes_to_node(g, node_id, entry.tags)
        
        # add the predecessors to the node
        if "is_a" in entry.tags:
            for predecessor in entry.tags["is_a"]:
                predecessor = str(predecessor)
                
                # add the predecessor node to the graph (don't worry about attributes, we'll add 
                # those when we parse the entry for the predecessor)
                if node_id not in g:
                    g.add_node(predecessor)
                
                g.add_edge(predecessor, node_id)    
   
    return g, alt_id_mapper, ref_id_mapper


def main():
    hpo_filename = "/Volumes/jm33/code/phenotype_ontology_filter/hpo_data/hp.obo"
    hpo_header, hpo = load_hpo_database(hpo_filename)
    g = create_hpo_graph(hpo_header, hpo)

if __name__ == '__main__':
    main()


