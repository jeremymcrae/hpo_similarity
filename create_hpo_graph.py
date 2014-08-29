""" loads a HPO OBO file as a graph network
"""

import networkx as nx

import obo


class loadHPONetwork(object):
    """ load a HPO obo file into a graph network for subsequent analysis
    """
    
    def __init__(self, hpo_path):
        """ load the hpoe file, and process it into a network
        """
        header, hpo = self.load_hpo_database(hpo_path)
        self.create_hpo_graph(header, hpo)
    
    def load_hpo_database(self, hpo_path):
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

    def add_hpo_attributes_to_node(self, node_id, obo_tags):
        """ add hpo attributes to a graph node
        
        Args:
            graph: networkx graph object
            node_id: string ID for a graph node
            obo_tags: tags for an obo entry
            
        Returns:
            nothing, updates the graph node within this function
        """
        
        hpo_keys = ['comment', 'subset', 'xref', 'synonym', 'name', \
                    'created_by', 'creation_date', 'def', 'alt_id']
        
        for key in hpo_keys:
            if key in obo_tags:
                if key not in self.g.node[node_id]:
                    self.g.node[node_id][key] = []
                
                value = str(obo_tags[key][0])
                self.g.node[node_id][key].append(value)

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
                self.alt_id_mapper[alt_id] = node_id
                
                if node_id not in self.ref_id_mapper:
                    self.ref_id_mapper[node_id] = []
                
                self.ref_id_mapper[node_id].append(alt_id)

    def check_if_obsolete(self, obo_tags):
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

    def create_hpo_graph(self, hpo_header, hpo_list):
        """builds a networkx graph from obo parsed data
        
        Args:
            hpo_header: header for obo file
            hpo_list: list of obo entries
        
        Returns:
            nothing
        """
        
        self.g = nx.DiGraph()
        
        # add the hpo header values as attributes for the graph
        for header_id in hpo_header:
            self.g.graph[header_id] = hpo_header[header_id]
        
        # track alternate HPO IDs (since we use HPO IDs as node IDs)
        self.alt_id_mapper = {}
        self.ref_id_mapper = {}
        
        for entry in hpo_list:
            tags = entry.tags
            # ignore obsolete HPO entries
            if self.check_if_obsolete(tags):
                continue
            
            node_id = str(tags["id"][0])
            self.g.add_node(node_id)
            
            # make sure we can convert between HPO ID and their alternate IDs
            self.track_alt_ids(tags, node_id)
            
            # include the attribute data for the node
            self.add_hpo_attributes_to_node(node_id, tags)
            
            # add the predecessors to the node
            if "is_a" in tags:
                for predecessor in tags["is_a"]:
                    predecessor = str(predecessor)
                    self.g.add_edge(predecessor, node_id)
    
    def get_graph(self):
        """ passes the graph object for external functions
        """
        
        try:
            return self.g
        except AttributeError:
            return None
    
    def get_alt_ids(self):
        """ passes the alt ID dictionary for external functions
        """
        
        try:
            return self.alt_id_mapper
        except AttributeError:
            return None

