""" some graph analyses of HPO terms and their usage in patients and gene hits
"""

# dependencies: 
#    go-parser (http://bazaar.launchpad.net/~ntamas/+junk/go-parser/files)
#    networkx (http://networkx.github.io/, pip install --user networkx)
#    matplotlib (pip install --user matplotlib)
#    HPO (http://compbio.charite.de/hudson/job/hpo/ or http://www.human-phenotype-ontology.org/)
#    DDG2P

import os
import networkx as nx
import matplotlib as plt

import obo
import load_files

# load the HPO database as graph database
    # convert terms between database versions, if necessary

# load clinical reporting file, find frequently identified genes, check their HPO terms from DDG2P, 
# load the patients HPO terms, check for matches between top genes and 

USER_PATH = "/nfs/users/nfs_j/jm33/"
HPO_PATH = os.path.join(USER_PATH, "apps", "hpo_filter", "hpo_data", "hp.obo")
CANDIDATE_VARIANTS_PATH = os.path.join(USER_PATH, "clinical_reporting.txt")

ddd_freeze = "/nfs/ddd0/Data/datafreeze/1139trios_20131030/"
DDG2P_PATH = os.path.join(ddd_freeze, "DDG2P_with_genomic_coordinates_20131107.tsv")
PHENOTYPES_PATH = os.path.join(ddd_freeze, "phenotypes.shared.20131129.txt")
ALTERNATE_IDS_PATH = os.path.join(ddd_freeze, "person_sanger_decipher.private.txt")

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
                self.g.node[node_id] = obo_tags[key]

    def track_alt_ids(self, node_id, obo_tags):
        """ track alternate HPO IDs, to map between alternate and canonical
        
        Args:
            node_id: string ID for a graph node
            obo_tags: tags for an obo entry
        
        Returns:
            nothing, updates the dictionaries within this function
        """
        
        # make sure we can convert between the alternate IDs and their HPO ID
        if "alt_id" in obo_tags:
            for alt_id in obo_tags["alt_id"]:
                alt_id  = str(alt_id)
                self.alt_mapper[alt_id] = node_id
                
                if node_id not in ref_mapper:
                    self.ref_mapper[node_id] = []
                
                self.ref_mapper[node_id].append(alt_id)

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


def main():
    # build a graph of DDG2P terms, so we can trace paths between terms
    hpo = loadHPONetwork(HPO_PATH)
    graph = hpo.get_graph()
    
    ddg2p_genes = load_files.load_ddg2p(DDG2P_PATH)
    proband_hpo_terms = load_files.load_participants_hpo_terms(PHENOTYPES_PATH, ALTERNATE_IDS_PATH)
    genes_index, probands_index = load_candidate_genes(CANDIDATE_VARIANTS_PATH)
    
    # nx.draw(g)
    # plt.savefig("test.png")

if __name__ == '__main__':
    main()


