""" find matching HPO terms for genes to HPO terms in probands
"""

import matplotlib as plt
plt.use("Agg")
import networkx as nx

class checkHPOMatches(object):
    """ class to find matching HPO terms for genes to HPO terms in probands
    """
    
    def __init__(self, family_hpos, graph, genes_index):
        """ initiate the class with some standard values
        
        Args:
            family_hpos: family list, with proband HPO terms and parents HPO terms
            genes_index: dict of genes, listing proband and inheritance tuples
            graph: graph of hpo terms
        """
        
        self.family_hpos = family_hpos
        self.graph = graph
        self.genes_index = genes_index
        
        self.cached_subterms = {}
    
    def find_matches(self, obligate_terms):
        """ finds whether probands have an obligate HPO term for a gene
        
        Args:
            obligate_terms: obligate HPO terms indexed by gene name
        
        Returns:
            probands with hpo terms that match terms for certain genes as a dict 
            of gene lists indexed by proband ID
        """
        
        hpo_matches = {}
        for gene in obligate_terms:
            # don't bother to check genes that don't occur in the probands
            if gene not in self.genes_index:
                continue
            
            obligate_hpos = obligate_terms[gene]
            probands = self.genes_index[gene]
            
            for proband in probands:
                has_obligate = False
                proband = proband[0]
                
                # pull out the terms used for the proband
                family_terms = self.family_hpos[proband]
                proband_terms = family_terms.get_child_hpo()
                
                for obligate_term in obligate_terms[gene]:
                    subterms = self.get_subterms(self.graph, obligate_term)
                    
                    for proband_term in proband_terms:
                        if proband_term in subterms:
                            has_obligate = True
                            break
                    
                    if has_obligate:
                        break
                
                if has_obligate:
                    subgraph = graph.subgraph(subterms)
                    self.plot_subgraph(subgraph, obligate_term, proband_term)
                    
                    if proband not in hpo_matches:
                        hpo_matches[proband] = []
                    hpo_matches[proband].append(gene)
                
                # if set(subterms) != set(alt_subterms):
                #     self.plot_compare_sets(graph, subterms, alt_subterms)
        
        return hpo_matches
    
    def get_subterms(self, graph, top_term):
        """ finds the set of subterms that descend from a top level HPO term
        """
        
        if top_term in self.cached_subterms:
            subterms = self.cached_subterms[top_term]
        else:
            subterms = nx.dfs_successors(graph, top_term)
            successor_list = subterms.values()
            subterms = set([item for sublist in successor_list for item in sublist])
            subterms.add(top_term)
            self.cached_subterms[top_term] = subterms
        
        return subterms
    
    def plot_compare_sets(self, graph, subterms, alt_subterms):
        """ plots a subgraph to compare nodes from different lists
        """
        
        # find which nodes are in both sets
        subterms = set(subterms)
        alt_subterms = set(alt_subterms)
        intersection = subterms & alt_subterms
        union = subterms | alt_subterms
        
        # find the terms in each set that are not in both
        diff_subterms = subterms - intersection
        diff_alt_subterms = alt_subterms - intersection
       
        subgraph = graph.subgraph(union)
        subnodes = subgraph.nodes()
        
        cols = []
        sizes = []
        for node in subnodes:
            if node in intersection:
                color = "black"
                size = 5
            elif node in diff_subterms:
                color = "blue"
                size = 50
            elif node in diff_alt_subterms:
                color = "red"
                size = 50
            else:
                sys.exit("node not in a defined set!")
            colors.append(color)
            sizes.append(size)
        
        nx.draw(subgraph, with_labels=False, width=0.5, node_color=cols, node_size=sizes, alpha=0.5)
        plt.pyplot.savefig("test.pdf")
    
    def plot_subgraph(self, graph, top_term, found_term):
        """ plots a subgraph to compare nodes from different lists
        """
        
        # find which nodes are in the subgraph
        nodes = graph.nodes()
        
        cols = []
        sizes = []
        for node in nodes:
            if node == top_term:
                color = "blue"
                size = 200
            elif node == found_term:
                color = "green"
                size = 200
            else:
                color = "red"
                size = 50
            cols.append(color)
            sizes.append(size)
        
        labels = {}
        labels[top_term] = top_term
        labels[found_term] = found_term
        
        # show the edges for the path from the top level term to the found term
        edges = graph.edges()
        path_between_nodes = nx.shortest_path(graph, top_term, found_term)
        edges_to_highlight = []
        for pos in path_between_nodes[:-1]:
            edges.highlight.append(path_between_nodes[pos], path_between_nodes[pos + 1])
        
        pos = nx.spectral_layout(graph)
        nx.draw_networkx_nodes(graph, pos=pos, nodelist=nodes, with_labels=False, \
                               width=0.01, node_color=cols, node_size=sizes, \
                               alpha=0.2)
        nx.draw_networkx_edges(graph, pos, edge_list=edges, width=0.01)
        nx.draw_networkx_egdes(graph, pos, edge_list=edges_to_highlight, width=0.5, egde_color="red")
        nx.draw_networkx_labels(graph, pos=pos, labels=labels, font_size=10, \
                                font_color="red")
        plt.pyplot.savefig("test.pdf")
        sys.exit("exited after plotting graph")

