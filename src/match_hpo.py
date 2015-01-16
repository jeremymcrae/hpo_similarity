""" find matching HPO terms for genes to HPO terms in probands
"""

import sys
import matplotlib as plt
plt.use("Agg")
import networkx as nx

class CheckHPOMatches(object):
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
        
        self.cached_terms = {}
        self.plotted_graphs = set()
         
        self.node_positions = {}
    
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
            
            probands = self.genes_index[gene]
            
            for proband in probands:
                has_obligate = False
                proband = proband[0]
                
                # pull out the terms used for the proband
                family_terms = self.family_hpos[proband]
                proband_terms = family_terms.get_child_hpo()
                
                for obligate_term in obligate_terms[gene]:
                    subterms = self.get_subterms(obligate_term)
                    
                    for proband_term in proband_terms:
                        if proband_term in subterms:
                            has_obligate = True
                            break
                    
                    if has_obligate:
                        break
                
                if has_obligate:
                    # subgraph = self.graph.subgraph(subterms)
                    # self.plot_subgraph(subgraph, obligate_term, proband_term)
                    
                    if proband not in hpo_matches:
                        hpo_matches[proband] = []
                    hpo_matches[proband].append(gene)
                
                # if set(subterms) != set(alt_subterms):
                #     self.plot_compare_sets(graph, subterms, alt_subterms)
        
        return hpo_matches
    
    def get_subterms(self, top_term):
        """ finds the set of subterms that descend from a top level HPO term
        """
        
        if top_term not in self.cached_terms:
            subterms = nx.descendants(self.graph, top_term)
            subterms.add(top_term)
            self.cached_terms[top_term] = subterms
        
        subterms = self.cached_terms[top_term]
        
        return subterms
    
    def plot_subgraph(self, g, top_term, found_term):
        """ plots a subgraph showing routes between different nodes
        """
       
        plotname = "hpo_graphs." + top_term + "-" + found_term + ".pdf"
        plotname = plotname.replace(":", "_")
        
        # don't replot the graph if it was created earlier for another proband
        if plotname in self.plotted_graphs:
            return
        
        # find which nodes are in the subgraph
        nodes = g.nodes()
        nodes_path = nx.shortest_path(g, top_term, found_term)
        
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
        
        # get some labels for the end points
        end_labels = {}
        end_labels[top_term] = top_term
        end_labels[found_term] = found_term
        
        # get the labels for the intermediate nodes
        intermediate_labels = {}
        for node in nodes_path:
            if node != top_term and node != found_term:
                intermediate_labels[node] = node
        
        # show the edges for the path from the top level term to the found term
        path_edges = []
        for pos in range(len(nodes_path))[:-1]:
            path_edges.append((nodes_path[pos], nodes_path[pos + 1]))
        
        pos = self.find_node_positions(g, top_term)
        nx.draw_networkx_nodes(g, pos=pos, nodelist=nodes, with_labels=False, \
                               node_color=cols, node_size=sizes, alpha=0.1)
        nx.draw_networkx_edges(g, pos, width=0.01, alpha=0.1)
        nx.draw_networkx_nodes(g, pos=pos, nodelist=nodes_path, with_labels=False)
        nx.draw_networkx_edges(g, pos, edgelist=path_edges, width=2, edge_color="blue")
        nx.draw_networkx_labels(g, pos=pos, labels=end_labels, font_size=10, \
                                font_color="black")
        nx.draw_networkx_labels(g, pos=pos, labels=intermediate_labels, \
                                font_size=6, font_color="black")
        
        plt.pyplot.savefig(plotname)
        plt.pyplot.close()
        
        # only plot each graph once per run
        self.plotted_graphs.add(plotname)
    
    def find_node_positions(self, graph, top_term):
        """ positions the nodes in the graph, and caches to repeat the node
        """
        
        # cache the graph node positioning
        if top_term not in self.node_positions:
            self.node_positions[top_term] = nx.spring_layout(graph)
        
        return self.node_positions[top_term]
        

