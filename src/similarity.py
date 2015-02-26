""" determines similarity scores between proband HPO terms and the DDG2P
genes predicted to be relevant to them.
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import math
import networkx

class CalculateSimilarity(object):
    """ calculate graph similarity scores
    """
    
    def __init__(self, family_hpos, hpo_graph, alt_node_ids):
        """
        
        Args:
            family_hpos: hpo term object, to get terms for probands, or parents
            graph: graph of hpo terms, as networkx object
        """
        
        self.graph = hpo_graph
        self.family_hpos = family_hpos
        self.alt_node_ids = alt_node_ids
        
        self.graph_depth = networkx.eccentricity(self.graph, v="HP:0000001")
        
        self.descendant_cache = {}
        self.ancestor_cache = {}
        
        self.hpo_counts = {}
        self.total_freq = 0
    
    def get_similarity_scores(self, ddg2p, probands):
        """ determines similarity scores for the probands
        
        Args:
            probands: dict of HPO tags for probands, indexed by proband ID
            score_type: similarity score to return (eg max_IC, path_length)
        
        Returns:
            scores: dict of similarity scores, indexed by proband ID
        """
        
        scores = {}
        for proband in probands:
            scores[proband] = self.get_proband_score(ddg2p, proband, probands[proband])
        
        return scores
    
    def get_proband_score(self, ddg2p, proband_ID, genes):
        """ calculate the similarity scores for the genes found for a proband
        
        Args:
            proband_ID: identifier string for an individual
            genes: list of (gene, inheritance) tuples for a proband
            score_type: the type of similarity measure to return
        """
        
        similarity_scores = {}
        for gene in genes:
            name = gene[0]
            inheritance = gene[1]
            
            if "," in inheritance:
                inheritances = inheritance.split(",")
                gene_terms = []
                for inh in inheritances:
                    gene_terms += self.get_gene_hpo_terms(ddg2p, name, inh)
            else:
                gene_terms = self.get_gene_hpo_terms(ddg2p, name, inheritance)
            
            proband_terms = self.family_hpos[proband_ID].get_child_hpo()
            
            max_values = []
            for proband_term in proband_terms:
                proband_term = self.fix_alternate_id(proband_term)
                # for each proband HPO term, we look for the highest value
                # across all the HPO terms for the gene
                max_value = None
                for gene_term in gene_terms:
                    if gene_term == "":
                        continue
                    gene_term = self.fix_alternate_id(gene_term)
                    max_value = self.update_value(proband_term, gene_term, max_value)
                
                # only include the IC value if we actually found one (avoids
                # genes lacking HPO terms)
                if max_value is not None:
                    max_values.append(max_value)
            
            similarity = "NA"
            if len(max_values) > 0:
                similarity = "{0:.2f}".format(sum(max_values)/len(max_values))
            
            if name not in similarity_scores:
                similarity_scores[name] = {}
            
            similarity_scores[name][inheritance] = similarity
        
        return similarity_scores
    
    def fix_alternate_id(self, term):
        """ converts HPO terms using alternate IDs to the standard term
        
        some of the HPO terms recorded for the probands, or in the DDG2P
        database are alternate IDs for the HPO terms. If so, swap over to the
        standard HPO term, as these are the node names in the HPO graph.
        """
        
        if self.graph.has_node(term):
            pass
        elif term in self.alt_node_ids:
            term = self.alt_node_ids[term]
        else:
            raise KeyError
        
        return term
    
    def get_gene_hpo_terms(self, ddg2p, gene_name, inheritance):
        """ pulls out the hpo terms for a DDG2P gene
        
        Sometimes the DDG2P database does not have the inheritance mode
        listed for a gene for which we have found variants for a proband. This
        is because occasionally the DDG2P genes have modes such as "Both",
        which means either "Monoallelic or "Biallelic", and we have reported
        the variants under one of those specific modes. We need to convert the
        inheritance back in these cases, as well as for "Monoallelic on the
        X-chrom, which is under "X-linked dominant".
        
        Args:
            gene_name: string gene name
            inheritance: the inheritance mode listed for the gene
        """
        
        if gene_name not in ddg2p:
            gene_terms = []
        elif inheritance not in ddg2p[gene_name]:
            if inheritance == "X-linked dominant":
                gene_terms = ddg2p[gene_name]["Monoallelic"]
            elif "Both" in ddg2p[gene_name]:
                if inheritance == "Monoallelic" or inheritance == "Biallelic":
                    gene_terms = ddg2p[gene_name]["Both"]
            else:
                print("missing another type")
                gene_terms = []
        elif inheritance in ddg2p[gene_name]:
            gene_terms = ddg2p[gene_name][inheritance]
        
        return gene_terms
    
    def tally_hpo_terms(self, hpo_terms, source="ddg2p"):
        """ tallies each HPO term across the DDG2P genes
        
        Args:
            hpo_terms: hpo term dictionary
        """
        
        assert source in ["ddg2p", "child_hpo"]
        
        for item in hpo_terms:
            if source == "ddg2p":
                for mode in hpo_terms[item]:
                    for term in hpo_terms[item][mode]:
                        self.add_hpo(term)
            elif source == "child_hpo":
                child_terms = hpo_terms[item].get_child_hpo()
                for term in child_terms:
                    self.add_hpo(term)
    
    def add_hpo(self, term):
        """ increments the count for an HPO term
        
        This increments a) the count for the specific term, and b) the total
        count of all terms.
        
        Args:
            term: string for HPO term
        """
        
        if term not in self.hpo_counts:
            # don't use terms which cannot be placed on the graph
            if not self.graph.has_node(term):
                return
            
            if term not in self.hpo_counts:
                self.hpo_counts[term] = 0
            
            self.hpo_counts[term] += 1
        
        self.total_freq += 1
    
    def get_subterms(self, top_term):
        """ finds the set of subterms that descend from a top level HPO term
        
        Args:
            top_term: hpo term to find descendants of
        
        Returns:
            set of descendant HPO terms
        """
        
        if top_term not in self.descendant_cache:
            self.descendant_cache[top_term] = networkx.descendants(self.graph, top_term)
        
        return self.descendant_cache[top_term]
    
    def get_ancestors(self, bottom_term):
        """ finds the set of subterms that are ancestors of a HPO term
        
        Args:
            bottom_term: hpo term to find ancestors of
        
        Returns:
            set of ancestor HPO terms
        """
        
        if bottom_term not in self.ancestor_cache:
            subterms = networkx.ancestors(self.graph, bottom_term)
            subterms.add(bottom_term)
            self.ancestor_cache[bottom_term] = subterms
        
        return self.ancestor_cache[bottom_term]
    
    def find_common_ancestors(self, term_1, term_2):
        """ finds the common ancestors of two hpo terms
        
        Args:
            term_1: hpo term, eg HP:0000002
            term_2: hpo term, eg HP:0000003
        
        Returns:
            a list of all the common ancestors for the two terms
        """
        
        # ignore terms that are obsolete (ie are not in the graph)
        if term_1 not in self.graph or term_2 not in self.graph:
            return set()
        
        return set(self.get_ancestors(term_1)) & set(self.get_ancestors(term_2))


class ICSimilarity(CalculateSimilarity):
    """ calculate similarity by IC score
    """
    
    counts_cache = {}
    ic_cache = {}
    
    def find_most_informative_term(self, terms=None):
        """ find the most infomative HPO term in the probands
        """
        
        informative_term = ""
        max_IC = 0
        
        if terms is None:
            terms = self.hpo_counts
        
        for term in terms:
            ic = self.calculate_information_content(term)
            if ic > max_IC:
                informative_term = term
                max_IC = ic
        
        return informative_term
    
    def calculate_information_content(self, term):
        """ calculates the information content for an hpo term
        
        For discussion of information content and similarity scores, see:
        Van der Aalst et al., (2007) Data & Knowledge Engineering 61:137-152
        
        Args:
            term: hpo term, eg HP:0000001
        
        Returns:
            the information content value for a single hpo term
        """
        
        if term not in self.ic_cache:
            term_count = self.get_term_count(term)
            
            if term not in self.graph:
                return 0
            
            # cache the IC, so we don't have to recalculate for the term
            self.ic_cache[term] = -math.log(term_count/self.total_freq)
        
        return self.ic_cache[term]
    
    def get_term_count(self, term):
        """ Count how many times a term (or its subterms) was used.
        
        Args:
            term: hpo term, eg HP:0000001
        
        Returns:
            the number of times a term (or its subterms) was used.
        """
        
        if term not in self.counts_cache:
            if term not in self.graph:
                return 0
            
            descendants = self.get_subterms(term)
            
            count = 0
            if term in self.hpo_counts:
                count += self.hpo_counts[term]
            for subterm in descendants:
                if subterm in self.hpo_counts:
                    count += self.hpo_counts[subterm]
            
            self.counts_cache[term] = count
        
        return self.counts_cache[term]
    
    def update_value(self, term_1, term_2, max_ic):
        """ finds the maximum information content value between two hpo terms
        
        Calculates simlarity by Resnick's most informative ancestor,
        
        Resnik P. Using information content to evaluate semantic similarity in
        a taxonomy, in: Proceedings of the 14th International Joint Conference
        on Artificial Intelligence, 1995.
        
        Kohler S et al. (2009). Clinical diagnostics in human genetics with
        semantic similarity searches in ontologies. American Journal of Human
        Genetics, 85:457-464.
        
        Args:
            term_1: hpo term, eg HP:0000001
            term_2: hpo term, eg HP:0000001
            max_ic: current maximum IC, or None
        
        Returns:
            the maximum IC value between the two terms
        """
        
        common_terms = self.find_common_ancestors(term_1, term_2)
        for term in common_terms:
            ic = self.calculate_information_content(term)
            if max_ic is None or ic > max_ic:
                max_ic = ic
        
        return max_ic


class ICDistanceSimilarity(ICSimilarity):
    """ calculate similarity by IC distance measure
    """
    
    def update_value(self, term_1, term_2, max_distance):
        """ finds the minimum distance value between two hpo terms
        
        Calculates similarity by Jian and Conrath's distance measure
        
        Jiang J, Conrath D. Semantic similarity based on corpus statistics and
        lexical taxonomy, in: Proceedings of the 10th International Conference
        on Research on Computational Linguistics, 1997.
        
        Args:
            term_1: hpo term, eg HP:0000001
            term_2: hpo term, eg HP:0000001
            max_ic: current maximum IC, or None
        
        Returns:
            the maximum IC value between the two terms
        """
        
        common_terms = self.find_common_ancestors(term_1, term_2)
        
        # if there are no common terms
        if len(common_terms) == 0:
            return max_distance
        
        term_1_ic = self.calculate_information_content(term_1)
        term_2_ic = self.calculate_information_content(term_2)
        for term in common_terms:
            ic = self.calculate_information_content(term)
            distance = term_1_ic + term_2_ic - 2 * ic
            if max_distance is None or -distance > max_distance:
                max_distance = -distance
        
        return max_distance


class PathLengthSimilarity(ICSimilarity):
    """ calculate similarity score by path length
    """
    
    path_cache = {}
    
    def get_shortest_path(self, term_1, term_2):
        """ finds the shortest path between two terms (and caches the result)
        
        Args:
            term_1: HPO ID for graph node
            term_2: HPO ID for graph node
        
        Returns:
            list of nodes for path
        """
        
        term_1 = self.fix_alternate_id(term_1)
        term_2 = self.fix_alternate_id(term_2)
        
        if (term_1, term_2) not in self.path_cache:
            try:
                path = networkx.shortest_path(self.graph, term_1, term_2)
            except networkx.exception.NetworkXNoPath:
                path = self.get_path_between_nondescendants(term_1, term_2)
                
            self.path_cache[(term_1, term_2)] = path
            self.path_cache[(term_2, term_1)] = path
        
        return self.path_cache[(term_1, term_2)]
    
    def get_path_between_nondescendants(self, term_1, term_2):
        """ gets the shortest path between terms not from the same branch
        
        Args:
            term_1: HPO ID for graph node
            term_2: HPO ID for graph node
        
        Returns:
            shortest_path: list of nodes for path
        """
        
        common_terms = self.find_common_ancestors(term_1, term_2)
        ancestor = self.find_closest_ancestor(term_1, common_terms)
        
        # get the paths from the two terms to their closest ancestor
        path_1 = self.get_shortest_path(ancestor, term_1)[::-1]
        path_2 = self.get_shortest_path(ancestor, term_2)
        
        return path_1[:-1] + path_2
    
    def find_closest_ancestor(self, node, ancestors):
        """ finds the closest ancestor of a term from a list of ancestor terms
        
        Args:
            node: node ID to search from
            ancestors: list of ancestor node IDs
        
        Returns:
            ancestor_to_use: closest ancestor node ID
        """
        
         # find the path to the closest common ancestor
        shortest = None
        ancestor_to_use = ""
        for ancestor in ancestors:
            length = len(self.get_shortest_path(ancestor, node))
            if shortest is None or length < shortest:
                shortest = length
                ancestor_to_use = ancestor
        
        return ancestor_to_use
    
    def update_value(self, term_1, term_2, max_length):
        """ finds the maximum length between two hpo terms
        
        This implements the Leacock and Chodorow pathfinding measure, as
        described in Pedersen et al., (2007) Journal of Biomedical Informatics,
        40:288-299. Note that as a directed acyclic graph, we can't find the
        path between terms on two different branches, we can only natively
        find paths from ancestors to descendants. We work around this by
        finding the closest ancestor for the terms, and then getting the paths
        to that.
        
        Leacock C, Chodorow M. Combining local context and WordNet similarity
        for word sense identification. In: Fellbaum C, editor. WordNet: An
        electronic lexical database. Cambridge, MA: MIT Press; 1998. p. 265-83
        
        Args:
            term_1: hpo term, eg HP:0000001
            term_2: hpo term, eg HP:0000001
            max_length: current maximum path length
        
        Returns:
            the maximum path length between the two terms
        """
        
        # find the shortest path length between the nodes
        path = self.get_path_between_nondescendants(term_1, term_2)
        
        path_length = float(len(path))
        value = -math.log(path_length/(2 * self.graph_depth))
        
        if max_length is None or value > max_length:
            max_length = value
        
        return max_length

class JaccardSimilarity(CalculateSimilarity):
    """ calculate similarity by Jaccard score
    """
    
    def get_proband_score(self, proband_ID, genes):
        """ calculate the similarity scores for the genes found for a proband
        
        Args:
            proband_ID: identifier string for an individual
            genes: list of (gene, inheritance) tuples for a proband
            score_type: the type of similarity measure to return
        """
        
        similarity_scores = {}
        for gene in genes:
            name = gene[0]
            inheritance = gene[1]
            
            if "," in inheritance:
                inheritances = inheritance.split(",")
                gene_terms = []
                for inh in inheritances:
                    gene_terms += self.get_gene_hpo_terms(name, inh)
            else:
                gene_terms = self.get_gene_hpo_terms(name, inheritance)
            
            proband_terms = self.family_hpos[proband_ID].get_child_hpo()
            
            intersection = set(gene_terms) & set(proband_terms)
            union = set(gene_terms) | set(proband_terms)
            jaccard_similarity = float(len(intersection))/float(len(union))
            
            if name not in similarity_scores:
                similarity_scores[name] = {}
            
            similarity_scores[name][inheritance] = "{0:.2f}".format(jaccard_similarity)
        
        return similarity_scores
        
        
    
    
