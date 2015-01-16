""" determines similarity scores between proband HPO terms and the DDG2P
genes predicted to be relevant to them.
"""

import math
import networkx as nx
    

class CalculateSimilarity(object):
    """ calculate graph similarity scores
    """
    
    def __init__(self, family_hpos, ddg2p_genes, hpo_graph, alt_node_ids, proband_phenotypes):
        """
        
        Args:
            family_hpos: hpo term object, to get terms for probands, or parents
            graph: graph of hpo terms, as networkx object
        """
        
        self.graph = hpo_graph
        self.ddg2p_genes = ddg2p_genes
        self.family_hpos = family_hpos
        self.alt_node_ids = alt_node_ids
        
        self.graph_depth = nx.eccentricity(self.graph, v="HP:0000001")
        
        self.ic_cache = {}
        self.descendant_cache = {}
        self.ancestor_cache = {}
        self.path_cache = {}
        
        self.test_tally_hpo_terms(proband_phenotypes)
    
    def get_similarity_scores(self, probands):
        """ determines similarity scores for the probands
        
        Args:
            probands: dict of HPO tags for probands, indexed by proband ID
            score_type: similarity score to return (eg max_IC, path_length)
        
        Returns:
            scores: dict of similarity scores, indexed by proband ID
        """
        
        scores = {}
        for proband in probands:
            scores[proband] = self.get_proband_score(proband, probands[proband])
        
        return scores
    
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
            
            max_values = []
            for proband_term in proband_terms:
                proband_term = self.fix_alternate_ID(proband_term)
                # for each proband HPO term, we look for the highest value
                # across all the HPO terms for the gene
                max_value = None
                for gene_term in gene_terms:
                    if gene_term == "":
                        continue
                    gene_term = self.fix_alternate_ID(gene_term)
                    max_value = self.update_value(proband_term, gene_term, max_value)
                
                # only include the IC value if we actually found one (avoids
                # genes lacking HPO terms)
                if max_value is not None:
                    max_values.append(max_value)
            
            similarity = "NA"
            if len(max_values) > 0:
                similarity = "{0:.2f}".format(sum(max_values)/len(max_values))
                # similarity = "{0:.2f}".format(max(max_values))
                # similarity = "{0}".format(self.get_h_index(max_values))
            
            if name not in similarity_scores:
                similarity_scores[name] = {}
            
            similarity_scores[name][inheritance] = similarity
        
        return similarity_scores
    
    def get_h_index(self, scores):
        """ get the h-index value from a list of scores
        
        Args:
            scores: list of float scores
        """
        
        h = 0
        
        for score in scores:
            if score > h:
                h += 1
        
        return h
    
    def geomean(self, numbers):
        """ get the geometric mean of a list of values
        """
        
        power = float(1)/len(numbers)
        
        total = 1
        for x in numbers:
            total *= x
        
        return total ** power
    
    def fix_alternate_ID(self, term):
        """ converts HPO terms using alternate IDs to the standard term
        
        some of the HPO terms recorded for the probands, or in the DDG2P
        database are alternate IDs for the HPO terms. If so, swap over to the
        standard HPO term, as these are the node names in the HPO graph.
        """
        
        if term in self.alt_node_ids:
            term = self.alt_node_ids[term]
        
        return term
    
    def get_gene_hpo_terms(self, gene_name, inheritance):
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
        
        if gene_name not in self.ddg2p_genes:
            gene_terms = []
        elif inheritance not in self.ddg2p_genes[gene_name]:
            if inheritance == "X-linked dominant":
                gene_terms = self.ddg2p_genes[gene_name]["Monoallelic"]
            elif "Both" in self.ddg2p_genes[gene_name]:
                if inheritance == "Monoallelic" or inheritance == "Biallelic":
                    gene_terms = self.ddg2p_genes[gene_name]["Both"]
            else:
                print("missing another type")
                gene_terms = []
        elif inheritance in self.ddg2p_genes[gene_name]:
            gene_terms = self.ddg2p_genes[gene_name][inheritance]
        
        return gene_terms
    
    def tally_hpo_terms(self):
        """ tallies each HPO term across the DDG2P genes
        
        Args:
            genes: hpo term dictionary
        """
        
        self.hpo_counts = {}
        self.total_freq = 0.0
        
        for gene in self.ddg2p_genes:
            for mode in self.ddg2p_genes[gene]:
                for term in self.ddg2p_genes[gene][mode]:
                    if term not in self.hpo_counts:
                        self.hpo_counts[term] = 0.0
                    
                    self.hpo_counts[term] += 1
                    self.total_freq += 1
    
    def test_tally_hpo_terms(self, patient_pheno):
        """ tallies each HPO term across the proband phenotypes
        
        Args:
            genes: hpo term dictionary
        """
        
        self.hpo_counts = {}
        self.total_freq = 0.0
        
        for line in patient_pheno:
            term = line[1]
            if term not in self.hpo_counts:
                self.hpo_counts[term] = 0.0
            
            self.hpo_counts[term] += 1
            self.total_freq += 1
    
    def get_subterms(self, top_term):
        """ finds the set of subterms that descend from a top level HPO term
        
        Args:
            top_term: hpo term to find descendants of
        
        Returns:
            set of descendant HPO terms
        """
        
        if top_term in self.descendant_cache:
            subterms = self.descendant_cache[top_term]
        else:
            subterms = nx.descendants(self.graph, top_term)
            self.descendant_cache[top_term] = subterms
        
        return subterms
    
    def get_ancestors(self, bottom_term):
        """ finds the set of subterms that are ancestors of a HPO term
        
        Args:
            bottom_term: hpo term to find ancestors of
        
        Returns:
            set of ancestor HPO terms
        """
        
        if bottom_term not in self.ancestor_cache:
            subterms = nx.ancestors(self.graph, bottom_term)
            subterms.add(bottom_term)
            self.ancestor_cache[bottom_term] = subterms
        
        subterms = self.ancestor_cache[bottom_term]
        
        return subterms
    
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
        
        term1_ancestors = set(self.get_ancestors(term_1))
        term2_ancestors = set(self.get_ancestors(term_2))
        
        return term1_ancestors & term2_ancestors
    

class ICSimilarity(CalculateSimilarity):
    """ calculate similarity by IC score
    """
    
    def find_most_informative_term(self):
        """ find the most infomative HPO term in the probands
        """
        
        informative_term = ""
        max_IC = 0
        for term in self.hpo_counts:
            ic = self.calculate_information_content(term)
            if ic > max_IC:
                informative_term = term
                max_IC = ic
        
        # print("max_IC value =", max_IC)
        
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
            if term not in self.graph:
                return 0
            
            descendants = self.get_subterms(term)
            
            term_count = 1.0
            if term in self.hpo_counts:
                term_count += self.hpo_counts[term]
            for subterm in descendants:
                if subterm in self.hpo_counts:
                    term_count += self.hpo_counts[subterm]
            
            # cache the IC, so we don't have to recalculate for the term
            self.ic_cache[term] = -math.log(term_count/self.total_freq)
        
        return self.ic_cache[term]
    
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


class PathLengthSimilarity(CalculateSimilarity):
    """ calculate similarity score by path length
    """
    
    def get_shortest_path(self, term_1, term_2):
        """ finds the shortest path between two terms (and caches the result)
        
        Args:
            term_1: HPO ID for graph node
            term_2: HPO ID for graph node
        
        Returns:
            shortest_path: list of nodes for path
        """
        
        if (term_1, term_2) not in self.path_cache:
            path = nx.shortest_path(self.graph, term_1, term_2)
            self.path_cache[(term_1, term_2)] = path
            self.path_cache[(term_2, term_1)] = path
       
        shortest_path = self.path_cache[(term_1, term_2)]
        
        return shortest_path
    
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
        
        common_terms = self.find_common_ancestors(term_1, term_2)
        
        if len(common_terms) == 0:
            return max_length
        
        ancestor = self.find_closest_ancestor(term_1, common_terms)
        
        # get the paths from the two terms to their closest ancestor
        path_1 = self.get_shortest_path(ancestor, term_1)
        path_2 = self.get_shortest_path(ancestor, term_2)
        
        # find the shortest path length between the nodes
        path_length = float(len(path_1 + path_2) - 1) # -1 as ancestor is included twice
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
        
        
    
    
