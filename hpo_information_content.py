""" determines similarity scores between proband HPO terms and the DDG2P 
genes predicted to be relevant to them.
"""

import math
import networkx as nx

class CalculateSimilarity(object):
    """ calculate graph similarity scores
    """
    
    def __init__(self, family_hpos, ddg2p_genes, hpo_graph, alt_node_ids):
        """
        
        Args:
            family_hpos: hpo term object, to get terms for probands, or parents 
            graph: graph of hpo terms, as networkx object
        """
        
        self.graph = hpo_graph
        self.ddg2p_genes = ddg2p_genes
        self.family_hpos = family_hpos
        self.alt_node_ids = alt_node_ids
        
        self.ic_cache = {}
        self.descendant_cache = {}
        self.ancestor_cache = {}
        
        self.tally_hpo_terms()
    
    def get_similarity_scores(self, probands):
        """ determines similarlity scores for the probands
        
        Args:
            probands: dict of HPO tags for probands, indexed by proband ID
        
        Returns:
            scores: dict of similarity scores, indexed by proband ID
        """
        
        scores = {}
        for proband in probands:
            score = self.get_ic_score_for_proband(proband, probands[proband])
            scores[proband] = score
        
        return scores
    
    def get_ic_score_for_proband(self, proband_ID, genes):
        """ calculate the similarity scores for the genes found for a proband
        
        Args:
            proband_ID: identifier string for an individual
            genes: list of (gene, inheritance) tuples for a proband
        """
        
        similarity_scores = {}
        for gene in genes:
            name = gene[0]
            inheritance = gene[1]
            
            gene_terms = self.get_gene_hpo_terms(name, inheritance)
            proband_terms = self.family_hpos[proband_ID].get_child_hpo()
            
            max_ics = []
            for proband_term in proband_terms:
                # for each proband HPO term, we look for the highest IC value
                # across all the HPO terms for the gene 
                max_ic = None
                for gene_term in gene_terms:
                    if gene_term == "":
                        continue
                    
                    # some of the HPO terms recorded for the probands, or in 
                    # the DDG2P database are alternate IDs for the HPO terms. 
                    # If so, swap over to the standard HPO term, as these are 
                    # the node names in the HPO graph.
                    if proband_term in self.alt_node_ids:
                        proband_term = self.alt_node_ids[proband_term]
                    if gene_term in self.alt_node_ids:
                        gene_term = self.alt_node_ids[gene_term]
                    
                    # find the max IC value for the term pair, and update the
                    # max_ic value if it is larger than the current value
                    max_ic = self.update_max_ic(proband_term, gene_term, max_ic)
                
                # only include the IC value if we actually found one (avoids
                # genes lacking HPO terms)
                if max_ic is not None:
                    max_ics.append(max_ic)
            
            similarity = "NA"
            if len(max_ics) > 0:
                similarity = "{0:.2f}".format(sum(max_ics)/len(max_ics))
            
            if name not in similarity_scores:
                similarity_scores[name] = {}
            
            similarity_scores[name][inheritance] = similarity
        
        return similarity_scores
    
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
        
        if inheritance not in self.ddg2p_genes[gene_name]:
            if inheritance == "X-linked dominant":
                gene_terms = self.ddg2p_genes[gene_name]["Monoallelic"]
            elif "Both" in self.ddg2p_genes[gene_name]:
                if inheritance == "Monoallelic" or inheritance == "Biallelic":
                    gene_terms = self.ddg2p_genes[gene_name]["Both"]
            else:
                print("missing another type")
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
        
        if bottom_term in self.ancestor_cache:
            subterms = self.ancestor_cache[bottom_term]
        else:
            subterms = nx.ancestors(self.graph, bottom_term)
            subterms.add(bottom_term)
            self.ancestor_cache[bottom_term] = subterms
        
        return subterms
    
    def find_common_ancestors(self, term_1, term_2):
        """ finds the common ancestors of two hpo terms
        
        Args:
            term_1: hpo term, eg HP:0000002
            term_2: hpo term, eg HP:0000003
        
        Returns:
            a list of all the common ancestors for the two terms
        """
        
        term1_ancestors = set(self.get_ancestors(term_1))
        term2_ancestors = set(self.get_ancestors(term_2))
        
        return term1_ancestors & term2_ancestors
    
    def update_max_ic(self, term_1, term_2, max_ic):
        """ finds the maximum information content value between two hpo terms
        
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
            if max_ic is None:
                max_ic = ic
            elif ic > max_ic:
                max_ic = ic
        
        return max_ic
    
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

