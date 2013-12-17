""" some graph analyses of HPO terms and their usage in patients and gene hits
"""

# dependencies: 
#    go-parser (http://bazaar.launchpad.net/~ntamas/+junk/go-parser/files)
#    networkx (http://networkx.github.io/, pip install --user networkx)
#    matplotlib (pip install --user matplotlib)
#    HPO (http://compbio.charite.de/hudson/job/hpo/ or http://www.human-phenotype-ontology.org/)
#    DDG2P

import networkx as nx
import os
import matplotlib as plt

import load_files
import create_hpo_graph as hpo

# load the HPO database as graph database
    # convert terms between database versions, if necessary

# load clinical reporting file, find frequently identified genes, check their HPO terms from DDG2P, 
# load the patients HPO terms, check for matches between top genes and 

USER_PATH = "/nfs/users/nfs_j/jm33/"
HPO_FOLDER = os.path.join(USER_PATH, "apps", "hpo_filter")
HPO_PATH = os.path.join(HPO_FOLDER, "hpo_data", "hp.obo")
OBLIGATE_GENES_PATH = os.path.join(HPO_FOLDER, "obligate_terms", "obligate_hpo_terms.frequent_genes.txt")
CANDIDATE_VARIANTS_PATH = os.path.join(USER_PATH, "clinical_reporting.txt")

ddd_freeze = "/nfs/ddd0/Data/datafreeze/1139trios_20131030/"
DDG2P_PATH = os.path.join(ddd_freeze, "DDG2P_with_genomic_coordinates_20131107.tsv")
PHENOTYPES_PATH = os.path.join(ddd_freeze, "phenotypes.shared.version.20131129.txt")
ALTERNATE_IDS_PATH = os.path.join(ddd_freeze, "person_sanger_decipher.private.txt")

def count_genes(genes_index):
    """ count the number of times each gene was found in the candidate genes
    
    Args:
        genes_index: dictionary of proband IDs and inheritance tuples indexed by gene
    
    Returns:
        list of genes, sorted by number of occurences
    """
    
    genes_count = []
    for gene in genes_index:
        genes_count.append((len(genes_index[gene]), gene))
    
    genes_count.sort()
    genes_count.reverse()
    
    return genes_count

def check_for_obligate_terms(family_hpos, genes_index, obligate_terms, graph):
    """ finds whether probands have an obligate HPO term given a gene that requires one
    
    Args:
        family_hpos: family list, with proband HPO terms and parents HPO terms
        genes_index: dict of genes, with proband and inheritance tuple lists as values
        obligate_terms: obligate HPO terms indexed by gene name
        graph: graph of hpo terms
    """
    
    for gene in obligate_terms:
        obligate_hpos = obligate_terms[gene]
        probands = genes_index[gene]
        has_obligate = False   
        # for proband in probands:
        for proband in ["DDDP111162"]:
            proband = proband[0]
            family_terms = family_hpos[proband]
            proband_terms = family_terms.get_child_hpo()
            # mother_terms = family_terms.get_maternal_hpo()
            # father_terms = family_terms.get_paternal_hpo()
            
            for obligate_term in obligate_terms[gene]:
                subterms = nx.dfs_successors(graph, obligate_term)
                # print subterms
                for proband_term in proband_terms:
                    if proband_term == obligate_term:
                        has_obligate = True
                        break
                    elif proband_term in subterms:
                        has_obligate = True
                        break
        # print subterms
        if has_obligate:
            print proband
            print proband_term, obligate_term, subterms
        
        subgraph = graph.subgraph(subterms.keys())
        nx.draw(subgraph)
        nx.show()
            

def main():
    # build a graph of DDG2P terms, so we can trace paths between terms
    hpo_file = hpo.loadHPONetwork(HPO_PATH)
    graph = hpo_file.get_graph()
    
    # load gene HPO terms, proband HPO terms, and HPO terms from candidate variants
    ddg2p_genes = load_files.load_ddg2p(DDG2P_PATH)
    family_hpo_terms = load_files.load_participants_hpo_terms(PHENOTYPES_PATH, ALTERNATE_IDS_PATH)
    genes_index, probands_index = load_files.load_candidate_genes(CANDIDATE_VARIANTS_PATH)
    obligate_hpo_terms = load_files.load_obligate_terms(OBLIGATE_GENES_PATH)
    
    # for gene in count_genes(genes_index):
    #     print str(gene[1]) + "\t" + str(gene[0])
    check_for_obligate_terms(family_hpo_terms, genes_index, obligate_hpo_terms, graph)
    
    # nx.draw(g)
    # plt.savefig("test.png")

if __name__ == '__main__':
    main()


