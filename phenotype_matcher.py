""" some graph analyses of HPO terms and their usage in patients and gene hits
"""

# dependencies: 
#    go-parser (http://bazaar.launchpad.net/~ntamas/+junk/go-parser/files)
#    networkx (http://networkx.github.io/, pip install --user networkx)
#    matplotlib (pip install --user matplotlib)
#    HPO (http://compbio.charite.de/hudson/job/hpo/ or http://www.human-phenotype-ontology.org/)
#    DDG2P

import matplotlib as plt
plt.use("Agg")
import networkx as nx
import os
import sys

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
ORGAN_TO_HPO_MAPPER_PATH = os.path.join(HPO_FOLDER, "obligate_terms", "ddg2p_organ_code_to_hpo_term.txt")
DDG2P_ORGAN_PATH = os.path.join(HPO_FOLDER, "obligate_terms", "ddg2p_organ_specificity.txt")
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

def plot_subgraph(graph, subterms, alt_subterms):
    """ plots a subgraph to compares node from different lists
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

def find_descendants(graph, start_node):
    """ recursively find all the descendants of a node
    """
    
    if len(graph.successors(start_node)) == 0:
        return [start_node]
    else:
        successors = graph.successors(start_node)
        downstream = []
        for node in successors:
            downstream += find_descendants(graph, node)
        return [start_node] + successors + downstream

def check_for_hpo_matches(family_hpos, genes_index, obligate_terms, graph):
    """ finds whether probands have an obligate HPO term for a gene
    
    Args:
        family_hpos: family list, with proband HPO terms and parents HPO terms
        genes_index: dict of genes, listing proband and inheritance tuples
        obligate_terms: obligate HPO terms indexed by gene name
        graph: graph of hpo terms
    
    Returns:
        probands with hpo terms that match terms for certain genes as a dict 
        of gene lists indexed by proband ID
    """
    
    hpo_matches = {}
    cached_subterms = {}
    for gene in obligate_terms:
        # don't bother to check genes that don't occur in the prbands
        if gene not in genes_index:
            continue
        
        obligate_hpos = obligate_terms[gene]
        probands = genes_index[gene]
        
        for proband in probands:
            has_obligate = False
            proband = proband[0]
            
            family_terms = family_hpos[proband]
            proband_terms = family_terms.get_child_hpo()
            # mother_terms = family_terms.get_maternal_hpo()
            # father_terms = family_terms.get_paternal_hpo()
            
            for obligate_term in obligate_terms[gene]:
                # subterms = nx.dfs_successors(graph, obligate_term)
                if obligate_term in cached_subterms:
                    subterms = cached_subterms[obligate_term]
                else:
                    subterms = find_descendants(graph, obligate_term)
                    cached_subterms[obligate_term] = subterms
                
                for proband_term in proband_terms:
                    if proband_term == obligate_term:
                        has_obligate = True
                        break
                    elif proband_term in subterms:
                        has_obligate = True
                        break
            
            if has_obligate:
                if proband not in hpo_matches:
                    hpo_matches[proband] = []
                hpo_matches[proband].append(gene)
            
            # if set(subterms) != set(alt_subterms):
            #     plot_subgraph(graph, subterms, alt_subterms)
    
    return hpo_matches

def annotate_clinical_report(path, matches, searched_genes, column_label):
    """
    """
    
    f = open(path)
    header = f.readline().strip().split("\t")
    proband_label = "proband"
    gene_label = "gene"
    proband_column = header.index(proband_label)
    gene_column = header.index(gene_label)
    
    header.append(column_label + "_searched")
    header.append(column_label + "_passed")
    header = "\t".join(header) + "\n"
    
    parsed_lines = [header]
    for line in f:
        if line == "\n":
            parsed_lines.append(line)
            continue
        
        line = line.strip().split("\t")
        proband = line[proband_column]
        gene = line[gene_column]
        
        searched = "not_checked"
        if gene in searched_genes:
            searched = "gene_checked"
        
        if searched == "not_checked":
            matched = "NA"
        else:
            matched = "FAIL"
        
        if proband in matches:
            if gene in matches[proband]:
                matched = "PASS"
        
        line.append(searched)
        line.append(matched)
        
        line = "\t".join(line) + "\n"
        parsed_lines.append(line)
    
    f.close()
    
    output = open(path, "w")
    output.writelines(parsed_lines)
    output.close()

def main():
    # build a graph of DDG2P terms, so we can trace paths between terms
    hpo_file = hpo.loadHPONetwork(HPO_PATH)
    graph = hpo_file.get_graph()
    
    # load gene HPO terms, proband HPO terms, and HPO terms from candidate variants
    ddg2p_genes = load_files.load_ddg2p(DDG2P_PATH)
    family_hpo_terms = load_files.load_participants_hpo_terms(PHENOTYPES_PATH, ALTERNATE_IDS_PATH)
    genes_index, probands_index = load_files.load_candidate_genes(CANDIDATE_VARIANTS_PATH)
    obligate_hpo_terms = load_files.load_obligate_terms(OBLIGATE_GENES_PATH)
    ddg2p_hpo_organ_terms = load_files.load_organ_terms(ORGAN_TO_HPO_MAPPER_PATH, DDG2P_ORGAN_PATH)
    
    # for gene in count_genes(genes_index):
    #     print str(gene[1]) + "\t" + str(gene[0])
    obligate_matches = check_for_hpo_matches(family_hpo_terms, genes_index, obligate_hpo_terms, graph,)
    ddg2p_organ_matches = check_for_hpo_matches(family_hpo_terms, genes_index, ddg2p_hpo_organ_terms, graph)
    
    annotate_clinical_report(CANDIDATE_VARIANTS_PATH, obligate_matches, obligate_hpo_terms, "obligate_terms")
    annotate_clinical_report(CANDIDATE_VARIANTS_PATH, ddg2p_organ_matches, ddg2p_hpo_organ_terms, "ddg2p_organ_terms")
    
    # nx.draw(g)
    # plt.savefig("test.png")

if __name__ == '__main__':
    main()


