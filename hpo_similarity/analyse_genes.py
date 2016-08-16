"""
Copyright (c) 2015 Wellcome Trust Sanger Institute

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from hpo_similarity.check_proband_terms import check_terms_in_graph
from hpo_similarity.test_similarity import test_similarity

def analyse_genes(hpo_graph, hpo_by_proband, probands_by_gene, output_path, iterations, score_type):
    """ tests genes to see if their probands share HPO terms more than by chance.
    
    Args:
        hpo_graph: ICSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        hpo_by_proband: dictionary of HPO terms per proband
        probands_by_gene: dictionary of genes, to the probands who have variants
            in those genes.
        output_path: path to file to write the results to, or sys.stdout object.
        iterations: number of iterations to run.
    """
    
    check_terms_in_graph(hpo_graph, hpo_by_proband)
    
    # Sometimes output_path is actually sys.stdout, other times it is a path.
    try:
        output = open(output_path, "w")
    except TypeError:
        output = output_path
    
    output.write("hgnc\thpo_similarity_p_value\n")
    
    for gene in sorted(probands_by_gene):
        probands = probands_by_gene[gene]
        
        p_value = None
        if len(probands) > 1:
            p_value = test_similarity(hpo_graph, hpo_by_proband, probands, iterations, score_type)
        
        if p_value is None:
            continue
        
        output.write("{0}\t{1}\n".format(gene, p_value))
    
    output.close()
