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

===================================================
Assess similarity of HPO terms in sets of probands.

If we have a set of probands who we know share some genetic feature in a gene,
we would like to know what is the probability of them sharing sharing their
Human Phenotype Ontology (HPO; standardised phenotypes) terms.

The HPO terms form a graph, so in order to estimate similarity, we use common
ancestors of sets of terms.
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import argparse

from hpo_similarity.load_files import load_participants_hpo_terms, load_genes
from hpo_similarity.ontology import open_ontology
from hpo_similarity.permute_probands import permute_probands
from hpo_similarity.analyse_genes import analyse_genes

def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description="Examines the likelihood of \
        obtaining similar HPO terms in probands with variants in the same gene.")
    parser.add_argument("--genes", dest="genes_path", required=True, \
        help="Path to JSON file listing probands per gene. See \
            data/example_genes.json for format.")
    parser.add_argument("--phenotypes", dest="phenotypes_path", required=True, \
        help="Path to JSON file listing phenotypes per proband. See \
            data/example_phenotypes.json for format.")
    parser.add_argument("--ontology", \
        help="Optional path to HPO ontology obo file, see http://human-phenotype-ontology.org. " \
              "By default, this uses the hpo file stored in hpo_similarity/data/hpo.obo")
    parser.add_argument("--output", default=sys.stdout, \
        help="path to output file, defaults to standard out.")
    parser.add_argument("--permute", action="store_true", default=False,
        help="whether to permute the probands across genes, in order to assess \
            method robustness.")
    parser.add_argument("--iterations", type=int, default=100000,
        help="whether to permute the probands across genes, in order to assess \
            method robustness.")
    
    # allow for using different similarity scoring metrics
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--resnik", default="False", action="store_true",
        help="whether to use Resnik's measure of similarity (the default).")
    group.add_argument("--simgic", "--simGIC", default="False", action="store_true",
        help="whether to use simGIC measure of similarity.")
    group.add_argument("--lin", default="False", action="store_true",
        help="whether to use Lin's measure of semantic similarity.")
    
    args = parser.parse_args()
    
    # figure out the score type
    if args.resnik:
        args.score_type = "resnik"
    elif args.simgic:
        args.score_type = "simGIC"
    elif args.lin:
        args.score_type = "lin"
    else:
        args.score_type = "resnik"
    
    return args

def main():
    
    options = get_options()
    
    # build a graph of HPO terms, so we can trace paths between terms
    graph, alt_ids, obsolete = open_ontology(options.ontology)
    
    # load HPO terms and probands for each gene
    print("loading HPO terms and probands by gene")
    hpo_by_proband = load_participants_hpo_terms(options.phenotypes_path, \
        alt_ids, obsolete)
    probands_by_gene = load_genes(options.genes_path)
    
    if options.permute:
        probands_by_gene = permute_probands(probands_by_gene)
    
    graph.tally_hpo_terms(hpo_by_proband)
    
    print("analysing similarity")
    try:
        analyse_genes(graph, hpo_by_proband, probands_by_gene, \
            options.output, options.iterations, options.score_type)
    except KeyboardInterrupt:
        sys.exit("HPO similarity exited.")

if __name__ == '__main__':
    main()
