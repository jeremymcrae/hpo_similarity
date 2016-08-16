"""
Copyright (c) 2016 Wellcome Trust Sanger Institute

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

def check_terms_in_graph(graph, hpo_by_proband):
    ''' check that all of the proband terms occur in the HPO ontology graph
    
    Args:
        graph: ICSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        hpo_by_proband: dictionary of HPO terms per proband
    
    Raises:
        ValueError if any term is not present in the ontology graph
    '''
    
    for proband in hpo_by_proband:
        for term in hpo_by_proband[proband]:
            if term not in graph:
                raise ValueError('{0} has a term ({1}) missing from the ontology. '
                    '{1} might be in a more recent ontology, see '
                    'http://human-phenotype-ontology.org. You can define '
                    'ontology files with --ontology.'.format(proband, term))
