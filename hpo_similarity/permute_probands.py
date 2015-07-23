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

import random

def permute_probands(probands):
    """ permute the probands by gene, so that every gene gets a random sample
    
    We occasionally want to permute the probands for each gene, so that each
    gene receives a new set of randomly sampled probands. This is so that we can
    demonstrate that the phenotype matching method gives P values that follow a
    null distribution. This shows that the inflation that occurs in the
    unmodified P value distribution is meaningful.
    
    Args:
        probands: dictionary of proband lists for each gene, for example
            {"ADNP": [DDD01, DDD02], "ANKRD1": ["DDD03", "DDD04"]}.
    
    Returns:
        dictionary where each gene now has randomly sampled probands.
    """
    
    all_probands = probands.values()
    all_probands = set.union(*[set(x) for x in all_probands])
    
    permuted = {}
    for gene in probands:
        # get a set of probands, except without those probands who are in the
        # current gene
        excluding_gene = all_probands - set(probands[gene])
        
        # randomly sample a new set of probands, as many as the current gene has
        sample = random.sample(excluding_gene, len(probands[gene]))
        permuted[gene] = sample
    
    return permuted
