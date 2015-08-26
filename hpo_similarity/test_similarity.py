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

from __future__ import division

import bisect
import random

from hpo_similarity.get_scores import get_proband_similarity

def test_similarity(hpo_graph, hpo_by_proband, probands, n_sims, score_type="resnik"):
    """ find if groups of probands per gene share HPO terms more than by chance.
    
    We simulate a distribution of similarity scores by randomly sampling groups
    of probands. I tried matching the number of sampled HPO terms to the numbers
    in the probands for the gene. For that, I gave each term the chance of being
    sampled as the rate at which it was observed in all the probands. However,
    these sampled terms gave abberant QQ plots, with excessive numbers of
    extremely signficant P values. I suspect this is due to underlying
    relationships between HPO terms.
    
    Args:
        hpo_graph: ICSimilarity object for the HPO term graph, with
            information on how many times each term has been used across all
            probands.
        hpo_by_proband: dictionary of HPO terms per proband
        probands: list of proband IDs.
        n_sims: number of simulations to run.
        score_type: type of similarity metric to use ["resnik", "lin", "simGIC"]
    
    Returns:
        The probability that the HPO terms used in the probands match as well as
        they do.
    """
    
    probands = [hpo_by_proband[x] for x in probands if x in hpo_by_proband]
    other_probands = [x for x in hpo_by_proband if x not in probands]
    
    # We can't test similarity from a single proband. We don't call this
    # function for genes with a single proband, however, sometimes only one of
    # the probands has HPO terms recorded. We cannot estimate the phenotypic
    # similarity between probands in this case, so return None instead.
    if len(probands) < 2:
        return None
    
    observed = get_proband_similarity(hpo_graph, probands, score_type)
    
    # get a distribution of scores for randomly sampled HPO terms
    distribution = []
    for x in range(n_sims):
        sampled = random.sample(other_probands, len(probands))
        simulated = [hpo_by_proband[n] for n in sampled]
        predicted = get_proband_similarity(hpo_graph, simulated, score_type)
        distribution.append(predicted)
    
    distribution = sorted(distribution)
    
    # figure out where in the distribution the observed value occurs
    pos = bisect.bisect_left(distribution, observed)
    sim_prob = (abs(pos - len(distribution)))/(1 + len(distribution))
    
    if sim_prob == 0:
        sim_prob = 1 / (1 + len(distribution))
    
    return sim_prob
