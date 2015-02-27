""" permutes probands through genes, so we can see whether phenotype matching is
overinflated.
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
    all_probands = set.union(*all_probands)
    
    permuted = {}
    for gene in probands:
        # get a set of probands, except without those probands who are in the
        # current gene
        excluding_gene = all_probands - set(probands[gene])
        
        # randomly sample a new set of probands, as many as the current gene has
        sample = random.sample(excluding_gene, len(probands[gene]))
        permuted[gene] = sample
    
    return permuted
