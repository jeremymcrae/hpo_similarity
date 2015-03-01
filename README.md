#### Similarity of Human-Phenotype-Ontology (HPO) terms in patient groups.

This code estimates the probability of a group of individuals sharing their HPO
terms. This is useful when we have a set of probands all of whom have similar
genetic variation in a single gene (for example they all have de novo mutations
in the same gene). We have standardised phenotype terms, some of which may be
shared amongst the individuals. We want to know how likely the group is to share
those terms.

To estimate this probability, we need three things:
- a way to quantify similarity of HPO terms for a pair of probands.
- a way to quantify similarity across more than two probands.
- a null distribution of similarity scores for those probands

###### Proband pair HPO similarity
For a pair of probands, we have two lists of HPO terms. The code calculates a
matrix of information content scores, for each term in proband A versus each
term in proband B. The information content is from the most informative common
ancestor of the two terms.

The IC scores are collapsed into a single score as the geometric mean of the IC
scores. I tried using the maximum IC rather than the geometric mean, but the QQ
plots were distorted, and the P values were excessively correlated with P values
from gene enrichment testing.

###### Proband group similarity
For a group of probands, I get a matrix of proband pair HPO similarity scores,
one score for each pairing. The scores are summed are across all the pairs to
get an overall similarity estimate for the group of probands.

###### Null distribution of similarity scores
A null distribution of similarity scores is simulated by randomly sampling
groups of probands. The similarity scores are calculated as above. The P value
is calculated as the proportion of simulated scores greater than the observed
probands' score.

Initially I thought we'd need to match the number of sampled HPO terms to the
numbers in the probands for the gene. For that, I sampled terms at the rate at
which they were observed in all the probands. However, the QQ plots from this
had an excess of extremely significant P values. I suspect this is due to
underlying relationships between HPO terms. Sampling probands seems better and
simpler.


##### Requirements
There are two python package requirements, both of which are pip installable:
- networkx
- matplotlib

This code depends on:
- a [python ontology parser](https://github.com/ntamas/gfam/blob/master/gfam/go/obo.py)
  written by Tamás Nepusz.
- the [obo file](http://purl.obolibrary.org/obo/hp.obo) from the
  [Human Phenotype Ontology Consortium](http://human-phenotype-ontology.org/).

##### Running the code
The similarity of HPO terms in groups of individuals can be analysed with:
```sh
python3 hpo similarity.py \
  --variants VARIANT_PATH \
  --phenotypes PHENOTYPES_PATH \
  --output OUTFILE
```
