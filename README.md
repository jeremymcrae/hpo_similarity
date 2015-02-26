### Test similarity of Human-Phenotype-Ontology (HPO) terms in groups of patients.

There are two python package requirements, both of which can be installed with pip:
- networkx
- matplotlib

This code depends on:
- a [python ontology parser](https://github.com/ntamas/gfam/blob/master/gfam/go/obo.py)
  written by Tamás Nepusz.
- the [obo file](http://purl.obolibrary.org/obo/hp.obo) from the
  [Human Phenotype Ontology Consortium](http://human-phenotype-ontology.org/).

The similarity of HPO terms in groups of individuals can be analysed with:
```sh
python3 hpo similarity.py \
  --variants VARIANT_PATH \
  --output OUTFILE
```
