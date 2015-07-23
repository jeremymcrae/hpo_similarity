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

import os
import subprocess
import json
import time

SIMILARITY_DIR = "/nfs/users/nfs_j/jm33/apps/hpo_similarity"
RECESSIVE_DIR = "/nfs/users/nfs_j/jm33/apps/recessiveStats"
DATA_DIR = os.path.join(SIMILARITY_DIR, "data")
ALL_GENES_PATH = os.path.join(RECESSIVE_DIR, "data-raw", "recessive_probands_by_gene.json")
PHENOTYPES_PATH = os.path.join(DATA_DIR, "phenotypes_by_proband.json")
SIMILARITY_CODE = os.path.join(SIMILARITY_DIR, "hpo_similarity.py")

# get the full list of genes and their probands
with open(ALL_GENES_PATH, "r") as handle:
    genes = json.load(handle)


iteration = 1
for_json = {}
for gene in sorted(genes):
    # don't include genes with insufficient probands for similarity analysis
    if len(genes[gene]) < 2:
        continue
    
    for_json[gene] = genes[gene]
    
    if len(for_json) > 1:
        gene_path = os.path.join(DATA_DIR, "recessive_probands.{0}.json".format(iteration))
        output_path = os.path.join(SIMILARITY_DIR, "recessive.hpo_similarity.{0}.txt".format(iteration))
        
        # write an input file for the hpo similarity to run on
        with open(gene_path, "w") as output:
            json.dump(for_json, output, indent=4, sort_keys=True)
        
        # define the similarity command
        preamble = ["bsub", "-q", "long", "-o", "similarity_hpo.bjob_output.txt", \
            "-R 'select[mem>{0}] rusage[mem={0}]' -M {0}".format(2000)]
        command = ["python", SIMILARITY_CODE, \
            "--genes", gene_path, \
            "--phenotypes", PHENOTYPES_PATH, \
            "--output", output_path]
        command = ["bash", "-c", "\""] + command + ["\""]
        
        # submit the command via bsub
        command = " ".join(preamble + command)
        subprocess.call(command, shell=True)
        
        for_json = {}
        iteration += 1
