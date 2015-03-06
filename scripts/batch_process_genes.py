""" submit genes for HPO similarity analysis in batches, to parallelise analysis
"""

import os
import subprocess
import json
import time

SIMILARITY_DIR = "/nfs/users/nfs_j/jm33/apps/hpo_similarity"
DATA_DIR = os.path.join(SIMILARITY_DIR, "data")
ALL_GENES_PATH = os.path.join(DATA_DIR, "recessive_probands_by_gene.json")
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
    
    if len(for_json) > 10:
        gene_path = os.path.join(DATA_DIR, "probands.{0}.json".format(iteration))
        output_path = os.path.join(SIMILARITY_DIR, "probands_results.{0}.txt".format(iteration))
        
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
