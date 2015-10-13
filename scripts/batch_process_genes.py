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
import argparse
import json
import time
import random
import tempfile

SIMILARITY_DIR = "/nfs/users/nfs_j/jm33/apps/hpo_similarity"
RECESSIVE_DIR = "/nfs/users/nfs_j/jm33/apps/recessiveStats"
DATA_DIR = os.path.join(SIMILARITY_DIR, "data")
# GENES_PATH = os.path.join(RECESSIVE_DIR, "data-raw", "recessive_probands_by_gene.silent.json")
GENES_PATH = os.path.join(DATA_DIR, "probands_by_gene.without_diagnosed.json")
PHENOTYPES_PATH = os.path.join(DATA_DIR, "phenotypes_by_proband.json")
SIMILARITY_CODE = os.path.join(SIMILARITY_DIR, "scripts", "proband_similarity.py")
LUSTRE_DIR = "/lustre/scratch113/projects/ddd/users/jm33/temp/tmp"

def get_bjobs():
    """ get a list of submitted jobs
    """
    
    command = ["bjobs", "-o", "\"JOBID", "USER", "STAT", "QUEUE", "JOB_NAME", "delimiter=';'\""]
    command = " ".join(command)
    output = subprocess.check_output(command, shell=True, stderr=open(os.devnull, "w"))
    
    bjobs = []
    for line in output.split("\n"):
        if line.startswith("JOBID") or line == "":
            continue
        
        line = line.strip().split(";")
        entry = {"jobid": line[0], "user":line[1], "stat":line[2], \
            "queue":line[3], "job_name":line[4]}
        bjobs.append(entry)
    
    return bjobs

def submit_bsub_job(command, job_id=None, dependent_id=None, memory=None, requeue_code=None, logfile=None, queue="normal", cpus=1):
    """ construct a bsub job submission command
    
    Args:
        command: list of strings that forma unix command
        job_id: string for job ID for submission
        dependent_id: job ID, or list of job IDs which the current command needs
            to have finished before the current command will start. Note that
            the list can be empty, in which case there are no dependencies.
        memory: minimum memory requirements (in megabytes)
    
    Returns:
        nothing
    """
    
    if job_id is None:
        job_id = get_random_string()
    
    job = "-J \"{0}\"".format(job_id)
    
    threads=""
    if cpus >1:
        threads="-n{0} -R 'span[hosts=1]'".format(cpus)
    
    mem = ""
    if memory is not None:
        mem = "-R 'select[mem>{0}] rusage[mem={0}]' -M {0}".format(memory)
    
    requeue = ""
    if requeue_code is not None:
        requeue = "-Q 'EXCLUDE({0})'".format(requeue_code)
    
    dependent = ""
    if dependent_id is not None:
        if type(dependent_id) == list:
            dependent_id = " && ".join(dependent_id)
        dependent = "-w '{0}'".format(dependent_id)
    
    log = "bjob_output.txt"
    if logfile is not None:
        log = logfile
    
    preamble = ["bsub", job, dependent, requeue, "-q", queue, "-o", log, mem, threads]
    command = ["bash", "-c", "\""] + command + ["\""]
    
    command = " ".join(preamble + command)
    subprocess.call(command, shell=True)

def is_number(string):
    """ check whether a string can be converted to a number
    
    Args:
        string: value as a string, could be a number
        
    Returns:
        True/False for whether the value can be converted to a number
    """
    
    try:
        number = float(string)
    except ValueError:
        return False
    
    return True

def get_random_string(prefix=None):
    """ make a random string, which we can use for bsub job IDs, so that
    different jobs do not have the same job IDs.
    """
    
    # set up a random string to associate with the run
    hash_string = "%8x" % random.getrandbits(32)
    hash_string = hash_string.strip()
    
    # done't allow the random strings to be equivalent to a number, since
    # the LSF cluster interprets those differently from letter-containing
    # strings
    while is_number(hash_string):
        hash_string = "%8x" % random.getrandbits(32)
        hash_string = hash_string.strip()
    
    if prefix is not None:
        hash_string = prefix + hash_string
    
    return hash_string

# get the full list of genes and their probands
with open(GENES_PATH, "r") as handle:
    genes = json.load(handle)

basename = os.path.basename(GENES_PATH)
basename, ext = os.path.splitext(basename)
# iteration = 1
for_json = {}
job_ids = []
result_paths = []
input_paths = []
for gene in sorted(genes):
    # don't include genes with insufficient probands for similarity analysis
    if len(genes[gene]) < 2:
        continue
    
    while len(get_bjobs()) > 200:
        time.sleep(30)
    
    for_json[gene] = genes[gene]
    
    if len(for_json) > 1:
        # gene_path = os.path.join(DATA_DIR, "{}.{}.json".format(basename, iteration))
        # output_path = os.path.join(SIMILARITY_DIR, "{}.{}.txt".format(basename, iteration))
        outfile = tempfile.NamedTemporaryFile(prefix=LUSTRE_DIR, suffix=".txt", delete=False)
        
        # write an input file for the hpo similarity to run on
        infile = tempfile.NamedTemporaryFile(prefix=LUSTRE_DIR, suffix=".json", delete=False)
        json.dump(for_json, infile, indent=4, sort_keys=True)
        infile.close()
        
        command = ["python", SIMILARITY_CODE, \
            "--genes", infile.name, \
            "--phenotypes", PHENOTYPES_PATH, \
            "--output", outfile.name,
            "--resnik"]
        job_id = get_random_string()
        submit_bsub_job(command, job_id, memory=2000, logfile="similarity_hpo.bjob_output.txt")
        time.sleep(0.5)
        
        for_json = {}
        # iteration += 1
        job_ids.append(job_id)
        result_paths.append(outfile.name)
        input_paths.append(infile.name)

while len([ x for x in get_bjobs() if x["job_name"] in job_ids ]) > 0:
    time.sleep(30)

# delete all the temporary input files
a = [ os.remove(x) for x in input_paths ]

output = open(os.path.join(SIMILARITY_DIR, "{}.txt".format(basename)), "w")
output.write("hgnc\thpo_similarity_p_value\n")
genes = []
for x in result_paths:
    with open(x) as handle:
        header = handle.readline()
        genes += handle.readlines()
    os.remove(x)

genes = sorted(genes)
output.writelines(genes)
output.close()
