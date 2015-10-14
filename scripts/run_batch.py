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

def get_options():
    
    parser = argparse.ArgumentParser(description="helper script to run similarity" \
        "testing on a LSF cluster.")
    parser.add_argument("--script", required=True, help="Path to hpo similarity script")
    parser.add_argument("--phenotypes", required=True, help="path to phenotypes by proband JSON")
    parser.add_argument("--genes", required=True, help="path to probands by gene JSON")
    parser.add_argument("--temp-dir", required=True, help="path to hold intermediate files")
    parser.add_argument("--out", required=True, help="path to final output file")
    
    args = parser.parse_args()
    
    return args

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

def split_genes(genes_path, temp_dir):
    """ split probands by gene into many separate files with fewer genes
    """
    
    # get the full list of genes and their probands
    with open(genes_path, "r") as handle:
        genes = json.load(handle)
    
    iteration = 1
    for_json = {}
    for gene in sorted(genes):
        # don't include genes with insufficient probands for similarity analysis
        if len(genes[gene]) < 2:
            continue
        
        for_json[gene] = genes[gene]
        
        if len(for_json) > 1:
            # write an input file for the hpo similarity to run on
            path = os.path.join(temp_dir, "tmp.{}.txt".format(iteration))
            with open(path, "w") as output:
                json.dump(for_json, output, indent=4, sort_keys=True)
            
            for_json = {}
            iteration += 1
    
    return iteration - 1

def main():
    
    args = get_options()
    
    temp_dir = tempfile.mkdtemp(prefix=args.temp_dir)
    
    count = split_genes(args.genes, temp_dir)
    
    # set up run parameters
    job_name = "hpoSauce"
    job_id = "{0}[1-{1}]%50".format(job_name, count)
    
    infile = os.path.join(temp_dir, "tmp.\$LSB_JOBINDEX\.txt")
    outfile = os.path.join(temp_dir, "tmp.\$LSB_JOBINDEX\.output")
    
    command = ["python", args.script, \
        "--genes", infile, \
        "--phenotypes", args.phenotypes, \
        "--output", outfile,
        "--resnik"]
    submit_bsub_job(command, job_id, memory=2000, logfile="hpo_similarity.bjob")
    
    # merge the array output after the array finishes
    merge_id = "merge1_" + job_name
    command = ["head", "-n", "1", os.path.join(temp_dir, "tmp.1.output"), ">", args.out, \
        "; tail", "-q", "-n", "+2", os.path.join(temp_dir, "tmp.*.output"), "|", "sort", ">>", args.out]
    submit_bsub_job(command, merge_id, dependent_id=job_id, logfile="hpo_similarity.bjob")
    time.sleep(2)
    
    # submit a cleanup job to the cluster
    submit_bsub_job(["rm", "-r", temp_dir], job_id="cleanup", \
        dependent_id=merge_id, logfile="hpo_similarity.bjob")

if __name__ == '__main__':
    main()
