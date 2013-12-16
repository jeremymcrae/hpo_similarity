""" simple functions to load and parse data for HPO terms analyses
"""

def load_ddg2p(ddg2p_path):
    """ load a DDG2P gene file, so we can extract HPO terms for each gene
    
    Args:
        ddg2p_path: path to DDG2P file
    
    Returns:
        dictionary of hpo terms indexed by gene name, then inheritance type
    """
    
    f = open(ddg2p_path)
    
    # allow for gene files with different column names and positions
    header = f.readline().strip().split("\t")
    if "type" in header:
        gene_label = "gene"
        confirmed_status_label = "type"
        inheritance_label = "mode"
        hpo_label = "hpo_ids"
    else:
        raise ValueError("The gene file lacks expected header column names")
    
    # get the positions of the columns in the list of header labels
    gene_column = header.index(gene_label)
    confirmed_status_column = header.index(confirmed_status_label)
    inheritance_column = header.index(inheritance_label)
    hpo_column = header.index(hpo_label)
    
    genes = {}
    for line in f:
        print line, hpo_column
        line = line.strip().split("\t")
        
        gene = line[gene_column]
        inheritance = line[inheritance_column]
        hpo_terms = line[hpo_column]
        
        if gene not in genes:
            genes[gene] = {}
        
        if inheritance not in genes[gene]:
            genes[gene][inheritance] = set()
        
        # split the hpo terms, then update the gene/inheritance set
        hpo_terms = hpo_terms.split(";")
        hpo_terms = [x.strip() for x in hpo_terms]
        genes[gene][inheritance].update(hpo_terms)
    
    return genes

class familyHPO(object):
    """small class to handle HPO terms for family members of a trio
    """
    
    def __init__(self, child_hpo, maternal_hpo, paternal_hpo):
        """initiate the class
        
        Args:
            child_hpo: HPO string for the proband, eg "HP:0005487|HP:0001363"
            maternal_hpo: string of HPO terms for the mother
            paternal_hpo: string of HPO terms for the father
        """
        
        self.child_hpo = self.format_hpo(child_hpo)
        self.maternal_hpo = self.format_hpo(maternal_hpo)
        self.paternal_hpo = self.format_hpo(paternal_hpo)
    
    def format_hpo(hpo_terms):
        """ formats a string of hpo terms to a list
        
        Args:
            hpo_terms: string of hpo terms joined with "|" 
        
        Returns:
            list of hpo terms, or None
        """
        
        # account for no hpo terms recorded for a person
        if hpo_terms.strip() == ".":
            return None
        
        # account for multiple hpo terms for an individual
        if "|" in hpo_terms:
            return hpo_terms.strip().split("|")
        
        return [hpo_terms.strip()]
    
    def get_child_hpo(self):
        return self.child_hpo
    
    def get_maternal_hpo(self):
        return self.maternal_hpo
    
    def get_paternal_hpo(self):
        return self.paternal_hpo

def load_participants_hpo_terms(pheno_path, alt_id_path):
    """ loads patient data, and obtains
    """
    
    # loads the decipher to DDD ID mapping file
    alt_ids = {}
    f = open(alt_id_path)
    for line in f:
        ref_id = line[0]
        alt_id = line[1]
        if ":" in alt_id:
            alt_id = alt_id.split(':"')[0]
        alt_ids[alt_id] = ref_id
    
    # load the phenotype data for each participant
    f = open(pheno_path)
    
    # allow for gene files with different column names and positions
    header = f.readline().strip().split("\t")
    if "type" in header:
        proband_label = "decipher_id"
        child_hpo_label = "child_hpo"
        maternal_hpo_label = "maternal_hpo"
        paternal_hpo_label = "paternal_hpo"
    else:
        raise ValueError("The gene file lacks expected header column names")
    
    # get the positions of the columns in the list of header labels
    proband_column = header.index(proband_label)
    child_hpo_column = header.index(child_hpo_label)
    maternal_hpo_column = header.index(maternal_hpo_label)
    paternal_hpo_column = header.index(paternal_hpo_label)
    
    participant_hpo = {}
    for line in f:
        proband_id = line[proband_column]
        child_hpo = line[child_hpo_column]
        maternal_hpo = line[maternal_hpo_column]
        paternal_hpo = line[paternal_hpo_column]
        
        # swap the proband across to the DDD ID if it exists
        if proband_id in alt_ids:
            proband_id = alt_ids[proband_id]
        
        participant_hpo[proband_id] = familyHPO(child_hpo, maternal_hpo, paternal_hpo)
    
    return participant_hpo

def load_candidate_genes(candidate_genes_path):
    """ loads a list of candidate genes for the participants
    """
    
    f = open(candidate_genes_path)
    
     # allow for gene files with different column names and positions
    header = f.readline().strip().split("\t")
    if "type" in header:
        proband_label = "proband"
        gene_label = "gene_hpo"
        inheritance_label = "inheritance_hpo"
    else:
        raise ValueError("The gene file lacks expected header column names")
    
    proband_column = header.index(proband_label)
    gene_column = header.index(gene_label)
    inheritance_column = header.index(inheritance_label)
    
    genes_index = {}
    probands_index = {}
    
    for line in f:
        # ignore blank lines
        if line == "\n":
            continue
        line = line.strip().split("\t")
        proband_ID = line[proband_column]
        gene = line[gene_column]
        inheritance = line[inheritance_column]
        
        if gene not in genes_index:
            genes_index[gene] = set()
            
        if proband not in probands_index:
            probands_index = set()
        
        genes_index[gene].update((proband_ID, inheritance))
        probands_index[proband_ID].update((gene, inheritance))
    
    return genes_index, probands_index
        
    
    






    
