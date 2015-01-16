""" class to append hpo search results back into the clinical reporting file
"""


class Reporting(object):
    """ small class to insert hpo findings into the clinical reporting file
    """
    
    def __init__(self, report_path, output_path):
        """ checks the file header
        
        Args:
            report_path: path to clinical reporting file
            output_path: path to write results to
        """
        
        self.report_path = report_path
        self.output_path = output_path
        
        self.check_header()
    
    def check_header(self):
        """ finds the important columns from the reporting file header
        """
        
        with open(self.report_path, "r") as f:
            header = f.readline().strip().split("\t")
        
        self.proband_column = header.index("proband")
        self.gene_column = header.index("gene")
        self.inheritance_column = header.index("inheritance")
    
    def start_modified_file(self, labels):
        """ starts the list of modified lines, by adding column labels
        
        Args:
            labels: column names to be added to the file header.
        
        Returns:
            Open file handle for output file.
        """
        
        handle = open(self.report_path)
        header = handle.readline().strip().split("\t")
        
        for label in labels:
            header.append(label)
        header = "\t".join(header) + "\n"
        
        self.parsed_lines = [header]
        
        return handle
    
    def get_insert_values(self, data, datatype, line):
        """ finds values to add to the clinical reporting line.
        
        Args:
            data: dictionary of values for a given proband, gene and inheritance
            datatype: one of ("matches", "scores")
            line: list of values for a clinical filtering variant line
        
        Returns:
            finds the value for a specified proband, gene and inheritance state.
        """
        
        proband = line[self.proband_column]
        gene = line[self.gene_column]
        inh = line[self.inheritance_column]
        
        if datatype == "matches":
            searched_genes = data[1]
            matches = data[0]
            
            values = ["not_checked", "NA"]
            if gene in searched_genes:
                values[0] = "gene_checked"
            
            if values[0] != "not_checked":
                values[1] = "FAIL"
            
            if proband in matches:
                if gene in matches[proband]:
                    values[1] = "PASS"
        elif datatype == "scores":
            values = [data[proband][gene][inh]]
        
        return values
    
    def modify_file(self, datatype, labels, data):
        """ modifies the clinical reporting file, depending on the result types
        
        Args:
            datatype: "matches" or "scores"
            labels: list of new column header names
            data: dictionary of values, indexed by proband_id, gene and
                inheritance.
        """
        
        f = self.start_modified_file(labels)
        for line in f:
            if line == "\n":
                self.parsed_lines.append(line)
                continue
            
            line = line.strip().split("\t")
            proband = line[self.proband_column]
            gene = line[self.gene_column]
            inh = line[self.inheritance_column]
            
            insert = self.get_insert_values(data, datatype, line)
            
            line += insert
            line = "\t".join(line) + "\n"
            self.parsed_lines.append(line)
        f.close()
        
        output = open(self.output_path, "w")
        output.writelines(self.parsed_lines)
        output.close()
        
        # we initially use the reporting file as input, and make a new file
        # from that. After the first set of data has been inserted, read from
        # the modified file.
        self.report_path = self.output_path
    
    def add_matches_to_report(self, matches, searched_genes, label):
        """ annotates candidate variants using proband matches indexed by gene
        
        Args:
            matches: dict of probands that were matched for each gene
            searched_genes: list of searched genes, so we can add whether a
                gene was not found because it was not searched for
            column_label: label to use in the file header
        """
        
        label = [label + s for s in ["_searched", "_passed"]]
        data = [matches, searched_genes]
        self.modify_file("matches", label, data)
        
        self.report_path = self.output_path
    
    def add_scores_to_report(self, scores, labels):
        """ annotates candidate variants using scores for each variant
        
        Args:
            scores: dict of scores, indexed by proband, gene, inheritance mode
            labels: column header labels for the scores
        """
        
        self.modify_file("scores", labels, scores)
