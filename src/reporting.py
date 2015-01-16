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
        
        f = open(self.report_path, "r")
        header = f.readline().strip().split("\t")
        f.close()
        
        proband_label = "proband"
        gene_label = "gene"
        inheritance_label = "inheritance"
        
        self.proband_column = header.index(proband_label)
        self.gene_column = header.index(gene_label)
        self.inheritance_column = header.index(inheritance_label)
    
    def start_modified_file(self):
        """ starts the list of modified lines, by adding column labels
        """
        
        f = open(self.report_path)
        header = f.readline().strip().split("\t")
        
        for label in self.column_labels:
            header.append(label)
        header = "\t".join(header) + "\n"
        
        self.parsed_lines = [header]
        
        return f
    
    def get_insert_values(self):
        """ finds values to add to the clinical reporting line
        """
        
        if self.type == "matches":
            searched_genes = self.data[1]
            matches = self.data[0]
            
            values = ["not_checked", "NA"]
            if self.gene in searched_genes:
                values[0] = "gene_checked"
            
            if values[0] != "not_checked":
                values[1] = "FAIL"
            
            if self.proband in matches:
                if self.gene in matches[self.proband]:
                    values[1] = "PASS"
        elif self.type == "scores":
            values = [self.data[self.proband][self.gene][self.inheritance]]
        
        return values
    
    def modify_file(self):
        """ modifies the clinical reporting file, depending on the result types
        """
        
        f = self.start_modified_file()
        for line in f:
            if line == "\n":
                self.parsed_lines.append(line)
                continue
            
            line = line.strip().split("\t")
            self.proband = line[self.proband_column]
            self.gene = line[self.gene_column]
            self.inheritance = line[self.inheritance_column]
            
            insert = self.get_insert_values()
            
            line += insert
            line = "\t".join(line) + "\n"
            self.parsed_lines.append(line)
        f.close()
        
        output = open(self.output_path, "w")
        output.writelines(self.parsed_lines)
        output.close()
        
        # we intially use the reporting file as input, and make a new file
        # from that. After the first set of data has been inserted, read from
        # the modified file.
        self.report_path = self.output_path
    
    def add_matches_to_report(self, matches, searched_genes, label):
        """ annotates candidate variants using proband matches indexed by gene
        
        Args:
            matches: fict of probands that were matched for each gene
            searched_genes: list of searched genes, so we can add whether a
                gene was not found because it was not searched for
            column_label: label to use in the file header
        """
        
        self.type = "matches"
        self.column_labels = [label + s for s in ["_searched", "_passed"]]
        self.data = [matches, searched_genes]
        self.modify_file()
        
        self.report_path = self.output_path
    
    def add_scores_to_report(self, scores, label):
        """ annotates candidate variants using scores for each variant
        
        Args:
            scores: dict of scores, indexed by proband, gene, inheritance mode
            label: column header label for the scores
        """
        
        self.type = "scores"
        self.column_labels = [label]
        self.data = scores
        self.modify_file()


