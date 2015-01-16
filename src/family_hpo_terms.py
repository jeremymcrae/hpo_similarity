""" class to hold HPO terms for members of a trio
"""

class FamilyHPO(object):
    """ small class to handle HPO terms for family members of a trio
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
    
    def format_hpo(self, hpo_terms):
        """ formats a string of hpo terms to a list
        
        Args:
            hpo_terms: string of hpo terms joined with "|"
        
        Returns:
            list of hpo terms, or None
        """
        
        hpo_terms = hpo_terms.strip()
        
        # account for no hpo terms recorded for a person
        if hpo_terms == ".":
            return None
        
        # account for multiple hpo terms for an individual
        if "|" in hpo_terms:
            hpo_terms = hpo_terms.split("|")
        else:
            hpo_terms = [hpo_terms]
        
        return hpo_terms
    
    def get_child_hpo(self):
        return self.child_hpo
    
    def get_maternal_hpo(self):
        return self.maternal_hpo
    
    def get_paternal_hpo(self):
        return self.paternal_hpo
