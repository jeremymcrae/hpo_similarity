""" class to hold HPO terms for members of a trio
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

class FamilyHPO(object):
    """ small class to handle HPO terms for family members of a trio
    """
    
    # these should be set externally to the class, so that all members share the
    # lists
    alt_ids = {}
    obsolete_ids = set()
    
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
    
    def format_hpo(self, terms):
        """ formats a string of hpo terms to a list
        
        Args:
            terms: string of hpo terms joined with "|"
        
        Returns:
            list of hpo terms, or None
        """
        
        terms = terms.strip()
        
        # account for no hpo terms recorded for a person
        if terms == ".":
            return None
        
        # account for multiple hpo terms for an individual
        if "|" in terms:
            terms = terms.split("|")
        else:
            terms = [terms]
        
        # strip out the obsolete terms, currently there are two probands (out of
        # >4000) who each have an obsolete term, so it's not worth converting the
        # obsolete terms to a more appropriate term
        terms = [term for term in terms if term not in self.obsolete_ids]
        
        # convert each term to it's standard HPO ID if the term is in the HPO IDs,
        # otherwise just assume it is a standard HPO Id already.
        terms = [self.alt_ids[term] if term in self.alt_ids else term for term in terms]
        
        return terms
    
    def get_child_hpo(self):
        return self.child_hpo
    
    def get_maternal_hpo(self):
        return self.maternal_hpo
    
    def get_paternal_hpo(self):
        return self.paternal_hpo
