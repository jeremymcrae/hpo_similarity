""" get a weighted random selection, from:
http://stackoverflow.com/questions/3679694/a-weighted-version-of-random-choice
"""

import random
import bisect

class WeightedChoice(object):
    """ class for weighted choices, rather than a function, so we don't have to 
    keep resumming the weights.
    """
    
    def __init__(self, choices):
        """ set up a list of cumulative probabilities
        
        Args:
            choices: list of (name, probability) tuples
        """
        self.choices = choices
        
        # make a list of the cumulative probabilities
        cum_prob = 0
        self.cum_probs = []
        for pair in self.choices:
            cum_prob += pair[1]
            self.cum_probs.append(cum_prob)
        
    def choice(self):
        """ chooses a random element using a set of probability weights
        
        Returns:
            the name of the randomly selected element (e.g. position)
        """
        
        # figure out where in the list a random probability would fall
        r = random.uniform(0, self.cum_probs[-1])
        pos = bisect.bisect_left(self.cum_probs, r)
        
        return self.choices[pos][0]
    
    def get_positions(self):
        """ gets the allowed CDS positions/names for the weighted choices
        """
        
        # we might want to reload this every simulation, just generate it once
        if hasattr(self, "sites"):
            return self.sites
        else:
            self.sites = set([])
            for (name, prob) in self.choices:
                self.sites.add(name)
        
        return self.sites
    
    def choose_amongst_positions(self, sites):
        """ restrict sampling to a specific set of positions
        """
        
        cum_prob = 0
        cum_probs = []
        positions = []
        for (position, probability) in self.choices:
            if position not in sites:
                continue
                
            cum_prob += probability
            cum_probs.append(cum_prob)
            positions.append(position)
        
        # figure out where in the list a random probability would fall
        r = random.uniform(0, cum_probs[-1])
        pos = bisect.bisect_left(cum_probs, r)
        
        return positions[pos]
            
            
        
        
