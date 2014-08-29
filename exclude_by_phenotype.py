

class ScreenPhenotypes(Object):
	"""
	"""
	
	def __init__(self, phenotypes, family_hpos, graph, terms_to_screen):
		"""
		"""
		
		self.phenotypes = phenotypes
		self.family_hpos = family_hpos
		self.graph = graph
		self.terms_to_screen = terms_to_screen
	
	def exclude_by_phenotype(self):
		
		for proband in probands:
			for gene in proband:
				for term in gene:
					if term in self.terms_to_screen:
						
						
						
