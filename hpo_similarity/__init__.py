from pkg_resources import get_distribution

__version__ = get_distribution('hpo_similarity').version

from hpo_similarity.ontology import open_ontology
