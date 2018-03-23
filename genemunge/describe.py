import os
import pandas
from cytoolz import identity

from . import convert
from . import search

datapath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
gtexpath = os.path.join(datapath, 'gtex')

class Describer(object):

    def __init__(self, identifier='symbol'):
        """
        Create an object to grab the information that describes a gene.

        Args:
            identifier (optional; str): the type of gene identifer you will use
                e.g., 'symbol', 'ensembl_gene_id'

        Returns:
            Describer

        """
        if identifier is not 'ensembl_gene_id':
            self.get_ensembl = convert.IDConverter(identifier, 'ensembl_gene_id').convert
        else:
            self.get_ensembl = convert.clean_ensembl_id
        self.get_name = convert.IDConverter('ensembl_gene_id', 'name').convert
        self.get_symbol = convert.IDConverter('ensembl_gene_id', 'symbol').convert
        self.searcher = search.Searcher()

    def get_tissue_expression(self, identifier):
        """

        """
        pass


    def plot_tissue_expression(self, idenifier):
        """

        """
        pass

