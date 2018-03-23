import os
import pandas

from . import convert
from . import search


datapath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
gtexpath = os.path.join(datapath, 'gtex')


class Describer(object):

    __stats__ = ['mean', 'median', 'std', 'lower_quartile', 'upper_quartile']

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
        self.tissue_stats = pandas.HDFStore(os.path.join(gtexpath, 'tissue_stats.h5'))

    def close(self):
        """

        """
        self.tissue_stats.close()

    def get_tissue_expression(self, identifier):
        """

        """
        i = self.get_ensembl(identifier)
        return pandas.concat(
                {k: self.tissue_stats[k].loc[i] for k in self.__stats__}, axis=1)

    def plot_tissue_expression(self, idenifier):
        """

        """
        pass

