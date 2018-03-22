import os
import pandas
import numpy

from . import convert

class Normalizer(object):

    def __init__(self, identifier='symbol'):
        """
        Tools to normalize expression data and transform into TPM.

        Args:
            None

        Returns:
            None

        """
        # read the gene lengths
        p = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'gtex')
        gene_info = pandas.read_csv(os.path.join(p, 'gene_info.csv'), sep='\t')
        gene_info.set_index('gene_id', inplace=True)
        self.gene_lengths = gene_info['bp_length']
        # clean the ensemble gene ids
        self.gene_lengths.index = list(map(convert.clean_ensembl_id, self.gene_lengths.index))
        self.gene_lengths = self.gene_lengths[~self.gene_lengths.index.duplicated(keep='first')]

    def tpm_from_rpkm(self, data, gene_list=None):
        """
        Transform data from RPKM to TPM.

        Args:
            data (pandas.DataFrame ~ (num_samples, num_genes))
            gene_list (optional; List[str]): a list of gene ids

        Returns:
            pandas.DataFrame

        """
        if gene_list is not None:
            subset = data[gene_list]
        else:
            subset = data
        return subset.divide(subset.sum(axis=1), axis='rows')

    def tpm_from_counts(self, data, gene_list=None):
        """
        Transform data from counts to TPM.

        Args:
            data (pandas.DataFrame ~ (num_samples, num_genes))
            gene_list (optional; List[str]): a list of gene ids

        Returns:
            pandas.DataFrame

        """
        if gene_list is not None:
            common_genes = list(set(gene_list) & set(self.gene_lengths.index))
        else:
            common_genes = list(self.gene_lengths.index)
        subset = data[common_genes].divide(self.gene_lengths[common_genes], axis='columns')
        return subset.divide(subset.sum(axis=1), axis='rows')


