import os
import pandas
import numpy

from . import convert


def deduplicate(data):
    """
    Adds the values from any duplicated genes.

    Args:
        data (pandas.DataFrame ~ (num_samples, num_genes))

    Returns:
        pandas.DataFrame

    """
    return data.groupby(data.columns, axis=1).sum()


def impute(data, scale=0.5):
    """
    Replace any zeros in each row with a fraction of the smallest non-zero
    value in the corresponding row.

    Args:
        data (pandas.DataFrame ~ (num_samples, num_genes))
        scale (optional; float)

    Returns:
        imputed data (pandas.DataFrame ~ (num_samples, num_genes))

    """
    v = scale * data[data > 0].min(axis=1)
    return data.fillna(0).T.replace(to_replace=0, value=v).T


def remove_unwanted_variation(Y, hk_genes):
    """
    Perform the 2-step Remove Unwanted Variation (RUV-2) algorithm.

    "Correcting gene expression data when neither the unwanted variation nor the
    factor of interest are observed."
    Biostatistics 17.1 (2015): 16-28.
    Jacob, Laurent, Johann A. Gagnon-Bartsch, and Terence P. Speed.

    Args:
        data (pandas.DataFrame ~ (num_samples, num_genes))
        hk_genes (List[str]): list of housekeeping genes

    Returns:
        batch corrected data (pandas.DataFrame ~ (num_samples, num_genes))

    """
    H = Y[hk_genes]
    U, L, V = numpy.linalg.svd(H)
    W = U * L
    B = numpy.dot(W, numpy.dot(numpy.linalg.inv(numpy.dot(W.T, W)), numpy.dot(W.T, Y)))
    return Y - B


class Normalizer(object):

    def __init__(self, identifier='symbol'):
        """
        Tools to normalize expression data and transform into TPM.

        Args:
            identifer (str)

        Returns:
            Normalizer

        """
        # read the gene lengths
        p = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'gtex')
        gene_info = pandas.read_csv(os.path.join(p, 'gene_info.csv'), sep='\t')
        gene_info.set_index('gene_id', inplace=True)
        self.gene_lengths = gene_info['bp_length']
        # clean the ensemble gene ids
        self.gene_lengths.index = convert.clean_ensembl_ids(self.gene_lengths.index)
        self.gene_lengths = self.gene_lengths[~self.gene_lengths.index.duplicated(keep='first')]
        # convert the gene ids
        if identifier is not 'ensembl_gene_id':
            c = convert.IDConverter('ensembl_gene_id', identifier)
            self.gene_lengths.index = c.convert_list(list(self.gene_lengths.index))
        # drop any NaN and duplicate ids
        self.gene_lengths = self.gene_lengths[~self.gene_lengths.index.isnull()]
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
        return 10**6 * subset.divide(subset.sum(axis=1), axis='index')

    def tpm_from_counts(self, data, gene_list=None):
        """
        Transform data from counts to TPM.
        Any genes not in the gene_lengths index is removed,
            as the gene length is not known.

        Args:
            data (pandas.DataFrame ~ (num_samples, num_genes))
            gene_list (optional; List[str]): a list of gene ids

        Returns:
            pandas.DataFrame

        """
        if gene_list is not None:
            common_genes = [gene for gene in gene_list if gene in self.gene_lengths.index]
        else:
            common_genes = [gene for gene in data.columns if gene in self.gene_lengths.index]
        subset = data[common_genes].divide(self.gene_lengths[common_genes], axis='columns')
        return 10**6 * subset.divide(subset.sum(axis=1), axis='rows')

    def tpm_from_subset(self, data, gene_list=None):
        """
        Renormalize a subset of genes already in TPM.

        Args:
            data (pandas.DataFrame ~ (num_samples, num_genes))
            gene_list (optional; List[str]): a list of gene ids

        Returns:
            pandas.DataFrame

        """
        return self.tpm_from_rpkm(data, gene_list)
