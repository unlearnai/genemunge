import os
import pandas
import numpy
import pickle
from pathlib import Path

from . import convert

def do_nothing(data):
    """
    A function that does nothing.

    Args:
        Anything

    Returns:
        Anything

    """
    return data


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


class Normalizer(object):

    def __init__(self, identifier='symbol'):
        """
        Tools to normalize expression data and transform into TPM.

        Args:
            identifier (str)

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

    def clr_from_tpm(self, data, gene_list=None, imputer=do_nothing):
        """
        Compute the centered log ratio transform of data in TPM format.

        Args:
            data (pandas.DataFrame ~ (num_samples, num_genes))
            gene_list (optional; List[str]): a list of gene ids
            imputer (optional; callable)

        Returns:
            pandas.DataFrame

        """
        imputed = self.tpm_from_subset(imputer(data), gene_list)
        log_transformed = numpy.log(imputed)
        return log_transformed.subtract(log_transformed.mean(axis=1), axis=0)

    def tpm_from_clr(self, data, gene_list=None):
        """
        Compute data in TPM format from centered log ratio transformed data.

        Args:
            data (pandas.DataFrame ~ (num_samples, num_genes))
            gene_list (optional; List[str]): a list of gene ids

        Returns:
            pandas.DataFrame

        """
        return self.tpm_from_rpkm(numpy.exp(data), gene_list)


class RemoveUnwantedVariation(object):

    def __init__(self, alpha=None):
        """
        Perform the 2-step Remove Unwanted Variation (RUV-2) algorithm
        defined in:

        "Correcting gene expression data when neither the unwanted variation nor the
        factor of interest are observed."
        Biostatistics 17.1 (2015): 16-28.
        Laurent Jacob, Johann A. Gagnon-Bartsch, and Terence P. Speed.

        The algorithm is modified slightly so that batch correction can be
        applied out-of-sample.

        Args:
            alpha (optional; numpy array ~ (num_singular_values, num_genes))

        Returns:
            RemoveUnwantedVariation

        """
        self.alpha = alpha

    def _nonzero_svd(self, matrix):
        """
        Compute the singular value decomposition of a matrix and get rid
        of any singular vectors with a zero singular value.

        Args:
            matrix (numpy array)

        Returns:
            U, L, V where M = U L V^{T}

        """
        U, L, V = numpy.linalg.svd(matrix, full_matrices=False)
        mask = ~numpy.isclose(L, 0)
        return U[:, mask], L[mask], V[mask, :]

    def fit(self, data, hk_genes, ridge=True, nu=None):
        """
        Perform a singular value decomposition of the housekeeping genes to
        fit the transform.

        Suppose that we measure data on the expression of N genes in M samples
        and store these (after CLR transformation) in a matrix Y \in R^{M, N}.
        We consider a linear model Y = X B + W A + noise where
            X \in R^{M, Q} are some unobserved, but biologically interesting, factors
            B \in R^{Q, N} describes how the genes are coupled to the interesting factors
            W \in R^{M, K} are some unobserved and uninteresting factors
            A \in R^{K, N} describes how the genes are coupled to the uninteresting factors

        We assume that there are some housekeeping genes Y_c for which we are
        sure that B_c = 0. That is, the housekeeping genes are not coupled to
        any biologically interesting factors. Therefore, we have Y_c = W A_c + noise.
        Let Y_c = U L V^{T} be the singular value decomposition of Y_c. Then,
        we can estiamte W = U L.

        Now, if we fix W and assume that X B = 0 for all genes then we can
        estimate A = (W W^{T})^{-1} W^{T} Y. This matrix stores K patterns of
        variation that are usually not biologically interesting.

        Args:
            data (pandas.DataFrame ~ (num_samples, num_genes)): clr transformed
                expression data
            hk_genes (List[str]): list of housekeeping genes
            ridge (bool): True to use an L2 penalty in regressing A.
            nu (float): A coefficient for the L2 penalty. If None, uses a heuristic
                from the RUV-2 reference.

        Returns:
            None

        """
        # restrict to available housekeeping genes
        hk_genes_data = [gene for gene in hk_genes if gene in data.columns]
        # solve for W ~ (num_samples, num_singular_values)
        housekeeping = data[hk_genes]
        U, L, V = self._nonzero_svd(housekeeping)
        W = U * L
        # solve for alpha ~ (num_singular_values, num_genes)
        if ridge:
            if nu is None:
                sigma1 = numpy.linalg.eigvals(numpy.dot(W.T, W)).max()
                nu = 0.001 * sigma1
        else:
            nu = 0
        penalty_term = nu*numpy.eye(W.shape[1])
        self.alpha = numpy.dot(numpy.linalg.inv(numpy.dot(W.T, W) + penalty_term),
                               numpy.dot(W.T, data))
        # store inverse of inner products J ~ (num_singular_values, num_singular_values)
        self.J = numpy.linalg.inv(numpy.dot(self.alpha, self.alpha.T))

    def transform(self, data):
        """
        Perform the 2-step Remove Unwanted Variation (RUV-2) algorithm.

        The `fit` method estimates the matrix
            A \in R^{K, N} which describes how the genes are coupled to the
            uninteresting factors

        We can estimate the activity of these factors from a new dataset \tilde{Y}
        by computing \tilde{W} = \tilde{Y} A^{T} (A A^{T})^{-1} using the right
        pseudoinverse of A.

        Finally, we can subtract \tilde{W} A from the data by computing
        \tilde{Y} - \tilde{Y} A^{T} (A A^{T})^{-1} A.

        Essentially, we are removing the components of the data that project
        onto the pre-defined axes of uninteresting variation.

        Args:
            data (pandas.DataFrame ~ (num_samples, num_genes)): clr transformed
                expression data
            hk_genes (List[str]): list of housekeeping genes

        Returns:
            batch corrected data (pandas.DataFrame ~ (num_samples, num_genes))

        """
        delta = numpy.dot(numpy.dot(numpy.dot(data, self.alpha.T), self.J), self.alpha)
        return data - delta

    def fit_transform(self, data, hk_genes):
        """
        Perform the 2-step Remove Unwanted Variation (RUV-2) algorithm.

        Args:
            data (pandas.DataFrame ~ (num_samples, num_genes)): clr transformed
                expression data
            hk_genes (List[str]): list of housekeeping genes

        Returns:
            batch corrected data (pandas.DataFrame ~ (num_samples, num_genes))

        """
        self.fit(data, hk_genes)
        return self.transform(data)

    def save(self, filename, overwrite_existing=False):
        """
        Save the RUV object to filename.

        Args:
            filename (string): absolute path to save file
            overwrite_existing (bool): whether or not to overwrite existing file

        Returns:
            None

        """
        path = Path(filename)
        assert overwrite_existing or not path.exists(), \
            "Must allow overwriting existing files"
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls, filename):
        """
        Create an RUV from a saved object.

        Args:
            filename (str)

        Returns:
            RemoveUnwantedVaraition

        """
        with open(filename, 'rb') as f:
            return pickle.load(f)
