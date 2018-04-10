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
    """
    Tools to change units of expression data, primarily to convert to TPM.

    Attributes:
        gene_lengths (DataFrame): bp lengths for genes.

    """
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
    """
    The RUV-2 algorithm.

    Attributes:
        hk_genes (List[str]): a list of housekeeping gene names used in fitting.
        means (pandas.Series): the means of each gene from the training data.
        U (optional; numpy array ~ (num_training_samples, num_factors)):
            left eigenvectors from SVD of housekeeping genes in training set
        L (optional; numpy array ~ (num_factors,))
            eigenvalues from SVD of housekeeping genes in training set
        Vt (optional; numpy_array ~ (num_factors, num_hk_genes)
            right eigenvectors from SVD of housekeeping genes in training set

    """
    def __init__(self, center=True, hk_genes=None, means=None, U=None, L=None, Vt=None):
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
            center (optional; bool): whether to center the gene means in the fit.
            hk_genes (optional; List[str]): list of housekeeping genes
            means (optional; numpy array ~ (num_genes,))
            U (optional; numpy array ~ (num_training_samples, num_factors))
            L (optional; numpy array ~ (num_factors,))
            Vt (optional; numpy_array ~ (num_factors, num_hk_genes))

        Returns:
            RemoveUnwantedVariation

        """
        self.center = center
        self.hk_genes =None
        self.means = None
        self.U = None
        self.L = None
        self.Vt = None

    def _is_fit(self):
        """
        Check if the batch effect transformation has been fit.

        Args:
            None

        Returns:
            bool

        """
        return (self.hk_genes is not None) or \
               (self.means is not None) or \
               (self.U is not None) or \
               (self.L is not None) or \
               (self.Vt is not None)

    def _cutoff_svd(self, matrix, variance_cutoff=1):
        """
        Compute the singular value decomposition of a matrix and get rid
        of any singular vectors below a cumulative variance threshold.

        Args:
            matrix (numpy array): the data
            variance_cutoff (float): retains only elements of L that contribute
                to the cumulative fractional variance up to the cutoff.

        Returns:
            U, L, Vt where M = U L V^{T}

        """
        U, L, Vt = numpy.linalg.svd(matrix, full_matrices=False)
        # trim eigenvalues close to 0, exploit the fact that L is ordered
        L = L[:(~numpy.isclose(L, 0)).sum()]
        cumul_variance_fracs = numpy.cumsum(L**2) / numpy.sum(L**2)
        L_cutoff = min(len(L), 1+numpy.searchsorted(cumul_variance_fracs, variance_cutoff))
        return U[:, :L_cutoff], L[:L_cutoff], Vt[:L_cutoff, :]

    def fit(self, data, hk_genes, variance_cutoff=0.9):
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
        we can estiamte W = U L.  Additionally, A_c = V^{T}.

        Now, if we fix W and assume that X B = 0 for all genes then we can
        estimate A = W^+ Y = (W W^{T})^{-1} W^{T} Y.
        This matrix stores K patterns of variation that are
        usually not biologically interesting.

        Args:
            data (pandas.DataFrame ~ (num_samples, num_genes)): clr transformed
                expression data
            hk_genes (List[str]): list of housekeeping genes
            variance_cutoff (float): the cumulative variance cutoff on SVD
                eigenvalues of Y_c (the variance fraction of the factors).

        Returns:
            None

        """
        self.means = data.mean(axis=0)
        # restrict to available housekeeping genes
        self.hk_genes = [gene for gene in hk_genes if gene in data.columns]
        # center the data along genes
        if self.center:
            housekeeping = data[self.hk_genes] - self.means[self.hk_genes]
        else:
            housekeeping = data[self.hk_genes]
        self.U, self.L, self.Vt = self._cutoff_svd(housekeeping, variance_cutoff)

    def _delta(self, W, data_centered, penalty):
        """
        Compute the corrections for RUV2.

        Args:
            W (numpy array ~ (num_samples, num_factors))
            data_centered (pandas.DataFrame ~ (num_samples, num_genes))
            penalty (float)

        Returns:
            delta (numpy array ~ (num_samples, num_genes))

        """
        penalty_term = penalty * numpy.eye(W.shape[1])
        J = numpy.linalg.inv(penalty_term + numpy.dot(W.T, W))
        return numpy.dot(W, numpy.dot(J, numpy.dot(W.T, data_centered)))

    def transform(self, data, penalty=0):
        """
        Perform the 2-step Remove Unwanted Variation (RUV-2) algorithm.

        The `fit` method estimates the matrix
            A \in R^{K, N} which describes how the genes are coupled to the
            uninteresting factors

        We can estimate the activity of these factors from a new dataset \tilde{Y}
        by using the housekeeping genes on this new dataset and computing
        \tilde{W} = \tilde{Y}_c A_c^{+}.  Since A_c = V^{T} from the SVD,
        the right pseudoinverse A_c^{+} = A_c^{T}.

        Finally, we can subtract \tilde{W} A from the data,
        \tilde{Y} - \tilde{W} A.

        Essentially, we are removing the components of the data that project
        onto the pre-defined axes of uninteresting variation.

        Args:
            data (pandas.DataFrame ~ (num_samples, num_genes)): clr transformed
                expression data
            penalty (float): regularization on the regression step

        Returns:
            batch corrected data (pandas.DataFrame ~ (num_samples, num_genes))

        """
        assert self._is_fit(), "RUV has not been fit!"
        if self.center:
            data_trans = data - self.means
        else:
            data_trans = data
        W = numpy.dot(data_trans[self.hk_genes], self.Vt.T)
        return data - self._delta(W, data_trans, penalty)

    def fit_transform(self, data, hk_genes, penalty=0, variance_cutoff=0.9):
        """
        Perform the 2-step Remove Unwanted Variation (RUV-2) algorithm.

        Args:
            data (pandas.DataFrame ~ (num_samples, num_genes)): clr transformed
                expression data
            hk_genes (List[str]): list of housekeeping genes
            penalty (float): regularization on the regression step
            variance_cutoff (float): the cumulative variance cutoff on SVD
                eigenvalues of Y_c.

        Returns:
            batch corrected data (pandas.DataFrame ~ (num_samples, num_genes))

        """
        self.fit(data, hk_genes, variance_cutoff)
        return self.transform(data, penalty)

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
