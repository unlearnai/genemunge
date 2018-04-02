import numpy as np
import pandas as pd
from collections import namedtuple

from genemunge import normalize

import pytest

np.random.seed(137)

ExpressionData = namedtuple("ExpressionData", ["counts", "tpm", "rpkm"])

@pytest.fixture
def expression_data():
    """Create some example data."""
    num_samples = 100
    num_genes = 1000
    max_read_count = 1234
    counts = pd.DataFrame(np.round(max_read_count*np.random.rand(num_samples, num_genes)))

    # create a Normalizer object and get gene lengths
    norm = normalize.Normalizer(identifier='symbol')
    gene_lengths = norm.gene_lengths[:num_genes]
    counts.columns = gene_lengths.index
    # TPM
    tpk = counts.divide(gene_lengths/1e3, axis='columns')
    tpm = tpk.divide(tpk.sum(axis=1)/1e6, axis='index')
    # RPKM
    cpm = counts.divide(counts.sum(axis=1)/1e6, axis='index')
    rpkm = cpm.divide(gene_lengths/1e3, axis='columns')

    return ExpressionData(counts, tpm, rpkm)


def test_deduplicate():
    x = np.random.rand(10, 5)
    df = pd.DataFrame(x, columns=['a', 'a', 'b', 'c', 'b'])

    df_dedup = normalize.deduplicate(df)
    assert np.allclose(x[:, [0, 1]].sum(axis=1), df_dedup.values[:,0])
    assert np.allclose(x[:, [2, 4]].sum(axis=1), df_dedup.values[:,1])
    assert np.allclose(x[:, 3], df_dedup.values[:,2])


def test_impute(expression_data):
    scale = 0.5
    counts = expression_data.counts
    imputed_data = normalize.impute(counts, scale)

    zero_mask = counts != 0
    assert np.allclose((imputed_data*zero_mask).values, counts.values)
    rowwise_impute_expected = ~zero_mask.all(axis=1) * scale \
                                * counts[zero_mask].min(axis=1)
    assert np.allclose((~zero_mask).multiply(rowwise_impute_expected, axis='index'),
                       imputed_data * ~zero_mask)


def test_normalizer_tpm_from_rpkm(expression_data):
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)

    rpkm = expression_data.rpkm
    tpm = expression_data.tpm
    tpm_calc = norm.tpm_from_rpkm(rpkm)
    assert (tpm.columns == tpm_calc.columns).all()
    assert (tpm.index == tpm_calc.index).all()
    assert np.allclose(tpm.values, tpm_calc.values)


def test_normalizer_tpm_from_counts(expression_data):
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)

    counts = expression_data.counts
    tpm = expression_data.tpm
    tpm_calc = norm.tpm_from_counts(counts)
    assert (tpm.columns == tpm_calc.columns).all()
    assert (tpm.index == tpm_calc.index).all()
    assert np.allclose(tpm.values, tpm_calc.values)


def test_normalizer_tpm_from_subset(expression_data):
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)

    tpm = expression_data.tpm
    tpm_fullset_calc = norm.tpm_from_subset(tpm)
    assert np.allclose(tpm.values, tpm_fullset_calc.values)

    tpm_subset_calc = norm.tpm_from_subset(tpm, tpm.columns[:100])
    tpm_subset_norm = tpm_subset_calc.sum(axis=1).values
    assert np.allclose(tpm_subset_norm - 1e6, np.zeros_like(tpm_subset_norm))


def test_clr_functions(expression_data):
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)

    tpm = normalize.impute(expression_data.tpm)
    clr = norm.clr_from_tpm(tpm)
    tpm_from_clr = norm.tpm_from_clr(clr)

    assert np.allclose(tpm, tpm_from_clr)


def test_remove_unwanted_variation():
    num_samples = 100
    num_genes = 1000
    num_hk_genes = 100
    num_hidden_factors = 10

    W = np.random.randn(num_samples, num_hidden_factors)
    alpha = np.random.randn(num_hidden_factors, num_genes)

    Y = pd.DataFrame(np.dot(W, alpha))

    ruv = normalize.RemoveUnwantedVariation()
    Y_tilde = ruv.fit_transform(Y, np.arange(num_genes))

    assert False





if __name__ == "__main__":
    pytest.main([__file__])
