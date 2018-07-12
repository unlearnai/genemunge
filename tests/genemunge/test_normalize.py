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
    """Check the deduplication of some data."""
    x = np.random.rand(10, 5)
    df = pd.DataFrame(x, columns=['a', 'a', 'b', 'c', 'b'])

    df_dedup = normalize.deduplicate(df)
    assert np.allclose(x[:, [0, 1]].sum(axis=1), df_dedup.values[:,0])
    assert np.allclose(x[:, [2, 4]].sum(axis=1), df_dedup.values[:,1])
    assert np.allclose(x[:, 3], df_dedup.values[:,2])


def test_impute(expression_data):
    """Check the imputation of some expression data."""
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
    """Test the RPKM -> TPM conversion for some expression data."""
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)

    rpkm = expression_data.rpkm
    tpm = expression_data.tpm
    tpm_calc = norm.tpm_from_rpkm(rpkm, gene_list=rpkm.columns)
    assert (tpm.columns == tpm_calc.columns).all()
    assert (tpm.index == tpm_calc.index).all()
    assert np.allclose(tpm.values, tpm_calc.values)


def test_normalizer_tpm_from_rpkm_allgenes(expression_data):
    """Test the RPKM -> TPM conversion for some expression data.
    Include all genes."""
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)

    rpkm = expression_data.rpkm
    tpm = expression_data.tpm
    tpm_calc = norm.tpm_from_rpkm(rpkm)
    assert (tpm.index == tpm_calc.index).all()
    assert np.allclose(tpm.values, tpm_calc[expression_data.tpm.columns].values)


def test_normalizer_tpm_from_counts(expression_data):
    """Test the counts -> TPM conversion for some expression data."""
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)

    counts = expression_data.counts
    tpm = expression_data.tpm
    tpm_calc = norm.tpm_from_counts(counts, counts.columns)

    assert (tpm.columns == tpm_calc.columns).all()
    assert (tpm.index == tpm_calc.index).all()
    assert np.allclose(tpm.values, tpm_calc.values)


def test_normalizer_tpm_from_counts_allgenes(expression_data):
    """Test the counts -> TPM conversion for some expression data.
    Include all genes."""
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)

    counts = expression_data.counts
    tpm = expression_data.tpm
    tpm_calc = norm.tpm_from_counts(counts)

    assert (tpm.index == tpm_calc.index).all()
    assert np.allclose(tpm.values, tpm_calc[expression_data.tpm.columns].values)


def test_normalizer_tpm_from_subset(expression_data):
    """Test the TPM -> TPM subset conversion for some expression data."""
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)

    tpm = expression_data.tpm
    tpm_allgenes_calc = norm.tpm_from_subset(tpm)
    assert np.allclose(tpm.values, tpm_allgenes_calc[expression_data.tpm.columns].values)

    tpm = expression_data.tpm
    tpm_fullset_calc = norm.tpm_from_subset(tpm, tpm.columns)
    assert np.allclose(tpm.values, tpm_fullset_calc.values)

    tpm_subset_calc = norm.tpm_from_subset(tpm, tpm.columns[:100])
    tpm_subset_norm = tpm_subset_calc.sum(axis=1).values
    assert np.allclose(tpm_subset_norm - 1e6, np.zeros_like(tpm_subset_norm))


def test_clr_functions(expression_data):
    """Test the TPM -> CLR and CLR -> TPM transforms for some expression data."""
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)

    tpm = normalize.impute(expression_data.tpm)
    clr = norm.clr_from_tpm(tpm, gene_list=tpm.columns)
    tpm_from_clr = norm.tpm_from_clr(clr, gene_list=clr.columns)

    assert np.allclose(tpm, tpm_from_clr)


def test_alr_functions(expression_data):
    """Test the TPM -> ALR  transform for some expression data."""
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)

    tpm = normalize.impute(expression_data.tpm)

    all_genes = list(tpm.columns)
    reference_genes = all_genes[:1]
    genes_to_keep = all_genes[1:]
    alr = norm.alr_from_tpm(tpm, reference_genes, gene_list=all_genes)
    alr_genes = list(alr.columns)

    assert (genes_to_keep == alr_genes)


def test_alr_functions_dirimpute(expression_data):
    """Test the TPM -> ALR  transform for some expression data.
    Directly impute in the alr calculation."""
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)

    all_genes = list(expression_data.tpm.columns)
    reference_genes = all_genes[:1]
    genes_to_keep = all_genes[1:]
    alr = norm.alr_from_tpm(expression_data.tpm, reference_genes,
                            gene_list=all_genes, imputer=normalize.impute)
    alr_genes = list(alr.columns)

    assert (genes_to_keep == alr_genes)


def test_gtex_reindex(expression_data):
    """Test that all of the transforms reindex to GTEx properly."""
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)
    tpm = norm.tpm_from_counts(expression_data.counts)
    clr = norm.clr_from_tpm(tpm, imputer=normalize.impute)
    assert (tpm.columns == clr.columns).all()
    assert (tpm.columns == norm.gene_lengths.index).all()


def test_zscore_from_clr(expression_data):
    """Test the z-score transformation on CLR data."""
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)

    tpm = normalize.impute(expression_data.tpm)
    clr = norm.clr_from_tpm(tpm, gene_list=tpm.columns)

    tissues = pd.Series('Liver', index=clr.index)
    zscore = norm.z_score_from_clr(clr, tissues)


def test_ordinalize(expression_data):
    """Test the ordinalize transformation on CLR data."""
    identifier = 'symbol'
    norm = normalize.Normalizer(identifier=identifier)

    tpm = normalize.impute(expression_data.tpm)
    clr = norm.clr_from_tpm(tpm, gene_list=tpm.columns)

    cutoffs = [0.37]
    min_value = 5
    ords = norm.ordinalize(clr, cutoffs, min_value=min_value)
    assert ((clr <= cutoffs[0]) == (ords == min_value)).all().all()
    assert ((clr > cutoffs[0]) == (ords == 1+min_value)).all().all()


def test_remove_unwanted_variation_noX():
    """Test the RUV2 implementation for data with no X."""
    num_samples = 100
    num_genes = 1000
    num_hidden_factors = 10

    W = np.random.randn(num_samples, num_hidden_factors)
    alpha = np.random.randn(num_hidden_factors, num_genes)

    Y = pd.DataFrame(np.dot(W, alpha))
    ruv = normalize.RemoveUnwantedVariation(center=True)
    Y_tilde = ruv.fit_transform(Y, np.arange(num_genes), variance_cutoff=1)

    assert np.allclose(Y.mean(axis=0), Y_tilde.mean(axis=0))
    assert np.allclose(Y_tilde - Y_tilde.mean(axis=0), 0)


def test_remove_unwanted_variation():
    """Test that the RUV2 implementation correctly infers the number
    of irrelevant factors on a constructed example."""
    num_samples = 200
    num_genes = 1000
    num_hk_genes = 100
    num_X_factors = 20
    num_W_factors = 10

    X = np.random.randn(num_samples, num_X_factors)
    beta = np.random.randn(num_X_factors, num_genes)

    W = np.random.randn(num_samples, num_W_factors)
    alpha = np.random.randn(num_W_factors, num_genes)

    Yc = np.dot(W, alpha)[:, :num_hk_genes]
    Yr = np.dot(X, beta)[:, num_hk_genes:] + np.dot(W, alpha)[:, num_hk_genes:]

    Y = pd.DataFrame(np.hstack([Yc, Yr]))
    ruv = normalize.RemoveUnwantedVariation(center=False)
    Y_tilde = ruv.fit_transform(Y, np.arange(num_hk_genes), variance_cutoff=1)

    # check the W factor count estimate
    assert ruv.L.shape[0] == num_W_factors


def test_remove_unwanted_variation_component_cut():
    """Test that the RUV2 implementation correctly trims the number of
    irrelevant factors on a constructed example."""
    num_samples = 200
    num_genes = 1000
    num_hk_genes = 100
    num_X_factors = 20
    num_W_factors = 10
    num_components = 3

    X = np.random.randn(num_samples, num_X_factors)
    beta = np.random.randn(num_X_factors, num_genes)

    W = np.random.randn(num_samples, num_W_factors)
    alpha = np.random.randn(num_W_factors, num_genes)

    Yc = np.dot(W, alpha)[:, :num_hk_genes]
    Yr = np.dot(X, beta)[:, num_hk_genes:] + np.dot(W, alpha)[:, num_hk_genes:]

    Y = pd.DataFrame(np.hstack([Yc, Yr]))
    ruv = normalize.RemoveUnwantedVariation(center=False)
    Y_tilde = ruv.fit_transform(Y, np.arange(num_hk_genes), variance_cutoff=1,
                                num_components=num_components)

    # check the W factor count estimate
    assert ruv.L.shape[0] == num_components


if __name__ == "__main__":
    pytest.main([__file__])
