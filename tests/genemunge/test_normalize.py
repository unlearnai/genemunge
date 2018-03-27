import numpy
import pandas

from genemunge import normalize

import pytest


num_samples = 100
num_genes = 1000
max_read_count = 1234
expression_data = pd.DataFrame(np.round(max_read_count*np.random.rand(num_samples, num_genes)))


def test_deduplicate():
    x = np.random.rand(10, 5)
    df = pd.DataFrame(x, columns=['a', 'a', 'b', 'c', 'b'])

    df_dedup = normalize.deduplicate(df)
    assert np.allclose(x[:, [0, 1]].sum(axis=1), df_dedup.values[:,0])
    assert np.allclose(x[:, [2, 4]].sum(axis=1), df_dedup.values[:,1])
    assert np.allclose(x[:, 3].sum(axis=1), df_dedup.values[:,2])


def test_impute():
    scale = 0.5
    imputed_data = impute(expression_data, scale)

    zero_mask = expression_data != 0
    assert np.allclose((imputed_data*zero_mask).values, expression_data.values)
    rowwise_impute_expected = ~zero_mask.all(axis=1) * scale \
                                * expression_data[zero_mask].min(axis=1)
    assert np.allclose((~zero_mask).multiply(rowwise_impute_expected, axis='index'),
                       imputed_data * ~zero_mask)


def test_normalizer():




if __name__ == "__main__":
    pytest.main([__file__])
