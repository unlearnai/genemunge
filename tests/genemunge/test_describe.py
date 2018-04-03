import tempfile

from genemunge import describe

import pytest


def test_describe_construct():
    """Try to construct a Describer object."""
    desc = describe.Describer(identifier='symbol')
    desc.close()


def test_describe_get_tissue_expression():
    """Try to get the tissue expression statistics for a given gene."""
    gene_name = 'TP53'
    desc = describe.Describer(identifier='symbol')
    expression = desc.get_tissue_expression(gene_name)


def test_describe_get_tissue_expression_unknown():
    """Check that an unknown gene raises an error from get_tissue_expression."""
    gene_name = 'foo'
    desc = describe.Describer(identifier='symbol')
    try:
        expression = desc.get_tissue_expression(gene_name)
    except KeyError:
        return
    assert False


def test_describe_plot_tissue_expression():
    """Try to plot the GTEx tissue expression for a given gene."""
    gene_name = 'TP53'
    desc = describe.Describer(identifier='symbol')

    with tempfile.NamedTemporaryFile() as tf:
        desc.plot_tissue_expression(gene_name, sortby='std', show=False, filename=tf.name)


def test_describe_get_gene_info():
    """Try to get name and GO info for a given gene."""
    gene_name = 'TP53'
    desc = describe.Describer(identifier='symbol')
    gene_info = desc.get_gene_info(gene_name)


if __name__ == "__main__":
    pytest.main([__file__])
