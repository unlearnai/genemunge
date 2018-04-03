from genemunge import convert

import pytest


def test_clean_ensembl_id():
    """Try to clean an Ensembl ID."""
    ensembl_id = 'ENSG00000002822.15'
    expected_clean = 'ENSG00000002822'
    assert convert.clean_ensembl_id(ensembl_id) == expected_clean


def test_clean_ensembl_ids():
    """Try to clean a list of Ensembl IDs."""
    ensembl_ids = ['foo.bar', 'bar.baz', 'baz.bop']
    assert convert.clean_ensembl_ids(ensembl_ids) == ['FOO', 'BAR', 'BAZ']


def test_converter_construct():
    """Try to construct an IDConverter object."""
    converter = convert.IDConverter('symbol', 'name')


def test_converter_convert():
    """Try to convert an example Ensembl ID to a gene symbol."""
    gene_id = 'ENSG00000000003.14'
    gene_symbol = 'TSPAN6'

    cleaned_id = convert.clean_ensembl_id(gene_id)
    converter = convert.IDConverter('ensembl_gene_id', 'symbol')
    assert converter.convert(cleaned_id) == gene_symbol


def test_converter_all_targets():
    """Try to convert an example Ensembl ID to all allowed identifier types."""
    conversion_targets = convert.IDConverter.potential_ids
    gene_id = 'ENSG00000000003.14'
    cleaned_id = convert.clean_ensembl_id(gene_id)

    for target_id in conversion_targets:
        converter = convert.IDConverter('ensembl_gene_id', target_id)
        converter.convert(cleaned_id)


def test_converter_convert_list():
    """Try to convert a list of Ensembl IDs to gene symbols."""
    gene_ids = ['ENSG00000000003.14', 'ENSG00000000005.5', 'ENSG00000000419.12']
    gene_symbols = ['TSPAN6', 'TNMD', 'DPM1']

    cleaned_ids = convert.clean_ensembl_ids(gene_ids)
    converter = convert.IDConverter('ensembl_gene_id', 'symbol')
    assert converter.convert_list(cleaned_ids) == gene_symbols


if __name__ == "__main__":
    pytest.main([__file__])
