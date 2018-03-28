import os
import numpy as np
import pandas as pd

from genemunge import convert

import pytest


def test_clean_ensembl_id():
    ensembl_id = 'ENSG00000002822.15'
    expected_clean = 'ENSG00000002822'
    assert convert.clean_ensembl_id(ensembl_id) == expected_clean


def test_clean_ensembl_ids():
    ensembl_ids = ['foo.bar', 'bar.baz', 'baz.bop']
    assert convert.clean_ensembl_ids(ensembl_ids) == ['FOO', 'BAR', 'BAZ']


def test_converter_construct():
    converter = convert.IDConverter('symbol', 'name')


def test_converter_convert():
    gene_id = 'ENSG00000000003.14'
    gene_symbol = 'TSPAN6'

    cleaned_id = convert.clean_ensembl_id(gene_id)
    converter = convert.IDConverter('ensembl_gene_id', 'symbol')
    assert converter.convert(cleaned_id) == gene_symbol


def test_converter_convert_list():
    gene_ids = ['ENSG00000000003.14', 'ENSG00000000005.5', 'ENSG00000000419.12']
    gene_symbols = ['TSPAN6', 'TNMD', 'DPM1']

    cleaned_ids = convert.clean_ensembl_ids(gene_ids)
    converter = convert.IDConverter('ensembl_gene_id', 'symbol')
    assert converter.convert_list(cleaned_ids) == gene_symbols


if __name__ == "__main__":
    pytest.main([__file__])
