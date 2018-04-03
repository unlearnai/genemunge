from genemunge import search

import pytest


biological_process_id = 'GO:0008150'
biological_process_namespace = 'biological_process'


def test_searcher_traverse():
    """Try to traverse the GO for a given GO ID."""
    searcher = search.Searcher()
    descendants = searcher.traverse(biological_process_id)


def test_searcher_select_namespace():
    """Try to select all IDs in a given namespace."""
    searcher = search.Searcher()
    ids = searcher.select_namespace(biological_process_namespace)


def test_searcher_traverse_vs_namespace():
    """Check that traversing the GO with a namespace ID returns the same result
    as selecting the IDs in the namesapce."""
    searcher = search.Searcher()
    descendants = searcher.traverse(biological_process_id)
    ids = searcher.select_namespace(biological_process_namespace)
    assert sorted(descendants) == sorted(ids)


def test_searcher_get_genes():
    """Try to get genes associated with a given GO ID."""
    searcher = search.Searcher()
    bp_genes = searcher.get_genes([biological_process_id])


def test_searcher_get_genes_allbp():
    """Try to get all genes associated with IDs in a namespace."""
    searcher = search.Searcher()
    ids = searcher.select_namespace(biological_process_namespace)
    all_bp_genes = searcher.get_genes(ids)


def test_searcher_get_housekeeping_genes():
    """Try to get the list of housekeeping genes.  Check a known HK gene."""
    searcher = search.Searcher()

    hk_genes = searcher.get_housekeeping_genes()
    # check ATF1 in housekeeping genes
    assert "ENSG00000123268" in hk_genes


def test_searcher_get_transcription_factors():
    """Try to get the list of transcription factors.  Check a known TF."""
    searcher = search.Searcher()

    tfactors = searcher.get_transcription_factors()
    # check ATF1 in transcription factors
    assert "ENSG00000123268" in tfactors


if __name__ == "__main__":
    pytest.main([__file__])
