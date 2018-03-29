from genemunge import search

import pytest


biological_process_id = 'GO:0008150'
biological_process_namespace = 'biological_process'


def test_searcher_traverse():
    searcher = search.Searcher()
    descendants = searcher.traverse(biological_process_id)


def test_searcher_select_namespace():
    searcher = search.Searcher()
    ids = searcher.select_namespace(biological_process_namespace)


def test_searcher_traverse_vs_namespace():
    searcher = search.Searcher()
    descendants = searcher.traverse(biological_process_id)
    ids = searcher.select_namespace(biological_process_namespace)
    assert sorted(descendants) == sorted(ids)


def test_searcher_get_genes():
    searcher = search.Searcher()
    bp_genes = searcher.get_genes([biological_process_id])


def test_searcher_get_genes_allbp():
    searcher = search.Searcher()
    ids = searcher.select_namespace(biological_process_namespace)
    all_bp_genes = searcher.get_genes(ids)


def test_searcher_get_housekeeping_genes():
    searcher = search.Searcher()

    hk_genes = searcher.get_housekeeping_genes()
    # check ATF1 in housekeeping genes
    assert "ENSG00000123268" in hk_genes


def test_searcher_get_transcription_factors():
    searcher = search.Searcher()

    tfactors = searcher.get_transcription_factors()
    # check ATF1 in transcription factors
    assert "ENSG00000123268" in tfactors



if __name__ == "__main__":
    pytest.main([__file__])
