from genemunge.data import parse_go

import pytest


example_group = [
    'id: GO:0000083',
    'name: regulation of transcription involved in G1/S transition of mitotic cell cycle',
    'namespace: biological_process',
    ('def: "Any process that regulates transcription such that the target genes'
     ' are involved in the transition between G1 and S phase of the mitotic cell cycle."'
     ' [GOC:mtg_cell_cycle]'),
    'xref: Reactome:REACT_683 "G1/S-Specific Transcription, Homo sapiens"',
    'is_a: GO:0006355 ! regulation of transcription, DNA-templated',
    'is_a: GO:1903047 ! mitotic cell cycle process',
    'relationship: part_of GO:0000082 ! G1/S transition of mitotic cell cycle',
    '',
    '[Term]']


def test_begins_with_pattern():
    """Check that the begins_with_pattern function correctly matches
    the beginning of a string on examples."""
    string = example_group[0]
    pattern = parse_go.go_id
    assert parse_go.begins_with_pattern(string, pattern)

    string = "foo " + example_group[0]
    assert not parse_go.begins_with_pattern(string, pattern)


def test_first_match():
    """Check that the first_match function correctly finds
    the first match to a pattern."""
    pattern = parse_go.name
    assert parse_go.first_match(example_group, pattern) == example_group[1]


def test_all_matches():
    """Check that the all_matches function correctly finds
    all matches to a pattern."""
    pattern = parse_go.is_a
    assert parse_go.all_matches(example_group, pattern) == example_group[5:7]


def test_get_id():
    """Check that a GO ID can be retrieved from example parsed GO data."""
    assert parse_go.get_id(example_group) == 'GO:0000083'


def test_get_name():
    """Check that a gene name can be retrieved from example parsed GO data."""
    assert parse_go.get_name(example_group) == \
        'regulation of transcription involved in G1/S transition of mitotic cell cycle'


def test_get_namespace():
    """Check that a GO namesapce can be retrieved from example parsed GO data."""
    assert parse_go.get_namespace(example_group) == 'biological_process'


def test_get_definition():
    """Check that a definition can be retrieved from example parsed GO data."""
    assert parse_go.get_definition(example_group) == \
        ('Any process that regulates transcription such that the target genes'
         ' are involved in the transition between G1 and S phase of the mitotic cell cycle.')


def test_get_parents():
    """Check that the GO ID parents can be retrieved from example parsed GO data."""
    assert parse_go.get_parents(example_group) == \
        ['GO:0006355', 'GO:1903047', 'GO:0000082']


if __name__ == "__main__":
    pytest.main([__file__])
