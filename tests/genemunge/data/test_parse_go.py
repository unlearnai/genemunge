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
    string = example_group[0]
    pattern = parse_go.go_id
    assert parse_go.begins_with_pattern(string, pattern)

    string = "foo " + example_group[0]
    assert not parse_go.begins_with_pattern(string, pattern)


def test_first_match():
    pattern = parse_go.name
    assert parse_go.first_match(example_group, pattern) == example_group[1]


def test_all_matches():
    pattern = parse_go.is_a
    assert parse_go.all_matches(example_group, pattern) == example_group[5:7]


def test_get_id():
    assert parse_go.get_id(example_group) == 'GO:0000083'


def test_get_name():
    assert parse_go.get_name(example_group) == \
        'regulation of transcription involved in G1/S transition of mitotic cell cycle'


def test_get_namespace():
    assert parse_go.get_namespace(example_group) == 'biological_process'


def test_get_definition():
    assert parse_go.get_definition(example_group) == \
        ('Any process that regulates transcription such that the target genes'
         ' are involved in the transition between G1 and S phase of the mitotic cell cycle.')


def test_get_parents():
    assert parse_go.get_parents(example_group) == \
        ['GO:0006355', 'GO:1903047', 'GO:0000082']


if __name__ == "__main__":
    pytest.main([__file__])
