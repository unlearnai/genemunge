import os, re

filepath = os.path.dirname(os.path.abspath(__file__))
gofile = os.path.join(filepath, "go-basic.obo")

id_pattern = 'GO:[0-9]{7}'
go_id = re.compile('id: ' + id_pattern)
is_a = re.compile('is_a: ' + id_pattern)
part_of = re.compile('part_of ' + id_pattern)
name = re.compile('name: .*')
definition = re.compile('def: .*')
namespace = re.compile('namespace: .*')
obsolete = re.compile('is_obsolete: true')

def has_pattern(string, pattern):
    x = re.match(pattern, string)
    if x is None:
        return False
    else:
        return True

def first_match(iterable, pattern):
    for string in iterable:
        match = re.search(pattern, string)
        if match is not None:
            return match.string
    return None

def all_matches(iterable, pattern):
    result = []
    for string in iterable:
        match = re.search(pattern, string)
        if match is not None:
            result.append(match.string)
    return result

def get_id(group):
    return first_match(group, go_id).strip().split(':',1)[1].strip()

def get_name(group):
    return first_match(group, name).strip().split(':',1)[1].strip()

def get_namespace(group):
    return first_match(group, namespace).strip().split(':',1)[1].strip()

def get_definition(group):
    tmp = first_match(group, definition).strip().split(':',1)[1].strip()
    match = re.search(r'"([^"]*)"', tmp)
    return match.groups()[0].strip('"')

def get_parents(group):
    parents = all_matches(group, is_a) + all_matches(group, part_of)
    return [re.search(id_pattern, p).group() for p in parents]

def parse_group(group, dictionary):
    dictionary[get_id(group)] = {
            'name': get_name(group),
            'namespace': get_namespace(group),
            'def': get_definition(group),
            'parents': get_parents(group),
            'genes': []
            }

def make_godict(gofile):
    # id: name, namespace, def, parents, children, genes
    # connections (parent/child): 'is_a' or 'part_of'
    # ignore if 'is_obsolete: true'

    # read in the ontology file
    with open(gofile, "r") as go:
        unparsed = go.readlines()

    # find the indices marking the beginning of each term
    indices = [i for i, x in enumerate(unparsed) if has_pattern(x, "id:")]

    # group the terms
    grouped = [unparsed[indices[i]: indices[i+1]] for i in range(len(indices)-1)]

    # get rid of obselete terms
    not_obsolete = [g for g in grouped if first_match(g, obsolete) is None]

    # get rid of any term that don't have ids
    has_id = [g for g in not_obsolete if first_match(g, go_id) is not None]

    # create the go dictionary
    godict = {}
    for group in has_id:
        parse_group(group, godict)

    return godict

godict = make_godict(gofile
                     )
