import os
import re
import gzip
import itertools
import json


FILEPATH = os.path.dirname(os.path.abspath(__file__))
GOFILE = os.path.join(FILEPATH, "go-basic.obo")
ANNOTATIONFILE = os.path.join(FILEPATH, "goa_human.gaf.gz")
OUTPUTFILE = os.path.join(FILEPATH, 'go.json')


id_pattern = 'GO:[0-9]{7}'
go_id = re.compile('id: ' + id_pattern)
is_a = re.compile('is_a: ' + id_pattern)
part_of = re.compile('part_of ' + id_pattern)
name = re.compile('name: .*')
definition = re.compile('def: .*')
namespace = re.compile('namespace: .*')
obsolete = re.compile('is_obsolete: true')


def begins_with_pattern(string, pattern):
    """
    Check if a string begins with a pattern.

    Args:
        string (str)
        pattern (str or regular expression)

    Returns:
        bool

    """
    x = re.match(pattern, string)
    if x is None:
        return False
    return True


def first_match(iterable, pattern):
    """
    Get the first item that matches a pattern.
    Return None if there are no matches.

    Args:
        iterable (Iterable[str])
        pattern (str or regular expression)

    Returns:
        str or None

    """
    for string in iterable:
        match = re.search(pattern, string)
        if match is not None:
            return match.string
    return None


def all_matches(iterable, pattern):
    """
    Get all of the items that match a pattern.

    Args:
        iterable (Iterable[str])
        pattern (str or regular expression)

    Returns:
        List[str]

    """
    result = []
    for string in iterable:
        match = re.search(pattern, string)
        if match is not None:
            result.append(match.string)
    return result


def get_id(group):
    """
    Get the GO identifier from a list of GO term properties.
    Finds the first match to the id pattern.

    Args:
        group (List[str])

    Returns:
        str

    """
    return first_match(group, go_id).split(':',1)[1].strip()


def get_name(group):
    """
    Get the name of the GO category from a list of GO term properties.
    Finds the first match to the name pattern.

    Args:
        group (List[str])

    Returns:
        str

    """
    return first_match(group, name).split(':',1)[1].strip()


def get_namespace(group):
    """
    Get the namespace of the GO category from a list of GO term properties.
    Finds the first match to the namespace pattern.

    Args:
        group (List[str])

    Returns:
        str

    """
    return first_match(group, namespace).split(':',1)[1].strip()


def get_definition(group):
    """
    Get the definition of the GO category from a list of GO term properties.
    Finds the first match to the definition pattern and removes an annotation label.

    Args:
        group (List[str])

    Returns:
        str

    """
    tmp = first_match(group, definition).split(':',1)[1].strip()
    match = re.search(r'"([^"]*)"', tmp)
    return match.groups()[0].strip('"')


def get_parents(group):
    """
    Get all of the parents of a GO category from a list of GO term properties.
    We say that X is a parent of Y if (Y 'is_a' X) or if (Y is 'part_of' X).
    Finds matches to the is_a and part_of patterns, then matches GO ids.

    Args:
        group (List[str])

    Returns:
        parents (List[str])

    """
    parents = all_matches(group, is_a) + all_matches(group, part_of)
    return [re.search(id_pattern, p).group() for p in parents]


def parse_group(group, dictionary):
    """
    Parse a GO category and add it to a dictionary.

    Notes:
        Modifies dictionary in place!

    Args:
        group (List[str])
        dictionary (dict)

    Returns:
        None

    """
    dictionary[get_id(group)] = {
            'name': get_name(group),
            'namespace': get_namespace(group),
            'def': get_definition(group),
            'parents': get_parents(group),
            'children': [],
            'genes': {
                        'EXP': [],
                        'IDA': [],
                        'IPI': [],
                        'IMP': [],
                        'IGI': [],
                        'IEP': [],
                        'ISS': [],
                        'ISO': [],
                        'ISA': [],
                        'ISM': [],
                        'IGC': [],
                        'IBA': [],
                        'IBD': [],
                        'IKR': [],
                        'IRD': [],
                        'RCA': [],
                        'TAS': [],
                        'NAS': [],
                        'IC': [],
                        'ND': [],
                        'IEA': []
                    }
            }


def make_godict(gofile, force=False):
    """
    Parses the Gene Ontology file and creates a dictionary that is easier
    to work with. Saves the dictionary as a json file.

    Notes:

        uniprot id: column 1
        gene symbol: column 2
        GO Evidence codes: column 5

        Experiment:
            Inferred from Experiment (EXP)
            Inferred from Direct Assay (IDA)
            Inferred from Physical Interaction (IPI)
            Inferred from Mutant Phenotype (IMP)
            Inferred from Genetic Interaction (IGI)
            Inferred from Expression Pattern (IEP)

        Computational:
            Inferred from Sequence or structural Similarity (ISS)
            Inferred from Sequence Orthology (ISO)
            Inferred from Sequence Alignment (ISA)
            Inferred from Sequence Model (ISM)
            Inferred from Genomic Context (IGC)
            Inferred from Biological aspect of Ancestor (IBA)
            Inferred from Biological aspect of Descendant (IBD)
            Inferred from Key Residues (IKR)
            Inferred from Rapid Divergence(IRD)
            Inferred from Reviewed Computational Analysis (RCA)

        Literature:
            Traceable Author Statement (TAS)
            Non-traceable Author Statement (NAS)

        Other:
            Inferred by Curator (IC)
            No biological Data available (ND) evidence code
            Inferred from Electronic Annotation (IEA)

    Args:
        gofile (str): path to the gene ontology file
        force (optional; bool): overwrite the json file if true

    Returns:
        None

    """
    from genemunge import convert
    converter = convert.IDConverter('uniprot_ids', 'ensembl_gene_id')

    # check if the outputfile already exists
    if not force and os.path.exists(OUTPUTFILE):
        return

    # id: {name, namespace, def, parents, children, genes}
    # connections (parent/child): 'is_a' or 'part_of'
    # ignore if 'is_obsolete: true'

    # read in the ontology file
    with open(gofile, "r") as go:
        unparsed = [line.rstrip() for line in go]

    # find the indices marking the beginning of each term
    indices = [i for i, x in enumerate(unparsed) if begins_with_pattern(x, "id:")]

    # group the terms
    grouped = [unparsed[indices[i]: indices[i+1]] for i in range(len(indices)-1)]

    # get rid of obselete terms
    not_obsolete = [g for g in grouped if first_match(g, obsolete) is None]

    # get rid of any term that doesn't have ids
    has_id = [g for g in not_obsolete if first_match(g, go_id) is not None]

    # create the go dictionary
    godict = {}
    for group in has_id:
        parse_group(group, godict)

    # add the children terms
    for term in godict:
        parents = godict[term]['parents']
        for p in parents:
            if term not in godict[p]['children']:
                godict[p]['children'] += [term]

    # add the annotations
    with gzip.open(ANNOTATIONFILE ,'rb') as annotfile:
        for raw_line in annotfile:
            line = raw_line.decode('utf-8')
            if line[0] != '!': # comments
                parsed = line.strip().split('\t')

                database = parsed[0] # currently, this is always UniProtKB
                database_id = parsed[1]
                symbol = parsed[2] # ORF for unnamed
                qualifier = parsed[3]
                go_term = parsed[4]
                database_reference = parsed[5]
                evidence = parsed[6]

                # what to do about colocalizes_with and contributes_to?
                if 'NOT' not in qualifier:
                    try:
                        ensembl = converter.convert(database_id)
                        # add the identifier if it is not NaN
                        if ensembl == ensembl:
                            godict[go_term]['genes'][evidence] += [ensembl]
                    except KeyError:
                        # we have filtered out obsolete go terms
                        # therefore, we have to catch this exception
                        pass

    # write to the file
    with open(OUTPUTFILE, "w") as outfile:
        json.dump(godict, outfile)


if __name__ == "__main__":
    make_godict(GOFILE, force=True)
