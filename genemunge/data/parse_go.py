import os, re

filepath = os.path.dirname(os.path.abspath(__file__))
gofile = os.path.join(filepath, "go-basic.obo")

patter ='id:'

def has_pattern(string, pattern):
    x = re.match(pattern, string)
    if x is None:
        return False
    else:
        return True

def split_entry(entry):
    return [line.strip().split(':', 1) for line in entry]

def parse_entry(entry):
    pass

# id: name, namespace, def, parents, children
# connections (parent/child): 'is_a' or 'part_of' without 'not'
# ignore if 'is_obsolete: true'

# read in the ontology file
with open(gofile, "r") as go:
    unparsed = go.readlines()

# find the indices marking the beginning of each term
indices = [i for i, x in enumerate(unparsed) if has_pattern(x, "id:")]

# group the terms
grouped = [unparsed[indices[i]: indices[i+1]] for i in range(len(indices)-1)]

