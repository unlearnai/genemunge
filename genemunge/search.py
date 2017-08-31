import os, json, itertools

FILEPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
FILENAME = os.path.join(FILEPATH, 'go.json')

class Searcher(object):

    def __init__(self):
        with open(FILENAME, 'r') as infile:
            self.go = json.load(infile)

    def traverse(self, term, inclusive=True):
        if inclusive:
            relatives = [term]
        else:
            relatives = []

        new_relatives = self.go[term]['children']

        while True:
            update = list(set(relatives).union(new_relatives))
            if len(relatives) == len(update):
                break
            else:
                relatives = update
            new_relatives = itertools.chain.from_iterable([self.go[t]['children'] for t in new_relatives])

        return relatives

    def select_namespace(self, namespace):
        return [term for term in self.go if self.go[term]['namespace'] == namespace]

    def keyword_search(self, keywords, fields=['name', 'def'], exact=True):
        assert type(keywords) == list, \
        "keywords must be a list"
        exact_terms = [term for term in self.go if
         any(any(k in self.go[term][f] for f in fields) for k in keywords)]
        if exact:
            return exact_terms
        else:
            return list(set(itertools.chain.from_iterable([self.traverse(t) for t in exact_terms])))


