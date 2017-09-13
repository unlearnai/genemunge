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

    def _get_proteins_from_term(self, term, evidence_codes):
        proteins = []
        for code in self.go[term]["genes"]:
            if evidence_codes is None or code in evidence_codes:
                proteins += self.go[term]["genes"][code]
        return proteins

    def get_proteins(self, terms, evidence_codes = None):
        if evidence_codes is not None:
            assert type(evidence_codes) == list, \
            "evidence_codes must be None or a list of GO evidence codes"
        all_proteins = [self._get_proteins_from_term(t, evidence_codes) for t in terms]
        return list(set().union(*all_proteins))


