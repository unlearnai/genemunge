import os, json
from itertools.chain import from_iterable


FILEPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
FILENAME = os.path.join(FILEPATH, 'go.json')


class Searcher(object):

    def __init__(self):
        """
        Create a object to search through the gene ontology.

        Args:
            None

        Returns:
            Searcher

        """
        with open(FILENAME, 'r') as infile:
            self.go = json.load(infile)

    def traverse(self, term, inclusive=True):
        """
        Get all of the children of a given GO category.

        Args:
            term (str): GO id
            inclusive (optional; bool): include the given term in the final
                list of GO ids

        Returns:
            list of GO ids (List[str])

        """
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
            new_relatives = from_iterable([self.go[t]['children'] for t in new_relatives])

        return relatives

    def select_namespace(self, namespace):
        """
        Get all of the GO identifiers associated with a particular namespace.

        Args:
            namespace (str): 'biological_process' or 'molecular_function'

        Returns:
            list of go ids (List[str])


        """
        return [term for term in self.go if self.go[term]['namespace'] == namespace]

    def keyword_search(self, keywords, fields=['name', 'def'], exact=True):
        """
        Search for GO identifiers associated with some keywords.

        Args:
            keywords (List[str]): words to look for
            fields (List[str]): fields to look in
            extact (optional; bool): if true, then this function only returns
                the GO ids where there is a match. if false, then this function
                will also select all of the child GO terms.

        Returns:
            list of GO ids (List[str])

        """
        assert type(keywords) == list, \
        "keywords must be a list"
        exact_terms = [term for term in self.go if
         any(any(k in self.go[term][f] for f in fields) for k in keywords)]
        if exact:
            return exact_terms
        return list(set(from_iterable([self.traverse(t) for t in exact_terms])))

    def _get_proteins_from_term(self, term, evidence_codes):
        """
        Get all of the proteins (as UniprotKB ids) associated with a
        given GO identifier and some evidence codes.

        Args:
            term (str): a GO id
            evidence_codes (None or List[str]):

        Returns:
            proteins (List[str])

        """
        proteins = []
        for code in self.go[term]["genes"]:
            if evidence_codes is None or code in evidence_codes:
                proteins += self.go[term]["genes"][code]
        return proteins

    def get_proteins(self, terms, evidence_codes = None):
        """
        Get all of the proteins (as UniprotKB ids) associated with a list of
        GO idenifiers and some evidence codes.

        Args:
            terms (List[str]): a list of GO ids
            evidence_codes (None or List[str]):

        Returns:
            proteins (List[str])

        """
        if evidence_codes is not None:
            assert type(evidence_codes) == list, \
            "evidence_codes must be None or a list of GO evidence codes"
        all_proteins = [self._get_proteins_from_term(t, evidence_codes) for t in terms]
        return list(set().union(*all_proteins))


