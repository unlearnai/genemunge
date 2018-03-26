import os, json
from itertools import chain


FILEPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
GONAME = os.path.join(FILEPATH, 'go.json')
ATTRIBUTENAME = os.path.join(FILEPATH, 'gene_attributes.json')


class Searcher(object):

    def __init__(self):
        """
        Create a object to search through the gene ontology.

        Args:
            None

        Returns:
            Searcher

        """
        with open(GONAME, 'r') as infile:
            self.go = json.load(infile)
        with open(ATTRIBUTENAME, 'r') as infile:
            self.attributes = json.load(infile)

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
            new_relatives = chain.from_iterable([self.go[t]['children'] for t in new_relatives])

        return relatives

    def select_namespace(self, namespace):
        """
        Get all of the GO identifiers associated with a particular namespace.

        Args:
            namespace (str): 'biological_process' or 'molecular_function'

        Returns:
            list of go ids (List[str])


        """
        return sorted([term for term in self.go if self.go[term]['namespace'] == namespace])

    def _keyword_match(self, term, keyword, fields):
        """
        Check if any of the fields of a given term of the gene ontology
        are (or are not) associatedd with a given keyword.

        Args:
            term (str): a GO id
            keyword (str): keyword to look for
            fields (List[str]): fields to look in
        Returns:
            bool

        """
        return any(keyword in self.go[term][f] for f in fields)

    def keyword_search(self, keywords, fields=['name', 'def'], exact=True,
                       exclude_keywords=None, exclude_ids=None):
        """
        Search for GO identifiers associated with some keywords.

        Args:
            keywords (List[str]): words to look for
            fields (List[str]): fields to look in
            extact (optional; bool): if true, then this function only returns
                the GO ids where there is a match. if false, then this function
                will also select all of the child GO terms.
            exclude_keyworks (optional; List[str]): do NOT include GO categories
                that contain these keywords
            exclude_ids (optional; List[str]): do NOT include these GO categories

        Returns:
            list of GO ids (List[str])

        """
        assert type(keywords) == list, "keywords must be a list"
        matches = [term for term in self.go if
                   any(self._keyword_match(term, k, fields) for k in keywords)]
        anti_matches = [] if exclude_keywords is None else [term for term in self.go if
                   any(self._keyword_match(term, k, fields) for k in exclude_keywords)]
        anti_ids = [] if exclude_ids is None else exclude_ids
        # TODO: this exclusion works with 'exact' but doesn't exclude
        # things that are picked up while traversing the graph
        exact_terms = list(set(matches) - set(anti_matches) - set(anti_ids))
        if exact:
            return exact_terms
        return sorted(list(set(chain.from_iterable([self.traverse(t) for t in exact_terms]))))

    def _get_proteins_from_term(self, term, evidence_codes):
        """
        Get all of the genes associated with a
        given GO identifier and some evidence codes.

        Args:
            term (str): a GO id
            evidence_codes (None or List[str]):

        Returns:
            genes (List[str]): list of genes by UniprotKB id

        """
        proteins = []
        for code in self.go[term]["genes"]:
            if evidence_codes is None or code in evidence_codes:
                proteins += self.go[term]["genes"][code]
        return proteins

    def get_genes(self, terms, evidence_codes = None):
        """
        Get all of the genes associated with a list of
        GO idenifiers and some evidence codes.

        Args:
            terms (List[str]): a list of GO ids
            evidence_codes (None or List[str]):

        Returns:
            genes (List[str]): list of genes by ensembl_gene_id

        """
        if evidence_codes is not None:
            assert type(evidence_codes) == list, \
            "evidence_codes must be None or a list of GO evidence codes"
        all_proteins = [self._get_proteins_from_term(t, evidence_codes) for t in terms]
        return sorted(list(set().union(*all_proteins)))

    def get_housekeeping_genes(self):
        """
        Get a list of genes that are designated to be "housekeeping genes".

        Args:
            None

        Returns:
            genes (List[str]): list of genes by ensembl_gene_id

        """
        return sorted(self.attributes["housekeeping_genes"])

    def get_transcription_factors(self):
        """
        Get a list of genes that are designated to be transcription factors
        in humans.

        Args:
            None

        Returns:
            genes (List[str]): list of genes by ensembl_gene_id

        """
        return sorted(self.attributes["transcription_factors"])


