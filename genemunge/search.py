import os, json, pandas
from itertools import chain


FILEPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
GONAME = os.path.join(FILEPATH, 'go.json')
ATTRIBUTENAME = os.path.join(FILEPATH, 'gene_attributes.json')
GTEXPATH = os.path.join(FILEPATH, 'gtex')


class Searcher(object):
    """
    A utility for searching the Gene Ontology.

    Attributes:
        go (dict): the GO data.
        attributes (dict): gene attributes

    """
    def __init__(self):
        """
        Create a object to search through the Gene Ontology.

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
        descendants = set([term]) if inclusive else set([])
        children = [term]
        descendant_count = -1

        while descendant_count < len(descendants):
            descendant_count = len(descendants)
            children = set(sum([self.go[t]['children'] for t in children], []))
            children = children - set(descendants)
            descendants = descendants.union(children)

        return list(descendants)

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
        are (or are not) associated with a given keyword.

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
        codes = evidence_codes or self.go[term]["genes"].keys()
        for code in codes:
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

    def get_control_genes(self, cutoff):
        """
        Get a list of genes that are designated to be "housekeeping genes"
        and that have similar expression across tissues in GTEx.

        Each gene has been assigned a score (the 'Hellinger' distance) that
        describes the similarity in its expression across tissues. The score
        ranges between 0 (all tissues have the same expression) to 1 (at least
        one pair of tissues are easily distinguished from the expression of the
        gene).

        Args:
            cutoff (float \in [0,1]):

        Returns:
            genes (List[str]): list of genes by ensembl_gene_id

        """
        tissue_status_filename = os.path.join(GTEXPATH, 'tissue_stats.h5')
        hellinger = pandas.read_hdf(tissue_status_filename, 'hellinger')
        hk_genes = hellinger.loc[self.get_housekeeping_genes()]
        return list(hk_genes[hk_genes['hellinger'] < cutoff].index)

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
