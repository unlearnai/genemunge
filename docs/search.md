# Documentation for Search (search.py)

## class Searcher
A utility for searching the Gene Ontology.<br /><br />Attributes:<br />    go (dict): the GO data.<br />    attributes (dict): gene attributes
### \_\_init\_\_
```py

def __init__(self)

```



Create a object to search through the Gene Ontology.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;None<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;Searcher


### get\_genes
```py

def get_genes(self, terms, evidence_codes=None)

```



Get all of the genes associated with a list of<br />GO idenifiers and some evidence codes.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;terms (List[str]): a list of GO ids<br />&nbsp;&nbsp;&nbsp;&nbsp;evidence_codes (None or List[str]):<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;genes (List[str]): list of genes by ensembl_gene_id


### get\_housekeeping\_genes
```py

def get_housekeeping_genes(self)

```



Get a list of genes that are designated to be "housekeeping genes".<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;None<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;genes (List[str]): list of genes by ensembl_gene_id


### get\_transcription\_factors
```py

def get_transcription_factors(self)

```



Get a list of genes that are designated to be transcription factors<br />in humans.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;None<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;genes (List[str]): list of genes by ensembl_gene_id


### keyword\_search
```py

def keyword_search(self, keywords, fields=['name', 'def'], exact=True, exclude_keywords=None, exclude_ids=None)

```



Search for GO identifiers associated with some keywords.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;keywords (List[str]): words to look for<br />&nbsp;&nbsp;&nbsp;&nbsp;fields (List[str]): fields to look in<br />&nbsp;&nbsp;&nbsp;&nbsp;extact (optional; bool): if true, then this function only returns<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the GO ids where there is a match. if false, then this function<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;will also select all of the child GO terms.<br />&nbsp;&nbsp;&nbsp;&nbsp;exclude_keyworks (optional; List[str]): do NOT include GO categories<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;that contain these keywords<br />&nbsp;&nbsp;&nbsp;&nbsp;exclude_ids (optional; List[str]): do NOT include these GO categories<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;list of GO ids (List[str])


### select\_namespace
```py

def select_namespace(self, namespace)

```



Get all of the GO identifiers associated with a particular namespace.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;namespace (str): 'biological_process' or 'molecular_function'<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;list of go ids (List[str])


### traverse
```py

def traverse(self, term, inclusive=True)

```



Get all of the children of a given GO category.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;term (str): GO id<br />&nbsp;&nbsp;&nbsp;&nbsp;inclusive (optional; bool): include the given term in the final<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;list of GO ids<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;list of GO ids (List[str])




## class chain
chain(*iterables) --> chain object<br /><br />Return a chain object whose .__next__() method returns elements from the<br />first iterable until it is exhausted, then elements from the next<br />iterable, until all of the iterables are exhausted.

