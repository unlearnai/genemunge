# Documentation for Convert (convert.py)

## class IDConverter
Convert between gene identifiers.
### \_\_init\_\_
```py

def __init__(self, source_id: str, target_id: str)

```



Create IDConverter.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;source_id (str): the id type to convert<br />&nbsp;&nbsp;&nbsp;&nbsp;target_id (str): the desired id type<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;IDConverter


### convert
```py

def convert(self, identifier)

```



Convert a gene identifier.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;identifier (str): gene identifier to convert<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;str: converted gene identifier


### convert\_list
```py

def convert_list(self, ids: List) -> List

```



Convert a list of gene identifiers.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;ids (List[str]): list of gene identifiers to convert<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;List[str]: list of converted gene identifiers




## class List
list() -> new empty list<br />list(iterable) -> new list initialized from iterable's items


## functions

### clean\_ensembl\_id
```py

def clean_ensembl_id(identifier)

```



Formats an ensembl gene identifier to drop the version number.<br /><br />E.g., ENSG00000002822.15 -> ENSG00000002822<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;identifier (str)<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;identifier (str)


### clean\_ensembl\_ids
```py

def clean_ensembl_ids(identifiers)

```



Formats ensembl gene identifiers to drop the version number<br /><br />E.g., ENSG00000002822.15 -> ENSG00000002822<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;identifiers (List[str])<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;identifiers (List[str])

