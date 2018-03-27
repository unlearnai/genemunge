# Documentation for Convert (convert.py)

## class IDConverter
Convert between gene identifiers.
### \_\_init\_\_
```py

def __init__(self, source_id: str, target_id: str)

```



Create IDConverter.<br /><br />Args:<br />	source_id (str): the id type to convert<br />	target_id (str): the desired id type<br /><br />Returns:<br />	IDConverter


### convert
```py

def convert(self, identifier)

```



Convert a of gene identifier.<br /><br />Args:<br />	id (str): gene identifier to convert<br /><br />Returns:<br />	str: converted gene identifier


### convert\_list
```py

def convert_list(self, ids: List) -> List

```



Convert an list of gene identifiers.<br /><br />Args:<br />	ids (List[str]): list of gene identifiers to convert<br /><br />Returns:<br />	List[str]: list of converted gene identifiers




## class List
list() -> new empty list<br />list(iterable) -> new list initialized from iterable's items


## functions

### clean\_ensembl\_id
```py

def clean_ensembl_id(identifier)

```



Formats ensembl gene identifiers to drop the version number.<br /><br />E.g., ENSG00000002822.15 -> ENSG00000002822<br /><br />Args:<br />	identifier (str)<br /><br />Returns:<br />	identifier (str)


### clean\_ensembl\_ids
```py

def clean_ensembl_ids(identifiers)

```



Formats ensembl gene identifiers to drop the version number<br /><br />E.g., ENSG00000002822.15 -> ENSG00000002822<br /><br />Args:<br />	identifiers (List[str])<br /><br />Returns:<br />	identifier (List[str])

