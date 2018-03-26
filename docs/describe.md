# Documentation for Describe (describe.py)

## class Describer
### \_\_init\_\_
```py

def __init__(self, identifier='symbol')

```



Create an object to grab the information that describes a gene.<br /><br />Args:<br /> ~ identifier (optional; str): the type of gene identifer you will use<br /> ~  ~ e.g., 'symbol', 'ensembl_gene_id'<br /><br />Returns:<br /> ~ Describer


### close
```py

def close(self)

```



Close the HDF5 store.<br /><br />Args:<br /> ~ None<br /><br />Returns:<br /> ~ None


### get\_gene\_info
```py

def get_gene_info(self, identifier)

```



Get some information about a gene such as:<br /> ~ ensemble_gene_id<br /> ~ gene symbol<br /> ~ name<br /> ~ associated categories in the Gene Ontology<br /><br />Args:<br /> ~ identifier (str)<br /><br />Returns:<br /> ~ None


### get\_tissue\_expression
```py

def get_tissue_expression(self, identifier)

```



Get statistics describing the expression of a gene across tissues<br />in healthy people (from GTEx).<br /><br />Args:<br /> ~ identifier (str)<br /><br />Returns:<br /> ~ pandas.DataFrame


### plot\_tissue\_expression
```py

def plot_tissue_expression(self, identifier, sortby=None, filename=None)

```



Plot the expression of a gene across tissues in health people<br />(from GTEx).<br /><br />Args:<br /> ~ identifier (str)<br /> ~ sortby (optional; str): 'median', 'mean', 'std', 'lower_quartile'<br /> ~  ~ or 'upper_quartile'. if None, then tissues are alphabetical<br /> ~ filename (optional; str)<br /><br />Returns:<br /> ~ None



