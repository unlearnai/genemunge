# Documentation for Normalize (normalize.py)

## class Normalizer
### \_\_init\_\_
```py

def __init__(self, identifier='symbol')

```



Tools to normalize expression data and transform into TPM.<br /><br />Args:<br /> ~ None<br /><br />Returns:<br /> ~ None


### tpm\_from\_counts
```py

def tpm_from_counts(self, data, gene_list=None)

```



Transform data from counts to TPM.<br /><br />Args:<br /> ~ data (pandas.DataFrame ~ (num_samples, num_genes))<br /> ~ gene_list (optional; List[str]): a list of gene ids<br /><br />Returns:<br /> ~ pandas.DataFrame


### tpm\_from\_rpkm
```py

def tpm_from_rpkm(self, data, gene_list=None)

```



Transform data from RPKM to TPM.<br /><br />Args:<br /> ~ data (pandas.DataFrame ~ (num_samples, num_genes))<br /> ~ gene_list (optional; List[str]): a list of gene ids<br /><br />Returns:<br /> ~ pandas.DataFrame


### tpm\_from\_subset
```py

def tpm_from_subset(self, data, gene_list=None)

```



Renormalize a subset of genes.<br /><br />Args:<br /> ~ data (pandas.DataFrame ~ (num_samples, num_genes))<br /> ~ gene_list (optional; List[str]): a list of gene ids<br /><br />Returns:<br /> ~ pandas.DataFrame




## functions

### deduplicate
```py

def deduplicate(data)

```



Adds the values from any duplicated genes.<br /><br />Args:<br /> ~ data (pandas.DataFrame ~ (num_samples, num_genes))<br /><br />Returns:<br /> ~ pandas.DataFrame


### impute
```py

def impute(data, scale=0.5)

```



Replace any zeros in each row with a fraction of the smallest non-zero<br />value in the corresponding row.<br /><br />Args:<br /> ~ data (pandas.DataFrame ~ (num_samples, num_genes))<br /> ~ scale (optional; float)<br /><br />Returns:<br /> ~ imputed data (pandas.DataFrame ~ (num_samples, num_genes))

