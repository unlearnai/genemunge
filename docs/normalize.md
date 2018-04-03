# Documentation for Normalize (normalize.py)

## class RemoveUnwantedVariation
### \_\_init\_\_
```py

def __init__(self, alpha=None)

```



Perform the 2-step Remove Unwanted Variation (RUV-2) algorithm<br />defined in:<br /><br />"Correcting gene expression data when neither the unwanted variation nor the<br />factor of interest are observed."<br />Biostatistics 17.1 (2015): 16-28.<br />Laurent Jacob, Johann A. Gagnon-Bartsch, and Terence P. Speed.<br /><br />The algorithm is modified slightly so that batch correction can be<br />applied out-of-sample.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;alpha (optional; numpy array ~ (num_singular_values, num_genes))<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;RemoveUnwantedVariation


### fit
```py

def fit(self, data, hk_genes)

```



Perform a singular value decomposition of the housekeeping genes to<br />fit the transform.<br /><br />Suppose that we measure data on the expression of N genes in M samples<br />and store these (after CLR transformation) in a matrix Y \in R^{M, N}.<br />We consider a linear model Y = X B + W A + noise where<br />&nbsp;&nbsp;&nbsp;&nbsp;X \in R^{M, Q} are some unobserved, but biologically interesting, factors<br />&nbsp;&nbsp;&nbsp;&nbsp;B \in R^{Q, N} describes how the genes are coupled to the interesting factors<br />&nbsp;&nbsp;&nbsp;&nbsp;W \in R^{M, K} are some unobserved and uninteresting factors<br />&nbsp;&nbsp;&nbsp;&nbsp;A \in R^{K, N} describes how the genes are coupled to the uninteresting factors<br /><br />We assume that there are some housekeeping genes Y_c for which we are<br />sure that B_c = 0. That is, the housekeeping genes are not coupled to<br />any biologically interesting factors. Therefore, we have Y_c = W A_c + noise.<br />Let Y_c = U L V^{T} be the singular value decomposition of Y_c. Then,<br />we can estiamte W = U L.<br /><br />Now, if we fix W and assume that X B = 0 for all genes then we can<br />estimate A = (W W^{T})^{-1} W^{T} Y. This matrix stores K patterns of<br />variation that are usually not biologically interesting.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;data (pandas.DataFrame ~ (num_samples, num_genes)): clr transformed<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;expression data<br />&nbsp;&nbsp;&nbsp;&nbsp;hk_genes (List[str]): list of housekeeping genes<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;None


### fit\_transform
```py

def fit_transform(self, data, hk_genes)

```



Perform the 2-step Remove Unwanted Variation (RUV-2) algorithm.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;data (pandas.DataFrame ~ (num_samples, num_genes)): clr transformed<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;expression data<br />&nbsp;&nbsp;&nbsp;&nbsp;hk_genes (List[str]): list of housekeeping genes<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;batch corrected data (pandas.DataFrame ~ (num_samples, num_genes))


### save
```py

def save(self, filename, overwrite_existing=False)

```



Save the RUV object to filename.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;filename (string): absolute path to save file<br />&nbsp;&nbsp;&nbsp;&nbsp;overwrite_existing (bool): whether or not to overwrite existing file<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;None


### transform
```py

def transform(self, data)

```



Perform the 2-step Remove Unwanted Variation (RUV-2) algorithm.<br /><br />The `fit` method estimates the matrix<br />&nbsp;&nbsp;&nbsp;&nbsp;A \in R^{K, N} which describes how the genes are coupled to the<br />&nbsp;&nbsp;&nbsp;&nbsp;uninteresting factors<br /><br />We can estimate the activity of these factors from a new dataset&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ilde{Y}<br />by computing&nbsp;&nbsp;&nbsp;&nbsp;ilde{W} =&nbsp;&nbsp;&nbsp;&nbsp;   ilde{Y} A^{T} (A A^{T})^{-1} using the right<br />pseudoinverse of A.<br /><br />Finally, we can subtract&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ilde{W} A from the data by computing<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ilde{Y} -&nbsp;&nbsp;&nbsp;&nbsp;   ilde{Y} A^{T} (A A^{T})^{-1} A.<br /><br />Essentially, we are removing the components of the data that project<br />onto the pre-defined axes of uninteresting variation.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;data (pandas.DataFrame ~ (num_samples, num_genes)): clr transformed<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;expression data<br />&nbsp;&nbsp;&nbsp;&nbsp;hk_genes (List[str]): list of housekeeping genes<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;batch corrected data (pandas.DataFrame ~ (num_samples, num_genes))




## class Normalizer
### \_\_init\_\_
```py

def __init__(self, identifier='symbol')

```



Tools to normalize expression data and transform into TPM.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;identifer (str)<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;Normalizer


### clr\_from\_tpm
```py

def clr_from_tpm(self, data, gene_list=None, imputer=<function do_nothing>)

```



Compute the centered log ratio transform of data in TPM format.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;data (pandas.DataFrame ~ (num_samples, num_genes))<br />&nbsp;&nbsp;&nbsp;&nbsp;gene_list (optional; List[str]): a list of gene ids<br />&nbsp;&nbsp;&nbsp;&nbsp;imputer (optional; callable)<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;pandas.DataFrame


### tpm\_from\_clr
```py

def tpm_from_clr(self, data, gene_list=None)

```



Compute data in TPM format from centered log ratio transformed data.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;data (pandas.DataFrame ~ (num_samples, num_genes))<br />&nbsp;&nbsp;&nbsp;&nbsp;gene_list (optional; List[str]): a list of gene ids<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;pandas.DataFrame


### tpm\_from\_counts
```py

def tpm_from_counts(self, data, gene_list=None)

```



Transform data from counts to TPM.<br />Any genes not in the gene_lengths index is removed,<br />&nbsp;&nbsp;&nbsp;&nbsp;as the gene length is not known.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;data (pandas.DataFrame ~ (num_samples, num_genes))<br />&nbsp;&nbsp;&nbsp;&nbsp;gene_list (optional; List[str]): a list of gene ids<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;pandas.DataFrame


### tpm\_from\_rpkm
```py

def tpm_from_rpkm(self, data, gene_list=None)

```



Transform data from RPKM to TPM.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;data (pandas.DataFrame ~ (num_samples, num_genes))<br />&nbsp;&nbsp;&nbsp;&nbsp;gene_list (optional; List[str]): a list of gene ids<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;pandas.DataFrame


### tpm\_from\_subset
```py

def tpm_from_subset(self, data, gene_list=None)

```



Renormalize a subset of genes already in TPM.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;data (pandas.DataFrame ~ (num_samples, num_genes))<br />&nbsp;&nbsp;&nbsp;&nbsp;gene_list (optional; List[str]): a list of gene ids<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;pandas.DataFrame




## class Path
PurePath represents a filesystem path and offers operations which<br />don't imply any actual filesystem I/O.  Depending on your system,<br />instantiating a PurePath will return either a PurePosixPath or a<br />PureWindowsPath object.  You can also instantiate either of these classes<br />directly, regardless of your system.


## functions

### deduplicate
```py

def deduplicate(data)

```



Adds the values from any duplicated genes.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;data (pandas.DataFrame ~ (num_samples, num_genes))<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;pandas.DataFrame


### do\_nothing
```py

def do_nothing(data)

```



A function that does nothing.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;Anything<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;Anything


### impute
```py

def impute(data, scale=0.5)

```



Replace any zeros in each row with a fraction of the smallest non-zero<br />value in the corresponding row.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;data (pandas.DataFrame ~ (num_samples, num_genes))<br />&nbsp;&nbsp;&nbsp;&nbsp;scale (optional; float)<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;imputed data (pandas.DataFrame ~ (num_samples, num_genes))

