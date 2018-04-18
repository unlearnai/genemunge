# Documentation for Normalize (normalize.py)

## class RemoveUnwantedVariation
The RUV-2 algorithm.<br /><br />Attributes:<br />    hk_genes (List[str]): a list of housekeeping gene names used in fitting.<br />    means (pandas.Series): the means of each gene from the training data.<br />    U (optional; numpy array ~ (num_training_samples, num_factors)):<br />        left eigenvectors from SVD of housekeeping genes in training set<br />    L (optional; numpy array ~ (num_factors,))<br />        eigenvalues from SVD of housekeeping genes in training set<br />    Vt (optional; numpy_array ~ (num_factors, num_hk_genes)<br />        right eigenvectors from SVD of housekeeping genes in training set
### \_\_init\_\_
```py

def __init__(self, center=True, hk_genes=None, means=None, U=None, L=None, Vt=None)

```



Perform the 2-step Remove Unwanted Variation (RUV-2) algorithm<br />defined in:<br /><br />"Correcting gene expression data when neither the unwanted variation nor the<br />factor of interest are observed."<br />Biostatistics 17.1 (2015): 16-28.<br />Laurent Jacob, Johann A. Gagnon-Bartsch, and Terence P. Speed.<br /><br />The algorithm is modified slightly so that batch correction can be<br />applied out-of-sample.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;center (optional; bool): whether to center the gene means in the fit.<br />&nbsp;&nbsp;&nbsp;&nbsp;hk_genes (optional; List[str]): list of housekeeping genes<br />&nbsp;&nbsp;&nbsp;&nbsp;means (optional; numpy array ~ (num_genes,))<br />&nbsp;&nbsp;&nbsp;&nbsp;U (optional; numpy array ~ (num_training_samples, num_factors))<br />&nbsp;&nbsp;&nbsp;&nbsp;L (optional; numpy array ~ (num_factors,))<br />&nbsp;&nbsp;&nbsp;&nbsp;Vt (optional; numpy_array ~ (num_factors, num_hk_genes))<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;RemoveUnwantedVariation


### fit
```py

def fit(self, data, hk_genes, variance_cutoff=0.9)

```



Perform a singular value decomposition of the housekeeping genes to<br />fit the transform.<br /><br />Suppose that we measure data on the expression of N genes in M samples<br />and store these (after CLR transformation) in a matrix Y \in R^{M, N}.<br />We consider a linear model Y = X B + W A + noise where<br />&nbsp;&nbsp;&nbsp;&nbsp;X \in R^{M, Q} are some unobserved, but biologically interesting, factors<br />&nbsp;&nbsp;&nbsp;&nbsp;B \in R^{Q, N} describes how the genes are coupled to the interesting factors<br />&nbsp;&nbsp;&nbsp;&nbsp;W \in R^{M, K} are some unobserved and uninteresting factors<br />&nbsp;&nbsp;&nbsp;&nbsp;A \in R^{K, N} describes how the genes are coupled to the uninteresting factors<br /><br />We assume that there are some housekeeping genes Y_c for which we are<br />sure that B_c = 0. That is, the housekeeping genes are not coupled to<br />any biologically interesting factors. Therefore, we have Y_c = W A_c + noise.<br />Let Y_c = U L V^{T} be the singular value decomposition of Y_c. Then,<br />we can estiamte W = U L.  Additionally, A_c = V^{T}.<br /><br />Now, if we fix W and assume that X B = 0 for all genes then we can<br />estimate A = W^+ Y = (W W^{T})^{-1} W^{T} Y.<br />This matrix stores K patterns of variation that are<br />usually not biologically interesting.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;data (pandas.DataFrame ~ (num_samples, num_genes)): clr transformed<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;expression data<br />&nbsp;&nbsp;&nbsp;&nbsp;hk_genes (List[str]): list of housekeeping genes<br />&nbsp;&nbsp;&nbsp;&nbsp;variance_cutoff (float): the cumulative variance cutoff on SVD<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;eigenvalues of Y_c (the variance fraction of the factors).<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;None


### fit\_transform
```py

def fit_transform(self, data, hk_genes, penalty=0, variance_cutoff=0.9)

```



Perform the 2-step Remove Unwanted Variation (RUV-2) algorithm.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;data (pandas.DataFrame ~ (num_samples, num_genes)): clr transformed<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;expression data<br />&nbsp;&nbsp;&nbsp;&nbsp;hk_genes (List[str]): list of housekeeping genes<br />&nbsp;&nbsp;&nbsp;&nbsp;penalty (float): regularization on the regression step<br />&nbsp;&nbsp;&nbsp;&nbsp;variance_cutoff (float): the cumulative variance cutoff on SVD<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;eigenvalues of Y_c.<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;batch corrected data (pandas.DataFrame ~ (num_samples, num_genes))


### save
```py

def save(self, filename, overwrite_existing=False)

```



Save the RUV object to filename.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;filename (string): absolute path to save file<br />&nbsp;&nbsp;&nbsp;&nbsp;overwrite_existing (bool): whether or not to overwrite existing file<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;None


### transform
```py

def transform(self, data, penalty=0)

```



Perform the 2-step Remove Unwanted Variation (RUV-2) algorithm.<br /><br />The `fit` method estimates the matrix<br />&nbsp;&nbsp;&nbsp;&nbsp;A \in R^{K, N} which describes how the genes are coupled to the<br />&nbsp;&nbsp;&nbsp;&nbsp;uninteresting factors<br /><br />We can estimate the activity of these factors from a new dataset&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ilde{Y}<br />by using the housekeeping genes on this new dataset and computing<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ilde{W} =&nbsp;&nbsp;&nbsp;&nbsp;   ilde{Y}_c A_c^{+}.  Since A_c = V^{T} from the SVD,<br />the right pseudoinverse A_c^{+} = A_c^{T}.<br /><br />Finally, we can subtract&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ilde{W} A from the data,<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ilde{Y} -&nbsp;&nbsp;&nbsp;&nbsp;   ilde{W} A.<br /><br />Essentially, we are removing the components of the data that project<br />onto the pre-defined axes of uninteresting variation.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;data (pandas.DataFrame ~ (num_samples, num_genes)): clr transformed<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;expression data<br />&nbsp;&nbsp;&nbsp;&nbsp;penalty (float): regularization on the regression step<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;batch corrected data (pandas.DataFrame ~ (num_samples, num_genes))




## class Normalizer
Tools to change units of expression data, primarily to convert to TPM.<br /><br />Attributes:<br />    gene_lengths (DataFrame): bp lengths for genes.
### \_\_init\_\_
```py

def __init__(self, identifier='symbol')

```



Tools to normalize expression data and transform into TPM.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;identifier (str)<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;Normalizer


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

