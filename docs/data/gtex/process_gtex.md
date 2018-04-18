# Documentation for Process_Gtex (process_gtex.py)

## functions

### create\_tissue\_stats
```py

def create_tissue_stats()

```



Uses data from the GTEx project to estimate statistics of gene expression<br />across tissues.<br /><br />Expression is measured in Transcripts per Million (TPM).<br />Statistics are saved in csv format in the __file__ directory.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;None<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;None


### hellinger
```py

def hellinger(aves, stds)

```



Computes pairwise Hellinger distances between Gaussian distributions<br />from lists of the means and standard deviations.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;aves (numpy array): list of means (length n)<br />&nbsp;&nbsp;&nbsp;&nbsp;stds (numpy array): list of standard deviations (length n)<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;pairise hellinger distance (numpy array ~ (n, n))


### max\_hellinger
```py

def max_hellinger(aves, stds)

```



Computes maximum Hellinger distance from pairwise comparisions of Gaussian<br />distributions from lists of the means and standard deviations.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;aves (numpy array): list of means (length n)<br />&nbsp;&nbsp;&nbsp;&nbsp;stds (numpy array): list of standard deviations (length n)<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;maximum pairise hellinger distance (float)

