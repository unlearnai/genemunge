# Documentation for Downloads (downloads.py)

## functions

### download\_annotations
```py

def download_annotations(force=False)

```



Fetch the gene ontology annotations file.<br /><br />Args:<br /> ~ force (optional; bool): set to 'True' to overwrite existing files<br /><br />Returns:<br /> ~ None


### download\_everything
```py

def download_everything(force=False)

```



Download all of the files.<br /><br />Args:<br /> ~ force (optional; bool): set to 'True' to overwrite existing files<br /><br />Returns:<br /> ~ None


### download\_go
```py

def download_go(force=False)

```



Fetch the gene ontology file.<br /><br />Args:<br /> ~ force (optional; bool): set to 'True' to overwrite existing files<br /><br />Returns:<br /> ~ None


### download\_hgnc
```py

def download_hgnc(force=False)

```



Fetch the HGNC gene naming table.<br /><br />Args:<br /> ~ force (optional; bool): set to 'True' to overwrite existing files<br /><br />Returns:<br /> ~ None


### download\_housekeeping
```py

def download_housekeeping(force=False)

```



Fetch the list of housekeeping genes.<br /><br />Args:<br /> ~ force (optional; bool): set to 'True' to overwrite existing files<br /><br />Returns:<br /> ~ None


### download\_progress\_indicator
```py

def download_progress_indicator(count, blockSize, totalSize)

```



Print out a progress indicator to the sreen.<br /><br />Args:<br /> ~ count (int): number of bits downloaded<br /> ~ blockSize (int): number of bits per block<br /> ~ totalSize (int): number of bits in the file<br /><br />Returns:<br /> ~ None


### download\_transcription\_factors
```py

def download_transcription_factors(force=False)

```



Fetch the list of transcription factors.<br /><br />Args:<br /> ~ force (optional; bool): set to 'True' to overwrite existing files<br /><br />Returns:<br /> ~ None

