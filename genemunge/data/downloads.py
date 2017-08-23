import urllib, os, sys

FILEPATH = os.path.dirname(os.path.abspath(__file__))
LINK_TO_GO = "http://purl.obolibrary.org/obo/go/go-basic.obo"
GONAME = "go-basic.obo"
LINK_TO_HGNC = "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt"
HGNCNAME = "hgnc_complete_set.txt"
LINK_TO_ANNOTATIONS = "http://geneontology.org/gene-associations/goa_human.gaf.gz"
ANNOTATIONSNAME = "goa_human.gaf.gz"

def download_progress_indicator(count, blockSize, totalSize):
      percent = int(count*blockSize*100/totalSize)
      sys.stdout.write("\rProgress... %d%%" % percent)
      sys.stdout.flush()

def download_go(force=False):
    output_file = os.path.join(FILEPATH, GONAME)
    output_exists = os.path.exists(output_file)
    if not output_exists or force:
        print("downloading {}".format(GONAME))
        urllib.request.urlretrieve(LINK_TO_GO, output_file,
                               reporthook=download_progress_indicator)
        print("\ndone")

def download_annotations(force=False):
    output_file = os.path.join(FILEPATH, ANNOTATIONSNAME)
    output_exists = os.path.exists(output_file)
    if not output_exists or force:
        print("downloading {}".format(ANNOTATIONSNAME))
        urllib.request.urlretrieve(LINK_TO_ANNOTATIONS, output_file,
                               reporthook=download_progress_indicator)
        print("\ndone")

def download_hgnc(force=False):
    output_file = os.path.join(FILEPATH, HGNCNAME)
    output_exists = os.path.exists(output_file)
    if not output_exists or force:
        print("\ndownloading {}".format(HGNCNAME))
        urllib.request.urlretrieve(LINK_TO_HGNC, output_file,
                               reporthook=download_progress_indicator)
        print("\ndone")

if __name__ == "__main__":
    download_go()
    download_annotations()
    download_hgnc()
