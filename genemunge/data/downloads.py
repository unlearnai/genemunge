import urllib.request
import os
import sys


FILEPATH = os.path.dirname(os.path.abspath(__file__))
LINK_TO_GO = "http://purl.obolibrary.org/obo/go/go-basic.obo"
GONAME = "go-basic.obo"
LINK_TO_HGNC = "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt"
HGNCNAME = "hgnc_complete_set.txt"
LINK_TO_ANNOTATIONS = "http://geneontology.org/gene-associations/goa_human.gaf.gz"
ANNOTATIONSNAME = "goa_human.gaf.gz"
LINK_TO_HOUSEKEEPING = "http://www.tau.ac.il/~elieis/HKG/HK_genes.txt"
HOUSEKEEPINGNAME = "HK_genes.txt"
LINK_TO_TRANSCRIPTION_FACTORS = "http://www.tfcheckpoint.org/data/TFCheckpoint_download_180515.txt"
TFNAME = "TFCheckpoint_download_180515.txt"


def download_progress_indicator(count, blockSize, totalSize):
    """
    Print out a progress indicator to the sreen.

    Args:
        count (int): number of bits downloaded
        blockSize (int): number of bits per block
        totalSize (int): number of bits in the file

    Returns:
        None

    """
    percent = int(count*blockSize*100/totalSize)
    sys.stdout.write("\rProgress... %d%%" % percent)
    sys.stdout.flush()


def download_go(force=False):
    """
    Fetch the gene ontology file.

    Args:
        force (optional; bool): set to 'True' to overwrite existing files

    Returns:
        None

    """
    output_file = os.path.join(FILEPATH, GONAME)
    output_exists = os.path.exists(output_file)
    if not output_exists or force:
        print("downloading {}".format(GONAME))
        urllib.request.urlretrieve(LINK_TO_GO, output_file,
                               reporthook=download_progress_indicator)
        print("\ndone")


def download_annotations(force=False):
    """
    Fetch the gene ontology annotations file.

    Args:
        force (optional; bool): set to 'True' to overwrite existing files

    Returns:
        None

    """
    output_file = os.path.join(FILEPATH, ANNOTATIONSNAME)
    output_exists = os.path.exists(output_file)
    if not output_exists or force:
        print("downloading {}".format(ANNOTATIONSNAME))
        urllib.request.urlretrieve(LINK_TO_ANNOTATIONS, output_file,
                               reporthook=download_progress_indicator)
        print("\ndone")


def download_hgnc(force=False):
    """
    Fetch the HGNC gene naming table.

    Args:
        force (optional; bool): set to 'True' to overwrite existing files

    Returns:
        None

    """
    output_file = os.path.join(FILEPATH, HGNCNAME)
    output_exists = os.path.exists(output_file)
    if not output_exists or force:
        print("\ndownloading {}".format(HGNCNAME))
        urllib.request.urlretrieve(LINK_TO_HGNC, output_file,
                               reporthook=download_progress_indicator)
        print("\ndone")


def download_housekeeping(force=False):
    """
    Fetch the list of housekeeping genes.

    Args:
        force (optional; bool): set to 'True' to overwrite existing files

    Returns:
        None

    """
    output_file = os.path.join(FILEPATH, HOUSEKEEPINGNAME)
    output_exists = os.path.exists(output_file)
    if not output_exists or force:
        print("\ndownloading {}".format(HOUSEKEEPINGNAME))
        urllib.request.urlretrieve(LINK_TO_HOUSEKEEPING, output_file)
        print("done")


def download_transcription_factors(force=False):
    """
    Fetch the list of transcription factors.

    Args:
        force (optional; bool): set to 'True' to overwrite existing files

    Returns:
        None

    """
    output_file = os.path.join(FILEPATH, TFNAME)
    output_exists = os.path.exists(output_file)
    if not output_exists or force:
        print("\ndownloading {}".format(TFNAME))
        urllib.request.urlretrieve(LINK_TO_TRANSCRIPTION_FACTORS, output_file)
        print("done")


def download_everything(force=False):
    """
    Download all of the files.

    Args:
        force (optional; bool): set to 'True' to overwrite existing files

    Returns:
        None

    """
    download_go(force)
    download_annotations(force)
    download_hgnc(force)
    download_housekeeping(force)
    download_transcription_factors(force)


if __name__ == "__main__":
    download_everything(force=False)
