import os
import contextlib

filepath = os.path.dirname(os.path.abspath(__file__))


def remove_installed_data_files():
    """
    Removes some of the files that are downloaded during install.

    Args:
        None

    Returns:
        None

    """
    # remove the GO file
    with contextlib.suppress(FileNotFoundError):
        os.remove(os.path.join(filepath, 'go-basic.obo'))

    # remove the GO annotations
    with contextlib.suppress(FileNotFoundError):
        os.remove(os.path.join(filepath, 'goa_human.gaf.gz'))

    # remove the housekeeping genes file
    with contextlib.suppress(FileNotFoundError):
        os.remove(os.path.join(filepath, 'HK_genes.txt'))

    # remove the transcription factors genes file
    with contextlib.suppress(FileNotFoundError):
        os.remove(os.path.join(filepath, 'TFCheckpoint_download_180515.txt'))

    # remove the GTEx Rdata file
    with contextlib.suppress(FileNotFoundError):
        os.remove(os.path.join(filepath, 'gtex', 'rse_gene.Rdata'))

    # remove the GTEx expression file
    with contextlib.suppress(FileNotFoundError):
        os.remove(os.path.join(filepath, 'gtex', 'expression_data.csv'))