import os
import pandas
import json

def create_attributes_file():
    """
    Creates a json file with lists (as gene symbols) of housekeeping genes
    and transcription factors.

    Args:
        None

    Returns:
        None

    """
    filedir = os.path.dirname(os.path.abspath(__file__))
    hk_symbols = list(pandas.read_csv(os.path.join(filedir, "HK_genes.txt"),
                                      sep="\t", header=None)[0])
    tfs = pandas.read_csv(os.path.join(filedir, "TFCheckpoint_download_180515.txt"),
                          sep='\t')
    tf_symbols = list(tfs[tfs['TFClass_human'] == 'TFclass']['gene_symbol'])

    with open(os.path.join(filedir, 'gene_attributes.json'), 'w') as out:
        json.dump({'housekeeping_genes': hk_symbols,
                   'transcription_factors': tf_symbols}, out)