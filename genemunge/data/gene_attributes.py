import os
import pandas
import json

from .. import convert

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
    # read in the list of housekeeping genes
    hk_symbols = list(pandas.read_csv(os.path.join(filedir, "HK_genes.txt"),
                                      sep="\t", header=None)[0])
    hk_symbols = list(map(lambda x: x.strip(), hk_symbols))
    # read in the list of transcription factors
    tfs = pandas.read_csv(os.path.join(filedir, "TFCheckpoint_download_180515.txt"),
                          sep='\t')
    tf_symbols = list(tfs[tfs['TFClass_human'] == 'TFclass']['gene_symbol'])
    # convert gene symbols to uniprot ids to be consistent with GO
    converter = convert.IDConverter("symbol", "uniprot_ids")
    hk_uniprot = converter.convert_list(hk_symbols)
    tf_uniprot = converter.convert_list(tf_symbols)
    # write to a file
    with open(os.path.join(filedir, 'gene_attributes.json'), 'w') as out:
        json.dump({'housekeeping_genes': hk_uniprot,
                   'transcription_factors': tf_uniprot}, out)

