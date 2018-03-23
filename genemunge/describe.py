import os
import pandas

from . import convert

class Describer(object):

    def __init__(self, identifier='symbol'):
        """

        """
        c = convert.IDConverter('ensembl_gene_id', identifier)

