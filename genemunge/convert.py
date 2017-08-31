import os, pandas
from typing import List

FILEPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
FILENAME = os.path.join(FILEPATH, 'hgnc_complete_set.txt')

class IDConverter(object):
    """Convert between gene identifiers."""

    potential_ids = ['hgnc_id',
                     'symbol',
                     'name',
                     'locus_group',
                     'locus_type',
                     'status',
                     'location',
                     'location_sortable',
                     'alias_symbol',
                     'alias_name',
                     'prev_symbol',
                     'prev_name',
                     'gene_family',
                     'gene_family_id',
                     'date_approved_reserved',
                     'date_symbol_changed',
                     'date_name_changed',
                     'date_modified',
                     'entrez_id',
                     'ensembl_gene_id',
                     'vega_id',
                     'ucsc_id',
                     'ena',
                     'refseq_accession',
                     'ccds_id',
                     'uniprot_ids',
                     'pubmed_id',
                     'mgd_id',
                     'rgd_id',
                     'lsdb',
                     'cosmic',
                     'omim_id',
                     'mirbase',
                     'homeodb',
                     'snornabase',
                     'bioparadigms_slc',
                     'orphanet',
                     'pseudogene.org',
                     'horde_id',
                     'merops',
                     'imgt',
                     'iuphar',
                     'kznf_gene_catalog',
                     'mamit-trnadb',
                     'cd',
                     'lncrnadb',
                     'enzyme_id',
                     'intermediate_filament_db',
                     'rna_central_ids']

    def __init__(self, source_id: str, target_id: str):
        """
        Create IDConverter.

        Args:
            source_id (str): the id type to convert
            target_id (str): the desired id type

        Returns:
            IDConverter

        """
        assert source_id in self.potential_ids, \
        "unknown source_id type"
        assert target_id in self.potential_ids, \
        "unknown target_id type"

        self.source = source_id
        self.target = target_id
        self.conversion_table = pandas.read_table(FILENAME,
                        usecols=[source_id, target_id], dtype=str)
        self._clean_conversion_table()

    def _clean_conversion_table(self):
        # drop NaN from the source column
        self.conversion_table.dropna(subset=[self.source], inplace=True)
        # drop duplicates from the source column
        self.conversion_table.drop_duplicates(subset=[self.source], inplace=True, keep=False)
        # set the index to the source column
        self.conversion_table.set_index(self.source, inplace=True)

    def convert(self, ids: List) -> List:
        """
        Convert an list of gene identifiers.

        Args:
            ids (List[str]): list of gene identifiers to convert

        Returns:
            List[str]: list of converted gene identifiers
        """
        return list(self.conversion_table.loc[ids][self.target])
