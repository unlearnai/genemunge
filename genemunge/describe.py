import os
import pandas
import numpy
from matplotlib import pyplot as plt


from . import convert
from . import search


datapath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
gtexpath = os.path.join(datapath, 'gtex')


class Describer(object):

    __stats__ = ['mean', 'median', 'std', 'lower_quartile', 'upper_quartile']

    def __init__(self, identifier='symbol'):
        """
        Create an object to grab the information that describes a gene.

        Args:
            identifier (optional; str): the type of gene identifer you will use
                e.g., 'symbol', 'ensembl_gene_id'

        Returns:
            Describer

        """
        if identifier is not 'ensembl_gene_id':
            self.get_ensembl = convert.IDConverter(identifier, 'ensembl_gene_id').convert
        else:
            self.get_ensembl = convert.clean_ensembl_id
        self.get_name = convert.IDConverter('ensembl_gene_id', 'name').convert
        self.get_symbol = convert.IDConverter('ensembl_gene_id', 'symbol').convert
        self.searcher = search.Searcher()
        self.tissue_stats = pandas.HDFStore(os.path.join(gtexpath, 'tissue_stats.h5'))

    def close(self):
        """
        Close the HDF5 store.

        Args:
            None

        Returns:
            None

        """
        self.tissue_stats.close()

    def get_tissue_expression(self, identifier):
        """
        Get statistics describing the expression of a gene across tissues
        in healthy people (from GTEx).

        Args:
            identifier (str)

        Returns:
            pandas.DataFrame

        """
        i = self.get_ensembl(identifier)
        return pandas.concat(
                {k: self.tissue_stats[k].loc[i] for k in self.__stats__}, axis=1)

    def plot_tissue_expression(self, identifier, sortby=None, filename=None):
        """
        Plot the expression of a gene across tissues in health people
        (from GTEx).

        Args:
            identifier (str)
            sortby (optional; str): 'median', 'mean', 'std', 'lower_quartile'
                or 'upper_quartile'. if None, then tissues are alphabetical
            filename (optional; str)

        Returns:
            None

        """
        tissue_stats = self.get_tissue_expression(identifier)

        if sortby is not None:
            stats = tissue_stats.sort_values(sortby)
        else:
            stats = tissue_stats

        fig, ax = plt.subplots(figsize=(10, 4))

        n_box = len(stats)
        boxes = ax.boxplot([[-9, -4, 2, 4, 9],]*n_box)

        min_y, max_y = float('inf'), -float('inf')

        tissues = list(stats.index)
        for i, tissue in enumerate(tissues):
            mean = stats.loc[tissue]['mean']
            med = stats.loc[tissue]['median']
            std = stats.loc[tissue]['std']
            lq = stats.loc[tissue]['lower_quartile']
            uq = stats.loc[tissue]['upper_quartile']

            lower_whisker = mean - 1.5 * std
            upper_whisker = mean + 1.5 * std

            boxes['medians'][i].set_ydata([med, med])
            boxes['caps'][2*i].set_ydata([lower_whisker, lower_whisker])
            boxes['whiskers'][2*i].set_ydata([lower_whisker, lq])
            boxes['caps'][2*i+1].set_ydata([upper_whisker, upper_whisker])
            boxes['whiskers'][2*i + 1].set_ydata([uq, upper_whisker])
            boxes['boxes'][i].set_ydata([lq, lq, uq, uq, lq])

            min_y = min(min_y, lower_whisker)
            max_y = max(max_y, upper_whisker)

        ax.set_ylim([min_y - 0.1 * abs(min_y), max_y + 0.1 * abs(max_y)])
        plt.xticks(numpy.arange(len(tissues)) + 1, tissues, rotation='vertical')
        ax.set_title(identifier)
        ax.set_ylabel('TPM')
        plt.show(fig)
        if filename is not None:
            fig.savefig(filename)

    def _get_go_from_ensemble(self, ensembl):
        """
        Get a list of GO categories associated with the given ensembl id.

        Args:
            ensemble (str)

        Returns:
            GO identifiers (List[str])

        """
        terms = []
        for term in self.searcher.go:
            has_gene = any(ensembl in self.searcher.go[term]["genes"][code]
                        for code in self.searcher.go[term]["genes"])
            if has_gene:
                terms += [term]
        return list(set(terms))

    def get_gene_info(self, identifier):
        """
        Get some information about a gene such as:
            ensemble_gene_id
            gene symbol
            name
            associated categories in the Gene Ontology

        Args:
            identifier (str)

        Returns:
            None

        """
        gene_info = {}
        gene_info['ensembl'] = self.get_ensembl(identifier)
        gene_info['symbol'] = self.get_symbol(gene_info['ensembl'])
        gene_info['name'] = self.get_name(gene_info['ensembl'])
        gene_info['ontology'] = {i: self.searcher.go[i]['name']
            for i in self._get_go_from_ensemble(gene_info['ensembl'])}
        return gene_info
