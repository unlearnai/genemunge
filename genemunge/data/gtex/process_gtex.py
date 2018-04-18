import os
import pandas
import numpy

from ... import normalize
from ... import convert


def hellinger(aves, stds):
    """
    Computes pairwise Hellinger distances between Gaussian distributions
    from lists of the means and standard deviations.

    Args:
        aves (numpy array): list of means (length n)
        stds (numpy array): list of standard deviations (length n)

    Returns:
        pairise hellinger distance (numpy array ~ (n, n))

    """
    epsilon = 1e-6
    variances = (epsilon + stds) ** 2
    std_probs = numpy.multiply.outer(epsilon + stds, epsilon + stds)
    var_sums = numpy.add.outer(variances, variances)
    mean_diffs = numpy.subtract.outer(aves, aves)
    return 1 - numpy.exp(-mean_diffs ** 2 / (4 * var_sums)) * \
                numpy.sqrt(2 * std_probs / var_sums)


def max_hellinger(aves, stds):
    """
    Computes maximum Hellinger distance from pairwise comparisions of Gaussian
    distributions from lists of the means and standard deviations.

    Args:
        aves (numpy array): list of means (length n)
        stds (numpy array): list of standard deviations (length n)

    Returns:
        maximum pairise hellinger distance (float)

    """
    return numpy.max(hellinger(aves, stds))


def create_tissue_stats():
    """
    Uses data from the GTEx project to estimate statistics of gene expression
    across tissues.

    Expression is measured in Transcripts per Million (TPM).
    Statistics are saved in csv format in the __file__ directory.

    Args:
        None

    Returns:
        None

    """
    filepath = os.path.dirname(os.path.abspath(__file__))

    samples = pandas.read_csv(os.path.join(filepath, 'SRP012682.tsv'), sep='\t')
    expression = pandas.read_csv(os.path.join(filepath, 'expression_data.csv'), sep='\t')

    mean = pandas.DataFrame()
    median = pandas.DataFrame()
    std = pandas.DataFrame()
    lower_quartile = pandas.DataFrame()
    upper_quartile = pandas.DataFrame()
    fraction_zero = pandas.DataFrame()

    tissues = samples.groupby('smts').groups
    norm = normalize.Normalizer('ensembl_gene_id')
    for t in tissues:
        print('processing tissue: {}'.format(t))
        # select the sample ids corresponding to the tissue
        index = list(samples['run'].loc[tissues[t]].values)
        tissue_expression = (expression[index].T)
        # replace missing data with 0
        tissue_expression.fillna(0.0, inplace=True)
        # clean the ensemble ids and add up the duplicates
        tissue_expression.columns = convert.clean_ensembl_ids(tissue_expression.columns)
        tissue_expression = normalize.deduplicate(tissue_expression)
        # convert from counts to TPM
        tpm = norm.tpm_from_counts(tissue_expression)
        # compute some statistics
        mean = pandas.concat(
                [mean, pandas.DataFrame(tpm.mean(), columns=[t])], axis=1)
        median = pandas.concat(
                [median, pandas.DataFrame(tpm.median(), columns=[t])], axis=1)
        std = pandas.concat(
                [std, pandas.DataFrame(tpm.std(), columns=[t])], axis=1)
        lower_quartile = pandas.concat(
                [lower_quartile, pandas.DataFrame(
                        tpm.quantile(q=0.25, axis=0)).rename(columns={0.25: t})], axis=1)
        upper_quartile = pandas.concat(
                [upper_quartile, pandas.DataFrame(
                        tpm.quantile(q=0.75, axis=0)).rename(columns={0.75: t})], axis=1)
        fraction_zero = pandas.concat([
                fraction_zero, pandas.DataFrame(
                        (tpm == 0).mean().astype(float), columns=[t])], axis=1)

    # compute the maximum pairwise hellinger distance across tissues for each gene
    genes = list(mean.index)
    hellinger = []
    for gene in genes:
        aves = numpy.ravel(mean.loc[gene].as_matrix())
        stds = numpy.ravel(std.loc[gene].as_matrix())
        hellinger += [max_hellinger(aves, stds)]
    hellinger = pandas.DataFrame(hellinger, index=genes).rename(columns={0: 'hellinger'})

    with pandas.HDFStore(os.path.join(filepath, 'tissue_stats.h5'), 'w') as store:
        store.put('mean', mean.astype(numpy.float32))
        store.put('median', median.astype(numpy.float32))
        store.put('std', std.astype(numpy.float32))
        store.put('lower_quartile', lower_quartile.astype(numpy.float32))
        store.put('upper_quartile', upper_quartile.astype(numpy.float32))
        store.put('fraction_zero', fraction_zero.astype(numpy.float32))
        store.put('hellinger', hellinger.astype(numpy.float32))
