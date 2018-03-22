import os
import pandas

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
    gene_info = pandas.read_csv(os.path.join(filepath, 'gene_info.csv'), sep='\t')
    gene_info.set_index('gene_id', inplace=True)
    gene_lengths = gene_info['bp_length']

    tissues = samples.groupby('smts').groups

    for t in tissues:
        # select the sample ids corresponding to the tissue
        index = list(samples['run'].loc[tissues[t]].values)
        tissue_expression = (expression[index].T)
        # replace missing data with 0
        tissue_expression.fillna(0.0, inplace=True)
        # convert from counts to TPM
        tissue_expression = tissue_expression / gene_lengths
        total = tissue_expression.sum(axis=1)
        tpm = 10**6 * tissue_expression.divide(total, axis='rows')
        # compute some statistics
        mean = pandas.DataFrame(tpm.mean(), columns=["mean"])
        median = pandas.DataFrame(tpm.median(), columns=['median'])
        std = pandas.DataFrame(tpm.std(), columns=['std'])
        lower_quantile = pandas.DataFrame(
                tpm.quantile(q=0.25, axis=0)).rename(columns={0.25: 'lower_quartile'})
        upper_quantile = pandas.DataFrame(
                tpm.quantile(q=0.75, axis=0)).rename(columns={0.75: 'upper_quartile'})
        # save the statistics to a csv file
        stats = pandas.concat([mean, median, std, lower_quantile, upper_quantile], axis=1)
        stats.to_csv(os.path.join(filepath, t + '.csv'))
