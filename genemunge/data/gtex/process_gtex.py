import os
import pandas

filepath = os.path.dirname(os.path.abspath(__file__))

samples = pandas.read_csv(os.path.join(filepath, 'SRP012682.tsv'), sep='\t')
expression = pandas.read_csv(os.path.join(filepath, 'expression_data.csv'), sep='\t')
gene_info = pandas.read_csv(os.path.join(filepath, 'gene_info.csv'), sep='\t')
gene_info.set_index('gene_id', inplace=True)
gene_lengths = gene_info['bp_length']

tissues = samples.groupby('smts').groups

for t in tissues:
    print(t)
    index = list(samples['run'].loc[tissues[t]].values)
    tissue_expression = (expression[index].T)
    tissue_expression.fillna(0.0, inplace=True)
    tissue_expression = tissue_expression / gene_lengths
    total = tissue_expression.sum(axis=1)
    tpm = 10**6 * tissue_expression.divide(total, axis='rows')
    mean = pandas.DataFrame(tpm.mean(), columns=["mean"])
    median = pandas.DataFrame(tpm.median(), columns=['median'])
    std = pandas.DataFrame(tpm.std(), columns=['std'])
    lower_quantile = pandas.DataFrame(
            tpm.quantile(q=0.25, axis=0)).rename(columns={0.25: 'lower_quartile'})
    upper_quantile = pandas.DataFrame(
            tpm.quantile(q=0.75, axis=0)).rename(columns={0.75: 'upper_quartile'})
    stats = pandas.concat([mean, median, std, lower_quantile, upper_quantile], axis=1)


