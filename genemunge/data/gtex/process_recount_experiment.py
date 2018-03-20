import os
import numpy as np
import pandas

# seed the seed for determinism
np.random.seed(137)

from omics import data_processing as proc
from paysage import backends as be
from paysage import preprocess as pre


# 1: [directory name].tsv -- metadata
# 2: files_info.tsv -- information on file content
# 3: rse_gene.Rdata -- R format file for the experiment
# 4: sample_info.csv -- information about each sample (column)
# 6: gene_info.csv -- information about each gene (row)
# 7: expression_data.csv -- gene counts per sample

big_studies = ['SRP012682', 'TCGA']


def select_genes(gene_info, zero_threshold=1.0, exclude_housekeeping=False,
                 select_by_go=True, select_transcription_factors=False):
    """
    Select genes.

    Args:
        gene_info (pandas.DataFrame)
        zero_threshold (optional; float): keep genes where the fraction of
            gtex samples where they were zero is less than zero_threshold
        exclude_housekeeping (optional; bool)
        select_by_go (optional; bool): keep genes that belong to either the
            go biological_process or molecular_function namespaces
        select_transcription_factors (optional; bool): only keep transcription
            factors

    Returns:
        gene identifiers (List[str])

    """
    mask = gene_info["fraction_zero"] < zero_threshold
    if select_by_go:
        tmp = gene_info["biological_process"] | gene_info["molecular_function"]
        mask = mask & tmp
    if exclude_housekeeping:
        mask = mask & ~gene_info["housekeeping"]
    if select_transcription_factors:
        mask = mask & gene_info["transcription_factor"]
    return list(gene_info[mask].index)


def collect_studies(ignore=None):
    """
    Collect all of the study paths.

    Args:
        ignore (List[str]): study identifiers to ignore

    Returns:
        List[str]

    """
    file_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(file_dir, 'datasets')
    subdirs = [os.path.join(data_dir, f) for f in os.listdir(data_dir) if
           os.path.isdir(os.path.join(data_dir, f))]
    subdirs = sorted(subdirs)
    if ignore is not None:
        subdirs = [d for d in subdirs if os.path.basename(d) not in ignore
                   and "__" not in d]
    return subdirs


def tcga_path():
    """
    Collect all of the study paths.

    Args:
        None

    Returns:
        List[str]

    """
    file_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(file_dir, 'datasets')
    return [os.path.join(data_dir, 'TCGA')]


def process_recount_data(filename, gene_list, subdirs, chunksize=1000):
    """
    Processes the expression data.

    Args:
        filename (str): output file name
        subdirs (List[str]): list of directories with data
        chunksize (optional; int)

    Returns:
        None

    """
    # get a list of genes based on the gene ontology
    chunks = [gene_list[i:i + chunksize]
                for i in range(0, len(gene_list), chunksize)]

    # check if the output file already exists
    # if yes, remove the output file. we will overwrite it
    file_dir = os.path.dirname(os.path.abspath(__file__))
    output_file = os.path.join(file_dir, filename)
    if os.path.exists(output_file):
        os.remove(output_file)

    # collect all of the expression data into a single hdf5 file
    # using the common set of genes
    # the data have to be divided into chunks because there are too many
    # columns (i.e., genes) for the pandas hdf5
    with pandas.HDFStore(output_file) as store:
        number_of_subdirs = 0
        for subdir in subdirs:
            number_of_subdirs += 1
            print("{0}: {1}".format(number_of_subdirs, subdir))
            # read the expression data
            filename = os.path.join(subdir, "expression_data.csv")
            expression_data = proc.prepare_expression_data(filename, gene_list)
            # add the expression data to the output file
            for i in range(len(chunks)):
                store.append('expression/chunk_{}'.format(i),
                             expression_data[chunks[i]],
                             min_itemsize=200)


def impute_(dataframe):
    """
    Imputes any zeros in each row with half of the smallest observed value.

    Notes:
        Modifies the dataframe in place!

    Args:
        dataframe (pandas.DataFrame)

    Returns:
        None

    """
    impute_values = 0.5 * dataframe.replace(to_replace=0, value=np.nan).min(axis=1)
    for index, row in dataframe.iterrows():
        dataframe.loc[index] = row.replace(to_replace=0, value=impute_values[index])


def create_binarize_transform(gene_info, gene_list):
    """
    Create a transform to binarize expression data.
    Greater than GTEx median -> 1, Less than GTEx median -> 0

    Args:
        gene_info (pandas.DataFrame)
        gene_list (List(str))

    Returns:
        callable

    """
    median_expression = be.float_tensor(gene_info.loc[gene_list]["median"])
    def binarize(tensors):
        # tensors ~ List[tensor ~ (batch_size, num_genes]
        return [be.greater_equal(t, median_expression) for t in tensors]
    return binarize


def create_soft_binarize_transform(gene_info, gene_list):
    """
    Create a transform to binarize expression data.
    Greater than GTEx median -> 1, Less than GTEx median -> 0

    Args:
        gene_info (pandas.DataFrame)
        gene_list (List(str))

    Returns:
        callable

    """
    median_expression = be.float_tensor(gene_info.loc[gene_list]["median"])
    std_expression = be.float_tensor(gene_info.loc[gene_list]["std"])

    def soft_binarize(tensors):
        # tensors ~ List[tensor ~ (batch_size, num_genes]
        return [be.expit(be.divide(std_expression, t - median_expression)) for t in tensors]
    return soft_binarize


# TODO: RPKM transform


if __name__ == "__main__":
    filedir = os.path.dirname(os.path.abspath(__file__))
    dataset_dir = os.path.join(filedir, 'processed_datasets')
    os.makedirs(dataset_dir, exist_ok=True)

    gene_info = pandas.read_csv(os.path.join(filedir, 'gene_list.csv'), index_col=0)
    ot = select_genes(gene_info, zero_threshold=1.0, exclude_housekeeping=False,
                        select_by_go=True, select_transcription_factors=True)
    hotz = select_genes(gene_info, zero_threshold=0.05, exclude_housekeeping=True,
                        select_by_go=True, select_transcription_factors=True)

    config = [
            {"data_filename": os.path.join(dataset_dir, 'noTCGA_HOTZ.h5'),
             "shuffled_filename": os.path.join(dataset_dir, 'noTCGA_HOTZ_shuffled_split.h5'),
             "gene_list": hotz,
             "subdirs": collect_studies(ignore='TCGA')
            },
            {"data_filename": os.path.join(dataset_dir, 'onlyTCGA_HOTZ.h5'),
             "shuffled_filename": None,
             "gene_list": hotz,
             "subdirs": tcga_path()
            }
            ]

    transforms = {
        "clr": pre.partial(pre.centered_log_ratio, pseudocount=1e6),
        "binary": create_binarize_transform(gene_info, hotz),
        "soft-binary": create_binarize_transform(gene_info, hotz)
        }

    for c in config:
        process_recount_data(c["data_filename"], c["gene_list"], c["subdirs"])
        if c["shuffled_filename"] is not None:
            proc.shuffle_split_expression_data(c["data_filename"], c["shuffled_filename"])
            base_name = os.path.basename(c["shuffled_filename"]).split('.')[0]
            proc.transform_data(os.path.join(c["shuffled_filename"]), base_name, transforms)
        else:
            base_name = os.path.basename(c["data_filename"]).split('.')[0]
            proc.transform_data(os.path.join(c["data_filename"]), base_name, transforms,
                                keys=['expression'])
