# get the path to this script
args <- commandArgs(trailingOnly = TRUE)
data_path = args[1]

## Load recount R package
library('recount')
id <- 'SRP012682'

# download the data if necessary
data_file <- file.path(data_path, 'rse_gene.Rdata')
if (!file.exists(data_file)) {
  download_study(id, type="rse-gene", download=TRUE, outdir=file.path(data_path))
}

# download the file information if necessary
info_file <- file.path(data_path, 'files_info.tsv')
if (!file.exists(info_file)) {
  download_study(id, type="files-info", download=TRUE, outdir=file.path(data_path))
}

# download the phenotype data if necessary
pheno_file <- file.path(data_path, paste(id,'tsv', sep='.'))
if (!file.exists(pheno_file)) {
  download_study(id, type="phenotype", download=TRUE, outdir=file.path(data_path))
}

# load the data
print(data_file)
load(data_file)

# export the data in csv
gene_file <- file.path(data_path, "gene_info.csv")
sample_file <- file.path(data_path, "sample_info.csv")
expression_file <- file.path(data_path, "expression_data.csv")

if (!file.exists(gene_file)){
  write.table(rowData(rse_gene), gene_file, sep='\t')
}
if (!file.exists(sample_file)) {
  write.table(colData(rse_gene), sample_file, sep='\t')
}
if (!file.exists(expression_file)) {
  write.table(assays(rse_gene)$counts, expression_file, sep='\t')
}
