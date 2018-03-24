# genemunge
Tools for munging genomic data such as:
 - Converting between different types of gene identifiers
 - Searching for terms in the Gene Ontology (GO) associated with a keyword
 - Looking up housekeeping genes and transcription factors
 - Getting a list of GO terms associated with a given gene
 - Looking up how a gene is expressed across tissues
 - Normalizing a matrix of gene expression data by converting to TPM

```

## GO evidence codes

Experiment:
 - Inferred from Experiment (EXP)
 - Inferred from Direct Assay (IDA)
 - Inferred from Physical Interaction (IPI)
 - Inferred from Mutant Phenotype (IMP)
 - Inferred from Genetic Interaction (IGI)
 - Inferred from Expression Pattern (IEP)

Computational:
 - Inferred from Sequence or structural Similarity (ISS)
 - Inferred from Sequence Orthology (ISO)
 - Inferred from Sequence Alignment (ISA)
 - Inferred from Sequence Model (ISM)
 - Inferred from Genomic Context (IGC)
 - Inferred from Biological aspect of Ancestor (IBA)
 - Inferred from Biological aspect of Descendant (IBD)
 - Inferred from Key Residues (IKR)
 - Inferred from Rapid Divergence(IRD)
 - Inferred from Reviewed Computational Analysis (RCA)

Literature:
 - Traceable Author Statement (TAS)
 - Non-traceable Author Statement (NAS)

Other:
 - Inferred by Curator (IC)
 - No biological Data available (ND) evidence code
 - Inferred from Electronic Annotation (IEA)
 
 ## Common gene id types
 
 `['symbol','name','entrez_id','ensembl_gene_id','refseq_accession','uniprot_ids']`
 
```
## Install
This package requires `R` and the `recount` package from `bioconductor`.
```
source("https://bioconductor.org/biocLite.R")
biocLite("recount")
```
 
## Citations
 
 
## Similar Tools
 
 
