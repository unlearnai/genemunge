genemunge
=========

Tools for munging genomic data such as: - Converting between different
types of gene identifiers - Searching for terms in the Gene Ontology
(GO) associated with a keyword - Looking up housekeeping genes and
transcription factors - Getting a list of GO terms associated with a
given gene - Looking up how a gene is expressed across tissues -
Normalizing a matrix of gene expression data by converting to TPM

Unlearn.AI
----------

When we’re not developing super awesome open source packages like
``genemunge``, we help biopharma partners use unsupervised deep learning
to extract insights from their omics data. Learn more at
`unlearn.health <http://unlearn.health?utm_source=github&utm_medium=web&utm_campaign=genemunge>`__.

Install
-------

This library is accompanied by the following data sources: - The `Gene
Ontology <http://geneontology.org/>`__. The current version used here is
the 2018-03-27 release. -
`recount2 <https://jhubiostatistics.shinyapps.io/recount/>`__ data for
GTEx. - `HGNC <https://www.genenames.org/>`__ gene symbols. - A list of
`transcription factors <http://www.tfcheckpoint.org/>`__. - A list of
`housekeeping genes <https://www.tau.ac.il/~elieis/HKG/>`__.

Installing this package through ``pip`` (``pip install genemunge`` from PyPI,
``pip install .`` from GitHub) will use the static data that accompanies this repository.

If you wish to use the latest data from the above sources, you may
install in "develop" mode from GitHub with ``pip -e install .``. Notably, this will
download and process the recount2 GTEx data, requiring ``R`` and the
``recount`` package from ``bioconductor``:

::

    source("https://bioconductor.org/biocLite.R")
    biocLite("recount")

Citations
---------

Please cite the following papers if you make use of genemunge for a
publication.

This package:

Gene Ontology: Ashburner et al. Gene ontology: tool for the unification
of biology (2000) Nat Genet 25(1):25-9 GO Consortium, Nucleic Acids
Res., 2017

recount2: Collado-Torres L, Nellore A, Kammers K, Ellis SE, Taub MA,
Hansen KD, Jaffe AE, Langmead B, Leek JT. Reproducible RNA-seq analysis
using recount2. Nature Biotechnology, 2017.

HGNC: Gray KA, Yates B, Seal RL, Wright MW, Bruford EA. genenames.org:
the HGNC resources in 2015. Nucleic Acids Res. 2015 Jan;43(Database
issue):D1079-85.

Transcription factors: TFcheckpoint: a curated compendium of specific
DNA-binding RNA polymerase II transcription factors Konika Chawla;
Sushil Tripathi; Liv Thommesen; Astrid Laegreid; Martin Kuiper
Bioinformatics 2013.

Housekeeping genes: E. Eisenberg and E.Y. Levanon, Trends in Genetics
29, (2013)

Similar tools
-------------

If you know of similar tools that would be helpful references for users,
please contribute an attribution to them here.

1. `goatools <https://github.com/tanghaibao/goatools>`__
2. `goenrich <https://github.com/jdrudolph/goenrich>`__

GO evidence codes
-----------------

::

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

Common gene id types
--------------------

``['symbol','name','entrez_id','ensembl_gene_id','refseq_accession','uniprot_ids']``
