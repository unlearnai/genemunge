# genemunge
Tools for munging genomic data

Example usage: getting proteins associated with the keywords ["immune"]
```
import genemunge

# set up a searcher
searcher = genemunge.search.Searcher()

# get the gene ontology terms associated with the keywords
# setting 'exact=False' will traverse the ontology graph and
# get all of the child terms as well
associated_terms = searcher.keyword_search(["immune"], exact=False)

# get the proteins associated with these terms as uniprot_ids
# set 'evidence_codes=None' to get all of the proteins
uniprot_ids = searcher.get_proteins(associated_terms, evidence_codes=None)

# convert the uniprot_ids to gene symbols
converter = genemunge.convert.IDConverter("uniprot_ids", "symbol")
gene_symbols = converter.convert(uniprot_ids)
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
