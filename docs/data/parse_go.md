# Documentation for Parse_Go (parse_go.py)

## functions

### all\_matches
```py

def all_matches(iterable, pattern)

```



Get all of the items that match a pattern.<br /><br />Args:<br /> ~ iterable (Iterable[str])<br /> ~ pattern (str or regular expression)<br /><br />Returns:<br /> ~ List[str]


### first\_match
```py

def first_match(iterable, pattern)

```



Get the first item that matches a pattern. Return None if there are no<br />matches.<br /><br />Args:<br /> ~ iterable (Iterable[str])<br /> ~ pattern (str or regular expression)<br /><br />Returns:<br /> ~ str or None


### get\_definition
```py

def get_definition(group)

```



Get the definition of the GO category.<br /><br />Args:<br /> ~ group (List[str])<br /><br />Returns:<br /> ~ str


### get\_id
```py

def get_id(group)

```



Get the GO identifier.<br /><br />Args:<br /> ~ group (List[str])<br /><br />Returns:<br /> ~ str


### get\_name
```py

def get_name(group)

```



Get the name of the GO category.<br /><br />Args:<br /> ~ group (List[str])<br /><br />Returns:<br /> ~ str


### get\_namespace
```py

def get_namespace(group)

```



Get the namespace of the GO category.<br /><br />Args:<br /> ~ group (List[str])<br /><br />Returns:<br /> ~ str


### get\_parents
```py

def get_parents(group)

```



Get all of the parents of a GO category.<br />We say that X is a parent of Y if (Y 'is_a' X) or if (Y is 'part_of' X).<br /><br />Args:<br /> ~ group (List[str])<br /><br />Returns:<br /> ~ parents (List[str])


### has\_pattern
```py

def has_pattern(string, pattern)

```



Check if a string contains a pattern:<br /><br />Args:<br /> ~ string (str)<br /> ~ pattern (str or regular expression)<br /><br />Returns:<br /> ~ bool


### make\_godict
```py

def make_godict(gofile, force=False)

```



Parses the Gene Ontology file and creates a dictionary that is easier<br />to work with. Saves the dictionary as a json file.<br /><br />Notes:<br /><br /> ~ uniprot id: column 1<br /> ~ gene symbol: column 2<br /> ~ GO Evidence codes: column 5<br /><br /> ~ Experiment:<br /> ~  ~ Inferred from Experiment (EXP)<br /> ~  ~ Inferred from Direct Assay (IDA)<br /> ~  ~ Inferred from Physical Interaction (IPI)<br /> ~  ~ Inferred from Mutant Phenotype (IMP)<br /> ~  ~ Inferred from Genetic Interaction (IGI)<br /> ~  ~ Inferred from Expression Pattern (IEP)<br /><br /> ~ Computational:<br /> ~  ~ Inferred from Sequence or structural Similarity (ISS)<br /> ~  ~ Inferred from Sequence Orthology (ISO)<br /> ~  ~ Inferred from Sequence Alignment (ISA)<br /> ~  ~ Inferred from Sequence Model (ISM)<br /> ~  ~ Inferred from Genomic Context (IGC)<br /> ~  ~ Inferred from Biological aspect of Ancestor (IBA)<br /> ~  ~ Inferred from Biological aspect of Descendant (IBD)<br /> ~  ~ Inferred from Key Residues (IKR)<br /> ~  ~ Inferred from Rapid Divergence(IRD)<br /> ~  ~ Inferred from Reviewed Computational Analysis (RCA)<br /><br /> ~ Literature:<br /> ~  ~ Traceable Author Statement (TAS)<br /> ~  ~ Non-traceable Author Statement (NAS)<br /><br /> ~ Other:<br /> ~  ~ Inferred by Curator (IC)<br /> ~  ~ No biological Data available (ND) evidence code<br /> ~  ~ Inferred from Electronic Annotation (IEA)<br /><br />Args:<br /> ~ gofile (str): path to the gene ontology file<br /> ~ force (optional; bool): overwrite the json file if true<br /><br />Returns:<br /> ~ bool


### parse\_group
```py

def parse_group(group, dictionary)

```



Parse a GO category and add it to a dictionary.<br /><br />Notes:<br /> ~ Modifies dictionary in place!<br /><br />Args:<br /> ~ group (List[str])<br /> ~ dictionary (dict)<br /><br />Returns:<br /> ~ None

