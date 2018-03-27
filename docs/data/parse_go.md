# Documentation for Parse_Go (parse_go.py)

## functions

### all\_matches
```py

def all_matches(iterable, pattern)

```



Get all of the items that match a pattern.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;iterable (Iterable[str])<br />&nbsp;&nbsp;&nbsp;&nbsp;pattern (str or regular expression)<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;List[str]


### first\_match
```py

def first_match(iterable, pattern)

```



Get the first item that matches a pattern. Return None if there are no<br />matches.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;iterable (Iterable[str])<br />&nbsp;&nbsp;&nbsp;&nbsp;pattern (str or regular expression)<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;str or None


### get\_definition
```py

def get_definition(group)

```



Get the definition of the GO category.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;group (List[str])<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;str


### get\_id
```py

def get_id(group)

```



Get the GO identifier.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;group (List[str])<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;str


### get\_name
```py

def get_name(group)

```



Get the name of the GO category.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;group (List[str])<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;str


### get\_namespace
```py

def get_namespace(group)

```



Get the namespace of the GO category.<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;group (List[str])<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;str


### get\_parents
```py

def get_parents(group)

```



Get all of the parents of a GO category.<br />We say that X is a parent of Y if (Y 'is_a' X) or if (Y is 'part_of' X).<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;group (List[str])<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;parents (List[str])


### has\_pattern
```py

def has_pattern(string, pattern)

```



Check if a string contains a pattern:<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;string (str)<br />&nbsp;&nbsp;&nbsp;&nbsp;pattern (str or regular expression)<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;bool


### make\_godict
```py

def make_godict(gofile, force=False)

```



Parses the Gene Ontology file and creates a dictionary that is easier<br />to work with. Saves the dictionary as a json file.<br /><br />Notes:<br /><br />&nbsp;&nbsp;&nbsp;&nbsp;uniprot id: column 1<br />&nbsp;&nbsp;&nbsp;&nbsp;gene symbol: column 2<br />&nbsp;&nbsp;&nbsp;&nbsp;GO Evidence codes: column 5<br /><br />&nbsp;&nbsp;&nbsp;&nbsp;Experiment:<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Experiment (EXP)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Direct Assay (IDA)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Physical Interaction (IPI)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Mutant Phenotype (IMP)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Genetic Interaction (IGI)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Expression Pattern (IEP)<br /><br />&nbsp;&nbsp;&nbsp;&nbsp;Computational:<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Sequence or structural Similarity (ISS)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Sequence Orthology (ISO)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Sequence Alignment (ISA)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Sequence Model (ISM)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Genomic Context (IGC)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Biological aspect of Ancestor (IBA)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Biological aspect of Descendant (IBD)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Key Residues (IKR)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Rapid Divergence(IRD)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Reviewed Computational Analysis (RCA)<br /><br />&nbsp;&nbsp;&nbsp;&nbsp;Literature:<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Traceable Author Statement (TAS)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Non-traceable Author Statement (NAS)<br /><br />&nbsp;&nbsp;&nbsp;&nbsp;Other:<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred by Curator (IC)<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;No biological Data available (ND) evidence code<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Inferred from Electronic Annotation (IEA)<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;gofile (str): path to the gene ontology file<br />&nbsp;&nbsp;&nbsp;&nbsp;force (optional; bool): overwrite the json file if true<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;bool


### parse\_group
```py

def parse_group(group, dictionary)

```



Parse a GO category and add it to a dictionary.<br /><br />Notes:<br />&nbsp;&nbsp;&nbsp;&nbsp;Modifies dictionary in place!<br /><br />Args:<br />&nbsp;&nbsp;&nbsp;&nbsp;group (List[str])<br />&nbsp;&nbsp;&nbsp;&nbsp;dictionary (dict)<br /><br />Returns:<br />&nbsp;&nbsp;&nbsp;&nbsp;None

