# EDTI: Exogeneous drug target identification tool

##About
The availability of complete genome sequences of pathogenic bacteria and their protein complements in public domain has made it possible to determine potential drug targets in these pathogens using computer-based in-silico techniques. Intersection of two datasets, namely 
[i] a pathogen's subtractive proteome dataset with the host proteome, and 
[ii] the pathogen's minimal essential protein dataset, should represent a set of proteins whose manipulation may reasonably be expected to interfere with the pathogen's survival without adversely affecting the host. 

These proteins could thus act as potential targets for drugs acting against the particular pathogen.

##Citation

## PERL implementation
Active perl is used to develop the tool. 
###System requirements
Currently supports only Windows XP/7/8.

###Dependacies
Following modules are need to be installed
* Tkx
* GD
* DBI
* Graph

### Changes and bug fixes
* PPI data integration
* First version stable: June 12, 2015
* Wed, Dec 02, 2015  2:41:40 PM : Broadspectrum analysis, BLAST known target_databases added;
* Fri, Dec 04, 2015  4:20:21 PM : File path problem solved; Param seting s for Downstream anal added;
* Sat, Dec 05, 2015 9:29:32 PM : Broad spectrum analysis Settings; select pathogens based on taxonomy based selection;

##Futures supported/proposed
- [ ] Single project should contain all outputs
  - [x] Import all inputs to create a single project
  - [x] All outputs in the same project
  - [ ] Import a runned project to rerun a uncomplete task
- [x] Putative Bactrial Drug targets
 - [x] Single-copy genes (CD-hit) + Non-host proteins(BLAST) + Essntial proteins (BLAST vs DEG)
 - [x] Single-copy genes (CD-hit) + Non-host proteins (BLAST) + Hub genes via Protein-protein interaction analysis (STRING)
- [ ] Downstream annotation of Targets.
 - [x] Broadspectrum analysis
 - [x] Compare against Known target databases
 - [ ] Ontology prediction
 - [ ] Pathway prediction
 - [ ] Sub-cellular localization prediction
- [ ] Utility to create inputs
 - [ ] create PPI networks from STRING dataset
 - [ ] Create Drug traget database
 - [ ] create Broadspectrum analysis database
- [ ] Create help doc
 - [x] MS Word Documnet
 - [ ] PDF Doc
 - [ ] HTML version



