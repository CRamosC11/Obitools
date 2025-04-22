## Metabarcoding analysis using obitools 4.4.0
This document outlines the pipeline employed for the bioinformatic analysis of the trnL loop region of sedaDNA sequences obtained from Basa de la Mora (Pyrenees, Spain) as part of the CORREDORAS project. We provide the code necessary to replicate the analysis, including the use of OBITools version 4.4.0 for cleaning and filtering raw sequences, as well as the methodology used to construct a reference database for sequence comparison. We have mainly followed Obitools 4 tutorial (https://obitools4.metabarcoding.org/docs/cookbook/illumina/)

### Basic commands
module avail 2>&1 | grep obitools 
module unload obitools/3.0.1-beta24
module load obitools/4.4.0
Set a directory cd /home/scc/cramos

### Check the quality of the sequences with FastQC
```module load fastqc
