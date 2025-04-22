## Metabarcoding analysis using obitools 4.4.0
This document outlines the pipeline employed for the bioinformatic analysis of the trnL loop region of sedaDNA sequences obtained from Basa de la Mora (Pyrenees, Spain) as part of the CORREDORAS project. We provide the code necessary to replicate the analysis, including the use of OBITools version 4.4.0 for cleaning and filtering raw sequences, as well as the methodology used to construct a reference database for sequence comparison. We have mainly followed Obitools 4 tutorial (https://obitools4.metabarcoding.org/docs/cookbook/illumina/)

### Basic commands
`module avail 2>&1 | grep obitools´
`module unload obitools/3.0.1-beta24´
`module load obitools/4.4.0´
Set a directory `cd /home/scc/cramos´

### Check the quality of the sequences with FastQC
`module load fastqc´
Unzip folders with forward and reverse sequences
`gunzip -k 250207_A00902_A_L002_BFTV-2_R1.fastq.gz gunzip -k 250207_A00902_A_L002_BFTV-2_R2.fastq.gz´
Fastqc generates a zip file within the directory with images and a report of the analysis
## Quality report
In this report we are analysing the quality of the forward sequences. Reverse sequences don´t show any anomalies or quality problems therefore we just show the graphs of the forward sequences file. 
### Quality score per base of all sequences
This graph shows that most sequences have a high quality score. However, the last nucleotides, located between positions 140 and 150 bp, tend to have lower quality. This indicates that, across all sequences analyzed, the majority exhibit a decrease in nucleotide quality towards the end. The blue line represents the mean quality score, with values above 20 considered to indicate good quality. Overall, the sequence quality is quite good.
![Captura de pantalla 2025-04-22 163402](https://github.com/user-attachments/assets/842c9f3e-1f7f-4cbd-88fb-031840effdc5)
### Quality per sequences
What we wanna see is tight pattern. No more than 2 distributions within the data. The pattern observed in the quality of our dataset is good.
![Captura de pantalla 2025-04-22 163412](https://github.com/user-attachments/assets/0515c7a3-8b4c-44f2-bace-a70591ef47ea)
### Per base sequences contest
The pattern of the lines vary according the type of data you have. Working with different species, genus, families...mean that the pattern you are going to obtained is not homogeneous and parallels which means that all the sequences have the same base at the same position. In our case, we may have different organisms with a lot of variability among their bases which explains that we have all these different nucleotids contest among the different positions of the sequences. Also, in some cases, this situation could be explained by contamination or problems in the quality of the sequences but as we saw before, the overall quality is high. 
![Captura de pantalla 2025-04-22 163422](https://github.com/user-attachments/assets/6845ee5d-f1ec-4a3c-8f0c-801a42b5c054)

### Per sequence GC content
It presents the distribution of the content of GC over all the sequences. The blue line shows a theorical distribution while red line shows the GC content of our data. In this case, our data doesn´t fit with the theorical distribution because the GC content per sequence shows a normal distribution with some peaks which could mean contamination. A second option is that, If we observed an increase in the GC% , could mean that there is ADN from some organisims that has higher content of G and C (for example bacteria or even human DNA). Also some primers or indexs have high GC content which could be affecting the data.
![Captura de pantalla 2025-04-22 163430](https://github.com/user-attachments/assets/f9623f29-c2d1-47a4-8757-701f12ff0fac)

### Sequence Duplication Levels
Percent of sequences remainin if deduplicated is 15.33 which means that 15% of the sequwnces are non unique.
![Captura de pantalla 2025-04-22 163439](https://github.com/user-attachments/assets/3b3181e9-f565-41c3-9483-09db22700202)



