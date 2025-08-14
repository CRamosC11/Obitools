## trnL Metabarcoding analysis using obitools 4.4.0
This document outlines the pipeline employed for the bioinformatic analysis of the trnL loop region of sedaDNA sequences obtained from Basa de la Mora (Pyrenees, Spain) as part of the CORREDORAS project. We provide the code necessary to replicate the analysis, including the use of OBITools version 4.4.0 for cleaning and filtering raw sequences, as well as the methodology used to construct a reference database for sequence comparison. We have mainly followed Obitools 4 tutorial (https://obitools4.metabarcoding.org/docs/cookbook/illumina/)

### Basic commands
     module avail 2>&1 | grep obitools

     module unload obitools/3.0.1-beta24

     module load obitools/4.4.0

Set a directory

     cd /home/scc/cramos/trnL
 
### Check the quality of the sequences with FastQC
     module load fastqc

Unzip folders with forward and reverse sequences.

     gunzip -k 250207_A00902_A_L002_BFTV-2_R1.fastq.gz gunzip -k 250207_A00902_A_L002_BFTV-2_R2.fastq.gz

     fastqc 250207_A00902_A_L002_BFTV-2_R1.fastq.gz
     
     250207_A00902_A_L002_BFTV-2_R2.fastq.gz

Fastqc generates a zip file within the directory with images and a report of the analysis. 
Link of both reports in this repository: BSM24_FastQC_TrnL_R2_fastqc and BSM24_FastQC_TrnL_R1_fastqc
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
Overrepresented sequences: FastQC indicates that the possible source of the 2 sequences more overrepresented is TruSeq Adapter which could explain the high % of GC content that this file has. Reverse sequences present a overrepresented sequence (GGGGGGGGG....) that also could explain the high % of GC content of the file but I don´t know how these sequences come from or how they were originated.
![Captura de pantalla 2025-04-22 163439](https://github.com/user-attachments/assets/3b3181e9-f565-41c3-9483-09db22700202)


## Filtering and cleaning our dataset
### Importing forward and reverse sequences with `obiparing`

      obipairing --min-identity=0.8 --min-overlap=10 -F forward.fastq.gz -R reverse.fastqgz  > results/consensus.fastq

-min-overlap: minimum overlap between both the reads to consider the alignment (default: 20).
--min-identity: minimum identity between overlapped regions of the reads to consider the alignment (default: 0.900000).
These options allow to discard sequences with low alignment quality. A low alignment quality corresponds to paired-end reads overlapping over less than 10 base pairs, or to paired-end reads exhibiting an alignment of less than 80% of identity.

We can observed information about the first sequence extracted 

      obicount results/consensus.fastq 

Sequences number = **22.421.828**

    
      obihead -n 1 results/consensus.fastq
      
### Exclude unpaired reads
    
    obigrep -p 'annotations.mode != "join"' results/consensus.fastq > results/assembled.fastq  

Seqeunces number = **21.342.441** 

### Assign each sequence record to the corresponding sample and marker combination
Each sequence record is assigned to its corresponding sample, primer and tag. We need to compare with a csv file created with and specific structure. This strucutre only works for obitools 4.4.0. 

      obimultiplex -s ngsfile.csv -u results/unidentified_new.fastq results/assembled.fastq > results/assembled_assigned.fastq

Sequences number = **20.136.204**

Within the file ngsfile we can set up the number of error allowed between the sequence primer and the primer information contained in ngsfile. Default error is 1. The match between the primers and their corresponding sites in the obtained sequences can exhibit at most two mismatches.

### Reads dereplication

      obiuniq -m sample results/assembled_assigned.fastq > results/assembled_assigned_uniq.fasta

Sequences number = **256.282**

Command -m sample option stores the frequency of the sequence in each sample. To keep only these two attributes "count" and "merged_sample" in the sequence definition.

### Obigrep to delete ambiguous sequences

There are some sequences whose bases are in some cases ambiguous. We filter those sequencies that could lead to mistakes in the identification.

      obigrep -p 'sequence =~ "^[actg]+$"' results/assembled_assigned_uniq.fasta > results/clean_sequences.fasta

Sequences number = **148.315**

      obiannotate -k count -k merged_sample results/clean_sequences.fasta > results/assembled_assigned_simple.fasta
      
### Dataset denoising 

Having all sequences assigned to their respective samples does not mean that all these sequences are biologically meaningful. Some of these sequences can correspond to PCR/sequencing errors, or chimeras.

Parameter -r is the threshold ratio between counts (rare/abundant counts) of two sequence records so that the less abundant one is a variant of the more abundant (default: 1.00).

      obiclean -s sample -r 0.1 --detect-chimera -H results/assembled_assigned_simple.fasta > results/cleaned_chimeras_0.1.fasta

Sequences number = **87.544**

### Filter by abundance

If we are interested on filtering the file by the number of times sequences appear in the database, we can use `obigrep`

> According to Alsos et al., 2015, they filtered those sequences that occur at least 10 reads per PCR. The obigrep command below keeps only sequences that occur at least ten times in the data set.

      obigrep -c 10  results/cleaned_chimeras_0.1.fasta > results/no_singleton_0.1.fasta

Sequences number = **45.272**

We want to filter the sequences according to their length. The primer used to sequences allows us to obtain information about a  specific region of the trnL gen. We use g and h primer whcih focus on the trnL loop and the length of this region is inferior to 150 bp and this length also could vary among different organisms. We consider filtering sequences according differnet ranges: 
* length between 10 and 150 base pairs: 44.608

      obigrep -l 10 -L 150 results/no_singleton_0.1.fasta > results/length_10/length_10.fasta

Sequences number =  **9.796** 

# Distribution of counts per AVS
It could be interesting to explore the AVS (amplicon variants sequences) in order to filter them (in case we consider it neccessary) by the number of times each sequence appears in the dataset. 

      obicsv -k count results/length_10/length_10_0.1.fasta | tail -n +2 | sort -n | uniq -c | awk '{print $2,$1}' > results/length_10/counts_per_sequences_distribution.csv

With the csv, we plot the distribution of the count attribute in the data set in Rstudio (code attach in this repository). 

![Rplot](https://github.com/user-attachments/assets/28055f72-369a-4af0-9780-fdf6d72fb570)

The y-axis represents the ‘count’ attribute, which is the number of occurrences of a sequence in the dataset. The x-axis represents the number of sequences that occur that many times. A priori, we can see majority of the sequences occur between 2 and 74 times. We "zoom it".

![Rplot01](https://github.com/user-attachments/assets/02de6909-c87b-4abc-9167-cf845e5b496f)


We calculate the percent of times the sequences variants occur. Majority of sequences (36%) occur just 2 times and the 75,6% of sequence variants occur between 2 and 10 times. 

![Captura de pantalla 2025-05-08 120953](https://github.com/user-attachments/assets/44cd577a-fdfc-4f39-9a89-a398b4c88de6)

Calculating the median of the occurence of the sequences, the value is 3. In this case, there are many sequences occurring few times, therefore, i will not grep the variants by their counts. 


## Sequences taxonomic assignment 
Once the dataset is curated, the next step in a classical diet metabarcoding analysis is to assign the barcodes a taxon name (species, genus, etc.), in order to retrieve the list of taxa detected in each sample.

      obitag -t ncbitaxo.tgz -R database/database_EMBL.fasta results/length_10/length_10.fasta > results/length_10/taxo_seq_10.fasta

* before starting with obitag, we need to construct the database (in our case we have construct 2 db using EMBL and PhyloAlps data). In this script we are showing results with EMBL database. Construction of the database commands is shown in this repository. 

Exporting the results in a tabular format

      obiannotate  --delete-tag=obiclean_head --delete-tag=obiclean_headcount --delete-tag=obiclean_internalcount --delete-tag=obiclean_samplecount --delete-tag=obiclean_singletoncount results/length_10/taxo_seq_10.fasta > results/length_10/taxo_seq_red_10.fasta


## MOTUs and OCCURRENCE files

MOTU: molecular operational taxonomic units

To create the CSV metadata file describing the MOTUs attributes, you can use obicsv with the --auto option. 

      obicsv --auto -i -s results/length_10/taxo_seq_red_10.fasta > results/length_10/MOTUS_trnL.csv

The `obimatrix` command creates the CSV file using the obiclean_weight attribute to report the abundances of the MOTUs.

      obimatrix --map obiclean_weight results/length_10/taxo_seq_red_10.fasta > results/length_10/occurrence_trnL_EMBL.csv



      
 ## PhyloAlps database

      obitag -t ncbitaxo.tgz -R PhyloAlps/phyloalps_new.fasta results/length_10/length_10.fasta > results_phyloAlps/length_10/Phylo_taxo_seq_10.fasta

      obiannotate  --delete-tag=obiclean_head --delete-tag=obiclean_headcount --delete-tag=obiclean_internalcount --delete-tag=obiclean_samplecount --delete-tag=obiclean_singletoncount results_phyloAlps/length_10/Phylo_taxo_seq_10.fasta > results_phyloAlps/length_10/Phylo_taxo_seq_anno.fasta

      obicsv --auto -i -s results_phyloAlps/length_10/Phylo_taxo_seq_anno.fasta > results_phyloAlps/length_10/Phylo_MOTUS_trnL.csv

To create the CSV metadata file describing the MOTUs attributes, you can use obicsv with the --auto option. 

      obimatrix --map obiclean_weight results_phyloAlps/length_10/Phylo_taxo_seq_anno.fasta > results_phyloAlps/length_10/Phylo_occurrences_trnL.csv

