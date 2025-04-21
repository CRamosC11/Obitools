  **Obitools pipeline**
  
This document outlines the pipeline employed for the bioinformatic analysis of the trnL loop region of sedaDNA sequences obtained from Basa de la Mora (Pyrenees, Spain).
The research is part of the CORREDORAS project. We provide the code necessary to replicate the analysis, including the use of OBITools version 3.0.1 for cleaning and filtering raw sequences, as well as the methodology used to construct a reference database for sequence comparison.


`module load obitools/3.0.1-beta24` 

**Set a directory**
`cd /home/scc/cramos`

### **Check the Quality of the sequences**
`module load fastqc` 
Unzip folders with sequences

`gunzip -k 250207_A00902_A_L002_BFTV-2_R1.fastq.gz`
`gunzip -k 250207_A00902_A_L002_BFTV-2_R2.fastq.gz`

### **Funcion fastqc para ver la calidad**
`fastqc 250207_A00902_A_L002_BFTV-2_R1.fastq`
`fastqc 250207_A00902_A_L002_BFTV-2_R2.fastq`

Fastqc generates a zip file within the directory with images and a report of the analysis 

### **REPORT ANALISYS**
### **Quality score per base of all sequences.** 
This graph shows that most sequences have a high quality score.
However, the last nucleotides, located between positions 140 and 150 bp, tend to have lower quality. This indicates that, across all sequences analyzed, the majority exhibit a decrease in nucleotide quality towards the end. The blue line represents the mean quality score, with values above 20 considered to indicate good quality. Overall, the sequence quality is quite good.
![Captura de pantalla 2025-03-31 153237](https://github.com/user-attachments/assets/c2fb0b83-d693-4aed-8cf7-8b6f5f1b6d89)

### **Quality per sequences**
What we wanna see is tight pattern. No more than 2 distributions within the data. 
![Captura de pantalla 2025-03-31 154049](https://github.com/user-attachments/assets/415febb0-fc1b-4050-9b79-999d4e20911a)

### **Per base sequences contest**
According to tutorials what we must see are paralells lines which means that all the sequences you are analising present the same base at the same position. In our case, we may have different organisms with a lot of variability among their bases which could explain that we have all this different nucleotids among the different positions of the sequences. Also this could be produced by contamination or problems in the quality of the sequences but as we saw before, the overall quality is high.
![Captura de pantalla 2025-03-31 154203](https://github.com/user-attachments/assets/a51dce35-f4b0-4d1d-b812-7f9af083a8ac)

### **Per sequence GC content**
It presents the distribution of the content of GC over all the sequences. The blue line shows a theorical distribution while red line shows the GC content of our data. In this case, our data doesn´t fit with the theorical distribution because the GC content per sequence shows a normal distribution with some peaks which could mean contamination. A second option is that, If we observed an increase in the GC% , could mean that there is ADN from some organisims that has higher content of G and C (for example bacteria or even human DNA). Also some primers or indexs have high GC content which could be affecting the data. 
![Captura de pantalla 2025-03-31 154426](https://github.com/user-attachments/assets/9b03afa5-5202-41ea-8e66-2f999b4f40d3)

### **Sequence Duplication Levels**
Percent of sequences remainin if deduplicated is 15.33 which means that 15% of the sequwnces are non unique.
![Captura de pantalla 2025-03-31 160314](https://github.com/user-attachments/assets/c7946d52-0804-407a-bfa6-99687ffab790)
### **Overrepresented sequences**
FastQC indicates that the possible source of the 2 sequences more overrepresented is TruSeq Adapter which could explain the high % of GC content that this file has.  Reverse sequences present a overrepresented sequence (GGGGGGGGG....) that also could explain the high % of GC content of the file but I don´t know how these sequences come from or how they were originated.

Forward sequences FastQC analysis
["C:\Users\CRISTINA\NextCloud\CRC\1_data\BSM\2024\BSM24_FastQC_TrnL_R1_fastqc.html"](url)
Reverse sequences FastQC analysis
["C:\Users\CRISTINA\NextCloud\CRC\1_data\BSM\2024\BSM24_FastQC_TrnL_R2_fastqc.html"](url)


### **Import forward and reverse sequences using obi import**

`obi import --quality-sanger /data/scc/edna/data/BFTV-2_Pyrenees/250207_A00902_A_L002_BFTV-2_R1.fastq.gz trnL/forward`
`obi import --quality-sanger /data/scc/edna/data/BFTV-2_Pyrenees/250207_A00902_A_L002_BFTV-2_R2.fastq.gz trnL/reverse`

We could check the number of sequences we have:
`obi count trnL/forward`
`obi count trnL/reverse`
In our case we have 22.421.828 forward sequences and 22.421.828 reverse

### **Import ngsfilter file with information about tags and primers of the sequences**
`obi import --ngsfilter ngsfile_trnL.txt trnL/ngsfile`
(before importing to the directory, the file has to be imported manually)

**Different commands to get general information about the directory and archives**
`obi ls trnL`
`obi ls -l trnL`
`obi ls trnL/forward`
`obi ls trnL/ngsfile/sample`
`obi clean_dms trnL` to clean the DMS 

### **Aligned sequences. Recover the full sequences from the partial forward and reverse reads
• Aligns the two reads of a pair-end library sequenced optimazing the overlap.
• If the two reads overlap, it returns the consensus sequence together with its quality. It calculates the alignment score and indicates the quality of both reads´ assembly.
• Otherwise, it concatenates sequence merging the forward read and the reversed-complemented reverse read**

`obi alignpairedend -R trnL/reverse trnL/forward trnL/aligned_reads`


There are different parameters to filter the information you wanna get:
-R= indicates which is the reverse read
-m N= minimun number of bases that overlap in the consensus sequences
-M N= maximun number of mistakes allowed in the consensus sequence
-q= filters the low quality sequences.

`e.g: obi alignpairedend **-M 2** -R trnL/reverse trnL/forward trnL/aligned_reads`

### **Check the Average alignment score**

**SCORE NORM=A real value computed based on the alignment score divided by the alignment length.** [https://pythonhosted.org/OBITools/attributes/score_norm.html](url)
`obi stats -a score_norm trnL/aligned_reads`
Result= 0.948971

**SCORE=A real value computed based on the alignment of two paired-end reads** 
`obi stats -a score trnL/aligned_reads`
Result= mean score = 99.965519

The alignment score indicates the similarity among forward and reverse reads in the overlap.
The alignment score is based on the **number of bases that match** on the overlap, **the mismatches** (errors) between read 1 and read 2 and the **Phred quality score**

Image below shows a chromatogram of a  DNA fragment sequence. The numbers indicated Phred score which is the quality of each nucleotide. Based on  this score and the matches/mismatches of the forward and reverse reads, we obtain the alignment score. 
![Captura de pantalla 2025-03-06 164217](https://github.com/user-attachments/assets/d9658504-5f0b-454f-b8c2-36dd95b904c4)

To print the results of the alignment ->`aligned_reads.fastq`. The consensus sequences resulting from the pairing are stored in the FASTQ file with several annotations allowing for an easy filtering of reads based on the slignment quality. Recommendations: filtering out any consensus sequence with a **score < 40 or norm_score <0.9**

### **Removing unaligned sequences**
`obi grep -a mode:alignment trnL/aligned_reads trnL/good_sequences`

To remove  sequnces based on the quality score: 
`obi grep -p "sequence"[`score_norm´]>0.9" trnL/aligned_reads trnL/good_quality_sequences`

(collecting all the commands)
*`obi grep -p 'sequence[`score_norm´]>0.9' trnL/aligned_reads | obi grep -a mode:alignment > trnL/good_quality_aligned_sequences` 

According to #Alsos et al., (2015), > sequences having a low alignment quality score (threshold set at 40) were filtered out buit they use illuminapairedend function. 
> Only merged sequences with a high alignment quality score were retained (>40; this corresponds to a pair of reads that can align perfectly on at least 10 bp at each read end) #Pansú et al., 2015

 
### **Assign each sequence record to the corresponding PCR combination**
This step assigns each sequence to a PCR ID by the unique tag combination in the nsgfile. Primer and tag sequence wil be **trimmed off**
`obi ngsfilter -t trnL/ngsfile -u trnL/unid_sequences trnL/good_quality_sequences trnL/identified_sequences`
 
-t= Used to specify the file containing the samples description (with tags, primers, sample names,...)
-u= Filename used to store the sequences unassigned to any sample
unidentified_sequences would contain those sequences with error within the barcodes, sequences with incompleted primers or mutations or even contamination or artefacts.
-e=Used to specify the number of errors allowed for matching primers [default = 2]. According to Boyer et al (2016), the maximum number of mismatches allowed between primers and sequences can be set up using the -e option.  According to [https://pythonhosted.org/OBITools/scripts/ngsfilter.html] the default number of mismatches allowed is 2 (Obitools 1). In this case, there are some articles like  Alsos et al (2015) and Bienert et al (2012) that allowed a maximum of  3 mismatches with primers and no mismatches with tags. *By default *no mismatches are allowed in tags identifying samples** 

Using  the option -e
`obi ngsfilter -t trnL/ngsfile -u trnL/unidentified_sequences -e 3 trnL/good_quality_sequences trnL/ident_sequences`

with the -e option the number of sequences recovered is (IDENT_SEQUENCES): 18.664.681 and (UNIDENTIFIED_SEQUENCES): 2.101.698
witout the -e option the number of sequences recovered is 

To know the number of sequences that where correctly assigned to each sample:
`obi ls trnL`

**Check the score norm and score**
`obi stats -a score_norm trnL/ident_sequences`
0.996079 
`obi stats -a score trnL/ident_sequences`
105,32991 
After the filtering the mean score of the remain sequences has increased.

**Filters low quality sequences, improving the overall dataset quality by removing unreliable reads**
`obi grep -p "sequence['score']>=50 and sequence['score_norm']>=0.4" trnL/ident_sequences trnL/identified_sequences_filtered`
Remaining sequences: 18.586.167
Mean score: 105,599879

### **Dereplicate sequences**
Keep only one sequence and add the count to the header.
`obi uniq -m sample trnL/identified_sequences_filtered trnL/dereplicated_sequences_filtered`
Result: 136.474 sequences

### **Denoise sequences**
Once we have a uniq sequence of each sample we keep only useful tags. We are filtering by 3 parameters: COUNT, MERGE_sample and length. 

`obi annotate -k COUNT -k MERGED_sample -k seq_length trnL/dereplicated_sequences_filtered trnL/cleaned_metadata_sequences_filtered`

The useless metadata is cleaned and keep only COUNT and merged_sample. COUNT parameter indicates the number of times a unique sequence appears in the original dataset (after dereplication). Indicates the abundance of the sequences and could be good to calculate relative abundances (alpha/beta diversity...). MERGE_sample is a list of sample names where the sequence was observed. We could track which samples contain certain sequence.

To extract a fasta document with the information of the sequences, their count and the sample names where the sequence was observed: 
-n 10: print just 10 lines
 `obi head trnL/cleaned_metadata_sequences_filtered clean_metadata.fasta -n 10 --fasta-output`

**Filter sequence by length**
Keep those sequences with a lenght higher or equal to 105 bp. In the tutorial, the filtered length was 80 but Konstanz pdf indicates 105.([https://git.metabarcoding.org/obitools/obitools3/-/wikis/Wolf-tutorial-with-the-OBITools3])

`obi grep -p "len(sequence)>=80" trnL/cleaned_metadata_sequences_filtered trnL/denoised_sequences_filtered`

 or #`obi grep -p "len(sequence)>=80" and sequence ["COUNT"] >= 10 trnL/cleaned_metadata_sequences_filtered trnL/denoised_sequences_filtered`
COUNT 10 is the command used to delate sequences with low abundance rate which could be errors or contamination.

> The sequences were filtered using the obigrep program according to their length as follows: sequences shorter than 50 bp for the MamP007 and 10 bp for the gh primers were removed; and sequences that had an occurrence lower than 100 over the entire dataset were removed #Giguet-Covex et al., 2014.

**Clean the sequences from PCR/sequencing errors**
(Wolf tutorial -r=0.05. Sequences that represent less that 5% will be consider as mistakes, here 1%)

`obi clean -s MERGED_sample -r 0.01 -H trnL/denoised_sequences_filtered trnL/cleaned_sequences_filtered_r01`
result: **11803 sequences** 

#Alsos et al., 2015  > occurring as at least 10 reads per PCR repeat were kept

 I would like to see if i command COUNT>= 10 to delate sequences with low abundance rate which could be errors or contamination and the number of sequences that remains
`obi grep -p 'sequence["COUNT"] >= 10' trnL/cleaned_sequences_filtered_r01 trnL/high_abund_data`
Result: **805 sequences**

Extract result (just 10 lines)
 `obi head trnL/high_abund_data FinalData.fasta -n 10 --fasta-output`

 `obi export trnL/high_abund_data --output FORMAT=fasta > FINAL.fasta`
 `obi export trnL/high_abund_data -f fasta > FINAL.fasta`
 


***Problemas con Obiconvert**. No reconoce los archivos descargados en EMBL. Obitools 3 no reconoce obiconvert, obiconvert es de obitools 2 y 1. 
para ver que versiones hay disponibles en mobaXterm `module avail 2>&1 | grep obitools`
Para cargar una version diferente primero hay que eliminar la version 
`module unload obitools/3.0.1-beta24` y
`module load obitools /1.2.13` 
o
`obitools/4.0.3`
Para ver las versiones disponibles de obitools
`module avail 2>&1 | grep obitools`


### OBITOOLS 4 to built the reference database
**Download the sequences**
Before starting with the obipcr you need to download the target sequences from EMBL, genbank, PhyloAlps...The files could be storaged in a different directory. I storaged the data here cd /data/scc/edna/YESSS/metabarcoding/Tele02/EMBL
`wget -nH --cut-dirs=6 -A 'STD_*.dat.gz' -R 'STD_HUM*.dat.gz','STD_ENV*.dat.gz' -m -np ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/std/`

**Perform a in silico PCR amplification** 

In this case I downloaded only the sequences that are include in 3 files (PLN 1, 2 and 3) to make the analyses easier and faster
`obipcr -e 3 -l 10 -L 250 --embl --forward GGGCAATCCTGAGCCAA --reverse CCATTGAGTCTCTGCACCTATC --no-order /data/archiv/SCC/edna/AnnaC/AC153IL/Database/STD_PLN_1.dat.gz /data/archiv/SCC/edna/AnnaC/AC153IL/Database/STD_PLN_2.dat.gz /data/archiv/SCC/edna/AnnaC/AC153IL/Database/STD_PLN_3.dat.gz > embl_pcr.fasta`
The main difference with obitools 3 and obitools 4 is you don´t need to download all the sequences and import them to your directory. 
The final path is: 
cd /home/scc/cramos/trnL
`obipcr -e 3 -l 10 -L 160 -forward GGGCAATCCTGAGCCAA --reverse CCATTGAGTCTCTGCACCTATC --no-order /data/scc/edna/YESSS/metabarcoding/Tele02/EMBL* > COMPLETE_ecopcr.fasta`

### Clean the database
We choose to apply these different steps of filtering to clean up the sequences obtained with obipcr.
Keep the sequences with a taxid and a taxonomic description to family, genus and species ranks with obi grep.
Remove redundant sequences (dereplicate)
Ensure that the dereplicated sequences have a taxid (taxon identifier) at the family level
Ensure that sequences each have a unique identification ID with obiannotate
Index the database

**Keep annotated sequences**
To use the -t taxonomy option on all OBITools commands, you can either enter the path to the taxonomy taxdump file downloaded online with curl `curl http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`

Download the NCBI taxonomy to "/home/scc/cramos/trnL/taxonomy/" directory
`wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`
I am using obitools 4.0.3 wich is not the version documented. The only tutorial on internet ([https://obitools4.metabarcoding.org/docs/cookbook/reference_db/#keep-annotated-sequences]) is developed for obitools 4.4.0. Therefore, there is one step we need to do before going on with the next commands. We need to unzip taxdump.tar.dz into a new file, in ouur case "/home/scc/cramos/trnL/ncbitaxdump/"  where we can stored files like "citations.dmp", "nodes.dmp", names.dmp"....

The obigrep program allows to filter sequences, to keep only those with a taxid and a sufficient taxonomic description.
`obigrep -t ncbitaxdump -A taxid --require-rank species --require-rank genus --require-rank family COMPLETE_ecopcr.fasta > COMPLETE_grep_pcr.fasta`

**Dereplicate sequences **
`obiuniq -c taxid COMPLETE_grep_pcr.fasta > COMPLETE_uniq.fasta`

**Ensure that the dereplicated sequences have a taxid at the family level** 
Some sequences lose taxonomic information at the dereplication stage if certain versions of the sequence did not have this information beforehand. So we apply a second filter of this type.
`obigrep -t ncbitaxdump --require-rank=family COMPLETE_uniq.fasta > COMPLETE_Family_uniq.fasta`

**Ensure that sequences each have a unique identifier**
Index the database
`obirefidx -t ncbitaxdump COMPLETE_Family_uniq.fasta > COMPLETE_index_sequences.fasta`

**Assigning the sequences a taxon**
FinalData correspond with the filtered and cleaned sequences 
`obitag -t ncbitaxdump -R COMPLETE_index_sequences.fasta FinalData.fasta > assigned_sequences.fasta`

It is now possible to extract the useful information for our ecological analysis from our sequence file. The results of this extraction consists of two CSV files, one describing the occurrence of each sequence variant in the different samples, and one for the metadata describing each sequence variant, which can at this stage of the analysis be considered as a Molecular Taxonomic Unit, i.e. MOTU.

**The MOTU occurrence table**
In the results file wolf_final.fasta, two attributes inform us about the distribution of MOTU abundances across samples (which here correspond to individual PCR): the merge_sample attribute and the obiclean_weight attribute.

The merge_sample attribute was set by obiuniq during the initial reads dereplication procedure. It contains the observed number of reads for each sequence variant in the different samples. The obiclean_weight attribute is the number of reads assigned to each sequence variant after the obiclean denoising (or clustering) step. The number of reads shown in this attribute takes into account not only the number of reads observed for this variant, but also the number of reads observed for the erroneous sequences clustered to this estimated genuine sequence. According to obiclean , obiclean_weight is a better estimate of the true sequence occurrence than the merge_sample attribute.

The obimatrix command creates the CSV file representing any map attribute of a OBITools4 sequence file. By default, it dumps the merge_sample attribute, but you can specify any other map attribute. Here we decided to use the obiclean_weight attribute, as we prefer to report the abundances of the MOTUs.
`obimatrix --map obiclean_weight results/COMPLETE_index_sequences.fasta > final_assigned_sequences.csv`

To create the CSV metadata file describing the MOTUs attributes, you can use obicsv with the --auto option. This will create a CSV file from the wolf_final.fasta file and automatically determine which columns to include based on their contents from the first sequence records of the input dataset. In the example below, the -i and -s options are used to include the sequence identifier and the sequence itself in the output CSV file. The result can be viewed with csvlook:

`obicsv --auto -i -s results/wolf_final.fasta > trnL_final_motus.csv`

