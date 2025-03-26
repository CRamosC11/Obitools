# Obitools
This folder shows the commands used in CORREDORAS project to analysed sedaDNA sequences from BASA DE LA MORA (Pyrenees, Spain). We show the obitools commands used as  well as the way we built a database to compare with.
 **Load obitools package**
`module` load obitools/3.0.1-beta24

**Set a directory**
`cd` /data/scc/edna/YESSS/metabarcoding/Tele02

### **Check the Quality of the sequences**
`module` load fastqc
Unzip in case of need. 
`gunzip -k 241105_A00902_B_L002_AQHH-55_R2.fastq.gz`

### **Funcion fastqc para ver la calidad**
`fastqc` 241105_A00902_B_L002_AQHH-55_R2.fastq
Fastqc generates a zip file within the directory with images and a report of the analysis 

### **REPORT ANALISYS**
### **Quality score per base of all sequences.** 
This graph shows most of the sequences have a high quality score. The last nucleotids (110-150) have worse quality. This means that among all the sequencews analysed, majority of them have low quality nucleotids at the end of the sequences.
The blue line is the mean of the quality score. Above 20 is considered good quality.
![Captura de pantalla 2025-02-19 113727](https://github.com/user-attachments/assets/e8c04d93-8274-4aee-9f22-c5a158f558e2)

### **Quality per sequences**
What we wanna see is tight pattern. No more than 2 distributions within the data. 
![Captura de pantalla 2025-02-19 114633](https://github.com/user-attachments/assets/7dfbc644-d248-4e65-9f4d-cd6a673eee4c)

### **Per base sequences contest**
According to tutorials what we must see are paralells lines which means that all the sequences you are analising present the same base at the same position. In our case, we may have different organisms with a lot of variability among their bases which could explain that we have all this different nucleotids among the different positions of the sequences. Also this could be produced by contamination or problems in the quality of the sequences but as we saw before, the overall quality is high.
![Captura de pantalla 2025-02-19 121328](https://github.com/user-attachments/assets/0e30b34a-beca-4f59-8dd3-5a5f4aa33db2)

### **Per sequence GC content**
It presents the distribution of the content of GC over all the sequences.The blue line shows a theorical distribution while red line shows the GC content of our data. In this case, our data doesn´t fit with the theorical distribution which could mean contamination. This GC content per sequence graph, shows a bimodal distribution, indicating the presence of two distinct sequence populations with different GC percentages.  In metabarcoding studies, a bimodal GC content distribution may result from multiple organisms with significantly different GC compositions. If some organisms have high GC content while others have low GC content, two peaks will appear in the distribution. 
It could be expected to find anwers as the nature of our data: are these sequences ancient DNA or modern DNA from comtamination? Unfortunatly, it is not possible to check if we have fossil sequences or not by metabarcoding.  There are some tools like metaDamage or mapDamage that infer if ouir sequences present certain specific pattern that sedaDNA sequences do (this tools are only possible when using shotgun). Metabarcoding does not permit the detection of deamination because we are selecting a fragment of the whole DNA.

Fossil or ancient DNA tends to fragment and undergo deamination processes, which can alter base composition. However, in terms of GC content, degraded DNA is usually biased toward lower GC content, since guanine and cytosine bases are more prone to oxidative damage. We could see a graph pattern where the begging and the end of the DNA present more T and A.
![Captura de pantalla 2025-02-19 122550](https://github.com/user-attachments/assets/3e8092fb-9d37-4668-9fec-8ad48bb9e327)

If we observed an increase in the GC% , it could mean that there is ADN from another organisim which has a higher content of G and C (for example bacteria or even human DNA). Also some primers or indexs have high GC content which could be affecting the data. 

### **Import forward and reverse sequences using obi import**
```
obi import --quality-sanger /data/scc/edna/YESSS/AQHH-20240911b/241105_A00902_B_L002_AQHH-55_R2.fastq.gz DMS_Tele02/reads1
obi import --quality-sanger /data/scc/edna/YESSS/AQHH-20240911b/241105_A00902_B_L002_AQHH-55_R1.fastq.gz DMS_Tele02/reads2

```
### **Import ngsfilter file with information about tags and primers of the sequences**
`obi import --ngsfilter ngsfile.txt DMS_Tele02/ngsfile`
(before importing to the directory, the file has to be imported manually)

**Different commands to get general information about the directory and archives**
```
obi ls DMS_Name
obi ls -l DMS_Name
obi ls DMS_NAME/reads1
obi ls DMS_NAME/ngsfile/sample
obi clean_dms DMS_NAME to clean the DMS

```
### **Aligned sequences. Recover the full sequences from the partial forward and reverse reads
• Aligns the two reads of a pair-end library sequenced optimazing the overlap.
• If the two reads overlap, it returns the consensus sequence together with its quality. It calculates the alignment score and indicates the quality of both reads´ assembly.
• Otherwise, it concatenates sequence merging the forward read and the reversed-complemented reverse read**

`obi alignpairedend- R DMS_NAME/reads2 DMS_NAME/reads1 DMS_NAME/aligned_reads`

There are different parameters to filter the information you wanna get:
-R= indicates which is the reverse read
-m N= minimun number of bases that overlap in the consensus sequences
-M N= maximun number of mistakes allowed in the consensus sequence.
-q= filters the low quality sequences.
`e.g: obi alignpairedend **-M 2** -R DMS_NAME/reads2....`

### **Check the Average alignment score**

`obi stats -a score_norm DMS_NAME_aligned_reads`

The alignment score indicates the similarity among forward and reverse reads in the overlap.
The alignment score is based on the **number of bases that match** on the overlap, **the mismatches** (errors) between read 1 and read 2 and the **Phred quality score**

Image below shows a chromatogram of a  DNA fragment sequence. The numbers indicated Phred score which is the quality of each nucleotide. Based on  this score and the matches/mismatches of the forward and reverse reads, we obtain the alignment score. 
![Captura de pantalla 2025-03-06 164217](https://github.com/user-attachments/assets/d9658504-5f0b-454f-b8c2-36dd95b904c4)

To print the results of the alignment ->`aligned_reads.fastq`. The consensus sequences resulting from the pairing are stored in the FASTQ file with several annotations allowing for an easy filtering of reads based on the slignment quality. Recommendations: filtering out any consensus sequence with a **score < 40 or norm_score <0.9**

### **Removing unaligned sequences**
`obi grep -a mode:alignment DMS_NAME/aligned_reads DMS_NAME/good_sequences`

To remove  sequnces based on the quality score: 
`obi grep -p "sequence"["score_norm"]>0.9 DMS_NAME/aligned_reads DMS_NAME/good_quality_sequences`
`obi grep -p 'sequence["score_norm"]>0.9' DMS_NAME/aligned_reads | obi grep -a mode:alignment > DMS_NAME/good_quality_aligned_sequences` (collecting all the commands)

According to @alsos et al., 2015,  sequences having a low alignment quality score (threshold set at 40) were filtered out buit tehy use illuminapairedend function. Maybe it is not the same value with alignpairedends. Check

 
### **Assign each sequence record to the corresponding PCR combination**
This step assigns each sequence to a PCR ID by the unique tag combination in the nsgfile. Primer abd tag sequence wil be **trimmed off**
`obi ngsfilter -t DMS_NAME/ngsfile -u DMS_NAME/unidentified_sequences felchen/good_sequences felchen/identified_sequences`
-t= indicates the filter file
-u= name of the place where all of these sequences taht don´t match with any tag or primer are found.
unidentified_sequences would contain those sequences with error within the barcodes, sequences with incompleted primers or mutations or even contamination or artefacts.
To know the number of sequences taht where correctly assigned to aech sample:
`obi stats -a sample felchen/identified_sequences`

**Check the score norm and score**
`obi stats -a score_norm DMS_NAME/identified_sequences
obi stats -a score DMS_NAME/identified_sequences`

**Filters low quality sequences, improving the overall dataset quality by removing unreliable reads**
`obi grep -p "sequence ["score"]>=40 and sequence ["score_norm"]>=0.4 DMS_NAME/identified_sequences DMS_NAME/identified_sequences_filtered`

### **Dereplicate sequences**
Keep only one sequence and add the count to the header
`obi uniq -m sample DMS_NAME/identified_sequences_filtered DMS_NAME/dereplicated_sequences_filtered`

### **Denoise sequences**
Once we have a uniq sequence of each sample :
**Keep only useful tags**
`obi annotate -k COUNT -k MERGED_sample -k seq_length DMS_NAME/dereplicated_sequences_filtered DMS_NAME/cleaned_metadata_sequences_filtered`

**Filter sequence by length**
Keep those sequences with a lenght higher or equal to 105 bp. In the tutorial ([https://git.metabarcoding.org/obitools/obitools3/-/wikis/Wolf-tutorial-with-the-OBITools3]) used length of 80
`obi grep -p "len(sequence)>=105" DMS_NAME/cleaned_metadata_sequences_filtered DMS_NAME/denoised_sequences_filtered`
or
`obi grep -p "len(sequence)>=105" and sequence ["COUNT"] >= 10 DMS_NAME/cleaned_metadata_sequences_filtered DMS_NAME/denoised_sequences_filtered`
COUNT 10 is a command to delate sequences with low abundance rate which could be errors or contamination.

**Clean the sequences from PCR/sequencing errors**
`obi clean -s MERGED_sample -r 0.01 -H DMS_NAME/denoised_sequences_filtered DMS_NAME/cleaned_sequences_filtered_r01`

(Wolf tutorial -r=0.05. Sequences that represent less that 5% will be consider as mistakes)


