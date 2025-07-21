# Metabarcoding analysis of 16S region of BSM24 sequences using Obitools 4.4.0 

`cd /home/scc/cramos/16S`

`module load obitools/4.4.0`

***FASTQC analysis***
`module load fastqc`

`gunzip -k 16S_sequences_foward.fastq.gz gunzip -k 16S_sequences_reverse.fastq.gz`

`fastqc 16S_sequences_foward.fastq.gz 16S_sequences_reverse.fastq.gz`

**Obitools analysis**

`obipairing --min-identity=0.9 --min-overlap=20 -F 16S_sequences_forward.fastq -R 16S_sequences_reverse.fastq  > results/consensus_2.fastq`

        --min-overlap: minimum overlap between both the reads to consider the alignment (default: 20)
        --min-identity: minimum identity between overlapped regions of the reads to consider the alignment (default: 0.900000)

`obicount results/consensus_2.fastq`
sequence number : 57.976.189

`obigrep -p 'annotations.mode != "join" results/consensus.fastq > results/assembled.fastq`

Sequences number: 54.684.072

-----`obigrep -p 'annotations.mode != "join" results/consensus_2.fastq > results/assembled_2.fastq`
Sequence number: 54.161.789


`obimultiplex -s 16S_new.txt -u results/unidentified_new.fastq results/assembled.fastq > results/assembled_assigned.fastq`

Sequences obtained:46.377.578

------`obimultiplex -s 16S:new.txt -u results/unidentified_2.fastq results/assembled_2.fastq > results/assembled_assigned_2.fastq`
Sequences obtained: 46.919.806


`obiuniq -m sample results/assembled_assigned.fastq > results/assembled_assigned_uniq.fasta`
Sequences number:652.984

-----`obiuniq -m sample results/assembled_assigned_2.fastq > results/assembled_assigned_uniq_2.fasta`
Sequences number:547.302

PONER PARRAFO PARA ELIMINAR LOS NUCLEOTIDOS QUE NO SEAN A,C,G ,T 
Mantener solo las secuencias que contienen exclusivamente A, C, T, G (mayúsculas o minúsculas)
---------obigrep -p 'sequence =~ "^[ACTGactg]+$"' results/assembled_assigned_uniq_2.fasta > results/clean_sequences.fasta
Sequences number: 303.501

-----------`obiannotate -k count -k merged_sample results/clean_sequences.fasta > results/assembled_assigned_simple.fasta`

`obiannotate -k count -k merged_sample results/assembled_assigned_uniq.fasta > results/assembled_assigned_simple.fasta`
`obiclean -s sample -r 0.1 --detect-chimera -H results/assembled_assigned_simple.fasta > results/cleaned_chimeras_0.1.fasta`

Sequences number: 362.060


---------`obiclean -s sample -r 0.1 --detect-chimera -H results/clean_sequences.fasta > results/cleaned_chimeras_0.1.fasta`
Sequences number:

`obigrep -p 'sequence.Count() == 1' results/cleaned_chimeras_0.1.fasta`

---------`obigrep -p 'sequence.Count() == 1' results/cleaned_chimeras_0.1_2.fasta`

`obigrep -c 2  results/cleaned_chimeras_0.1.fasta > results/no_singleton_0.1.fasta`

sequences number: 46.638

---------`obigrep -c 2  results/cleaned_chimeras_0.1_2.fasta > results/no_singleton_0.1_2.fasta`
sequences number: 30.206

`obigrep -l 40 results/no_singleton_0.1.fasta > results/length_10/length_40_0.1.fasta`
variants= 42.261
reads=36.535.108
If -l=40 and -L=140:
variants= 30.624

--------`obigrep -l 40 results/no_singleton_0.1_2.fasta > results/length_40/length_40_0.1_2.fasta`
variants:26.677
reads:36.903.488


According to Walker et al, 2023, Sequences were then assigned to the samples they came from (ngsfilter; up to two errors allowed), while sequences that were unaligned, contained ambiguous bases, or were outside the expected barcode length (< 40 or > 140 bp) were removed.





##Sequences taxonomic assignment

`obitag -t ncbitaxo.tgz -R database/database.fasta results/length_40/length_40_0.1.fasta > results/length_40/taxo_40.fasta`
---------- `obitag -t ncbitaxo.tgz -R database/database.fasta results/length_40/length_40_0.1_2.fasta > results/length_40/taxo_40_2.fasta`
AQIUUIIIII 26/06/2025

`obiannotate  --delete-tag=obiclean_head --delete-tag=obiclean_headcount --delete-tag=obiclean_internalcount --delete-tag=obiclean_samplecount --delete-tag=obiclean_singletoncount results/length_40/taxo_40.fasta > results/length_40/taxo_red_40.fasta`

`obicsv --auto -i -s results/length_40/taxo_red_40.fasta > results/length_40/16S_MOTUS.csv`

`obimatrix --map obiclean_weight results/length_40/16S_MOTUS.csv > results/length_40/16S_occurrency.csv`

## Building database. 
`cd data/scc/cramos/16S/database`

`wget -nH --cut-dirs=6 -A 'STD_*.dat.gz' -R 'STD_HUM*.dat.gz','STD_ENV*.dat.gz','STD_PLN*.dat.gz'  -m -np ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/std/`

We don´t select sequences from Humans, environment and plants to make the database smaller.

`obipcr -e 2 -l 30 -L 150 --forward CGAGAAGACCCTATGGAGCT --reverse CCGAGGTCRCCCCAACC --no-order database/ > pcr_16S.fasta`

`obigrep -t database/taxdump.tar.gz -A taxid --require-rank species --require-rank genus --require-rank family pcr_16S.fasta > db_16S.fasta` 

`obiuniq -c taxid db_16S.fasta > db_uniq.fasta`

`obirefidx -t database/taxdump.tar.gz db_uniq.fasta > db_indexed.fasta`

`obitag -t ncbitaxo.tgz -R db_indexed.fasta results/length_40/taxo_red_40.fasta > results/length_40/assembled_taxo.fasta`

`obiannotate  --delete-tag=obiclean_head --delete-tag=obiclean_headcount --delete-tag=obiclean_internalcount --delete-tag=obiclean_samplecount --delete-tag=obiclean_singletoncount results/length_40/assembled_taxo.fasta > results/length_40/final_taxa_16S.fasta`

`obimatrix --map obiclean_weight results/length_40/final_taxa_16S.fasta > results/length_40/16S_occurrency.csv`

`obicsv --auto -i -s results/length_40/final_taxa_16S.fasta > results/length_40/16S_motus.csv`

# Distribution of the counts per sequences
`obicsv -k count results/length_40/taxo_red_40.fasta | tail -n +2 | sort -n | uniq -c | awk '{print $2,$1}'> results/length_40/distribution_count_per_sequences.csv`

Rstudio analysis

![Rplot](https://github.com/user-attachments/assets/900947ec-8343-4e07-805c-92620f202c98)

Mayority of the sequences are repeated 2,3,4 times...More than 13 times is weird so we show a histogram "zoom it"

![Rplot01](https://github.com/user-attachments/assets/dfb2673a-fa89-466c-b0ed-a5f3c0ae4f8d)

Calculating the mean and meadian of the distribution 








