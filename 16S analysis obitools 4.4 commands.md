# Metabarcoding analysis of 16S region of BSM24 sequences using Obitools 4.4.0 

`cd /home/scc/cramos/16S`

`module load fastqc`

`gunzip -k 16S_sequences_foward.fastq.gz gunzip -k 16S_sequences_reverse.fastq.gz`

`fastqc 16S_sequences_foward.fastq.gz 16S_sequences_reverse.fastq.gz`

`obipairing --min-identity=0.8 --min-overlap=10 -F 16S_sequences_forward.fastq -R 16S_sequences_reverse.fastq  > results/consensus.fastq`

`obicount results/consensus.fastq`

sequences number: 57.976.189 

`obigrep -p 'annotations.mode != "join"' results/consensus.fastq > results/assembled.fastq`

Sequences number: 54.684.072

`obimultiplex -s ngsfile_16S.csv -u results/unidentified_new.fastq results/assembled.fastq > results/assembled_assigned.fastq`

Sequences obtained:46.377.578

`obiuniq -m sample results/assembled_assigned.fastq > results/assembled_assigned_uniq.fasta`

Sequences number:652.984

`obiannotate -k count -k merged_sample results/assembled_assigned_uniq.fasta > results/assembled_assigned_simple.fasta`

`obiclean -s sample -r 0.1 --detect-chimera -H results/assembled_assigned_simple.fasta > results/cleaned_chimeras_0.1.fasta`

Sequences number: 362.060

`obigrep -p 'sequence.Count() == 1' results/cleaned_chimeras_0.1.fasta`

`obigrep -c 2  results/cleaned_chimeras_0.1.fasta > results/no_singleton_0.1.fasta`

sequences number: 46.638

`obigrep -l 40 results/no_singleton_0.1.fasta > results/length_10/length_40_0.1.fasta`

According to Walker et al, 2023, Sequences were then assigned to the samples they came from (ngsfilter; up to two errors allowed), while sequences that were unaligned, contained ambiguous bases, or were outside the expected barcode length (< 40 or > 140 bp) were removed.

variants= 42.261
reads=36.535.108

If -l=40 and -L=140:
variants= 30.624

##Sequences taxonomic assignment

`obitag -t ncbitaxo.tgz -R database/database.fasta results/length_40/length_40_0.1.fasta > results/length_40/taxo_40.fasta`

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



