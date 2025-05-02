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


`obiuniq -m sample results/assembled_assigned.fastq > results/assembled_assigned_uniq.fasta`

Sequences number:


`obiannotate -k count -k merged_sample results/assembled_assigned_uniq.fasta > results/assembled_assigned_simple.fasta`


`obiclean -s sample -r 0.1 --detect-chimera -H results/assembled_assigned_simple.fasta > results/cleaned_chimeras_0.1.fasta`

`obigrep -p 'sequence.Count() == 1' results/cleaned_chimeras_0.1.fasta`


`obigrep -c 2  results/cleaned_chimeras_0.1.fasta > results/no_singleton_0.1.fasta`

`obigrep -l 10 -L 150 results/no_singleton_0.1.fasta > results/length_10/length_10_0.1.fasta`

variants= 



