# Metabarcoding analysis of 16S region of BSM24 sequences using Obitools 4.4.0 

`cd /home/scc/cramos/16S`

`obipairing --min-identity=0.8 --min-overlap=10 -F 16S_sequences_forward.fastq.gz -R 16S_sequences_forwardreverse.fastqgz  > results/consensus.fastq`

`obigrep -p 'annotations.mode != "join"' results/consensus.fastq > results/assembled.fastq`

Sequences number:

`obimultiplex -s ngsfile_16S_Obitools_4.csv -u results/unidentified_new.fastq results/assembled.fastq > results/assembled_assigned.fastq`


`obiuniq -m sample results/assembled_assigned.fastq > results/assembled_assigned_uniq.fasta`

Sequences number:


`obiannotate -k count -k merged_sample results/assembled_assigned_uniq.fasta > results/assembled_assigned_simple.fasta`


`obiclean -s sample -r 0.1 --detect-chimera -H results/assembled_assigned_simple.fasta > results/cleaned_chimeras_0.1.fasta`

`obigrep -p 'sequence.Count() == 1' results/cleaned_chimeras_0.1.fasta`

