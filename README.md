# Salmon-Quant
############ GOAL!!  Quantification of reads - How many reads come from each gene? #######################
### downloaded from
https://salmon.readthedocs.io/en/latest/index.html
### Check version 1.4.0
**#Description of Salmon**
 - A module for transcript quantification, it estimates the number of reads in each sample mapping to a reference transcript
- This count is an estimate of the abundance or expression level of these transcripts
- By the way Salmon is an alignment free method, read quantification does not require an input BAM file. 
- Salmon uses a selective mapping algorithm to align the reads directly to a set of target transcripts from a reference database

### EXPERIMENT !!  **using Salmon to quantify _Anopheles funestus_ reads (from RNASeq sequencing) that maps to transcripts of AnFun3 reference genome** #################

# step 1
## download the reference Transcriptome File
## download An. funestus genome
wget https://vectorbase.org/common/downloads/Current_Release/AfunestusFUMOZ/fasta/data/VectorBase-57_AfunestusFUMOZ_AnnotatedTranscripts.fasta -o /home/nattohz/Fun_RNASeq/Refseq
## download An. funestus  GFF file
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-53/gff3/anopheles_funestus/Anopheles_funestus.AfunF3.53.gff3.gz -o /home/nattohz/Fun_RNASeq/Refseq
##

### Study Design File
Download the study design file from dropbox (in .txt file format)

# step 2
### Create Transcriptome Index
module load salmon
salmon index -t /home/nattohz/Fun_RNASeq/RefseqVectorBase-57_AfunestusFUMOZ_AnnotatedTranscripts.fasta \
    -i salmon_index -k 31
 ls /home/nattohz/Fun_RNASeq/Refseq/salmon_index
##### Check out the various index files created in the salmon_index folder created by salmon
#Step 3
#quality control
## Trim raw reads to remove illumina adapter sequences, NNNN and bad quality reads
module load trimmomatic
trimmomatic Siaya_Res_R1.fastq.gz Siaya_Res_R2.fastq.gz \
 Siaya_Res_trimmed_forward_paired.fastq.gz Siaya_Res_trimmed_forward_unpaired.fastq.gz \
 Siaya_Res_trimmed_reverse_paired.fastq.gzSiaya_Res_trimmed_reverse_unpaired.fastq.gz \
 ILLUMINACLIP:/home/nattohz/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
 LEADING:3 TRAILING:3 MINLEN:35
 
 # if there are more than two pairs to be trimmed, a simple script will do this as follows
 # create a script text editor with nano and past the script below, edit the sample names to coincide with appropriate names
 nano Trim_raw.sh
 
 #!/usr/bin/env bash
 for r1 in *_R1.fastq.gz
do
  echo $r1
  sample=$(basename $r1)
  sample=${sample%_R1.fastq.gz}
  echo "Processing sample: "$sample
 module load trimmomatic
trimmomatic PE ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
 ${sample}_trimmed_forward_paired.fq.gz ${sample}_trimmed_forward_unpaired.fq.gz \
 ${sample}_trimmed_reverse_paired.fq.gz ${sample}_trimmed_reverse_unpaired.fq.gz \
 ILLUMINACLIP:/home/nattohz/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
 LEADING:3 TRAILING:3 MINLEN:35
 
 ### once trimming is done 
 ## Confirm QC is acceptabe
 ## with RNASeq data I will worry much about adapters and NNNNN , wont mind so much about overrepresentation 
 module load fastqc
 fastqc *_paired.fq.gz
 module load multiqc
 multiqc .
 ### check combined records or you can do individually by opening the html file on firefox
 
 #Step 4
## Transcript Quantificaiton with salmon
## ensure you Now we used the trimmed and paired reads to estimate transcript abundance:

gunzip /home/nattohz/Fun_RNASeq/Refseq/Anopheles_funestus.AfunF3.53.gff3.gz
gffread /home/nattohz/Fun_RNASeq/Refseq/Anopheles_funestus.AfunF3.53.gff3 -T -o /home/nattohz/Fun_RNASeq/Refseq/Anopheles_funestus.AfunF3.53.gtf
/home/nattohz/Fun_RNASeq/Quantified
salmon quant \
 --geneMap /home/nattohz/Fun_RNASeq/Refseq/Anopheles_funestus.AfunF3.53.gtf \
 --threads 2 -l A \
 -i /home/nattohz/Fun_RNASeq/Refseq/salmon_index/VectorBase-57_AfunestusFUMOZ_AnnotatedTranscripts.fasta \
 -1 ${sample}_trimmed_forward_paired.fq.gz \
 -2 ${sample}_trimmed_reverse_paired.fq.gz -o /home/nattohz/Fun_RNASeq/Quantified
 
done
 
#############  DOWNSTREAM ANALYSIS IN R ############### as follows














