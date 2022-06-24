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# asa-RNA-seq


# Directories 
**config**
 - config

**workflow**
  - envs
    - fastqc 
  - report 
  - rules
  - scripts 
  - Snakefile

**resources** 
  - reference
     - chromosome.fa
     - annotation.gff3
  - adapters
  - fastq
  - samples.tsv
 
**results** 

  - logs
  - counts
    - kallisto
    - FeatureCounts
    - rsem
  - bam 
     - sorted (*or bam_sorted*
  - sam 
     - {sample}.sam
  - trimmed
- fastq
  - trimmed
     - {sample}_1.fastq.gz
     - {sample}_2.fastq.gz
