rule merge:
    input:
        r1= expand(lambda wildcards: samples.at[wildcards.sample, 'fq1'], sample = IDS) if config["skip_trimming"]=='True' else expand("../results/fastq/trimmed/{sample}_r1_P.fastq.gz", sample = IDS),
        r2= expand(lambda wildcards: samples.at[wildcards.sample, 'fq2'], sample = IDS) if config["skip_trimming"]=='True' else expand("../results/fastq/trimmed/{sample}_r2_P.fastq.gz", sample = IDS),
    output:
        r1 = "../results/fastq/merged_r1.fq",
        r2 = "../results/fastq/merged_r2.fq" # has to be fq file for trinity
    log:
        "../results/logs/merge.log"
    conda:
        "../envs/trinity.yaml"
    script:
        "../scripts/merge.py"

rule trinity:  # denovo assembly
    input:
        r1="../results/fastq/merged_r1.fq",
        r2="../results/fastq/merged_r2.fq"
    output:
        "../results/denovo_assembly/{sample}/contigs.fasta"
    conda:
        "../envs/trinity.yaml"
    log:
        "../results/logs/spades/{sample}.log"
    threads: 4
    shell:
        "Trinity --seqType fq --left {input.r1}  --right {input.r2} --CPU {threads} --output {output} 2> {log}"


#        r1=(lambda wildcards: samples.at[wildcards.sample, 'fq1']) if config["skip_trimming"]=='True' else "../results/fastq/trimmed/{sample}_r1_P.fastq.gz",
#        r2=(lambda wildcards: samples.at[wildcards.sample, 'fq2']) if config["skip_trimming"]=='True' else "
