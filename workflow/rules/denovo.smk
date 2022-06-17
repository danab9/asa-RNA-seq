rule merge:
    input:
        r1= expand(lambda wildcards: samples.at[wildcards.sample, 'fq1'], sample = IDS) if config["skip_trimming"]=='True' else expand("../results/fastq/trimmed/{sample}_r1_P.fastq.gz", sample = IDS),
        r2= expand(lambda wildcards: samples.at[wildcards.sample, 'fq2'], sample = IDS) if config["skip_trimming"]=='True' else expand("../results/fastq/trimmed/{sample}_r2_P.fastq.gz", sample = IDS),
    output:
        r1 = "../results/fastq/merged/merged_r1.fq",
        r2 = "../results/fastq/merged/merged_r2.fq" # has to be fq file for trinity
    log:
        "../results/logs/denovo/merge.log"
    conda:
        "../envs/trinity.yaml"
    script:
        "../scripts/merge.py"

rule trinity:  
    input:
        r1="../results/fastq/merged/merged_r1.fq",
        r2="../results/fastq/merged/merged_r2.fq"
    output:
        "../results/denovo_assembly/trinity.fasta"
    conda:
        "../envs/trinity.yaml"
    log:
        "../results/logs/denovo/trinity.log"
    threads: 10
    shell:
        "Trinity  --max_memory 10G --seqType fq --left {input.r1}  --right {input.r2} --CPU {threads} --output {output} &> {log}"

rule kallisto_index:
    input:
        fasta = "../results/denovo_assembly/trinity.fasta"
    output:
        index= directory("../results/denovo_assembly/index") #todo, I am not sure if this a directory, "-i index"used in the manual;
    conda:
        "../envs/kallisto.yaml"
    log:
        "../results/logs/denovo/kallisto_index.log"
    threads: 1
    shell:
        "kallisto index -i {output.index} {input.fasta} &> {log}"

rule kallisto_quant:
    input:
        index = "../results/denovo_assembly/index", #todo, I am not sure if this a directory, "-i index"used in the manual;
        r1= lambda wildcards: samples.at[wildcards.sample, 'fq1'] if config["skip_trimming"]=='True' else "../results/fastq/trimmed/{sample}_r1_P.fastq.gz",
        r2= lambda wildcards: samples.at[wildcards.sample, 'fq2'] if config["skip_trimming"]=='True' else "../results/fastq/trimmed/{sample}_r2_P.fastq.gz"
    output:
        dir = directory("../results/denovo_assembly/{sample}")
    conda:
        "../envs/kallisto.yaml"
    log:
        "../results/logs/denovo/kallisto_quant_{sample}.log"
    threads: 4
    shell:
        "kallisto quant -i {input.index} -o {ouptut.dir} -t {threads} {input.r1} {input.r2} &> {log}"
