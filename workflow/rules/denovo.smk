# rule merge:
#     input:
#         r1= expand(lambda wildcards: samples.at[wildcards.sample, 'fq1'], sample = IDS) if config["skip_trimming"]=='True' else expand("../results/fastq/trimmed/{sample}_r1_P.fastq.gz", sample = IDS),
#         r2= expand(lambda wildcards: samples.at[wildcards.sample, 'fq2'], sample = IDS) if config["skip_trimming"]=='True' else expand("../results/fastq/trimmed/{sample}_r2_P.fastq.gz", sample = IDS),
#     output:
#         r1 = "../results/fastq/merged/merged_r1.fq",
#         r2 = "../results/fastq/merged/merged_r2.fq" # has to be fq file for trinity
#     log:
#         "../results/logs/denovo/merge.log"
#     conda:
#         "../envs/trinity.yaml"
#     script:
#         "../scripts/merge.py"


rule merge:
    input:
        r1= expand(lambda wildcards: samples.at[wildcards.sample, 'fq1'], sample = IDS) if config["skip_trimming"]=='True' else expand("../results/fastq/trimmed/{sample}_r1_P.fastq.gz", sample = IDS),
        r2= expand(lambda wildcards: samples.at[wildcards.sample, 'fq2'], sample = IDS) if config["skip_trimming"]=='True' else expand("../results/fastq/trimmed/{sample}_r2_P.fastq.gz", sample = IDS)
    output:
        r1 = "../results/fastq/merged/merged_r1.fq",
        r2 = "../results/fastq/merged/merged_r2.fq" # has to be fq file for trinity
    log:
        "../results/logs/denovo/merge_r1.log"
    shell:
        """
        zcat {input.r1} > {output.r1}  # zcat is cat for gzipped files
        zcat {input.r2} > {output.r2} 
        """


rule trinity:
    input:
        r1="../results/fastq/merged/merged_r1.fq",
        r2="../results/fastq/merged/merged_r2.fq"
    output:
        dir = directory("../results/denovo_assembly/trinity"),
        fasta = "../results/denovo_assembly/trinity.Trinity.fasta"
    conda:
        "../envs/trinity.yaml"
    log:
        "../results/logs/denovo/trinity.log"
    threads: 10
    shell:
        "Trinity  --max_memory 10G --seqType fq --left {input.r1}  --right {input.r2} --CPU {threads} --output {output.dir} &> {log}"

rule kallisto_index:
    # Documentation: https://pachterlab.github.io/kallisto/manual, example: https://github.com/EnvGen/snakemake-workflows/blob/master/bio/ngs/rules/quantification/kallisto.rules
    input:
        fasta = "../results/denovo_assembly/trinity.Trinity.fasta"
    output:
        index = "../results/denovo_assembly/index.kaix",
    conda:
        "../envs/kallisto.yaml"
    log:
        "../results/logs/denovo/kallisto_index.log"
    threads: 1
    shell:
        "kallisto index -i {output.index} {input.fasta} &> {log}"

rule kallisto_quant:
    input:
        index = "../results/denovo_assembly/index.kaix",
        r1= lambda wildcards: samples.at[wildcards.sample, 'fq1'] if config["skip_trimming"]=='True' else "../results/fastq/trimmed/{sample}_r1_P.fastq.gz",
        r2= lambda wildcards: samples.at[wildcards.sample, 'fq2'] if config["skip_trimming"]=='True' else "../results/fastq/trimmed/{sample}_r2_P.fastq.gz",
    output:
        dir = directory("../results/counts/kallisto/{sample}"),
        file = "../results/counts/kallisto/{sample}/abundance.tsv",
    conda:
        "../envs/kallisto.yaml"
    log:
        "../results/logs/denovo/kallisto_quant_{sample}.log"
    threads: 4
    shell:
        """ 
        kallisto quant -i {input.index} -o {output.dir} -t {threads} {input.r1} {input.r2} &> {log}
        """

rule rsem_prepare_reference:
    input:
        "../results/denovo_assembly/trinity.Trinity.fasta"
    output:
        "../results/rsem_refseq/rsem_refseq"
    conda:
        "../envs/rsem.yaml"
    shell:
        "mkdir -p ../results/rsem_refseq/ && rsem-prepare-reference {input} {output}"

rule rsem:
   # Documentation: https://github.com/bli25/RSEM_tutorial#-single-sample-analysis,
   # http://deweylab.github.io/RSEM/README.html
   # Wrapper: https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/rsem/calculate-expression.html
   # Example: https://github.com/dohlee/snakemake-star-rsem/blob/master/rules/rsem.smk
   # (conda install -c bioconda rsem)
    input:
        fasta = "../results/denovo_assembly/trinity.Trinity.fasta",
        r1= lambda wildcards: samples.at[wildcards.sample, 'fq1'] if config["skip_trimming"]=='True' else "../results/fastq/trimmed/{sample}_r1_P.fastq.gz",
        r2= lambda wildcards: samples.at[wildcards.sample, 'fq2'] if config["skip_trimming"]=='True' else "../results/fastq/trimmed/{sample}_r2_P.fastq.gz"
    output:
        counts = "../results/counts/rsem/{sample}/abundance.tsv",
        dir = directory("../results/counts/rsem/{sample}")
    conda:
        "../envs/rsem.yaml"
    log:
        "../results/logs/denovo/rsem_quant_{sample}.log"
    params:
        trinity = directory("../results/denovo_assembly/trinity"),
    threads: 4
    shell:
        """
            rsem-calculate-expression -p {threads} --paired-end \
    					--bowtie2 
    					--estimate-rspd \
    					--append-names \
    					--output-genome-bam \
    					{input.r1} {input.r2} \
    					{params.trinity} 
      """
   # maybe we first need to "prepare" the reference.
   # needs a last input folder, something like exp/LPS_6h
    #probably needs a bowtie path --bowtie2-path software/bowtie2-2.2.6 \