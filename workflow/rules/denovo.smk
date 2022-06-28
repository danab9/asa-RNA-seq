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
    params: bootstraps=config["kallisto"]["bootstraps"]
    shell:
        """ 
        kallisto quant -i {input.index} -b {params.bootstraps} -o {output.dir} -t {threads} {input.r1} {input.r2} &> {log}
        """

rule rsem_prepare_reference:
    input:
        "../results/denovo_assembly/trinity.Trinity.fasta"
    output:
        multiext(
            "../results/rsem_refseq/rsem_refseq",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".grp",
            ".idx.fa",
            ".n2g.idx.fa",
            ".rev.1.bt2",
            ".rev.2.bt2",
            ".seq",
            ".ti",
            ".transcripts.fa"
        )
    conda:
        "../envs/rsem.yaml"
    log:
        "../results/logs/rsem/prepare_reference.log"
    shell:
        "mkdir -p ../results/rsem_refseq/ && rsem-prepare-reference --bowtie2 {input} ../results/rsem_refseq/rsem_refseq &> {log}"

rule rsem:
   # Documentation: https://github.com/bli25/RSEM_tutorial#-single-sample-analysis,
   # http://deweylab.github.io/RSEM/README.html
   # Wrapper: https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/rsem/calculate-expression.html
   # Example: https://github.com/dohlee/snakemake-star-rsem/blob/master/rules/rsem.smk
   # (conda install -c bioconda rsem)
    input:
        r1= lambda wildcards: samples.at[wildcards.sample, 'fq1'] if config["skip_trimming"]=='True' else "../results/fastq/trimmed/{sample}_r1_P.fastq.gz",
        r2= lambda wildcards: samples.at[wildcards.sample, 'fq2'] if config["skip_trimming"]=='True' else "../results/fastq/trimmed/{sample}_r2_P.fastq.gz",
        idx_ref= multiext(
            "../results/rsem_refseq/rsem_refseq",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".grp",
            ".idx.fa",
            ".n2g.idx.fa",
            ".rev.1.bt2",
            ".rev.2.bt2",
            ".seq",
            ".ti",
            ".transcripts.fa"
        )
    output:
        counts = "../results/counts/rsem/{sample}.genes.results"
        # there are more files as output
    conda:
        "../envs/rsem.yaml"
    log:
        "../results/logs/denovo/rsem_quant_{sample}.log"
    # params:
    #     trinity = directory("../results/denovo_assembly/trinity"),
    threads: 8
    shell:
        # some params I didn't include and idk if they're important:
        #                       --estimate-rspd
        #     					--append-names
        #     					--output-genome-bam
        """mkdir -p ../results/counts/rsem/ && rsem-calculate-expression -p {threads} \
        --paired-end {input.r1} {input.r2} --bowtie2 ../results/rsem_refseq/rsem_refseq \
        ../results/counts/rsem/{wildcards.sample} &> {log}"""
