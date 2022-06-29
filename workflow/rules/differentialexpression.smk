# SOURCES
# Sleuth tutorial: https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html
# DESeq2: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# Deseq inspiration: https://github.com/snakemake-workflows/rna-seq-star-deseq2/blob/master/workflow/scripts/deseq2.R



rule deseq:
    input:
        counts = expand("results/counts/featurecounts/{sample}.txt", sample=IDS),
        samples= config["samples"]
    output:
        pca = 'results/differential_expression/deseq/pca.svg',
        table = 'results/differential_expression/deseq/table.csv',
    params:
        design = config["design"]
    log:
        "results/logs/deseq.log"
    conda:
        "../envs/deseq.yaml"
    threads: 1
    script:
        "../scripts/deseq.R"


rule sleuth:
    input:
        counts = expand("results/counts/kallisto/{sample}/abundance.tsv", sample=IDS),
        samples = config["samples"]
    output:
        pca = 'results/differential_expression/sleuth/pca.png',
        table = 'results/differential_expression/sleuth/table.csv'
    conda:
        "../envs/sleuth.yaml"
    log:
        "results/logs/differential_expression/sleuth.log"
    threads: 1
    script:
        "../scripts/sleuth.R"