# SOURCES
# Sleuth tutorial: https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html
# DESeq2: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# merging feature counts: https://divingintogeneticsandgenomics.rbind.io/post/merge-featurecount-table-from-rnaseq/

# For both approaches return:
# - Table resulting from the differential analysis containing the p/q values, fold-changes
# - A PCA plot

#
# rule count_matrix:
#     input:
#         counts = expand("../results/counts/featurecounts/{sample}.txt", sample=IDS),
#     output:
#         "../results/counts/count_matrix.txt",
#     # log:
#     #     "logs/count-matrix.log",
#     # params:
#     #     samples=units["sample_name"].tolist(),
#     #     strand=get_strandedness(units),
#     # conda:
#     #     "../envs/pandas.yaml"
#     script: #../results/counts/featurecounts/*.txt
#         """
#         ls -1 {input.counts} | parallel 'cat {{}} | sed '1d' | cut -f7 {{}} > {{/.}}_clean.txt'
#         ls -1 {input.counts} | head -1 | xargs cut -f1 > genes.txt
#         paste genes.txt *featureCount_clean.txt > {output}
#         """
#


rule deseq:
    input:
        counts = expand("../results/counts/featurecounts/{sample}.txt", sample=IDS),
        samples= config["samples"]
    output:
        pca = '../results/differential_expression/deseq/pca.png',
        table = '../results/differential_expression/deseq/table.csv',
    params:
        design = config["design"]
    log:
        "../results/logs/deseq.log"
    conda:
        "../envs/deseq.yaml"
    threads: 1
    script:
        "../scripts/deseq.R"


rule sleuth:
    input:
        counts = expand("../results/counts/kallisto/{sample}/abundance.tsv", sample=IDS),
        samples = config["samples"]
    output:
        pca = '../results/differential_expression/sleuth/pca.png',
        table = '../results/differential_expression/sleuth/table.csv',
    conda:
        "../envs/env.yaml"
    threads: 1
    script:
        "../scripts/sleuth.R"