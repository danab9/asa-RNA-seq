# Sleuth tutorial: https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html
# DESeq2: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#
# For both approaches return:
# - Table resulting from the differential analysis containing the p/q values, fold-changes
# - A PCA plot

rule deseq:
    input:
        counts = expand("../results/counts/featurecounts/{sample}.txt", sample=IDS),
    output:
        pca = '../results/differential_expression/deseq/pca.png',
        table = '../results/differential_expression/deseq/table.csv',
    conda:
        "../envs/env.yaml"
    threads: 1
    script:
        "../scripts/deseq.R"