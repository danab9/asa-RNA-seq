log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("tidyverse")
library("dplyr")
library("biomaRt")



samples_file <- read.csv(snakemake@input[["samples"]], sep='\t')
sample_id <- samples_file$sample

kal_dirs <- file.path("../results/counts/kallisto", sample_id)

# load auxilary table and add paths to kallisto directories
s2c <- dplyr::select(samples_file, sample=sample, condition=condition)
s2c$path <- kal_dirs



# construct sleuth object
so <- sleuth_prep(s2c, extra_bootstrap=TRUE)
# fit full model
so <- sleuth_fit(so, ~condition, 'full')
# reduced model fit
so <- sleuth_fit(so, ~1, 'reduced')
# test
so <- sleuth_lrt(so, 'reduced', 'full')
#print(models(so))

# examine test results
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
write.csv(sleuth_significant, snakemake@output[["table"]])  # save to file

# for pca:
# 1. collect gene names
# mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",  # change to transcript level
  #dataset = "hsapiens_gene_ensembl",
  #host = 'ensembl.org')
#t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
#    "external_gene_name"), mart = mart)
#t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
 # ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
# 1.1. add gene names into the sleuth table
so <- sleuth_prep(s2c)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

jpeg(file=snakemake@output[["pca"]])
plot_pca(so, color_by = 'condition')
dev.off()




