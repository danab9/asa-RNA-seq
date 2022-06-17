import gzip
from Bio import SeqIO, bgzf
# inspiration from https://bioinformatics.stackexchange.com/questions/892/how-do-you-write-a-gz-fastq-file-with-biopython

def merge(in_paths, out_path):
    out_content = " "
    for in_path in in_paths:
        with gzip.open(in_path, "rt") as f:
            fq = SeqIO.parse(f, "fastq")
            out_content = fq
        outfile = open(out_path, "w")
        outfile.write(out_content)
        outfile.close()

    # with bgzf.BgzfWriter(out_path, "wb") as outgz:
    #     SeqIO.write(sequences=out_content, handle=outgz, format="fastq")
    # I don't have to zip again

in_paths_r1 = snakemake.input["r1"]
in_paths_r2 = snakemake.input["r2"]
out_path_r1 = snakemake.output["r1"]
out_path_r2 = snakemake.output["r2"]
merge(in_paths_r1, out_path_r1)
merge(in_paths_r2, out_path_r2)
