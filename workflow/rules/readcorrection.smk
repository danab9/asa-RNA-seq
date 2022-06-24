rule rcorrector:
    input:
        r1= expand(lambda wildcards: samples.at[wildcards.sample, 'fq1'], sample = IDS) if config["skip_trimming"]=='True' else expand("../results/fastq/trimmed/{sample}_r1_P.fastq.gz", sample = IDS),
        r2= expand(lambda wildcards: samples.at[wildcards.sample, 'fq2'], sample = IDS) if config["skip_trimming"]=='True' else expand("../results/fastq/trimmed/{sample}_r2_P.fastq.gz", sample = IDS)
    output:
        dir = directory("../results/fastq/rcorrected")
    log:
        "../results/logs/rcorrector.log"
    conda:
        "../envs/rcorrector.yaml"
    shell:
        """
        rcorrector -1 {input.r1} -2 {input.r2} -od {output.dir}
        """
