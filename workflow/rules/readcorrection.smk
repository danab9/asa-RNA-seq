rule rcorrector:
    input:
        r1= expand(lambda wildcards: samples.at[wildcards.sample, 'fq1'], sample = IDS) if config["skip_trimming"]=='True' else expand("../results/fastq/trimmed/{sample}_r1_P.fastq.gz", sample = IDS),
        r2= expand(lambda wildcards: samples.at[wildcards.sample, 'fq2'], sample = IDS) if config["skip_trimming"]=='True' else expand("../results/fastq/trimmed/{sample}_r2_P.fastq.gz", sample = IDS)
    output:
        dir = directory("../results/fastq/rcorrected"),
        test = "test.txt"
    log:
        "../results/logs/rcorrector.log"
    params:
        r1= expand("../results/fastq/trimmed/{sample}_r1_P.fastq.gz,",sample=IDS), #they must be comma seperated
        r2= expand("../results/fastq/trimmed/{sample}_r2_P.fastq.gz,",sample=IDS)

    conda:
        "../envs/rcorrector.yaml"
    shell:
        """
        rcorrector -1 {params.r1} -2 {params.r2} -od {output.dir} > {log}
        """
        # or just rcorrecor -1 ?
