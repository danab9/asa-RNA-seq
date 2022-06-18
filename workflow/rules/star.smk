
rule star_genomeGenerate:
    input:
        fa = config['reference'],
        gff = config['annotation']
    output:
        "../results/reference/Genome"
    conda:
        "../envs/star.yaml"
    log:
        "../results/log/star/genome-generate.log"
    threads: 6
    shell:
        "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir  ../results/reference --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gff} &> {log}"


rule star_align:
    input:
        index="../results/reference/Genome",
        r1=(lambda wildcards: samples.at[wildcards.sample, 'fq1']),
        r2=(lambda wildcards: samples.at[wildcards.sample, 'fq2'])
    output:
        "../results/bam/{sample}/Aligned.sortedByCoord.out.bam"
    conda:
        "../envs/star.yaml"
    log:
        "../results/log/star/{sample}_alignment.log"
    threads: 6
    shell:
        "STAR --genomeDir ../results/reference/ --readFilesCommand zcat --runThreadN {threads} --readFilesIn {input.r1} {input.r2} --outFileNamePrefix ../results/bam/{wildcards.sample}/ --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard &> {log}"
