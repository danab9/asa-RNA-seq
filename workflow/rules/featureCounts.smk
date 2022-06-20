rule featureCounts:
    input:
        gff = config['annotation'],
        bam = "../results/bam/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "../results/counts/featurecounts/{sample}.txt"
    conda:
        "../envs/featurecounts.yaml"
    log:
        "../results/log/featureCounts/{sample}.log"
    threads: 8
    params:
        g = "-g 'Parent'" if config['annotation'].endswith('.gff3') else ''
    shell:
        "featureCounts -T {threads} {params.g} -a {input.gff} -o {output} {input.bam} &> {log}"

