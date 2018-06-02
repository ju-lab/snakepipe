configfile: 'pathConfig.yaml'
configfile: 'sampleConfig.yaml'

rule all:
    input:
        'done'

rule bwa_tumor:
    input:
        bwa = config['bwa'],
        samtools = config['samtools'],
        ref = config['reference'],
        fq1 = lambda wildcards: config['samples'][wildcards.sample]['tumor_fastq'] + '1.fq',
        fq2 = lambda wildcards: config['samples'][wildcards.sample]['tumor_fastq'] + '2.fq'

    output:
        bam=temp("temp_dna/{sample}_T.temp.bam"),
        sortedbam='temp_dna/{sample}_T.temp.sorted.bam'
    params:
        rg="@RG\tID:{sample}_T\tSM:{sample}_T\tPL:Illumina"
    log:
        "logs/{sample}_T.bwa.log"
    threads: 4
    shell:
        "{input.bwa} mem -t {threads} -R '{params.rg}' {input.ref} {input.fq1} {input.fq2} |"
        "{input.samtools} view -Sb - > {output.bam}; {input.samtools} sort -@ {threads} -o {output.sortedbam} {output.bam}"


rule bwa_normal:
    input:
        bwa = config['bwa'],
        samtools = config['samtools'],
        ref = config['reference'],
        fq1 = lambda wildcards: config['samples'][wildcards.sample]['normal_fastq'] + '1.fq',
        fq2 = lambda wildcards: config['samples'][wildcards.sample]['normal_fastq'] + '2.fq'

    output:
        bam=temp("temp_dna/{sample}_N.temp.bam"),
        sortedbam='temp_dna/{sample}_N.temp.sorted.bam'
    params:
        rg="@RG\tID:{sample}_N\tSM:{sample}_N\tPL:Illumina"
    log:
        "logs/{sample}_N.bwa.log"
    threads: 4
    shell:
        "{input.bwa} mem -t {threads} -R '{params.rg}' {input.ref} {input.fq1} {input.fq2} |"
        "{input.samtools} view -Sb - > {output.bam}; {input.samtools} sort -@ {threads} -o {output.sortedbam} {output.bam}"


rule markdup:
    input:
        java = config['java'], 
        picard = config['picard'],
        sortedbam = "temp_dna/{sample}_{tn}.temp.sorted.bam", 
    output:
        mdbam = "temp_dna/{sample}_{tn}.temp.sorted.md.bam"
    threads: 4
    log:
        "logs/{sample}_{tn}.md.log"
    shell:
        "{input.java} -XX:ParallelGCThreads={threads} -Xmx8g -jar {input.picard} MarkDuplicates "
        "REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true I={input.sortedbam} O={output.mdbam} "
        "M={output.mdbam}.metric VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp_{wildcards.sample}_{wildcards.tn} QUIET=true; "
        "samtools index {output.mdbam}"

rule combine:
        input:
            outfiles=expand("temp_dna/{sample}_{tn}.temp.sorted.md.bam", sample=config["samples"], tn=['T', 'N'])
        output:
            "done"
        shell:
            "echo {input.outfiles} > {output}"

'''

rule realign:
    input:
        java = config['java']
        gatk = config['gatk']

    output:
    threads:
    shell:
    log:

rule baserecal:
    input:
    output:
    threads:
    shell:
    log:

'''
    
