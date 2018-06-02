# Snakefile to run FASTQ -> Processed BAM -> basic variant calling (SNV, SV, and CNV)
# 2018.06.03 Jongsoo Yoon

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
        samtools = config['samtools'],
        sortedbam = "temp_dna/{sample}_{tn}.temp.sorted.bam", 
    output:
        mdbam = temp("temp_dna/{sample}_{tn}.temp.sorted.md.bam"),
        mdbai = temp("temp_dna/{sample}_{tn}.temp.sorted.md.bam.bai")

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        "logs/{sample}_{tn}.md.log"
    shell:
        "{input.java} -XX:ParallelGCThreads={threads} -Xmx8g -jar {input.picard} MarkDuplicates "
        "REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true I={input.sortedbam} O={output.mdbam} "
        "M={output.mdbam}.metric VALIDATION_STRINGENCY=LENIENT TMP_DIR=md_temp QUIET=true; "
        "{input.samtools} index {output.mdbam}"

rule realign:
    input:
        java = config['java'],
        gatk = config['gatk'],
        samtools = config['samtools'],
        ref = config['reference'], 
        knownindel = config['knownindel'], 
        bam = "temp_dna/{sample}_{tn}.temp.sorted.md.bam", 

    output:
        realignedbam = temp("temp_dna/{sample}_{tn}.temp.sorted.md.ir.bam"),
        realignedbai = temp("temp_dna/{sample}_{tn}.temp.sorted.md.ir.bam.bai")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    shell:
        "{input.java} -Xmx4g -jar {input.gatk} -T RealignerTargetCreator -R {input.ref} "
        " -I {input.bam} --known {input.knownindel} -o {output.realignedbam}.interval; "
        "{input.java} -Xmx4g -jar {input.gatk} -T IndelRealigner -R {input.ref} -I {input.bam} "
        "-targetIntervals {output.realignedbam}.interval -o {output.realignedbam}; "
        "{input.samtools} index {output.realignedbam}"


rule combine:
        input:
            outfiles=expand("dna_bam/{sample}_{tn}.bam", sample=config["samples"], tn=['T', 'N'])
        output:
            "done"
        shell:
            "echo {input.outfiles} > {output}"

rule baserecal:
    input:
        java = config['java'], 
        gatk = config['gatk'], 
        samtools = config['samtools'],
        ref = config['reference'], 
        dbsnp = config['dbsnp'], 
        knownindel = config['knownindel'], 
        bam = "temp_dna/{sample}_{tn}.temp.sorted.md.ir.bam"
    output:
        recalTable = temp('temp_dna/{sample}_{tn}.recaltable'),
        recalbam = "dna_bam/{sample}_{tn}.bam"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    shell:
        "{input.java} -Xmx4g -jar {input.gatk} -T BaseRecalibrator -R {input.ref} -I {input.bam} -knownSites {input.dbsnp} --knownSites {input.knownindel} -o {output.recalTable}; "
        "{input.java} -Xmx4g -jar {input.gatk} -T PrintReads -R {input.ref} -I {input.bam} -BQSR {output.recalTable} -o {recalibratedBam} -nct {threads}; "
        "{input.samtools} index {output.recalbam}"

'''
rule mpileup:
    input:
        bam = 'dna_bam/{sample}_{tn}.bam', 
        samtools = config['samtools'], 
        ref = config['reference']
    output:
        mpileup = "mpileup/{sample}_{tn}.mpileup.gz"
'''
