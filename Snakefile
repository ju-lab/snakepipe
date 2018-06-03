# Snakefile to run FASTQ -> Processed BAM -> basic somatic variant calling (SNV, SV, and CNV)
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
        sortedbam=temp('temp_dna/{sample}_T.temp.sorted.bam')
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
        sortedbam=temp('temp_dna/{sample}_N.temp.sorted.bam')
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
        java = config['java8'], 
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
        java = config['java8'],
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
        " -I {input.bam} --known {input.knownindel} -o {output.realignedbam}.intervals; "
        "{input.java} -Xmx4g -jar {input.gatk} -T IndelRealigner -R {input.ref} -I {input.bam} "
        "-targetIntervals {output.realignedbam}.intervals -o {output.realignedbam}; "
        "{input.samtools} index {output.realignedbam}"


rule baserecal:
    input:
        java = config['java8'], 
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
        "{input.java} -Xmx4g -jar {input.gatk} -T PrintReads -R {input.ref} -I {input.bam} -BQSR {output.recalTable} -o {output.recalbam} -nct {threads}; "
        "{input.samtools} index {output.recalbam}"

rule mpileup:
    input:
        bam = 'dna_bam/{sample}_{tn}.bam', 
        samtools = config['samtools'], 
        ref = config['reference']
    output:
        mpileup = "mpileup/{sample}_{tn}.mpileup"
    threads: 1
    shell:
        "{input.samtools} mpileup -Q 20 -q 20 -f {input.ref} -o {output} {input.bam}"

rule varscan:
    input:
        tumor_mpileup = 'mpileup/{sample}_T.mpileup', 
        normal_mpileup = 'mpileup/{sample}_N.mpileup',
        varscan = config['varscan'], 
        java = config['java8'], 
    output:
        snvvcf = 'data_processing/varscan/{sample}.varscan.snp.vcf', 
        indelvcf = 'data_processing/varscan/{sample}.varscan.indel.vcf'
    threads: 1
    params:
        basename = 'data_processing/varscan/{sample}.varscan'
    shell:
        "{input.java} -jar  {input.varscan} somatic {input.normal_mpileup} {input.tumor_mpileup} {params.basename} --output-vcf 1 --strand-filter 1"

rule manta:
    input:
        python2 = config['python2'], 
        manta = config['manta'], 
        normal_bam = "dna_bam/{sample}_N.bam", 
        tumor_bam = "dna_bam/{sample}_T.bam", 
        ref = config['reference']
    params:
        rundir = 'data_processing/manta/{sample}', 
        workflow = 'data_processing/manta/{sample}/runWorkflow.py'
    threads: 4
    output: 
        candidateSmallIndel = 'data_processing/manta/{sample}/results/variants/candidateSmallIndels.vcf.gz'
    shell:
        "{input.python2} {input.manta} --normalBam {input.normal_bam} --tumorBam {input.tumor_bam} "
        "--referenceFasta {input.ref} --runDir {params.rundir}; "
        "{input.python2} {params.workflow} -m local -j 4 --quiet"


rule mutect:
    input:
        java = config['java7'],
        mutect = config['mutect'],
        reference = config['reference'], 
        cosmic = config['cosmic'], 
        dbsnp = config['dbsnp_mutect'], 
        normal_bam = "dna_bam/{sample}_N.bam", 
        tumor_bam = "dna_bam/{sample}_T.bam", 

    output:
        vcf = 'data_processing/mutect/{sample}.mutect.vcf', 
        out = 'data_processing/mutect/{sample}.mutect.out'
    threads: 1
    shell:
        "{input.java} -jar {input.mutect} --analysis_type MuTect --only_passing_calls "
        "--reference_sequence {input.reference} --cosmic {input.cosmic} --dbsnp {input.dbsnp} "
        "--input_file:normal {input.normal_bam} --input_file:tumor {input.tumor_bam} "
        " --vcf {output.vcf} --out {output.out}"

'''


rule strelka:
    input:
        python2 = 
    output:
    threads: 1
    shell:

rule delly:
    input:
    output:
    threads:
    shell:




rule sequenza:
    input:
    output:
    threads:
    shell:

'''
rule combine:
        input:
            outfiles=expand("data_processing/mutect/{sample}.mutect.vcf", sample=config["samples"], tn=['T', 'N'])
        output:
            "done"
        shell:
            "touch {output}"

