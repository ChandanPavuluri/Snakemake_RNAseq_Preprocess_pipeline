files= [expand("fastq/{fastq}_R1_fastqc.html",fastq=config["fastq"]),
        expand("trimmed_{fastq}_R1.fastq.gz",fastq=config["fastq"]),
        expand("bam/{fastq}_Aligned.sortedByCoord.out.bam",fastq=config["fastq"])]

rule Done:
    input: files

rule fastqc:
    input:
        fastq1 = "/d/tmp/rechp/snakemakernaseq/fastq/{fastq}_R1.fastq.gz",
        fastq2 = "/d/tmp/rechp/snakemakernaseq/fastq/{fastq}_R2.fastq.gz"
    output: 
        fastqc1 = "fastq/{fastq}_R1_fastqc.html",
        fastqc2 = "fastq/{fastq}_R2_fastqc.html"
    conda: "envs/rnaseq.yaml"
    shell:
        """
        fastqc {input.fastq1} {input.fastq2}
        """    


rule Trim_adapters:
    input:
        fastq1 = "/d/tmp/rechp/snakemakernaseq/fastq/{fastq}_R1.fastq.gz",
        fastq2 = "/d/tmp/rechp/snakemakernaseq/fastq/{fastq}_R2.fastq.gz"
    output:
        trimmed_fastq_R1= "trimmed_{fastq}_R1.fastq.gz",
        trimmed_fastq_R2= "trimmed_{fastq}_R2.fastq.gz"
    conda: "envs/rnaseq.yaml"
    shell: 
        """
        fastp -i {input.fastq1} -I {input.fastq2} -o {output.trimmed_fastq_R1} -O {output.trimmed_fastq_R2} -h trimmed_fastq/{wildcards.fastq}_fastp.html -j trimmed_fastq/{wildcards.fastq}_fastp.json
        """
    
    
rule Alignment:
    input:
        trimmed_fastq_R1= "trimmed_{fastq}_R1.fastq.gz",
        trimmed_fastq_R2= "trimmed_{fastq}_R2.fastq.gz"
    output:
        bam_name="bam/{fastq}_Aligned.sortedByCoord.out.bam"
    conda: "envs/rnaseq.yaml"
    shell: 
        """
        STAR --runMode alignReads \
        --genomeDir /proj/edith/regeps/regep00/studies/COPDGene/analyses/rechp/RNAseq/GRCh38_STAR/ \
        --readFilesIn {input.trimmed_fastq_R1} {input.trimmed_fastq_R2} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix bam/{wildcards.fastq}_
        """
