files= [expand("rawfastqc/{fastq}_R1_fastqc.html",fastq=config["fastq"]),
        expand("bam/{fastq}_Aligned.sortedByCoord.out.bam",fastq=config["fastq"]),
        "cleaned_counts.txt",
        "multiqc_report.html"]

rule Done:
    input: files

rule fastqc:
    input:
        fastq1 = lambda wildcards: config["raw_fastq"]+"{fastq}"+config["left_suffix"],
        fastq2 = lambda wildcards: config["raw_fastq"]+"{fastq}"+config["right_suffix"]

    output: 
        fastqc1 = "rawfastqc/{fastq}_R1_fastqc.html",
        fastqc2 = "rawfastqc/{fastq}_R2_fastqc.html"
    
    params:
        outfolder= srcdir("rawfastqc/")

    conda: "envs/rnaseq.yaml"
    
    shell:
        """
        mkdir -p {params.outfolder}
        fastqc {input.fastq1} {input.fastq2} -o {params.outfolder}
        """    



rule Trim_adapters:
    input:
        fastq1 = lambda wildcards: config["raw_fastq"]+"{fastq}"+config["left_suffix"],
        fastq2 = lambda wildcards: config["raw_fastq"]+"{fastq}"+config["right_suffix"]

    output:
        trimmed_fastq_R1= "trimmed_fastq/trimmed_{fastq}_R1.fastq.gz",
        trimmed_fastq_R2= "trimmed_fastq/trimmed_{fastq}_R2.fastq.gz"

    conda: "envs/rnaseq.yaml"

    shell: 
        """
        fastp -i {input.fastq1} -I {input.fastq2} -o {output.trimmed_fastq_R1} -O {output.trimmed_fastq_R2} -h trimmed_fastq/{wildcards.fastq}_fastp.html -j trimmed_fastq/{wildcards.fastq}_fastp.json
        """

rule Alignment:
    input:
        trimmed_fastq_R1= "trimmed_fastq/trimmed_{fastq}_R1.fastq.gz",
        trimmed_fastq_R2= "trimmed_fastq/trimmed_{fastq}_R2.fastq.gz"

    output:
        bam_name="bam/{fastq}_Aligned.sortedByCoord.out.bam"

    params:
        prefix= "bam/{fastq}_",
        genome_dir= lambda wildcards: config["genome_dir"]
 
    conda: "envs/rnaseq.yaml"

    shell: 
        """
        STAR --runMode alignReads \
        --genomeDir {params.genome_dir} \
        --readFilesIn {input.trimmed_fastq_R1} {input.trimmed_fastq_R2} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix {params.prefix}
        """

rule Quantification:
    input: 
        bam_files= expand("bam/{fastq}_Aligned.sortedByCoord.out.bam",fastq=config["fastq"])

    output:
        counts="featurecounts/counts_out.txt"

    params:
        GTF= lambda wildcards: config["GTF"]

    conda: "envs/rnaseq.yaml"

    shell:
        """
        featureCounts -a {params.GTF} \
        -o {output.counts} \
         {input.bam_files}
        """

rule Clean_counts_table:
    input:
        counts="featurecounts/counts_out.txt"

    output:
        cleaned_counts= "cleaned_counts.txt"

    shell:
        """
        # removing 2nd column to 6th column
        awk '{{$2=$3=$4=$5=$6=""; print $0}}' {input.counts} > {output.cleaned_counts}
        # removing the first line
        sed -i '1d' {output.cleaned_counts}
        # replacing the alignedbam in column name
        sed -i 's/_Aligned.sortedByCoord.out.bam//g' {output.cleaned_counts}
        # replacing bam/ in the column name
        sed -i 's+bam/++g' {output.cleaned_counts}
        """
        

rule Quality_control:
    input:
        cleaned="cleaned_counts.txt"

    output:
        report= "multiqc_report.html"

    params:
        fastqc_dir= srcdir("rawfastqc"),
        fastp_dir= srcdir("trimmed_fastq"),
        bam_dir= srcdir("bam"),
        subreads_dir= srcdir("featurecounts")

    conda: "envs/rnaseq.yaml"
    
    shell:
        """
        multiqc \
            --dirs {params.fastqc_dir} {params.fastp_dir} {params.bam_dir} {params.subreads_dir}
        """
