# Snakemake_RNAseq_Preprocess_pipeline

# Chandan Pavuluri

Basic snakemake workflow for preprocessing the bulk RNAseq data from fastq files to counts table which can be used as an input for Differential expression analysis.
The workflow can be executed just by installing [snakemake](https://snakemake.readthedocs.io/en/stable/) and [conda](https://docs.conda.io/projects/conda/en/latest/index.html), required dependencies gets automatically installed after executing snakemake.
This workflow can be implemented easily just by cloning the repository and making the necessary changes to the config file without changing anything in the Snakefile.


![rulegraph](https://user-images.githubusercontent.com/77353407/167415070-e265c78d-3e4b-4b90-b0d5-e8ba12a98589.png)

Tools:

fatsqc: [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

Trim_adapters: [fastp](https://github.com/OpenGene/fastp)

Alignment: [STAR](https://github.com/alexdobin/STAR)

Quantification: [Subread/featureCounts](http://subread.sourceforge.net/)

Quality_control: [MultiQC](https://multiqc.info/)
