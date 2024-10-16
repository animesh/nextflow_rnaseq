# RNA-seq analysis pipeline

## Description
A nextflow pipeline for the analysis of _messenger_ RNA-seq data. It processes paired-end mRNA-seq fastq files and delivers raw count tables. This pipeline also outputs a QC report per fastq file and a `.bam` mapping file to use with a genome browser for instance.

## Installation
https://www.nextflow.io/blog/2021/setup-nextflow-on-windows.html
```
curl -s https://get.nextflow.io | bash
./nextflow run hello
./nextflow run nf-core/rnaseq -profile docker,test --outdir test
./nextflow run nf-core/rnaseq --max_memory '6.GB' --max_cpus 2  --input config/samples.csv --outdir results --genome GRCh37 -profile docker
rm 
## Custom genome

https://nf-co.re/docs/usage/reference_genomes

    3  wget http://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
    5  wget http://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz



rm -rf work/ result/ .nextflow*
./nextflow run nf-core/rnaseq --max_memory '6.GB' --max_cpus 2  --input config/samples.csv --outdir results --gtf Homo_sapiens.GRCh38.112.gtf.gz --fasta Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz -profile docker 

sudo touch nohup.out
sudo nohup dockerd &
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz
sudo rm -rf work/ result/ .nextflow*

wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
tar -xvf HBR_UHR_ERCC_ds_5pc.tar
 9  ls -ltrh test/salmon/deseq2_qc/
   10  ls -ltrh test/salmon/
   11  ls -ltrh test/salmon/salmon.merged.gene_counts_length_scaled.tsv
   12  less  test/salmon/salmon.merged.gene_counts_length_scaled.tsv
   13  ls
   14  ls work/
   15  find . -iname "sample*csv"
   16  find . -iname "*gtf*"
   17  find . -iname "*gtf*" | less
   18  find . -iname "*gtf" | less
   19  find . -iname "*genome" | less
   20  find . -iname "*genome*fa*" | less
   21  hsi
   22  history 
   23  wget http://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
   24  wget http://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
   25  rm -rf work/ result/ .nextflow*
   26  ls
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
tar -xvf HBR_UHR_ERCC_ds_5pc.tar
echo "sample,fastq_1,fastq_2,strandedness" > samples.csv
rm c?
ls -1 *1.fastq.gz | awk -F "_" '{print $1 $2}' > c0
ls -1 $PWD/*1.fastq.gz > c1
ls -1 $PWD/*2.fastq.gz > c2
printf 'auto\n%.0s' {1..`ls *1.fastq.gz`} > c3
paste -d "," c? >> samples.csv
cat samples.csv
curl -s https://get.nextflow.io | bash
./nextflow run hello
./nextflow run nf-core/rnaseq -profile docker,test --outdir test
./nextflow run nf-core/rnaseq --max_memory '63.GB' --max_cpus 16  --input samples.csv --outdir results --genome GRCh37 -profile docker


sudo ./nextflow run nf-core/rnaseq --max_memory '16.GB' --max_cpus 4  --input config/samples.csv --outdir test --gtf GCF_000001405.40_GRCh38.p14_genomic.gtf.gz --fasta GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz -profile dock
er -resume
```

## Input files
* __RNA-seq fastq files__ as listed in the `config/samples.tsv` file. Specify a sample name (e.g. "Sample_A") in the `sample` column and the paths to the forward read (`fq1`) and to the reverse read (`fq2`).
* __A genomic reference in FASTA format__. For instance, a fasta file containing the chromosome 1 of human genome (*Homo sapiens*).
* __A genome annotation file in the [GTF format](https://useast.ensembl.org/info/website/upload/gff.html)__. You can convert a GFF annotation file format into GTF with the [gffread program from Cufflinks](http://ccb.jhu.edu/software/stringtie/gff.shtml): `gffread my.gff3 -T -o my.gtf`. :warning: for featureCounts to work, the _feature_ in the GTF file should be `exon` while the _meta-feature_ has to be `transcript_id`. 


# Installation and usage (local machine)

## Installation

You will need a local copy of the GitHub `snakemake_rnaseq` repository on your machine. You can either:
- use git in the shell: `git clone git@github.com:lyaqing/snakemake_rnaseq.git`.
- click on ["Clone or download"](https://github.com/lyaqing/snakemake_rnaseq/archive/master.zip) and select `download`.
- Then navigate inside the `snakemake_rnaseq` folder using Shell commands.

## Usage 

## Configuration :pencil2:
You'll need to change a few things to accomodate this pipeline to your needs. Make sure you have changed the parameters in the `config/config.yaml` file that specifies where to find the sample data file, the genomic and transcriptomic reference fasta files to use and the parameters for certains rules etc.    
This file is used so the `Snakefile` does not need to be changed when locations or parameters need to be changed.

### :round_pushpin: Option 1: mamba/conda (easiest)
Mamba is an extremely fast and robust replacement for the Conda package manager which is highly recommended. The default conda solver is a bit slow and sometimes has issues with selecting the latest package releases. Therefore, we recommend to in any case use Mamba.
Using the mamba package manager, you need to create an environment where core softwares such as `Snakemake` will be installed.
1. Install the [Miniconda3 distribution (>= Python 3.7 version)](https://docs.conda.io/en/latest/miniconda.html) for your OS (Windows, Linux or Mac OS X).
2. Install Mamba into any other Conda-based Python distribution with `conda install -n base -c conda-forge mamba`.
3. Inside a Shell window (command line interface), create a virtual environment named `rnaseq` using the `workflow/envs/environment.yaml` file with the following command: `mamba env create --name rnaseq --file workflow/envs/environment.yaml`
4. Then, before you run the Snakemake pipeline, activate this virtual environment with `mamba activate rnaseq`.


### :whale: Option 2: Docker/Singularity
1. Install Docker desktop for your operating system.


## Dry run
- With conda: use the `snakemake -np` to perform a dry run that prints out the rules and commands.
- With Docker: use the `docker run ` 

## Real run
With conda: `snakemake --cores 10`


# References :green_book:

## Authors
- Yaqing Liu, yaqing.liu@outlook.com
- Duo Wang, 18801232285@163.com


## Pipeline dependencies
#### Workflow management
* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
#### Preprocess and quality control
* [fastp](https://github.com/OpenGene/fastp)
#### Alignment
* [STAR](https://github.com/alexdobin/STAR)
* [subread](http://subread.sourceforge.net/)
* [HISAT2](https://daehwankimlab.github.io/hisat2/)
#### Quantification
* [STAR](https://github.com/alexdobin/STAR)
* [StringTie](https://ccb.jhu.edu/software/stringtie/)
* [featureCounts (subread)](http://subread.sourceforge.net/)
* [htseq](https://htseq.readthedocs.io/en/latest/)
* [salmon](https://github.com/COMBINE-lab/salmon)
* [sailfish](https://sailfish.readthedocs.io/en/master/sailfish.html)
* [kallisto](https://github.com/pachterlab/kallisto)
* [RSME](https://github.com/deweylab/RSEM) with [Bowtie](https://bowtie-bio.sourceforge.net/index.shtml)
#### Normalization/DEGs analysis/...
Developing...

# Citation
If you use this software, please use the following citation:
Wang D., Liu Y., et al. A real-world multi-center RNA-seq benchmarking study using the Quartet
and MAQC reference materials. *bioRxiv* 2023.12.09.570956; doi: https://doi.org/10.1101/2023.12.09.570956


sample,fastq_1,fastq_2,strandedness
UHR2,UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz,UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz,auto


