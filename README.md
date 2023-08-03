# Differential Expression Analysis of RNA-Seq Data

Welcome to the Differential Expression Analysis of RNA-Seq Data repository! This repository contains scripts and documentation for identifying Differentially Expressed Genes (DEGs) from RNA-Seq data. The analysis pipeline includes preprocessing, quality control, alignment, read counting, and differential expression analysis using popular tools such as HISAT2, HTSeq-Count, and DESeq2.

## Table of Contents

1. [Introduction](#introduction)
2. [Dependencies](#dependencies)
3. [Data](#data)
4. [Analysis Workflow](#analysis-workflow)
5. [Usage](#usage)
6. [Results](#results)
7. [License](#license)
8. [Contributing](#contributing)

## Introduction

Differential Expression Analysis is a powerful method for identifying genes that are differentially expressed between different experimental conditions. This analysis provides insights into gene regulation and can lead to the discovery of genes that play a crucial role in specific biological processes.

## Dependencies

Before running the Differential Expression Analysis pipeline, ensure you have the following tools and software installed on your system:

- HISAT2 (version 2.1.0): A fast and sensitive alignment tool for mapping RNA-Seq reads to a reference genome.
- HTSeq-Count (version 2.0.3): A tool for counting aligned reads in specific genomic features.
- DESeq2 (version 1.41.6): A powerful R package for differential gene expression analysis.
- edgeR (version 3.14.0): A powerful R package for differential gene expression analysis.

Make sure to update the versions with the appropriate ones you are using.

## Data

The raw RNA-Seq data used in this analysis is not included in this repository. Ensure you have access to the data and place it in the appropriate input directory before running the pipeline.

## Analysis Workflow

The analysis pipeline can be divided into several key steps:

1. **Quality Control**: Assess the quality of raw sequencing data using tools like FastQC, Trimmomatic.
2. **Read Alignment**: Align the processed reads to a reference genome using HISAT2.
3. **Read Counting**: Count the number of reads mapped to each gene using HTSeq-Count.
4. **Differential Expression Analysis**: Use DESeq2 or edgeR to identify genes that are differentially expressed between experimental conditions.

## Usage

To identify DEGs from RNA-Seq data, follow these steps:

1. Clone this repository to your local machine: `git clone https://github.com/kashyapshilpi/ RNA-Seq-DEGs-analysis.git`
2. Place the raw RNA-Seq data in the appropriate input directory.
3. Install the required dependencies listed in the "Dependencies" section.
4. Modify the configuration file if necessary to specify parameters and options.
5. Execute the analysis script: `bash ngs.sh`, `Rscript deseq2.R`, `Rscript edger.R` 
6. The results, including lists of DEGs and visualization plots, will be generated in the output directory.

## Results

The results of the Differential Expression Analysis will be stored in the output directory. This will include lists of differentially expressed genes, statistical significance values, and visualization plots depicting expression changes.

## License

This project is licensed under the [MIT License](license.txt).

## Contributing

If you wish to contribute to this project, feel free to open issues, submit pull requests, or suggest improvements. We welcome your contributions!

Thank you for using this Differential Expression Analysis repository. If you have any questions or encounter any issues, please don't hesitate to contact us.

Happy analyzing!
