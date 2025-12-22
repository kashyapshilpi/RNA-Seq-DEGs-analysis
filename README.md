# Differential Expression Analysis of RNA-Seq Data

## Project Context

This repository represents an early RNA-seq differential expression
analysis pipeline developed during my MSc training in Bioinformatics.
The project helped me build a strong foundation in transcriptomic
data processing, statistical analysis using DESeq2/edgeR, and
reproducible scripting with R and Bash.

This repository is intended as a learning and reference implementation
and is not actively maintained.

---

## Overview

This repository contains scripts and documentation for identifying
Differentially Expressed Genes (DEGs) from RNA-Seq data. The analysis
pipeline includes preprocessing, quality control, alignment, read
counting, and differential expression analysis using widely adopted
tools such as HISAT2, HTSeq-Count, DESeq2, and edgeR.

---

## Table of Contents

1. Introduction  
2. Dependencies  
3. Data  
4. Analysis Workflow  
5. Usage  
6. Results  
7. License  
8. Contributing  

---

## Introduction

Differential expression analysis is a key approach for identifying genes
that show statistically significant expression changes between
experimental conditions. Such analyses provide insights into gene
regulation and underlying biological mechanisms.

---

## Dependencies

The following tools and packages are required:

- **HISAT2 (v2.1.0)** – RNA-seq read alignment  
- **HTSeq-Count (v2.0.3)** – Gene-level read counting  
- **DESeq2 (v1.41.6)** – Differential expression analysis  
- **edgeR (v3.14.0)** – Differential expression analysis  

Please adjust versions as needed based on your environment.

---

## Data

Raw RNA-seq data are not included in this repository. Users should place
their input data in the appropriate directory structure prior to
execution.

---

## Analysis Workflow

1. **Quality Control** – FastQC, Trimmomatic  
2. **Read Alignment** – HISAT2  
3. **Read Counting** – HTSeq-Count  
4. **Differential Expression Analysis** – DESeq2 / edgeR  

---

## Usage

```bash
git clone https://github.com/kashyapshilpi/RNA-Seq-DEGs-analysis.git


## Results

The results of the Differential Expression Analysis will be stored in the output directory. This will include lists of differentially expressed genes, statistical significance values, and visualization plots depicting expression changes.

## License

This project is licensed under the [MIT License](license.txt).

## Contributing

If you wish to contribute to this project, feel free to open issues, submit pull requests, or suggest improvements. We welcome your contributions!

Thank you for using this Differential Expression Analysis repository. If you have any questions or encounter any issues, please don't hesitate to contact us.

Happy analyzing!
