# RNA-seq Differential Expression Analysis

Paired-end RNA-seq workflow from raw reads to differentially expressed genes, with DESeq2 as the primary method and edgeR as a cross-check.

> Developed during my MSc in Bioinformatics. Kept as a reference implementation; not actively maintained.

## Workflow

```
FASTQ ─► FastQC ─► Trimmomatic ─► HISAT2 ─► samtools sort ─► HTSeq-count
      ─► count_matrix.tsv ─┬─► DESeq2 (apeglm-shrunken fold changes)
                           └─► edgeR (quasi-likelihood F-test)
```

## Repository contents

| File | Purpose |
|---|---|
| `01_rnaseq_counts.sh` | QC, trimming, alignment, gene counting, merged count matrix, MultiQC report |
| `02_deseq2.R` | DESeq2 differential expression, MA / PCA / distance heatmap / volcano plots |
| `03_edger.R` | edgeR differential expression, MDS / BCV / volcano plots |
| `samplesheet.example.csv` | Example sample sheet (`sample_id,condition`) |

## Requirements

- FastQC, MultiQC, Trimmomatic 0.39, HISAT2, samtools, HTSeq ≥ 2.0
- R ≥ 4.1 with `DESeq2`, `apeglm`, `edgeR`, `ggplot2`, `pheatmap`

Example install with conda:

```bash
conda create -n rnaseq -c bioconda -c conda-forge \
    fastqc multiqc trimmomatic hisat2 samtools htseq \
    bioconductor-deseq2 bioconductor-apeglm bioconductor-edger r-ggplot2 r-pheatmap
```

## Usage

**Step 1: counts.** List sample IDs in `samples.txt`, with reads at `data/raw/<ID>/<ID>_1.fastq.gz` and `_2.fastq.gz`. Set `STRAND` to match your library prep (`no`, `yes` or `reverse`; most Illumina TruSeq stranded kits are `reverse`):

```bash
GENOME_FA=ref/genome.fna GTF=ref/genes.gtf STRAND=reverse THREADS=16 \
    bash 01_rnaseq_counts.sh samples.txt
```

**Step 2: differential expression.** Fill in a sample sheet with `control` / `treatment` conditions, then:

```bash
Rscript 02_deseq2.R results/counts/count_matrix.tsv samplesheet.csv results/deseq2
Rscript 03_edger.R  results/counts/count_matrix.tsv samplesheet.csv results/edger
```

Thresholds default to adjusted p < 0.05 and |log2FC| > 1. Change them with environment variables, e.g. `PADJ=0.01 LFC=2 Rscript 02_deseq2.R ...`.

## Outputs

| File | Contents |
|---|---|
| `results/counts/count_matrix.tsv` | Raw gene counts, one column per sample |
| `results/deseq2/DEGs.tsv`, `DEGs_up.tsv`, `DEGs_down.tsv` | Significant genes |
| `results/deseq2/all_genes.tsv` | Full results table |
| `results/deseq2/deseq2_plots.pdf` | Dispersion, MA, PCA, sample-distance heatmap, volcano |
| `results/edger/...` | Same tables plus MDS, BCV and volcano plots |
| `results/multiqc/` | Combined QC report |

## Data

Raw data are not included. The original analysis used public RNA-seq data from the tea plant (*Camellia sinensis*, NCBI assembly GCF_004153795.1).

## License

MIT. See [licence.txt](licence.txt).
