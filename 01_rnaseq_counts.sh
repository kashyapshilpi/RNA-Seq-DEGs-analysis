#!/usr/bin/env bash
# RNA-seq preprocessing: QC -> trimming -> HISAT2 alignment -> HTSeq gene counts
# -> merged count matrix for DESeq2 / edgeR.
#
# Usage:
#   bash 01_rnaseq_counts.sh samples.txt
#
# samples.txt holds one sample ID per line. Paired reads are expected at
#   ${READS_DIR}/<ID>/<ID>_1.fastq.gz and ${READS_DIR}/<ID>/<ID>_2.fastq.gz
#
# Settings can be overridden with environment variables, e.g.
#   THREADS=16 STRAND=reverse bash 01_rnaseq_counts.sh samples.txt

set -euo pipefail

READS_DIR=${READS_DIR:-data/raw}
GENOME_FA=${GENOME_FA:-ref/genome.fna}
GTF=${GTF:-ref/genes.gtf}
OUT_DIR=${OUT_DIR:-results}
ADAPTERS=${ADAPTERS:-TruSeq3-PE.fa}
TRIMMOMATIC=${TRIMMOMATIC:-trimmomatic}  # or e.g. "java -jar /path/to/trimmomatic-0.39.jar"
THREADS=${THREADS:-8}
STRAND=${STRAND:-no}                     # library strandedness for HTSeq: no | yes | reverse

samples_file=${1:?Usage: bash 01_rnaseq_counts.sh samples.txt}
mapfile -t samples < <(grep -v '^[[:space:]]*$' "$samples_file")

mkdir -p "$OUT_DIR"/{qc_raw,trimmed,qc_trimmed,index,bam,counts}
INDEX="$OUT_DIR/index/genome"

# 1. QC of raw reads
fastqc -t "$THREADS" -o "$OUT_DIR/qc_raw" "$READS_DIR"/*/*.fastq.gz

# 2. Build the HISAT2 index (only needed once)
if [[ ! -f "${INDEX}.1.ht2" ]]; then
    hisat2-build -p "$THREADS" "$GENOME_FA" "$INDEX"
fi

for s in "${samples[@]}"; do
    echo "=== $s ==="
    r1="$READS_DIR/$s/${s}_1.fastq.gz"
    r2="$READS_DIR/$s/${s}_2.fastq.gz"
    t="$OUT_DIR/trimmed/$s"

    # 3. Adapter and quality trimming
    $TRIMMOMATIC PE -threads "$THREADS" "$r1" "$r2" \
        "${t}_R1.paired.fq.gz" "${t}_R1.unpaired.fq.gz" \
        "${t}_R2.paired.fq.gz" "${t}_R2.unpaired.fq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    fastqc -t "$THREADS" -o "$OUT_DIR/qc_trimmed" "${t}_R1.paired.fq.gz" "${t}_R2.paired.fq.gz"

    # 4. Align and coordinate-sort (HISAT2 reads gzipped FASTQ directly)
    hisat2 -p "$THREADS" -x "$INDEX" \
        -1 "${t}_R1.paired.fq.gz" -2 "${t}_R2.paired.fq.gz" \
        --summary-file "$OUT_DIR/bam/${s}.hisat2.log" \
        | samtools sort -@ "$THREADS" -o "$OUT_DIR/bam/${s}.bam" -
    samtools index "$OUT_DIR/bam/${s}.bam"

    # 5. Gene-level counts (requires HTSeq >= 2.0)
    htseq-count -r pos -s "$STRAND" "$OUT_DIR/bam/${s}.bam" "$GTF" \
        > "$OUT_DIR/counts/${s}.counts.txt"
done

# 6. Merge per-sample counts into one matrix, dropping HTSeq's __summary rows
files=()
for s in "${samples[@]}"; do files+=("$OUT_DIR/counts/${s}.counts.txt"); done
{
    printf "gene_id"; printf "\t%s" "${samples[@]}"; printf "\n"
    paste "${files[@]}" | awk -F'\t' '$1 !~ /^__/ {
        printf "%s", $1
        for (i = 2; i <= NF; i += 2) printf "\t%s", $i
        printf "\n"
    }'
} > "$OUT_DIR/counts/count_matrix.tsv"

# 7. One summary report across all samples
multiqc -f -o "$OUT_DIR/multiqc" "$OUT_DIR"

echo "Done. Count matrix: $OUT_DIR/counts/count_matrix.tsv"
