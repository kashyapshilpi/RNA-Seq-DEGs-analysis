process MERGE_COUNTS {
    // Why: Join per-sample counts into one gene x sample table for DESeq2 and edgeR.
    label 'process_low'
    conda 'conda-forge::python=3.12'

    input:
    path counts

    output:
    path 'count_matrix.tsv', emit: matrix

    script:
    """
    merge_counts.py ${counts} > count_matrix.tsv
    """

    stub:
    """
    merge_counts.py ${counts} > count_matrix.tsv
    """
}
