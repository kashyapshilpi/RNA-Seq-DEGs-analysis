process DESEQ2 {
    label 'process_medium'
    conda 'bioconda::bioconductor-deseq2 bioconda::bioconductor-apeglm conda-forge::r-ggplot2 conda-forge::r-pheatmap'

    input:
    path counts
    path samplesheet

    output:
    path 'deseq2/*', emit: results

    script:
    """
    PADJ=${params.padj} LFC=${params.lfc} deseq2.R \\
        ${counts} ${samplesheet} deseq2 ${params.contrast_reference} ${params.contrast_target}
    """

    stub:
    """
    mkdir deseq2
    touch deseq2/DEGs.tsv deseq2/all_genes.tsv deseq2/deseq2_plots.pdf
    """
}
