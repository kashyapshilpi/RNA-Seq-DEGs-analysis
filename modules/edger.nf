process EDGER {
    label 'process_medium'
    conda 'bioconda::bioconductor-edger'

    input:
    path counts
    path samplesheet

    output:
    path 'edger/*', emit: results

    script:
    """
    FDR=${params.padj} LFC=${params.lfc} edger.R \\
        ${counts} ${samplesheet} edger ${params.contrast_reference} ${params.contrast_target}
    """

    stub:
    """
    mkdir edger
    touch edger/DEGs.tsv edger/all_genes.tsv edger/edger_plots.pdf
    """
}
