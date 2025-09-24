process PLOT_PEAK_NUMBER {
    label 'process_medium'

    input:
    path(bed)

    output:
    path("*.png"), emit: png

    script:
    def bin         = params.bin

    """
    ${bin}/plot/plot_miRNA_target_number.py --input ./ --output_name miRNA_target_number
    echo "peak_number" > test.png
    echo "peak_number" > test.pdf

    """
}
