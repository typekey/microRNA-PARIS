process PLOT_PEAK_ANNOTATION_CHIPSEEKER {
    label 'process_medium'

    input:
    path(anno_chipseeker)

    output:
    path("*.pdf"), emit: pdf

    script:
    def bin         = params.bin

    """
    ${bin}/plot/plot_miRNA_target_number.py --input ./ --output_name miRNA_target_number

    labels=()
    tsv_files=()
    for tsv in *.tsv; do
        base_name=\$(basename "\$tsv" .tsv)
        mirna=\$(echo \${base_name} | cut -d'_' -f4 | cut -d'.' -f1)
        group=\$(echo \${base_name} | cut -d'_' -f2)
        label="\${mirna} (\${group})"
        labels+=("\${label}")
        tsv_files+=("\${tsv}")
    done;

    #echo \${labels[@]}
    #echo \${tsv_files[@]}

    plot_stack_bar.py \
        --file_list \${tsv_files[@]} \
        --labels "\${labels[@]}" \
        --plot_column 'peak_annotation' \
        --mode 'value' \
        --output_file 'all_peak_annotation_value.pdf' \

    plot_stack_bar.py \
        --file_list \${tsv_files[@]} \
        --labels "\${labels[@]}" \
        --plot_column 'peak_annotation' \
        --output_file 'all_peak_annotation_percent.pdf' \
        
    plot_stack_bar.py \
        --file_list ./*.anno.chipseeker.tsv \
        --file_list \${tsv_files[@]} \
        --labels "\${labels[@]}" \
        --plot_column 'peak_annotation2' \
        --mode 'value' \
        --output_file 'all_peak_annotation_detail_value.pdf' \

    plot_stack_bar.py \
        --file_list ./*.anno.chipseeker.tsv \
        --file_list \${tsv_files[@]} \
        --labels "\${labels[@]}" \
        --plot_column 'peak_annotation2' \
        --output_file 'all_peak_annotation_detail_percent.pdf' \

    # plot_pie.py \
    #    --input_file ./*.anno.chipseeker.tsv \
    #    --plot_column 'peak_annotation' \
    #    --output_file 'peak_annotation_pie.pdf' \

    """
}
