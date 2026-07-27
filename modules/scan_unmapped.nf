process SCAN_UNMAPPED {
    tag "${sample_id}:unmapped"
    container params.samtools_container
    cpus params.threads
    memory '4 GB'

    input:
    tuple val(sample_id), path(bam_or_cram), path(index), val(file_type)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}.unmapped.tsv"), emit: stats
    tuple val(sample_id), path("${sample_id}.unmapped.crash_coords.tsv"), emit: crash_coords

    script:
    def ref_flag = ref ? "-T ${ref}" : ''
    """
    set -euo pipefail
    export LC_ALL=C

    total_reads=\$(samtools view -c -@ ${params.threads} -f 4 ${ref_flag} ${bam_or_cram})
    yc_tags_seen=\$(samtools view -c -@ ${params.threads} -f 4 ${ref_flag} ${bam_or_cram} -e '[YC]')

    samtools view -@ ${params.threads} -f 4 ${ref_flag} ${bam_or_cram} \
      -e '[YC] =~ "^[0-9]+\\+\\+"' \
      | awk 'BEGIN{OFS="\\t"} {
          yc_tag=""
          for (i=12;i<=NF;i++) { if (\$i ~ /^YC:Z:/) { yc_tag=\$i; break } }
          print "${sample_id}", \$1, "unmapped", "NA", yc_tag
        }' > ${sample_id}.unmapped.crash_coords.tsv

    crash_total=\$(wc -l < ${sample_id}.unmapped.crash_coords.tsv)

    printf "sample_id\\tregion\\ttotal_reads\\tyc_tags_seen\\tcrash_pattern_reads\\n" > ${sample_id}.unmapped.tsv
    printf "%s\\t%s\\t%d\\t%d\\t%d\\n" "${sample_id}" "unmapped" "\$total_reads" "\$yc_tags_seen" "\$crash_total" >> ${sample_id}.unmapped.tsv
    """
}
