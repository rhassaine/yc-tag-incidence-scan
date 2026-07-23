process SCAN_REGION {
    tag "${sample_id}:${region.length() > 40 ? region.take(40)+'...' : region}"
    container params.samtools_container
    cpus params.threads
    memory '4 GB'

    input:
    tuple val(sample_id), path(bam_or_cram), path(index), val(file_type), val(ref), val(region)

    output:
    tuple val(sample_id), path("${sample_id}.chunk${task.index}.tsv"), emit: stats
    tuple val(sample_id), path("${sample_id}.chunk${task.index}.crash_coords.tsv"), emit: crash_coords

    script:
    def ref_flag = (file_type == 'cram') ? "-T ${ref}" : ''
    def out_name = "${sample_id}.chunk${task.index}.tsv"
    def coords_name = "${sample_id}.chunk${task.index}.crash_coords.tsv"
    """
    set -euo pipefail
    export LC_ALL=C

    # Cheap, count-only passes via samtools' own -e filter expression --
    # this evaluates the tag directly on the decoded record and never
    # formats SEQ/QUAL/CIGAR/other-tags into text at all, unlike piping
    # full SAM text through awk. total_reads and yc_tags_seen never need
    # to see read content, only a count, so -c is the cheapest possible
    # way to get them.
    total_reads=\$(samtools view -c -@ ${params.threads} ${ref_flag} ${bam_or_cram} ${region})
    yc_tags_seen=\$(samtools view -c -@ ${params.threads} ${ref_flag} ${bam_or_cram} ${region} -e '[YC]')

    # The crash-pattern check is the only one that needs record content
    # (to capture coordinates) -- and -e means only the rare matching
    # reads (not all reads) ever get formatted into text at all.
    samtools view -@ ${params.threads} ${ref_flag} ${bam_or_cram} ${region} \
      -e '[YC] =~ "^[0-9]+\\+\\+"' \
      | awk 'BEGIN{OFS="\\t"} {
          yc_tag=""
          for (i=12;i<=NF;i++) { if (\$i ~ /^YC:Z:/) { yc_tag=\$i; break } }
          print "${sample_id}", \$1, \$3, \$4, yc_tag
        }' > ${coords_name}

    crash_total=\$(wc -l < ${coords_name})

    printf "sample_id\\tregion\\ttotal_reads\\tyc_tags_seen\\tcrash_pattern_reads\\n" > ${out_name}
    printf "%s\\t%s\\t%d\\t%d\\t%d\\n" "${sample_id}" "${region}" "\$total_reads" "\$yc_tags_seen" "\$crash_total" >> ${out_name}
    """
}
