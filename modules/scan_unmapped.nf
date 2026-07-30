process SCAN_UNMAPPED {
    tag "${sample_id}:unmapped"
    container params.samtools_container
    cpus 8
    memory '16 GB'
    machineType 'n2-*'

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

    # Single decode pass instead of three separate samtools invocations
    # (total count, YC-tag count, crash-pattern extraction) each
    # re-decoding the same unmapped-read block from scratch -- this was
    # a real, confirmed bottleneck (a single samtools invocation observed
    # pegged at 196% CPU for over an hour on just one of the three
    # passes). One awk program now computes all three metrics from one
    # decode.
    cat > onepass.awk << 'AWK_EOF'
BEGIN{OFS="\t"}
{
  total++
  yc_tag=""
  for (i=12;i<=NF;i++) { if (\$i ~ /^YC:Z:/) { yc_tag=\$i; break } }
  if (yc_tag != "") {
    yc_seen++
    if (yc_tag ~ /^YC:Z:[0-9]+\\+\\+/) {
      crash++
      print sample, \$1, "unmapped", "NA", yc_tag >> crashfile
    }
  }
}
END {
  print total+0, yc_seen+0, crash+0 > countfile
}
AWK_EOF

    samtools view -@ ${task.cpus} -f 4 ${ref_flag} ${bam_or_cram} \
      | awk -v sample="${sample_id}" -v crashfile="${sample_id}.unmapped.crash_coords.tsv" -v countfile="counts.tmp" -f onepass.awk

    # crash_coords.tsv may not have been created at all if there were
    # zero matches (awk's >> only creates the file on first write) --
    # make sure it exists either way, matching the original script's
    # always-present-file behavior
    touch "${sample_id}.unmapped.crash_coords.tsv"

    read total_reads yc_tags_seen crash_total < counts.tmp

    printf "sample_id\tregion\ttotal_reads\tyc_tags_seen\tcrash_pattern_reads\n" > ${sample_id}.unmapped.tsv
    printf "%s\t%s\t%d\t%d\t%d\n" "${sample_id}" "unmapped" "\$total_reads" "\$yc_tags_seen" "\$crash_total" >> ${sample_id}.unmapped.tsv
    """
}
