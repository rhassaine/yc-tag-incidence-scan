process SCAN_REGION {
    tag "${sample_id}:${region.length() > 40 ? region.take(40)+'...' : region}"
    container params.samtools_container
    cpus params.threads
    memory '4 GB'

    input:
    tuple val(sample_id), path(bam_or_cram), path(index), val(file_type), val(region)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}.chunk${task.index}.tsv"), emit: stats
    tuple val(sample_id), path("${sample_id}.chunk${task.index}.crash_coords.tsv"), emit: crash_coords

    script:
    def ref_flag = ref ? "-T ${ref}" : ''
    def out_name = "${sample_id}.chunk${task.index}.tsv"
    def coords_name = "${sample_id}.chunk${task.index}.crash_coords.tsv"
    // Detect a true split-window (single token, "contig:digits-digits") as
    // opposed to a bare contig name or a batched multi-contig line -- same
    // detection already proven in transform_region.nf/filter_region.nf.
    def is_split_window = (region ==~ /^[^\s]+:[0-9]+-[0-9]+$/)
    def wstart = is_split_window ? region.tokenize(':')[1].tokenize('-')[0] : '1'
    def wend   = is_split_window ? region.tokenize(':')[1].tokenize('-')[1] : '999999999999'
    """
    set -euo pipefail
    export LC_ALL=C

    # Single decode pass instead of three separate samtools invocations
    # (total count, YC-tag count, crash-pattern extraction) each
    # re-decoding the same region from scratch -- same redundant-decode
    # fix already applied to SCAN_UNMAPPED. Also adds a POS restriction
    # (a no-op wide-bounds filter for non-split-window chunks): samtools
    # region queries return any read that OVERLAPS the interval, not
    # just reads that start inside it, so without this, a read spanning
    # a window boundary gets counted -- and can appear in
    # crash_coords.tsv -- twice, once by each neighboring chunk.
    # Confirmed as a real, measurable effect (a ~6,700-read discrepancy
    # against the transform pipeline's own accounting), not hypothetical.
    cat > onepass.awk << 'AWK_EOF'
BEGIN{OFS="\t"}
{
  if (\$4 < s || \$4 > e) next
  total++
  yc_tag=""
  for (i=12;i<=NF;i++) { if (\$i ~ /^YC:Z:/) { yc_tag=\$i; break } }
  if (yc_tag != "") {
    yc_seen++
    if (yc_tag ~ /^YC:Z:[0-9]+\\+\\+/) {
      crash++
      print sample, \$1, \$3, \$4, yc_tag >> crashfile
    }
  }
}
END {
  print total+0, yc_seen+0, crash+0 > countfile
}
AWK_EOF

    samtools view -@ ${params.threads} ${ref_flag} ${bam_or_cram} ${region} \
      | awk -v s="${wstart}" -v e="${wend}" -v sample="${sample_id}" -v crashfile="${coords_name}" -v countfile="counts.tmp" -f onepass.awk

    # crash_coords.tsv may not have been created at all if there were
    # zero matches (awk's >> only creates the file on first write) --
    # make sure it exists either way, matching the original script's
    # always-present-file behavior
    touch "${coords_name}"

    read total_reads yc_tags_seen crash_total < counts.tmp

    printf "sample_id\\tregion\\ttotal_reads\\tyc_tags_seen\\tcrash_pattern_reads\\n" > ${out_name}
    printf "%s\\t%s\\t%d\\t%d\\t%d\\n" "${sample_id}" "${region}" "\$total_reads" "\$yc_tags_seen" "\$crash_total" >> ${out_name}
    """
}
