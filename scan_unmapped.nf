process SCAN_UNMAPPED {
    tag "${sample_id}:unmapped"
    cpus params.threads
    memory '4 GB'

    input:
    tuple val(sample_id), path(bam_or_cram), path(index), val(file_type), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}.unmapped.tsv"), emit: stats

    script:
    def ref_flag = (file_type == 'cram') ? "-T ${ref}" : ''
    """
    set -euo pipefail
    export LC_ALL=C

    samtools view -@ ${params.threads} -f 4 ${ref_flag} ${bam_or_cram} | awk -v softclip_min=${params.softclip_min} '
    {
      total++
      is_softclip = 0
      if (\$6 ~ /^[0-9]+S/) {
        n = \$6 + 0
        if (n >= softclip_min) { is_softclip = 1; softclip_total++ }
      }
      for (i=12; i<=NF; i++) {
        if (\$i ~ /^YC:Z:/) {
          yc_total++
          if (\$i ~ /^YC:Z:[0-9]+\\+\\+/) {
            crash_total++
            if (is_softclip) crash_and_softclip++
          }
        }
      }
    }
    END {
      printf "sample_id\\tregion\\ttotal_reads\\tyc_tags_seen\\tsoftclip_reads\\tcrash_pattern_reads\\tcrash_and_softclip\\n"
      printf "%s\\t%s\\t%d\\t%d\\t%d\\t%d\\t%d\\n", "${sample_id}", "unmapped", total+0, yc_total+0, softclip_total+0, crash_total+0, crash_and_softclip+0
    }' > ${sample_id}.unmapped.tsv
    """
}
