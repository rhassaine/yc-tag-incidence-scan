process SCAN_UNMAPPED {
    tag "${sample_id}:unmapped"
    container params.samtools_container
    cpus params.threads
    memory '4 GB'

    input:
    tuple val(sample_id), path(bam_or_cram), path(index), val(file_type), val(ref)

    output:
    tuple val(sample_id), path("${sample_id}.unmapped.tsv"), emit: stats

    script:
    def ref_flag = (file_type == 'cram') ? "-T ${ref}" : ''
    """
    set -euo pipefail
    export LC_ALL=C

    samtools view -@ ${params.threads} -f 4 ${ref_flag} ${bam_or_cram} | awk '
    {
      total++
      for (i=12; i<=NF; i++) {
        if ($i ~ /^YC:Z:/) {
          yc_total++
          if ($i ~ /^YC:Z:[0-9]+\+\+/) {
            crash_total++
          }
        }
      }
    }
    END {
      printf "sample_id\tregion\ttotal_reads\tyc_tags_seen\tcrash_pattern_reads\n"
      printf "%s\t%s\t%d\t%d\t%d\n", "${sample_id}", "unmapped", total+0, yc_total+0, crash_total+0
    }' > ${sample_id}.unmapped.tsv
    """
}
