process SCAN_REGION {
    tag "${sample_id}:${region}"
    cpus params.threads
    memory '4 GB'

    input:
    tuple val(sample_id), path(bam_or_cram), path(index), val(file_type), path(ref), val(region)

    output:
    tuple val(sample_id), path("${sample_id}.${safe_region(region)}.tsv"), emit: stats

    script:
    def ref_flag = (file_type == 'cram') ? "-T ${ref}" : ''
    def out_name = "${sample_id}.${safe_region(region)}.tsv"
    """
    set -euo pipefail
    export LC_ALL=C

    samtools view -@ ${params.threads} ${ref_flag} ${bam_or_cram} ${region} | awk '
    {
      total++
      for (i=12; i<=NF; i++) {
        if (\$i ~ /^YC:Z:/) {
          yc_total++
          if (\$i ~ /^YC:Z:[0-9]+\\+\\+/) {
            crash_total++
          }
        }
      }
    }
    END {
      printf "sample_id\\tregion\\ttotal_reads\\tyc_tags_seen\\tcrash_pattern_reads\\n"
      printf "%s\\t%s\\t%d\\t%d\\t%d\\n", "${sample_id}", "${region}", total+0, yc_total+0, crash_total+0
    }' > ${out_name}
    """
}

// region strings like "chr1:1-20000000" aren't safe as bare filename
// fragments on their own (colons especially); sanitize for the output name.
def safe_region(String region) {
    region.replaceAll('[:]', '_')
}
