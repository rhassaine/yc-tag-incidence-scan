process MAKE_WINDOWS {
    tag "$sample_id"
    container params.samtools_container
    cpus 1
    memory '2 GB'

    input:
    tuple val(sample_id), path(bam_or_cram), path(index), val(file_type), path(ref)

    output:
    tuple val(sample_id), path(bam_or_cram), path(index), val(file_type), path(ref), path("${sample_id}.windows.txt"), emit: windows

    script:
    def window_size = params.window_size
    """
    set -euo pipefail

    # Pull contig name + length straight from this file's own header --
    # this is what makes chunking automatically cover every contig actually
    # present (main chromosomes, MT, decoys, unplaced contigs, etc.) rather
    # than relying on a hardcoded chromosome list.
    samtools view -H ${bam_or_cram} | awk '
      /^@SQ/ {
        sn=""; ln=""
        for (i=1;i<=NF;i++) {
          if (\$i ~ /^SN:/) sn = substr(\$i,4)
          if (\$i ~ /^LN:/) ln = substr(\$i,4)
        }
        if (sn!="" && ln!="") print sn, ln
      }
    ' > contigs.tsv

    # Split each contig into fixed-size windows (last window per contig may
    # be shorter). One region string per line, e.g. "chr1:1-20000000".
    awk -v win=${window_size} '
      {
        contig=\$1; len=\$2
        for (start=1; start<=len; start+=win) {
          end = start+win-1
          if (end>len) end=len
          print contig":"start"-"end
        }
      }
    ' contigs.tsv > ${sample_id}.windows.txt
    """
}
