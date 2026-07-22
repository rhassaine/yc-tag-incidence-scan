process MAKE_WINDOWS {
    tag "$sample_id"
    container params.samtools_container
    cpus 1
    memory '2 GB'

    input:
    tuple val(sample_id), path(bam_or_cram), path(index), val(file_type), val(ref)

    output:
    tuple val(sample_id), path(bam_or_cram), path(index), val(file_type), val(ref), path("${sample_id}.windows.txt"), emit: windows

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

    # Contigs larger than the window get split into "contig:start-end"
    # chunks, one per line, as before. Contigs at or below the window size
    # (e.g. ALT/HLA/decoy scaffolds -- a reference can have thousands of
    # these, each far smaller than any reasonable window) are NOT given
    # one task each -- that would blow up task count independent of
    # --window_size. Instead they're batched into groups of
    # small_contig_batch contig names per line; samtools view accepts
    # multiple region arguments in one invocation, so one line still
    # becomes one task, just scanning many small contigs together.
    awk -v win=${window_size} -v batch=${params.small_contig_batch} '
      {
        contig=\$1; len=\$2
        if (len > win) {
          for (start=1; start<=len; start+=win) {
            end = start+win-1
            if (end>len) end=len
            print contig":"start"-"end
          }
        } else {
          buf = buf (buf=="" ? "" : " ") contig
          count++
          if (count>=batch) { print buf; buf=""; count=0 }
        }
      }
      END {
        if (buf!="") print buf
      }
    ' contigs.tsv > ${sample_id}.windows.txt
    """
}
