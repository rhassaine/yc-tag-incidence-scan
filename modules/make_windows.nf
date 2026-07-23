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

    # Three-way split, not two:
    #  1) Contigs larger than the window -> split into "contig:start-end"
    #     chunks, one per line, each its own task.
    #  2) Contigs at/below the window but above small_contig_threshold --
    #     these are "real" chromosomes that just happen to be under
    #     window_size (e.g. chr16-22 at a 100Mb window) and each still
    #     carries a substantial amount of read data. Each gets its own
    #     dedicated task, unsplit, rather than being lumped in with tiny
    #     scaffolds -- bundling a real ~50-90Mb chromosome together with
    #     thousands of KB-scale ALT/HLA/decoy contigs created one task
    #     that ended up processing far more data than any other task in
    #     the run, becoming the long-pole straggler for the whole pipeline.
    #  3) Contigs at/below small_contig_threshold (genuinely tiny
    #     ALT/HLA/decoy scaffolds -- a reference can have thousands of
    #     these) -- batched together, small_contig_batch names per line,
    #     since giving each one its own task blows up task count
    #     independent of --window_size.
    awk -v win=${window_size} -v small_thresh=${params.small_contig_threshold} -v batch=${params.small_contig_batch} '
      {
        contig=\$1; len=\$2
        if (len > win) {
          for (start=1; start<=len; start+=win) {
            end = start+win-1
            if (end>len) end=len
            print contig":"start"-"end
          }
        } else if (len > small_thresh) {
          print contig
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
