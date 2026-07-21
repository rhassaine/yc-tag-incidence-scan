process AGGREGATE {
    publishDir params.outdir, mode: 'copy'
    cpus 1
    memory '2 GB'

    input:
    path(tsv_files)

    output:
    path("yc_incidence_summary.tsv")
    path("yc_incidence_per_sample.tsv")

    script:
    """
    set -euo pipefail
    export LC_ALL=C

    # Concatenate all chunk TSVs, keeping exactly one header line.
    head -n1 ${tsv_files[0]} > combined.tsv
    for f in ${tsv_files.join(' ')}; do
        tail -n +2 "\$f" >> combined.tsv
    done

    # Per-sample rollup
    awk -F'\\t' '
      NR==1 { next }
      {
        total[\$1]      += \$3
        yc[\$1]          += \$4
        softclip[\$1]    += \$5
        crash[\$1]       += \$6
        crashsoft[\$1]   += \$7
      }
      END {
        printf "sample_id\\ttotal_reads\\tyc_tags_seen\\tsoftclip_reads\\tcrash_pattern_reads\\tcrash_and_softclip\\n"
        for (s in total) {
          printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\n", s, total[s], yc[s], softclip[s], crash[s], crashsoft[s]
        }
      }
    ' combined.tsv | sort > yc_incidence_per_sample.tsv

    # Cohort-wide totals + observed rate + rule-of-three bound if zero found
    awk -F'\\t' '
      NR==1 { next }
      {
        total      += \$3
        yc         += \$4
        softclip   += \$5
        crash      += \$6
        crashsoft  += \$7
      }
      END {
        printf "metric\\tvalue\\n"
        printf "total_reads_scanned\\t%d\\n", total+0
        printf "total_yc_tags_seen\\t%d\\n", yc+0
        printf "total_softclip_reads\\t%d\\n", softclip+0
        printf "total_crash_pattern_reads\\t%d\\n", crash+0
        printf "total_crash_and_softclip\\t%d\\n", crashsoft+0
        if (yc > 0) {
          rate = crash / yc
          printf "observed_rate_per_yc_tag\\t%.10g\\n", rate
        }
        if (crash == 0 && yc > 0) {
          # Rule of three: with zero events in N trials, the ~95%% upper
          # confidence bound on the true rate is approximately 3/N.
          bound = 3.0 / yc
          printf "rule_of_three_upper_bound_95pct\\t%.10g\\n", bound
          printf "interpretation\\tzero observed in %d YC tags; true rate very likely below ~%.3g (95%% conf.)\\n", yc, bound
        } else if (crash > 0) {
          printf "interpretation\\tobserved %d crash-pattern reads out of %d YC tags scanned\\n", crash, yc
        }
      }
    ' combined.tsv > yc_incidence_summary.tsv

    cat yc_incidence_summary.tsv
    """
}
