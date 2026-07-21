# yc-tag-incidence-scan

Nextflow pipeline to measure the genome-wide incidence of a specific SBX YC
tag pattern — a non-zero head length paired with a truly empty duplex region
(e.g. `YC:Z:9++`) — across a cohort of BAM/CRAM files.

## Background

This pattern triggered a crash in REDUX's `getDuplexIndelIndices` (see
[hartwigmedical/hmftools](https://github.com/hartwigmedical/hmftools) commit
`2e156cc`, fixed in a follow-up patch). The pattern is valid per the SBX-D YC
tag spec (matches Example 5, "Only one read has bases"), but its real-world
incidence in native new-format sequencer output was unknown. This pipeline
scans a cohort empirically to answer that, rather than relying on
speculation about read-length/library-prep constraints.

## What it does

For each sample:
1. Derives genomic-interval chunks directly from the file's own `@SQ`
   header (`MAKE_WINDOWS`) — automatically covers every contig present
   (autosomes, X/Y, MT, decoys, unplaced), not a hardcoded chromosome list.
2. Scans each chunk (`SCAN_REGION`) and, separately, all unmapped reads
   (`SCAN_UNMAPPED`) — full tally, no early exit — for:
   - total reads
   - total YC tags seen
   - reads matching the crash pattern (`YC:Z:[0-9]+\+\+`)

Across all samples/chunks, `AGGREGATE` sums everything into one cohort-wide
total and reports:
- the observed rate (crash-pattern reads / total YC tags seen), and
- if zero crash-pattern reads are found, a rule-of-three ~95% upper
  confidence bound on the true rate (`3 / N`), so "we found none" becomes a
  quantitative statement rather than an absence of evidence.

## Usage

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --window_size 100000000 \
  --threads 4 \
  --outdir results
```

Resource allocation, executor, and GCP Batch configuration are expected to
be supplied externally (own execution wrapper/config) — `nextflow.config`
here only sets safe process-level defaults.

## Samplesheet format

```csv
sample_id,file_path,file_type,ref_path,index_path
sample01,/path/to/sample01.bam,bam,,
sample02,/path/to/sample02.cram,cram,/path/to/ref_genome.fasta,
sample03,/path/to/sample03.bam,bam,,/path/to/sample03.bam.bai
```

- `file_type`: `bam` or `cram` (case-insensitive)
- `ref_path`: required for CRAM (needs to match the exact reference the
  file was compressed against — a mismatched reference will cause CRAM
  decode errors, not silent corruption); leave blank for BAM
- `index_path`: optional. If blank, the pipeline derives the index by
  appending `.bai`/`.crai` to `file_path`. If provided, that exact path is
  used instead — recommended on cloud storage (e.g. GCS), where the
  implicit existence check on a constructed sibling path can be unreliable.

Each input file must already be indexed (`.bai` / `.crai`), whether via
the default sibling location or an explicit `index_path`.

## Outputs

Written to `--outdir` (default `results/`):
- `yc_incidence_summary.tsv` — cohort-wide totals, observed rate, and
  (if applicable) the rule-of-three bound
- `yc_incidence_per_sample.tsv` — same breakdown, one row per sample

## Notes / caveats

- This measures incidence in **native, non-migrated** new-format files.
  Don't point it at files that went through an old→new tag conversion —
  the question here is about real sequencer output.
- A zero-count result should be read as "the true rate is very likely
  below the reported bound," not "this can never happen" — statistical
  absence of evidence, not proof of impossibility.
