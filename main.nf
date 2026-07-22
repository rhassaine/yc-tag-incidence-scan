#!/usr/bin/env nextflow

/*
 * YC tag empty-duplex incidence scan
 *
 * For each sample (BAM/CRAM), splits the genome into intervals based on the
 * file's own @SQ header, scans each interval (plus unmapped reads separately)
 * for:
 *   - total reads
 *   - total YC tags seen
 *   - reads matching the crash-triggering pattern: YC:Z:<nonzero-head>++
 *
 * Aggregates all chunk-level counts into one genome-wide, cohort-wide total,
 * and reports the observed rate (+ a rule-of-three upper bound if zero found).
 *
 * Params (defaults set in nextflow.config, not here, to avoid drift):
 *   input, window_size, small_contig_batch, threads, outdir,
 *   samtools_container, base_container
 */

nextflow.enable.dsl=2

include { MAKE_WINDOWS  } from './modules/make_windows.nf'
include { SCAN_REGION   } from './modules/scan_region.nf'
include { SCAN_UNMAPPED } from './modules/scan_unmapped.nf'
include { AGGREGATE     } from './modules/aggregate.nf'

workflow {

    if( !params.input )
        error "Provide --input path/to/input.csv"

    // sample_id,file_path,file_type,ref_path,index_path
    // (ref_path may be blank for BAM; index_path is optional -- if blank,
    // falls back to deriving <file_path>.bai / .crai)
    ch_samples = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def ref = row.ref_path?.trim() ? file(row.ref_path, checkIfExists: false) : file('NO_FILE')
            def explicit_idx = row.index_path?.trim() ? file(row.index_path, checkIfExists: false) : null
            tuple(row.sample_id, file(row.file_path, checkIfExists: false), row.file_type.trim().toLowerCase(), ref, explicit_idx)
        }

    // use the explicit index if the samplesheet provided one; otherwise fall
    // back to deriving it by appending .bai/.crai (checkIfExists: false --
    // the implicit existence check on constructed GCS paths is known to be
    // unreliable; let staging/samtools surface a real error if the index is
    // genuinely missing, rather than failing here on a possibly-spurious check)
    ch_samples_with_index = ch_samples.map { sample_id, f, type, ref, explicit_idx ->
        def idx = explicit_idx ?: (type == 'cram' ? file("${f}.crai", checkIfExists: false) : file("${f}.bai", checkIfExists: false))
        tuple(sample_id, f, idx, type, ref)
    }

    // 1. generate genomic-interval windows per sample from its own header
    MAKE_WINDOWS(ch_samples_with_index)

    // 2. fan each sample's window list out into one task per interval
    ch_regions = MAKE_WINDOWS.out.windows
        .flatMap { sample_id, f, idx, type, ref, windows_file ->
            windows_file.readLines().collect { region ->
                tuple(sample_id, f, idx, type, ref, region)
            }
        }

    SCAN_REGION(ch_regions)

    // 3. separately scan unmapped reads per sample (one task per sample)
    SCAN_UNMAPPED(ch_samples_with_index)

    // 4. collect every chunk's stats (mapped-region chunks + unmapped chunk)
    //    and aggregate into one global summary
    ch_all_stats = SCAN_REGION.out.stats
        .mix(SCAN_UNMAPPED.out.stats)
        .map { sample_id, tsv -> tsv }
        .collect()

    AGGREGATE(ch_all_stats)
}
