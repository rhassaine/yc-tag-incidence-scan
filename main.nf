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
 *   - overlap between the above two
 *
 * Aggregates all chunk-level counts into one genome-wide, cohort-wide total,
 * and reports the observed rate (+ a rule-of-three upper bound if zero found).
 */

nextflow.enable.dsl=2

params.samplesheet   = null          // CSV: sample_id,file_path,file_type,ref_path
params.window_size   = 20_000_000    // bp per chunk; smaller = more parallelism
params.threads       = 4             // samtools threads per chunk task
params.outdir        = 'results'

include { MAKE_WINDOWS  } from './modules/make_windows.nf'
include { SCAN_REGION   } from './modules/scan_region.nf'
include { SCAN_UNMAPPED } from './modules/scan_unmapped.nf'
include { AGGREGATE     } from './modules/aggregate.nf'

workflow {

    if( !params.samplesheet )
        error "Provide --samplesheet path/to/samplesheet.csv"

    // sample_id,file_path,file_type,ref_path (ref_path may be blank for BAM)
    ch_samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def ref = row.ref_path?.trim() ? file(row.ref_path) : file('NO_FILE')
            tuple(row.sample_id, file(row.file_path), row.file_type.trim().toLowerCase(), ref)
        }

    // derive index file alongside each input, so Nextflow stages it too
    ch_samples_with_index = ch_samples.map { sample_id, f, type, ref ->
        def idx = type == 'cram' ? file("${f}.crai") : file("${f}.bai")
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
