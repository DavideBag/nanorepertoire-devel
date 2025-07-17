// This subworkflows is used to create the report of nanorepertoire analysis
// modules to include in this subworkflow

include { REPORT    }             from '../../../modules/local/report/main.nf'

workflow NANOREPORT {

    // with take we define the input channels

    take:
    report
    loop_tree
    clustersummaries
    cdr3histograms
    cdr3tables
    metadata



    main:

    REPORT(report, loop_tree, clustersummaries, cdr3histograms, cdr3tables, metadata)



    emit:
    analysis_report    = REPORT.out.report              // channel: [ val(meta), [ bai ] ]
    //versions           = REPORT.out.versions          // channel: [ versions.yml ]
}

