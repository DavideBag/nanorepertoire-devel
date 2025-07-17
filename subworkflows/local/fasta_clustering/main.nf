// this subworkflow use fasta file input to analyse and study cd loop of nanobodies
// modules to include in this subworkflow

include {CDHIT_CDHIT   } from '../../../modules/nf-core/cdhit/cdhit/main.nf'
include {READCDHIT     } from '../../../modules/local/readcdhit/main.nf'
include {GETCDR3       } from '../../../modules/local/getcdr3/main.nf'
//include {MAFFT         } from '../../../modules/nf-core/mafft/main.nf'

// main workflow
workflow FASTA_CLUSTERING {

    // with take we define the input channels
    take:
    input      // channel: [[id], [reads_forward, reads_reverse]]

    main:
    // define the channels that will be used to store the versions of the software

    ch_versions = Channel.empty()
    input_ch = input


    CDHIT_CDHIT(input_ch)
    ch_versions = ch_versions.mix(CDHIT_CDHIT.out.versions.first())

    READCDHIT(CDHIT_CDHIT.out.clusters)
    ch_versions = ch_versions.mix(READCDHIT.out.versions.first())

    GETCDR3(CDHIT_CDHIT.out.fasta)
    ch_versions = ch_versions.mix(GETCDR3.out.versions.first())


    //MAFFT([[ id:'test', single_end:false ], [GETCDR3.out.fasta.groupTuple(by:[1, 2])]], [[:], []], [[:], []], [[:], []], [[:], []], [[:], []], false)
    //ch_versions = ch_versions.mix(MAFFT.out.versions.first())




    // defining the output channels that the workflow will
    emit:
    // emitted channels
    cluster          = CDHIT_CDHIT.out.fasta              // channel: [[id], [ fastq_merged ]]
    clusteread       = READCDHIT.out.summaryonly           // channel: [[id], [ fastq_renamed]]
    cdrpredicted     = GETCDR3.out.fasta
    cdrtsv           = GETCDR3.out.tsvonly
    cdrhistograms    = GETCDR3.out.histonly
    cdrmeta          = GETCDR3.out.metaonly
    //mafftfasta       = MAFFT.out.fas

    versions      = ch_versions                     // channel: [ versions.yml ]
}

