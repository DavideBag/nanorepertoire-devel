// this subworkflow prepares the inputs from fastq files and does the translation, in the meanwhile performs a fastqc
// modules to include in this subworkflow

include {CUTADAPT       } from '../../../modules/local/cutadapt/main.nf'
include {FLASH          } from '../../../modules/nf-core/flash/main.nf'
include {RENAME         } from '../../../modules/local/rename/main.nf'
include {NANOTRANSLATE  } from '../../../modules/local/nanotranslate/main.nf'


///Users/bagordo/Desktop/all/all_bioinformatics/nf-core-nanorepertoire/data/*_{1,2}_dummy2.fastq
// main workflow
workflow FASTQ_TO_FASTA {

    // with take we define the input channels
    take:
    input          // channel: [[id], [reads_forward, reads_reverse]]
    adapterfile    // channel: "adapter.fasta"

    main:
    // define the channels that will be used to store the versions of the software

    ch_versions = Channel.empty()
    input_ch = input

    CUTADAPT(input_ch, adapterfile)
    ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())

    FLASH(CUTADAPT.out.reads)
    ch_versions = ch_versions.mix(FLASH.out.versions.first())

    RENAME(FLASH.out.merged, single_end = true)
    ch_versions = ch_versions.mix(RENAME.out.versions.first())

    NANOTRANSLATE(RENAME.out.renamed)
    ch_versions = ch_versions.mix(NANOTRANSLATE.out.versions.first())

    // defining the output channels that the workflow will
    emit:
    // emitted channels
    trimmed_fastq = CUTADAPT.out.reads              // channel: [[id], [reads_forward, reads_reverse]]
    merged        = FLASH.out.merged                // channel: [[id], [ fastq_merged ]]
    renamed       = RENAME.out.renamed              // channel: [[id], [ fastq_renamed]]
    translated    = NANOTRANSLATE.out.fasta         // channel: [[id], [ fasta]]
    translog      = NANOTRANSLATE.out.log           // channel: [[id], [ fastalog]]
    versions      = ch_versions                     // channel: [ versions.yml ]
}

