<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-nanorepertoire_logo_dark.png">
    <img alt="nf-core/nanorepertoire" src="docs/images/nf-core-nanorepertoire_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/nanorepertoire/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/nanorepertoire/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/nanorepertoire/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/nanorepertoire/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/nanorepertoire/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/nanorepertoire)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23nanorepertoire-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/nanorepertoire)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/nanorepertoire** is a bioinformatics pipeline that enables the analysis and characterization of nanobody repertoires from NGS data.
It identifies expanded nanobody clusters potentially involved in antigen binding and immune response, and generates a detailed, reproducible report on the repertoire features of antibodies and nanobodies.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/index.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Functionality Overview

A graphical view of the pipeline can be seen below

<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nanorepertoire_metromap_dark.svg">
    <img alt="nf-core/nanorepertoire" src="docs/images/nanorepertoire_metromap_light.svg">
  </picture>
</h1>

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

### Input samplesheet

To run the pipeline, you need to prepare a CSV file containing metadata and paths to your input FASTQ files.
The samplesheet must contain the following columns:

`sample,fastq_1,fastq_2,individual,boost,immunisation`

Each row represents a pair of FASTQ files (paired-end) associated with a given biological sample and metadata.
Here's an example:

`sample,fastq_1,fastq_2,individual,boost,immunisation`
`sample_01,sample_1.fastq.gz,sample_2.fastq.gz,ind1,1,unenriched`
`sample_02,sample1_1.fastq.gz,sample2_2.fastq.gz,ind2,2,RBD_enriched`

#### Column descriptions:

- **sample**: Unique identifier for the sample.
- **fastq_1**: Path or URL to the first read (R1) FASTQ file (paired-end).
- **fastq_2**: Path or URL to the second read (R2) FASTQ file (paired-end).
- **individual**: Identifier for the biological individual or subject.
- **boost**: Boost number (e.g. 1 for priming, 2 for first boost, etc.).
- **immunisation**: Type of immunisation or enrichment strategy (e.g. `unenriched`, `RBD_enriched`, etc.).

Make sure all paths are accessible (local or remote) and the CSV is comma-separated (not tab-separated).

Now, you can run the pipeline using:

```bash
nextflow run nf-core/nanorepertoire \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
   -- adapterfile 'path-to-adapter-file.fa'
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/nanorepertoire/usage) and the [parameter documentation](https://nf-co.re/nanorepertoire/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/nanorepertoire/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/nanorepertoire/output).

## Credits

nf-core/nanorepertoire was originally written by Francesco Lescai and Davide Bagordo.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#nanorepertoire` channel](https://nfcore.slack.com/channels/nanorepertoire) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/nanorepertoire for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
