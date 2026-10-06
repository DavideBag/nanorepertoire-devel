# nf-core/nanorepertoire: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc) - Raw read QC
- [Clustering summaries](#clustering-summaries) - Per-cluster tables produced from the CD-HIT output
- [CDR3 annotation](#cdr3-annotation) - Per-sample CDR3 calls from nanoCDR-X, with the outcome of every cluster representative
- [Aggregated statistics](#aggregated-statistics) - Tabular exports shared by the two reports
- [Repertoire report](#repertoire-report) - Interactive HTML report and its QC tables
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### Clustering summaries

<details markdown="1">
<summary>Output files</summary>

- `readcdhit/`
  - `*_clusters.summary`: one row per CD-HIT cluster, in CSV format.

</details>

`*_clusters.summary` is produced by `bin/readcdout.py` from the `.clstr` file written by CD-HIT and has the following schema:

| Column | Type | Description |
| :--- | :--- | :--- |
| `Representative` | string | Identifier of the sequence CD-HIT chose as the cluster representative. |
| `Count` | integer | Number of sequences in the cluster, representative included. |
| `Identity` | float or `NA` | Mean percent identity of the **non-representative** members to the representative, as reported by CD-HIT. `NA` for singleton clusters: with a single member there is nothing to compare against the representative, so the mean is undefined. It is not 100%. |
| `Identity_n` | integer | Number of members contributing to `Identity`, i.e. `Count - 1`. `0` for singletons. |

Because `Identity` is nullable, any downstream consumer must exclude `NA` rather than coerce it to a number: coercing to `0` places singletons at the origin of the identity histogram, and coercing to `100` asserts a perfect identity that was never measured. Values are bounded below by the clustering threshold (`--cdhit_identity`), since CD-HIT does not admit a sequence into a cluster below it.

### CDR3 annotation

<details markdown="1">
<summary>Output files</summary>

- `nanocdrx/`
  - `*_cdr3.tsv`: one row per cluster representative, with the CDR3 called by nanoCDR-X and its status.
  - `*_cdr3.fasta`: the unique CDR3 sequences, one record per first occurrence.
  - `*_cdr3.hist`: number of unique CDR3s per length, from 0 to 50 residues.
  - `*_cdr3_summary.tsv`: one row counting the representatives by outcome.

</details>

`bin/run_nanocdrx.py` runs nanoCDR-X on the cluster representatives written by CD-HIT, not on every sequence. Representatives shorter than 70 or longer than 150 residues (150 is the input length of the nanoCDR-X model), and representatives containing a residue outside the 20 standard amino acids, such as the `X` produced when a codon contains an undetermined base, are not sent to the model. They are kept in the outputs with their own status, so that every representative appears in `*_cdr3.tsv` exactly once. If nanoCDR-X fails, or does not return a result for every sequence it was given, the step stops with an error instead of producing partial output.

`*_cdr3.tsv` has the following schema:

| Column | Type | Description |
| :--- | :--- | :--- |
| `ID` | string | Identifier of the cluster representative, as in the CD-HIT output. |
| `CDR3` | string or `NA` | CDR3 called by nanoCDR-X. `NA` when no CDR3 was called or the representative was not sent to the model. |
| `sequence` | string | Amino-acid sequence of the representative. |
| `unique` | string | Status of the representative, one of the values below. Despite its name, the column holds a status, not a flag. |

| Status | Meaning |
| :--- | :--- |
| `unique` | CDR3 called; first occurrence of this CDR3 in the sample. Counted in `*_cdr3.hist` and written to `*_cdr3.fasta`. |
| `non-unique` | CDR3 called; the same CDR3 was already seen in the sample. |
| `cdr3-too-long` | CDR3 called but longer than 50 residues; kept out of the histogram. |
| `no-cdr3` | Sent to the model; no CDR3 called. |
| `excluded-length` | Not sent to the model: shorter than 70 or longer than 150 residues. |
| `excluded-nonstandard` | Not sent to the model: contains a residue outside the 20 standard amino acids. |

`*_cdr3_summary.tsv` holds one row per sample:

| Column | Description |
| :--- | :--- |
| `sample` | Sample identifier. |
| `input_sequences` | Cluster representatives read from the CD-HIT output. |
| `excluded_length`, `excluded_nonstandard` | Representatives not sent to the model, by reason. |
| `sent_to_model` | Representatives passed to nanoCDR-X: `input_sequences − excluded_length − excluded_nonstandard`. |
| `with_cdr3`, `no_cdr3` | Outcome of the model: `sent_to_model = with_cdr3 + no_cdr3`. |
| `unique`, `non_unique`, `cdr3_too_long` | Breakdown of `with_cdr3` by status. |

The annotation coverage of the model is `with_cdr3 / sent_to_model`; the share of representatives that could be passed to the model at all is `sent_to_model / input_sequences`.

### Aggregated statistics

<details markdown="1">
<summary>Output files</summary>

- `aggregate_stats/`
  - `analysis_report.html`: static scientific summary rendered from `assets/analysis_report.qmd`.
  - `nanobodies_report.RData`: serialised R workspace for downstream analysis.
  - `clustercounts.csv`, `cdrcounts.csv`, `cdrhists.csv`, `clusterbig.csv`, `fastaSeq.csv`: tabular exports consumed by the interactive report.
  - `cdr3_boost_overview_table.tsv`, `sampledata.tsv`: merged CDR3 and sample metadata tables.

</details>

`clustercounts.csv` holds one row per sample:

| Column | Description |
| :--- | :--- |
| `Sample` | Sample identifier. |
| `Clusters` | Number of clusters. |
| `Clusters_of_5`, `Clusters_of_100`, `Clusters_of_1000` | Number of clusters with at least 5, 100 and 1000 members. |
| `Total_sequences` | Total number of clustered sequences, `sum(Count)`. This is the denominator of every clonal fraction. |
| `Top1_count` | Size of the largest cluster. |
| `Top10_count` | Combined size of the ten largest clusters. |

`clusterbig.csv` is the concatenation of the per-sample `*_clusters.summary` files, with a `Sample` column and the three size-class flags used by the report; it therefore carries the `Identity`/`Identity_n` semantics described above.

### Repertoire report

<details markdown="1">
<summary>Output files</summary>

- `report/`
  - `nanorepertoire_report.html`: interactive HTML report.
  - `cdr3_boundary_qc.tsv`: per-sample summary of the CDR3 annotation QC.
  - `cdr3_boundary_offsets.tsv`: full distribution of the CDR3 boundary offsets.

</details>

Every figure of the report is described in the [README](../README.md#what-the-interactive-report-shows).

The diversity metrics table at the top of section 1 reports, for each sample, the quantities defined below. Metrics that are undefined for a sample are shown as `NA` rather than as a number. Throughout, `S` is the number of clusters of the sample, `c_i` the size of cluster `i`, `N = Σ c_i` the number of clustered sequences and `p_i = c_i / N` the clonal frequency.

| Metric | Definition | Notes |
| :--- | :--- | :--- |
| `Clusters` | `S` | Cluster richness. Not comparable between libraries sequenced at different depths. |
| `Clustered sequences` | `N` | Denominator of every fraction below. |
| `Clusters_per_1k_seqs` | `1000 · S / N` | Richness normalised by sequencing depth. This is the readout to use when comparing libraries of different size. |
| `Shannon_abundance` | `H' = −Σ p_i ln p_i` | Shannon index over clone abundances. `NA` when `N = 0`. |
| `Shannon_norm` | `H' / ln(S)` | Pielou's evenness. `NA` when `S = 1`, where `ln(S) = 0`. |
| `Clonality` | `1 − Shannon_norm` | `0` when every clone is equally abundant, approaching `1` as the library concentrates into one lineage. `NA` when `Shannon_norm` is. |
| `Simpson` | `1 − Σ p_i²` | Gini–Simpson index: probability that two sequences drawn at random belong to different clones. |
| `Gini` | Gini coefficient of `{c_i}` | `0` for equally sized clusters, approaching `1` as the sequences concentrate into one cluster. |
| `D50` | `min{k : Σ_{i≤k} c_(i) ≥ N/2}` over clusters sorted by decreasing size | Smallest number of clones covering half of the sequences. |
| `Top1_pct`, `Top10_pct` | `100 · c_(1) / N`, `100 · Σ_{i≤10} c_(i) / N` | Share of the library held by the largest one and ten clones, the standard AIRR clonal-fraction readout. |
| `Shannon_len` | `−Σ q_l ln q_l` over the CDR3 **length** distribution | Labelled "CDR3 length evenness (H')" in the report. It measures how evenly paratope lengths are spread and is **not** a repertoire diversity index. |
| `Pct_expanded` | `100 · Clusters_of_5 / S` | Proportion of clusters with at least 5 members. |
| `Pct_large` | `100 · Clusters_of_1000 / S` | Proportion of clusters with at least 1000 members. |

`cdr3_boundary_qc.tsv` compares the CDR3 called by nanoCDR-X with the motif-based definition (start after the framework-3 `T..Y.C` anchor, end before the framework-4 `WG.G` anchor, falling back to the C-terminal `TVSS` motif when the conserved tryptophan is substituted). Offsets are in residues, relative to that reference: a positive C-terminal offset means the called CDR3 extends into framework 4. `Compared` counts the sequences for which both anchors could be located and the called CDR3 was found in the sequence; the three `Skipped_*` columns account for the remainder.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
