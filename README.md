<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-nanorepertoire_logo_dark.png">
    <img alt="nf-core/nanorepertoire" src="docs/images/nf-core-nanorepertoire_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new/nf-core/nanorepertoire)
[![GitHub Actions CI Status](https://github.com/nf-core/nanorepertoire/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/nanorepertoire/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/nanorepertoire/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/nanorepertoire/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/nanorepertoire/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.17379842-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.17379842)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.4.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.4.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/nanorepertoire)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23nanorepertoire-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/nanorepertoire)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction


**lescailab/nanorepertoire** is a Nextflow-based bioinformatics pipeline designed to characterise *nanobody repertoires* from raw sequencing data (FASTQ files).
It performs quality control, read preprocessing, translation, clustering, and CDR3 extraction, producing comprehensive reports that describe nanobody diversity and repertoire composition.

The workflow is divided into three main subworkflows:

### 1. `fastq_to_fasta`
Processes raw sequencing data to produce translated FASTA sequences:
- **FastQC** – quality control of raw reads
- **Cutadapt** – adapter trimming
- **FLASH** – paired-end read merging
- **Rename** – standardized renaming of merged reads
- **Nanotranslate** – translation from nucleotide to amino acid sequences

### 2. `fasta_clustering`
Clusters and analyses translated nanobody sequences:
- **CD-HIT** – clustering of identical or highly similar sequences (configurable identity threshold)
- **ReadCDHIT** – extraction and summarization of cluster statistics
- **Nanocdr-x** – CDR3 extraction and analysis tool

### 3. `repertoire_report`
Generates the reports. Statistical aggregation and visualisation are separate steps, so new panels can be added to the interactive report without touching the analytical modules:
- **MultiQC** – aggregation of QC metrics across all initial samples
- **Aggregate Stats (R-based)** – statistical data engine for production of summary tables (`.csv`, `.tsv`), serialized objects (`.RData`), and static scientific summaries.
- **Nanorepertoire Report (Python-based)** – generation of a modern, interactive HTML dashboard based on the aggregated metrics.

## Parameters

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `--cdhit_identity` | `0.9` | Sequence identity threshold for CD-HIT clustering. Also reported verbatim in the HTML report. |
| `--cdhit_word_size` | `5` | Word size (`-n`) used by CD-HIT. Must be compatible with `--cdhit_identity`. |
| `--cluster_size_threshold` | `1` | Minimum sequences per cluster for inclusion in reports. |
| `--calculate_tree` | `false` | Enable/disable phylogenetic tree calculation. |
| `--adapterfile` | `null` | Custom adapter file for Cutadapt. |

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

1. Prepare a samplesheet (`samplesheet.csv`):
```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

2. Run the pipeline:

> [!WARNING]
> This pipeline currently requires **Nextflow version 25.04.0** due to syntax updates in newer Nextflow versions. You must set the `NXF_VER` environment variable when running the pipeline.
> We highly recommend using the `conda` profile for managing software dependencies.

```bash
NXF_VER=25.04.0 nextflow run lescailab/nanorepertoire \
  -profile conda \
  --input samplesheet.csv \
  --outdir results/
```

## Pipeline output

The pipeline produces:
- **Translated Sequences**: Amino acid sequences in FASTA format.
- **Cluster Stats**: Detailed statistics on nanobody clusters.
- **QC Reports**: `multiqc_report.html` for technical quality metrics.
- **Hybrid Reports**:
    - `analysis_report.html`: Static scientific summary (R/Quarto).
    - `nanorepertoire_report.html`: Interactive dashboard (Python/Plotly).
    - `nanobodies_report.RData`: Serialized R object for downstream analysis.
- **Annotation QC**: `cdr3_boundary_qc.tsv` and `cdr3_boundary_offsets.tsv`, comparing the called CDR3 boundaries with the IMGT definition of the CDR3 (CDR3-IMGT).

The file-by-file description of the output directories, including the schema of every table, is in [`docs/output.md`](docs/output.md).

## What the interactive report shows

`nanorepertoire_report.html` contains twelve panels. This section states, for each one, what is plotted, how it is computed, how it should be read, and what it cannot support. Symbols used throughout: `S` is the number of CD-HIT clusters of a sample, `c_i` the number of sequences in cluster `i`, `N = Σ c_i` the number of clustered sequences, and `p_i = c_i / N` the frequency of clone `i`.

### Summary cards

Five figures at the top of the report: the number of sequences for which the CDR3 caller returned a CDR3 (rows of `fastaSeq.csv` whose `CDR3` is not `NA`); the total number of clusters over all samples, annotated with the identity threshold actually used in the run; the number of **distinct** CDR3 sequences (`Unique_CDR3s` summed over samples); the mean CDR3 length with its standard deviation and median; and the number of libraries analysed. The first and third cards answer different questions and will not agree: one counts sequences, the other counts distinct paratopes.

### 1. Clonal cluster analysis

**Diversity metrics table.** One row per sample, computed by `bin/repertoire_metrics.py` from the cluster sizes in `clusterbig.csv`. It reports cluster richness `S`, the number of clustered sequences `N`, the depth-normalised richness `1000·S/N`, the Shannon index of the clone abundances `H' = −Σ p_i ln p_i`, clonality `1 − H'/ln(S)`, the Gini–Simpson index `1 − Σ p_i²`, the Gini coefficient of the cluster sizes, `D50` (the smallest number of clones covering half of the sequences), the fraction of the library held by the largest clone and by the ten largest clones, the Shannon entropy of the CDR3 *length* distribution, and the proportion of clusters with ≥ 5 and ≥ 1000 members.

*How to read it.* Clonality is 0 when every clone is equally abundant and approaches 1 as the library concentrates into a single lineage; antigen selection is expected to raise it. `D50` and the Top1/Top10 fractions say the same thing in units a wet-lab reader can act on ("the dominant clone is X% of the library"). Richness (`S`) rises with sequencing depth, so **only the depth-normalised richness should be compared between libraries** sequenced at different depths.

*Caveats.* Clonality and evenness are undefined for a single cluster (`ln(S) = 0`) and are reported as `NA` rather than as 0 or 1. All abundance metrics are sensitive to the clustering threshold: a lower threshold merges clones and mechanically increases clonality. "CDR3 length evenness (H')" measures how evenly paratope **lengths** are spread and must not be read as a diversity index; it is kept because it is informative about length selection, and it is labelled separately for exactly that reason.

**Clonal Expansion Profile** (stacked bar). Clusters of each sample split into four size classes: singletons and small clusters (< 5 members), 5–99, 100–999, ≥ 1000. Thresholds follow Deschaght et al. 2017. Read the *shape*: immune repertoires are strongly skewed, with the large majority of clusters below 5 members and a thin tail of dominant clones. The absolute height of the bars scales with sequencing depth and carries no biological meaning on its own.

**Cluster Abundance at Size Thresholds** (grouped bar, logarithmic y axis). The same information as cumulative counts (total, ≥ 5, ≥ 100, ≥ 1000), on a log scale so that four orders of magnitude are visible at once. The slope from one threshold to the next is the steepness of the clonal hierarchy: a shallow decay means abundance is spread over many lineages, a steep one means it is concentrated. A log axis cannot display a zero count, so a missing bar means "none", not "small".

**Top-10 Largest Clonotypes per Sample** (horizontal bar). The ten largest clusters of each sample, labelled with the identifier of their representative sequence and sorted by size. This is the candidate list for downstream expression and affinity screening: the representative identifiers map back to `clusterbig.csv` and to the sequences in `fastaSeq.csv`. The ratio between rank 1 and rank 10 quantifies how much the library is dominated by a single lineage. Rank order is stable, absolute sizes are not: they depend on depth and on PCR.

### 2. CDR3 diversity

**Unique CDR3 Sequences per Sample** (bar). Number of distinct CDR3 amino-acid sequences per sample, from `cdrcounts.csv`. It is a richness measure of the paratope space and, like every raw count, it grows with sequencing depth; comparisons between libraries need the depth-normalised metrics of section 1.

**CDR3 Length Distribution (Violin)** (with embedded box plot). Distribution of the length in residues of the unique CDR3s of each sample; the box marks the interquartile range, the line the median, the width the density. VHH CDR3s are longer than the CDR3s of conventional VH domains (a mean around 19 residues has been reported for camelid VHH, against roughly 12 for human VH), and lengths above ~16 residues are associated with the protruding, finger-like loops able to reach cleft and cavity epitopes. Compare *medians and spread* between samples: a shift indicates that selection acted on paratope geometry, its absence that it did not. Only lengths of 1–40 residues are drawn, and the length depends on where the boundaries were called, which is why the annotation QC below matters.

**CDR3 Length Frequency Profile** (filled line chart). The same distribution as a frequency profile over lengths 0–45, from `cdrhists.csv`, with one curve per sample. Better than the violin for spotting multiple modes and for comparing the shoulders of two samples directly. Counts are of unique CDR3s, not of sequences, so a highly expanded clone contributes one observation, not thousands.

**CDR3 annotation QC — summary table and CDR3 Boundary Offsets.** For every sequence with a called CDR3, `bin/cdr3_boundary_qc.py` re-derives the CDR3 as defined by IMGT (CDR3-IMGT, positions 105–117, between the conserved cysteine 104 and tryptophan 118, both excluded; Lefranc et al. 2003), locating the two anchors with sequence motifs — start immediately after the cysteine of the framework-3 anchor `T..Y.C` (the `YYC` motif), end immediately before the tryptophan of the framework-4 anchor `WG.G` (the `WGQ` motif), falling back to the C-terminal `TVSS` motif in reads where that tryptophan is substituted — and compares it residue by residue with the boundaries called by nanoCDR-X. The table reports the number of comparable sequences, the exact agreement, and the median offset at each terminus; the CDR3 Boundary Offsets panel shows the full offset distribution per sample and terminus, on a logarithmic y axis so that rare disagreements stay visible next to the dominant zero-offset bar.

*How to read it.* An offset of 0 means the two definitions agree exactly. A systematically **positive C-terminal offset** means the called CDR3 extends into framework 4, which is rich in alanine, threonine and valine and would inflate those residues in the composition panels and lengthen the length distributions. A negative N-terminal offset means the call reaches back into framework 3. Sequences in which neither anchor can be located are not comparable and are counted separately: they are not evidence of agreement. The reference locates the IMGT anchors with sequence motifs rather than by alignment, so it is an approximation, not ground truth, and this panel measures *consistency between two annotations*, not accuracy against a curated reference.

### 3. CDR3 amino-acid composition

**Amino Acid Frequency Heatmap.** Frequency of each of the twenty standard amino acids as a percentage of all residues of the unique CDR3s of a sample; rows are samples, columns are residues, colour is the percentage. **Side-by-Side Amino Acid Usage** (grouped bar) shows the same matrix per residue, which is the better view for reading off a difference between two samples at a single position.

*How to read them.* VHH paratopes are characteristically rich in tyrosine and serine, whose hydroxyl and aromatic side chains provide most of the hydrogen bonding and stacking contacts; non-canonical cysteines are of particular interest because they can form an interloop disulfide that rigidifies the paratope, so a change in cysteine frequency between conditions is a hypothesis worth testing. Read differences *between* samples on the same residue rather than the absolute ranking, which is largely conserved across VHH repertoires.

*Caveats.* Frequencies are computed over unique CDR3s, so they describe the sequence space rather than the abundance-weighted repertoire: a dominant clone counts once. Composition is conditional on the annotated boundaries, so it must be read together with the annotation QC above. Non-standard characters are excluded from both the numerator and the denominator.

### 4. Clonality

**Mean Member Identity to Cluster Representative** (overlaid histograms). One value per cluster: the mean percent identity, as reported by CD-HIT, of the non-representative members of a cluster to its representative. Singleton clusters are excluded, the quantity being undefined for them, and the number of informative clusters is stated in the caption; the x axis starts just below the clustering threshold.

*How to read it.* This is a descriptive measure of **intra-clonal homogeneity**: values near 100% mean the members of a cluster are near-identical to their representative, lower values mean the cluster groups more heterogeneous sequences. A shift of the distribution towards 100% between conditions means the surviving lineages are internally more uniform.

*Caveats — this panel is the easiest one to over-read.* The distribution is **truncated by construction**: CD-HIT admits no member below the identity threshold, so nothing can be observed to the left of it and the metric cannot quantify diversification beyond that window. Identity is measured against an empirically chosen cluster member, not against an inferred germline, so this is **not** a somatic hypermutation rate and is not equivalent to the V-gene mutation frequency of MiXCR or IMGT/V-QUEST; the analogy is conceptual only. Finally, without upstream deduplication, identity at 100% cannot be attributed to clonal expansion rather than to PCR duplication, and differences below 100% cannot be attributed to somatic mutation rather than to sequencing error.

**Repertoire Diversity Landscape** (bubble chart). Each sample is one bubble: the Shannon index of the clone abundances on the x axis, clonality (1 − Pielou evenness) on the y axis, and the number of unique CDR3s as bubble area. Samples move to the lower right as sequences spread evenly over many clones and to the upper left as they concentrate into a few dominant lineages, which is the direction antigen selection is expected to push a library. Both axes derive from the same abundance vector and are therefore not independent; the hover text carries the depth-normalised richness and the dominant-clone fraction, which is what makes a comparison between libraries of different size legitimate. Samples for which the metrics are undefined are omitted from the panel and shown as `NA` in the table.

### 5. Methods summary

A summary of the pre-processing, clustering and metric definitions, with the identity threshold and word size of the run substituted into the text, and the key references.

## Credits

**lescailab/nanorepertoire** was originally written by Francesco Lescai.
Key contributors: Davide Bagordo.

## Citations

If you use lescailab/nanorepertoire, please cite: [10.5281/zenodo.17379842](https://doi.org/10.5281/zenodo.17379842)

An extensive list of references for the tools used can be found in [`CITATIONS.md`](CITATIONS.md).
for proper execution.
