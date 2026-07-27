#!/usr/bin/env python3
import csv, json, math, argparse, sys, os
from collections import defaultdict, Counter
from datetime import datetime

# repertoire_metrics.py sits next to this script in the pipeline bin/ directory,
# which is on PATH but not necessarily on PYTHONPATH.
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from repertoire_metrics import abundance_metrics

# ─────────────────────────────────────────────────────────────────────────────
# 0. I/O & ARGUMENT PARSING
# ─────────────────────────────────────────────────────────────────────────────
def read_any_csv(path):
    """Read CSV or TSV and return list of dicts."""
    if not os.path.exists(path):
        return []
    ext = os.path.splitext(path)[1].lower()
    sep = '\t' if ext in ['.tsv', '.summary', '.tab'] else ','
    with open(path, newline='', encoding='utf-8') as f:
        # Check first line for separator if ambiguous
        first_line = f.readline()
        f.seek(0)
        if ext == '.csv' and '\t' in first_line and ',' not in first_line:
            sep = '\t'
        return list(csv.DictReader(f, delimiter=sep))

def read_and_merge(paths):
    """Read multiple CSVs/TSVs and concatenate them, injecting Sample name if missing."""
    merged = []
    if not paths:
        return merged
    for p in paths:
        data = read_any_csv(p)
        if not data: continue
        
        # Derive sample name from filename if "Sample" or "sample" column is missing
        sample_name = os.path.basename(p).split('.')[0].replace('_clusters', '').replace('_cdr3', '')
        for row in data:
            if 'Sample' not in row and 'sample' not in row:
                row['Sample'] = sample_name
            # Standardize case for Sample column
            if 'sample' in row and 'Sample' not in row:
                row['Sample'] = row['sample']
        merged.extend(data)
    return merged

parser = argparse.ArgumentParser(description="Nanorepertoire Technical Report Generator")
parser.add_argument("--clustercounts", nargs="*", help="One or more cluster summary CSVs")
parser.add_argument("--cdrcounts", nargs="*", help="One or more CDR3 count CSVs")
parser.add_argument("--cdrhists", nargs="*", help="One or more CDR3 size histogram CSVs")
parser.add_argument("--clusterbig", nargs="*", help="One or more large cluster detail CSVs")
parser.add_argument("--fastaseq", nargs="*", help="One or more unique CDR3 sequence CSVs/FASTAs")
parser.add_argument("--output", default="nanorepertoire_report.html", help="Output HTML file path")
parser.add_argument("--identity-threshold", type=float, default=90.0,
                    help="CD-HIT sequence identity threshold used for clustering, as a percentage "
                         "(default: 90.0)")
parser.add_argument("--word-size", type=int, default=5,
                    help="CD-HIT word size (-n) used for clustering (default: 5)")

args = parser.parse_args()

print("Reading and merging data...")
clustercounts = read_and_merge(args.clustercounts)
cdrcounts     = read_and_merge(args.cdrcounts)
cdrhists      = read_and_merge(args.cdrhists)
clusterbig    = read_and_merge(args.clusterbig)
fastaseq      = read_and_merge(args.fastaseq)
output_html   = args.output

# Clustering parameters, reported verbatim in the text of the report so that the
# document always describes the analysis that was actually run.
identity_threshold = args.identity_threshold
word_size          = args.word_size
identity_pct       = f"{identity_threshold:g}%"

# ─────────────────────────────────────────────────────────────────────────────
# 1. DATA PREPARATION
# ─────────────────────────────────────────────────────────────────────────────
def to_int(v, default=0):
    try: return int(float(v))
    except: return default

def to_float(v, default=0.0):
    try: return float(v)
    except: return default

AA_STANDARD = set("ACDEFGHIKLMNPQRSTVWY")

# expansion holds summary stats per sample: total, small, medium, large, etc.
expansion = {}

# If clustercounts is missing, we try to derive it from clusterbig
if not clustercounts and clusterbig:
    print("Deriving cluster summaries from cluster detail data...")
    temp_stats = defaultdict(lambda: {"total": 0, "of5": 0, "of100": 0, "of1000": 0})
    for r in clusterbig:
        s = r.get("Sample", "Unknown")
        count = to_int(r.get("Count", r.get("count", 0)))
        temp_stats[s]["total"] += 1
        if count >= 5:    temp_stats[s]["of5"] += 1
        if count >= 100:  temp_stats[s]["of100"] += 1
        if count >= 1000: temp_stats[s]["of1000"] += 1
    
    for s, st in temp_stats.items():
        # Populate expansion for charts
        expansion[s] = {
            "total":     st["total"],
            "large":     st["of1000"],
            "medium":    st["of100"] - st["of1000"],
            "small":     st["of5"]   - st["of100"],
            "singleton": st["total"] - st["of5"],
            "of5": st["of5"], "of100": st["of100"], "of1000": st["of1000"]
        }
        # Also populate clustercounts for global stats
        clustercounts.append({
            "Sample": s,
            "Clusters": st["total"],
            "Clusters_of_5": st["of5"],
            "Clusters_of_100": st["of100"],
            "Clusters_of_1000": st["of1000"]
        })
else:
    for r in clustercounts:
        s = r["Sample"]
        total   = to_int(r.get("Clusters", 0))
        of5     = to_int(r.get("Clusters_of_5", 0))
        of100   = to_int(r.get("Clusters_of_100", 0))
        of1000  = to_int(r.get("Clusters_of_1000", 0))
        expansion[s] = {
            "total": total,
            "large":     of1000,
            "medium":    of100  - of1000,
            "small":     of5    - of100,
            "singleton": total  - of5,
            "of5": of5, "of100": of100, "of1000": of1000
        }

samples = list(expansion.keys())
if not samples:
    samples = sorted(list(set([r.get("Sample", r.get("sample", "Unknown")) for r in clusterbig + cdrcounts + cdrhists + fastaseq])))

PALETTE = {
    samples[i]: c for i, c in enumerate(
        ["#4E79A7","#F28E2B","#E15759","#76B7B2",
         "#59A14F","#EDC948","#B07AA1","#FF9DA7"]
    ) if i < len(samples)
}

# CDR3 analysis per sample
cdr3_by_sample = defaultdict(list)
for r in fastaseq:
    if r.get("unique") != "unique":
        continue
    cdr3 = r.get("CDR3", "").strip()
    if cdr3 and cdr3 != "NA":
        cdr3_by_sample[r.get("Sample", r.get("sample", "Unknown"))].append(cdr3)

all_cdr3_lens = [len(c) for s in cdr3_by_sample.values() for c in s]
mean_len   = sum(all_cdr3_lens) / len(all_cdr3_lens) if all_cdr3_lens else 0
var_len    = sum((x - mean_len)**2 for x in all_cdr3_lens) / len(all_cdr3_lens) if all_cdr3_lens else 0
sd_len     = math.sqrt(var_len)
median_len = sorted(all_cdr3_lens)[len(all_cdr3_lens)//2] if all_cdr3_lens else 0

# Number of translated sequences for which the CDR3 caller returned a CDR3.
# Rows with no called CDR3 carry "NA" in the CDR3 column of fastaSeq.csv.
seqs_with_cdr3 = sum(
    1 for r in fastaseq if (r.get("CDR3") or "").strip() not in ("", "NA")
)
total_clusters = sum(to_int(r["Clusters"]) for r in clustercounts)
total_cdr3     = sum(to_int(r["Unique_CDR3s"]) for r in cdrcounts) or sum(len(set(v)) for v in cdr3_by_sample.values())

# Shannon diversity (on CDR3 length distribution)
def shannon(lens):
    c = Counter(lens)
    total = sum(c.values())
    return -sum((v/total)*math.log(v/total + 1e-12) for v in c.values()) if total > 0 else 0

# Amino acid frequency per sample
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
aa_freq_by_sample = {}
for s, cdr3s in cdr3_by_sample.items():
    all_aa = [aa for cdr3 in cdr3s for aa in cdr3.upper() if aa in AA_STANDARD]
    total = len(all_aa)
    aa_freq_by_sample[s] = {
        aa: (all_aa.count(aa)/total*100 if total > 0 else 0) for aa in AA_ORDER
    }

# Cluster sizes per sample, the input of every clonal abundance metric
counts_by_sample = defaultdict(list)
for r in clusterbig:
    counts_by_sample[r.get("Sample", "Unknown")].append(to_int(r.get("Count", r.get("count", 0))))

# Diversity metrics table
div_metrics = []
for r in clustercounts:
    s = r["Sample"]
    n_cdr3 = next((to_int(x["Unique_CDR3s"]) for x in cdrcounts if x["Sample"] == s), 0)
    lens   = [len(c) for c in cdr3_by_sample.get(s, [])]
    H_len  = round(shannon(lens), 3)
    of5    = to_int(r["Clusters_of_5"])
    total  = to_int(r["Clusters"])
    of1000 = to_int(r["Clusters_of_1000"])
    pct_exp   = round(of5/total*100, 1) if total > 0 else 0
    pct_large = round(of1000/total*100, 2) if total > 0 else 0
    row = {
        "Sample": s, "Unique_CDR3s": n_cdr3, "Total_Clusters": total,
        # Entropy of the CDR3 *length* distribution: a measure of how evenly the
        # paratope lengths are spread, not of repertoire diversity.
        "Shannon_len": H_len, "Pct_expanded": pct_exp, "Pct_large": pct_large
    }
    # Clonal abundance metrics, computed on the cluster sizes of the sample.
    row.update(abundance_metrics(counts_by_sample.get(s, [])))
    row["Total_Clusters"] = total or row["Clusters"]
    div_metrics.append(row)

def fmt_metric(value, digits=3, suffix=""):
    """Render a metric, or NA when it is not defined for this sample."""
    if value is None:
        return '<span style="color:#95a5a6">NA</span>'
    if isinstance(value, int):
        return f"{value:,}{suffix}"
    return f"{value:.{digits}f}{suffix}"

# CDR3 histogram data per sample
hist_by_sample = defaultdict(lambda: defaultdict(int))
for r in cdrhists:
    hist_by_sample[r["Sample"]][to_int(r["Size"])] = to_int(r["Count"])

# ─────────────────────────────────────────────────────────────────────────────
# 2. BUILD PLOTLY TRACES AS JSON
# ─────────────────────────────────────────────────────────────────────────────

def hex_alpha(hex_color, alpha_hex):
    return hex_color + alpha_hex  # e.g. "#4E79A780"

# ── Fig 1: Cluster stacked bar ─────────────────────────────────────────────
fig1_traces = []
size_classes = [("large","≥1000 members","#d62728"),
                ("medium","100–999 members","#ff7f0e"),
                ("small","5–99 members","#1f77b4"),
                ("singleton","Singletons (<5)","#aec7e8")]
for key, label, color in size_classes:
    fig1_traces.append({
        "x": samples,
        "y": [expansion[s][key] for s in samples],
        "name": label,
        "type": "bar",
        "marker": {"color": color},
        "hovertemplate": "<b>%{x}</b><br>" + label + ": %{y:,}<extra></extra>"
    })
fig1_layout = {
    "barmode": "stack",
    "title": {"text": "Cluster Size Distribution per Sample", "font": {"size": 16, "color": "#2c3e50"}},
    "xaxis": {"title": "Sample"},
    "yaxis": {"title": "Number of Clusters"},
    "legend": {"title": {"text": "Size Class"}},
    "plot_bgcolor": "#fafafa", "paper_bgcolor": "#ffffff"
}

# ── Fig 2: Cluster count grouped bar (log scale) ───────────────────────────
thresholds = [
    ("Clusters",       "Total"),
    ("Clusters_of_5",  "≥5"),
    ("Clusters_of_100","≥100"),
    ("Clusters_of_1000","≥1000")
]
colors2 = ["#003f5c","#58508d","#bc5090","#ff6361"]
fig2_traces = []
cc_map = {r["Sample"]: r for r in clustercounts}
for (col, label), color in zip(thresholds, colors2):
    fig2_traces.append({
        "x": samples,
        "y": [to_int(cc_map[s][col]) for s in samples],
        "name": label,
        "type": "bar",
        "marker": {"color": color},
        "hovertemplate": "<b>%{x}</b><br>" + label + ": %{y:,}<extra></extra>"
    })
fig2_layout = {
    "barmode": "group",
    "title": {"text": "Cluster Counts at Multiple Size Thresholds", "font": {"size": 16, "color": "#2c3e50"}},
    "xaxis": {"title": "Sample"},
    "yaxis": {"title": "Count (log scale)", "type": "log"},
    "legend": {"title": {"text": "Threshold"}},
    "plot_bgcolor": "#fafafa", "paper_bgcolor": "#ffffff"
}

# ── Fig 3: CDR3 length violin ──────────────────────────────────────────────
fig3_traces = []
for s in samples:
    lens = [len(c) for c in cdr3_by_sample.get(s, []) if 1 <= len(c) <= 40]
    col  = PALETTE.get(s, "#888888")
    fig3_traces.append({
        "type": "violin",
        "y": lens,
        "name": s,
        "box": {"visible": True},
        "meanline": {"visible": True},
        "fillcolor": col + "80",
        "line": {"color": col},
        "hovertemplate": f"<b>{s}</b><br>CDR3 Length: %{{y}} AA<extra></extra>"
    })
fig3_layout = {
    "title": {"text": "CDR3 Length Distribution per Sample (VHH)", "font": {"size": 16, "color": "#2c3e50"}},
    "yaxis": {"title": "CDR3 Length (AA)", "range": [0, 42]},
    "xaxis": {"title": "Sample"},
    "violingap": 0.1,
    "plot_bgcolor": "#fafafa", "paper_bgcolor": "#ffffff"
}

# ── Fig 4: CDR3 length frequency area chart ───────────────────────────────
fig4_traces = []
for s in samples:
    col    = PALETTE.get(s, "#888888")
    sizes  = list(range(0, 46))
    counts = [hist_by_sample[s].get(sz, 0) for sz in sizes]
    fig4_traces.append({
        "x": sizes, "y": counts,
        "name": s,
        "type": "scatter", "mode": "lines+markers",
        "line": {"color": col, "width": 2},
        "marker": {"color": col, "size": 4},
        "fill": "tozeroy",
        "fillcolor": col + "30",
        "hovertemplate": f"<b>{s}</b><br>Length: %{{x}} AA<br>Count: %{{y:,}}<extra></extra>"
    })
fig4_layout = {
    "title": {"text": "CDR3 Length Frequency Profile", "font": {"size": 16, "color": "#2c3e50"}},
    "xaxis": {"title": "CDR3 Length (AA)", "range": [0, 45]},
    "yaxis": {"title": "Unique CDR3s"},
    "hovermode": "x unified",
    "legend": {"title": {"text": "Sample"}},
    "plot_bgcolor": "#fafafa", "paper_bgcolor": "#ffffff"
}

# ── Fig 5: Unique CDR3s bar ────────────────────────────────────────────────
cdr_map = {r["Sample"]: to_int(r["Unique_CDR3s"]) for r in cdrcounts}
fig5_traces = [{
    "x": samples,
    "y": [cdr_map.get(s, 0) for s in samples],
    "type": "bar",
    "marker": {"color": [PALETTE.get(s,"#888") for s in samples]},
    "text": [f"{cdr_map.get(s,0):,}" for s in samples],
    "textposition": "outside",
    "hovertemplate": "<b>%{x}</b><br>Unique CDR3s: %{y:,}<extra></extra>"
}]
fig5_layout = {
    "title": {"text": "Unique CDR3 Sequences per Sample", "font": {"size": 16, "color": "#2c3e50"}},
    "xaxis": {"title": "Sample"},
    "yaxis": {"title": "Unique CDR3s"},
    "showlegend": False,
    "plot_bgcolor": "#fafafa", "paper_bgcolor": "#ffffff"
}

# ── Fig 6: AA heatmap ─────────────────────────────────────────────────────
z_matrix = [[aa_freq_by_sample.get(s, {}).get(aa, 0) for aa in AA_ORDER] for s in samples]
fig6_traces = [{
    "type": "heatmap",
    "z": z_matrix,
    "x": AA_ORDER,
    "y": samples,
    "colorscale": "Viridis",
    "hovertemplate": "<b>%{y}</b><br>%{x}: %{z:.2f}%<extra></extra>",
    "colorbar": {"title": {"text": "Freq (%)"}}
}]
fig6_layout = {
    "title": {"text": "CDR3 Amino Acid Composition Heatmap", "font": {"size": 16, "color": "#2c3e50"}},
    "xaxis": {"title": "Amino Acid"},
    "yaxis": {"title": "Sample"},
    "plot_bgcolor": "#fafafa", "paper_bgcolor": "#ffffff"
}

# ── Fig 7: AA grouped bar ─────────────────────────────────────────────────
fig7_traces = []
for s in samples:
    col  = PALETTE.get(s, "#888888")
    freq = aa_freq_by_sample.get(s, {})
    fig7_traces.append({
        "x": AA_ORDER,
        "y": [round(freq.get(aa, 0), 2) for aa in AA_ORDER],
        "name": s,
        "type": "bar",
        "marker": {"color": col},
        "hovertemplate": f"<b>{s}</b><br>%{{x}}: %{{y:.2f}}%<extra></extra>"
    })
fig7_layout = {
    "barmode": "group",
    "title": {"text": "CDR3 Amino Acid Usage (% of all CDR3 residues)", "font": {"size": 16, "color": "#2c3e50"}},
    "xaxis": {"title": "Amino Acid"},
    "yaxis": {"title": "Frequency (%)"},
    "legend": {"title": {"text": "Sample"}},
    "plot_bgcolor": "#fafafa", "paper_bgcolor": "#ffffff"
}

# ── Fig 8: Diversity bubble chart ─────────────────────────────────────────
# Samples for which the abundance metrics are defined (at least two clusters
# holding at least one sequence).
fig8_metrics = [m for m in div_metrics
                if m["Shannon_abundance"] is not None and m["Clonality"] is not None]
fig8_traces = [{
    "x": [m["Shannon_abundance"] for m in fig8_metrics],
    "y": [m["Clonality"] for m in fig8_metrics],
    "mode": "markers+text",
    "type": "scatter",
    "text": [m["Sample"] for m in fig8_metrics],
    "textposition": "top center",
    "marker": {
        "size": [math.sqrt(m["Unique_CDR3s"])/4 for m in fig8_metrics],
        "color": [PALETTE.get(m["Sample"],"#888") for m in fig8_metrics],
        "opacity": 0.8
    },
    "customdata": [[m["Unique_CDR3s"], m["Sample"], m["Clusters_per_1k_seqs"] or 0,
                    m["Top1_pct"] or 0] for m in fig8_metrics],
    "hovertemplate": (
        "<b>%{customdata[1]}</b><br>"
        "Shannon H' (abundance): %{x:.3f}<br>"
        "Clonality: %{y:.3f}<br>"
        "Clusters per 1k sequences: %{customdata[2]:.2f}<br>"
        "Largest clone: %{customdata[3]:.2f}% of sequences<br>"
        "Unique CDR3s: %{customdata[0]:,}<extra></extra>"
    )
}]
fig8_layout = {
    "title": {"text": "Repertoire Diversity Landscape", "font": {"size": 16, "color": "#2c3e50"}},
    "xaxis": {"title": "Shannon Diversity Index (clone abundances)"},
    "yaxis": {"title": "Clonality (1 − Pielou evenness)", "range": [0, 1]},
    "showlegend": False,
    "plot_bgcolor": "#fafafa", "paper_bgcolor": "#ffffff"
}

# ── Fig 9: Identity to representative histogram (intra-clonal homogeneity) ─
# One value per *cluster*, not per sequence: Identity is the mean identity of
# the non-representative members of that cluster. Singleton clusters carry NA
# (no member to compare with the representative) and are excluded here rather
# than being coerced to 0, which would create a spurious peak at the origin.
def parse_identity(value):
    """Return the mean identity of a cluster, or None when it is not defined."""
    if value is None:
        return None
    value = str(value).strip()
    if value in ("", "NA", "NaN", "nan", "None"):
        return None
    try:
        return float(value)
    except ValueError:
        return None

ident_by_sample = defaultdict(list)
for r in clusterbig:
    val = parse_identity(r.get("Identity"))
    if val is not None:
        ident_by_sample[r["Sample"]].append(val)

n_informative = sum(len(v) for v in ident_by_sample.values())

fig9_traces = []
for s in samples:
    col = PALETTE.get(s, "#888888")
    vals = ident_by_sample.get(s, [])
    fig9_traces.append({
        "x": vals,
        "type": "histogram",
        "name": f"{s} (n = {len(vals):,})",
        "nbinsx": 50,
        "marker": {"color": col + "BB", "line": {"color": col, "width": 0.5}},
        "opacity": 0.7,
        "hovertemplate": f"<b>{s}</b><br>Mean identity: %{{x:.1f}}%<br>Clusters: %{{y:,}}<extra></extra>"
    })
fig9_layout = {
    "barmode": "overlay",
    "title": {"text": "Mean Member Identity to Cluster Representative", "font": {"size": 16, "color": "#2c3e50"}},
    "xaxis": {"title": "Mean % Identity of Members to Representative",
              "range": [identity_threshold - 1, 100.5]},
    "yaxis": {"title": "Number of Clusters"},
    "legend": {"title": {"text": "Sample"}},
    "plot_bgcolor": "#fafafa", "paper_bgcolor": "#ffffff"
}

# ── Fig 10: Top-10 clusters per sample horizontal bar ─────────────────────
top10_by_sample = defaultdict(list)
for r in clusterbig:
    top10_by_sample[r["Sample"]].append((to_int(r.get("Count",0)), r.get("Representative","?")))

fig10_traces = []
for s in samples:
    top = sorted(top10_by_sample[s], reverse=True)[:10]
    if not top: continue
    col = PALETTE.get(s, "#888888")
    labels = [f"Rank {i+1} ({rep})" for i, (cnt, rep) in enumerate(top)]
    counts = [cnt for cnt, _ in top]
    fig10_traces.append({
        "x": counts,
        "y": labels,
        "name": s,
        "type": "bar",
        "orientation": "h",
        "marker": {"color": col},
        "hovertemplate": f"<b>{s}</b><br>%{{y}}<br>Members: %{{x:,}}<extra></extra>"
    })
fig10_layout = {
    "barmode": "group",
    "title": {"text": "Top-10 Largest Clonotypes per Sample", "font": {"size": 16, "color": "#2c3e50"}},
    "xaxis": {"title": "Cluster Size (sequences)"},
    "yaxis": {"title": "", "autorange": "reversed"},
    "legend": {"title": {"text": "Sample"}},
    "height": 520,
    "plot_bgcolor": "#fafafa", "paper_bgcolor": "#ffffff"
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. HELPER: JSON‐serialize a figure pair
# ─────────────────────────────────────────────────────────────────────────────
def make_fig(fig_id, traces, layout, desc="", title=""):
    return f"""
<div class="chart-card">
  {'<div class="chart-card-header"><h3>' + title + '</h3><p>' + desc + '</p></div>' if title else ''}
  <div class="chart-body" id="{fig_id}" style="height:460px;"></div>
</div>
<script>
  Plotly.newPlot('{fig_id}',
    {json.dumps(traces)},
    {{...{json.dumps(layout)}, margin:{{t:60,l:60,r:30,b:60}}}},
    {{responsive:true, displayModeBar:false}}
  );
</script>
"""

# ─────────────────────────────────────────────────────────────────────────────
# 4. METRICS TABLE
# ─────────────────────────────────────────────────────────────────────────────
def badge(val, thresholds=((25,"badge-high"),(10,"badge-med"),(0,"badge-low"))):
    for t, cls in thresholds:
        if val >= t:
            return f'<span class="{cls}">{val}%</span>'
    return f'<span class="badge-low">{val}%</span>'

table_rows = ""
for m in div_metrics:
    table_rows += (
        f"<tr><td><strong>{m['Sample']}</strong></td>"
        f"<td>{m['Unique_CDR3s']:,}</td>"
        f"<td>{m['Total_Clusters']:,}</td>"
        f"<td>{m['Total_sequences']:,}</td>"
        f"<td>{fmt_metric(m['Clusters_per_1k_seqs'], 2)}</td>"
        f"<td>{fmt_metric(m['Shannon_abundance'])}</td>"
        f"<td>{fmt_metric(m['Clonality'])}</td>"
        f"<td>{fmt_metric(m['Simpson'])}</td>"
        f"<td>{fmt_metric(m['Gini'])}</td>"
        f"<td>{fmt_metric(m['D50'])}</td>"
        f"<td>{fmt_metric(m['Top1_pct'], 2, '%')}</td>"
        f"<td>{fmt_metric(m['Top10_pct'], 2, '%')}</td>"
        f"<td>{m['Shannon_len']:.3f}</td>"
        f"<td>{badge(m['Pct_expanded'])}</td>"
        f"<td>{m['Pct_large']:.2f}%</td></tr>"
    )

metrics_table = f"""
<div style="overflow-x:auto">
<table class="metrics-table">
  <thead><tr>
    <th>Sample</th><th>Unique CDR3s</th><th>Clusters</th><th>Clustered sequences</th>
    <th title="Number of clusters per 1000 clustered sequences: richness normalised for sequencing depth">Clusters / 1k seqs</th>
    <th title="Shannon index computed on clone abundances">Shannon H' (abundance)</th>
    <th title="1 - Pielou evenness of the clone abundances">Clonality</th>
    <th title="Gini-Simpson index, 1 - sum(p²)">Simpson</th>
    <th title="Gini coefficient of the cluster size distribution">Gini</th>
    <th title="Smallest number of clones covering 50% of the sequences">D50</th>
    <th title="Share of the sequences carried by the largest clone">Top1</th>
    <th title="Share of the sequences carried by the ten largest clones">Top10</th>
    <th title="Shannon entropy of the CDR3 length distribution: evenness of lengths, not repertoire diversity">CDR3 length evenness (H')</th>
    <th>% Expanded (≥5)</th><th>% Large (≥1000)</th>
  </tr></thead>
  <tbody>{table_rows}</tbody>
</table>
</div>"""

# ─────────────────────────────────────────────────────────────────────────────
# 5. CSS
# ─────────────────────────────────────────────────────────────────────────────
CSS = """
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:'Inter',sans-serif;background:#f0f4f8;color:#2c3e50}

.report-header{
  background:linear-gradient(135deg,#1a1a2e 0%,#16213e 50%,#0f3460 100%);
  color:white;padding:48px 60px 36px;border-bottom:4px solid #e94560
}
.report-header h1{font-size:2.3em;font-weight:700;margin-bottom:6px}
.report-header .subtitle{font-size:1.05em;opacity:.75;font-weight:300}
.report-header .meta{margin-top:22px;display:flex;gap:32px;flex-wrap:wrap}
.report-header .meta-item{font-size:.85em;opacity:.7}
.report-header .meta-item strong{display:block;font-size:1.2em;opacity:1;color:#e94560}

.nav-bar{
  background:#16213e;position:sticky;top:0;z-index:1000;
  display:flex;padding:0 60px;box-shadow:0 2px 12px rgba(0,0,0,.3);overflow-x:auto
}
.nav-bar a{
  color:rgba(255,255,255,.7);text-decoration:none;padding:13px 15px;
  font-size:.84em;font-weight:500;white-space:nowrap;
  border-bottom:3px solid transparent;transition:.2s
}
.nav-bar a:hover{color:#e94560;border-bottom-color:#e94560}

.container{max-width:1380px;margin:0 auto;padding:40px 60px}

.kpi-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(195px,1fr));gap:18px;margin-bottom:38px}
.kpi-card{
  background:white;border-radius:16px;padding:26px 22px;
  box-shadow:0 2px 14px rgba(0,0,0,.06);border-left:5px solid #e94560;
  transition:transform .2s,box-shadow .2s;cursor:default
}
.kpi-card:hover{transform:translateY(-3px);box-shadow:0 6px 22px rgba(0,0,0,.12)}
.kpi-card:nth-child(2){border-left-color:#3498db}
.kpi-card:nth-child(3){border-left-color:#2ecc71}
.kpi-card:nth-child(4){border-left-color:#f39c12}
.kpi-card:nth-child(5){border-left-color:#9b59b6}
.kpi-label{font-size:.78em;text-transform:uppercase;letter-spacing:.8px;color:#7f8c8d;font-weight:600;margin-bottom:9px}
.kpi-value{font-size:2.15em;font-weight:700;color:#2c3e50;line-height:1}
.kpi-sub{font-size:.8em;color:#95a5a6;margin-top:6px}

.section{margin-bottom:54px}
.section-header{display:flex;align-items:center;gap:14px;margin-bottom:20px;padding-bottom:14px;border-bottom:2px solid #e8ecef}
.section-number{
  background:linear-gradient(135deg,#e94560,#c0392b);color:white;
  width:38px;height:38px;border-radius:10px;display:flex;align-items:center;
  justify-content:center;font-size:1em;font-weight:700;flex-shrink:0
}
.section-header h2{font-size:1.45em;font-weight:600;color:#2c3e50}
.section-intro{
  background:#eef2f7;border-left:4px solid #3498db;padding:14px 18px;
  border-radius:0 8px 8px 0;font-size:.9em;color:#34495e;margin-bottom:20px;line-height:1.65
}
.section-intro cite{font-style:italic;color:#95a5a6;font-size:.92em}
.nav-num {
  background: linear-gradient(135deg, #e94560, #c0392b);
  color: white; width: 22px; height: 22px; border-radius: 6px;
  display: inline-flex; align-items: center; justify-content: center;
  font-size: 0.75em; font-weight: 700; margin-right: 8px; vertical-align: middle;
}

.chart-card{background:white;border-radius:16px;box-shadow:0 2px 14px rgba(0,0,0,.06);overflow:hidden;margin-bottom:22px}
.chart-card-header{padding:18px 22px 0}
.chart-card-header h3{font-size:1em;font-weight:600;color:#2c3e50}
.chart-card-header p{font-size:.84em;color:#7f8c8d;margin-top:4px;line-height:1.5}
.chart-body{padding:8px}

.metrics-table{width:100%;border-collapse:collapse;font-size:.88em}
.metrics-table thead tr{background:#2c3e50;color:white}
.metrics-table th{padding:11px 15px;text-align:left;font-weight:600;font-size:.84em;letter-spacing:.4px}
.metrics-table td{padding:10px 15px;border-bottom:1px solid #f0f4f8}
.metrics-table tbody tr:nth-child(even){background:#f8f9fb}
.metrics-table tbody tr:hover{background:#eef2f7}
.badge-high{display:inline-block;padding:3px 10px;border-radius:20px;font-size:.78em;font-weight:600;background:#fde8e8;color:#c0392b}
.badge-med {display:inline-block;padding:3px 10px;border-radius:20px;font-size:.78em;font-weight:600;background:#fef6e4;color:#d68910}
.badge-low {display:inline-block;padding:3px 10px;border-radius:20px;font-size:.78em;font-weight:600;background:#e8f8f0;color:#27ae60}

.grid-2{display:grid;grid-template-columns:1fr 1fr;gap:22px}
.methods-grid{display:grid;grid-template-columns:1fr 1fr;gap:22px}
.methods-card{background:white;border-radius:16px;padding:26px;box-shadow:0 2px 14px rgba(0,0,0,.06)}
.methods-card h3{font-size:1em;font-weight:600;margin-bottom:10px;color:#2c3e50}
.methods-card p{color:#555;font-size:.87em;line-height:1.7;margin-bottom:8px}
.methods-card ul{padding-left:18px;color:#555;font-size:.85em;line-height:1.85}
.ref-list{font-size:.82em;color:#7f8c8d;list-style:none;line-height:1.9}
.ref-list li::before{content:"📄 "}

.footer{background:#2c3e50;color:rgba(255,255,255,.55);text-align:center;padding:26px 60px;font-size:.81em;margin-top:38px}

@media(max-width:900px){
  .container{padding:22px 18px}
  .report-header{padding:30px 18px}
  .report-nav{padding:0 18px}
  .grid-2,.methods-grid{grid-template-columns:1fr}
  .card-grid{grid-template-columns:repeat(2,1fr)}
}
"""

# ─────────────────────────────────────────────────────────────────────────────
# 6. ASSEMBLE HTML
# ─────────────────────────────────────────────────────────────────────────────
now = datetime.now().strftime("%B %d, %Y")

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width,initial-scale=1.0">
  <title>Nanobody Repertoire Report</title>
  <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>{CSS}</style>
</head>
<body>

<!-- HEADER ─────────────────────────────────────────── -->
<div class="report-header">
  <h1>Nanorepertoire</h1>
  <div class="subtitle">Nanobody Repertoire Analysis (VHH Sequencing · Clonotype &amp; CDR3 Analysis Report)</div>
  <div class="meta" style="margin-top:20px; border-top:1px solid rgba(255,255,255,0.1); padding-top:15px; display:flex; gap:30px; flex-wrap:wrap;">
    <div class="meta-item" style="font-size:0.85em; opacity:0.8; max-width:400px;">
      <strong>Pipeline Execution Metrics</strong> 
      <div style="font-size:0.9em; margin-top:4px; line-height: 1.4;">
        If you are interested in detailed Start/End Times, CPU Efficiency, and CO₂ metrics, please see the dedicated reports produced upon completion of the pipeline. We sincerely thank the authors of the <strong>nf-co2footprint plugin</strong> and the Green Algorithms project for enabling these measurements.
      </div>
    </div>
    <div class="meta-item" style="margin-left:auto; text-align:right; font-size:0.85em; opacity:0.8;">
      <strong>Report Generated</strong>
      <div style="margin-top:4px;">{now}</div>
    </div>
  </div>
</div>

<!-- NAVIGATION ─────────────────────────────────────── -->

<!-- NAVIGATION ─────────────────────────────────────── -->
<nav class="nav-bar">
  <a href="#summary">Summary</a>
  <a href="#clusters"><span class="nav-num">1</span>Clonal Clusters</a>
  <a href="#cdr3"><span class="nav-num">2</span>CDR3 Diversity</a>
  <a href="#aminoacids"><span class="nav-num">3</span>Amino Acids</a>
  <a href="#clonality"><span class="nav-num">4</span>Clonality</a>
  <a href="#methods"><span class="nav-num">5</span>Methods</a>
</nav>

<div class="container">

<!-- SECTION 0: KPI CARDS ══════════════════════════════════════════════════ -->
<div id="summary" class="section">
  <div class="section-header" style="border-bottom:none; margin-bottom:5px;">
    <h2>Summary</h2>
  </div>
  <div class="kpi-grid">
    <div class="kpi-card">
      <div class="kpi-label">Sequences with a called CDR3</div>
      <div class="kpi-value">{seqs_with_cdr3:,}</div>
      <div class="kpi-sub">Total across all samples</div>
    </div>
    <div class="kpi-card">
      <div class="kpi-label">Total Clusters (CD-HIT)</div>
      <div class="kpi-value">{total_clusters:,}</div>
      <div class="kpi-sub">{identity_pct} identity threshold</div>
    </div>
    <div class="kpi-card">
      <div class="kpi-label">Unique CDR3 Sequences</div>
      <div class="kpi-value">{total_cdr3:,}</div>
      <div class="kpi-sub">Distinct paratopes</div>
    </div>
    <div class="kpi-card">
      <div class="kpi-label">Mean CDR3 Length</div>
      <div class="kpi-value">{mean_len:.1f} AA</div>
      <div class="kpi-sub">SD: {sd_len:.1f} · Median: {median_len}</div>
    </div>
    <div class="kpi-card">
      <div class="kpi-label">Samples Analysed</div>
      <div class="kpi-value">{len(samples)}</div>
      <div class="kpi-sub">Independent libraries</div>
    </div>
  </div>
</div>

<!-- ══════════════════════════════════════════════════
     SECTION 1: CLONAL CLUSTERS
══════════════════════════════════════════════════ -->
<div id="clusters" class="section">
  <div class="section-header">
    <div class="section-number">1</div>
    <h2>Clonal Cluster Analysis</h2>
  </div>
  <div class="section-intro">
    Amino acid sequences were clustered using <strong>CD-HIT</strong> at {identity_pct} identity, following the landmark approach of Deschaght et al. (2017).
    <strong>Expanded clonotypes</strong> (≥5 members) represent B‑cell lineages that have undergone antigen‑driven somatic expansion.
    Large clusters (≥1000 members) are dominant clonal responses — prime candidates for high‑affinity nanobodies.
    The cluster size distribution mirrors the "clonotype expansion" readout of tools like MiXCR and IMGT/VQuest.
    <br><cite>Deschaght et al. 2017. Front. Immunol. doi:10.3389/fimmu.2017.00420 · Bolotin et al. 2015 (MiXCR). Nature Methods.</cite>
  </div>

  <div style="margin-bottom:20px">{metrics_table}</div>

  {make_fig("fig1", fig1_traces, fig1_layout,
    "Breakdown of clusters by size class per sample. Singletons represent rare or noisy sequences. "
    "Expanded clones (≥5 members, orange/red) indicate B-cell lineages selected by antigen exposure.",
    "Clonal Expansion Profile")}

  <div class="grid-2">
    {make_fig("fig2", fig2_traces, fig2_layout,
      "Cluster counts at four size thresholds (log scale). Reveals the hierarchical repertoire structure from singletons to dominant clones.",
      "Cluster Abundance at Size Thresholds")}
    {make_fig("fig10", fig10_traces, fig10_layout,
      "The 10 most abundant clusters per sample. Dominant clones in immunized animals likely represent the most affinity-matured antigen binders.",
      "Top-10 Largest Clonotypes per Sample")}
  </div>
</div>

<!-- ══════════════════════════════════════════════════
     SECTION 2: CDR3 DIVERSITY
══════════════════════════════════════════════════ -->
<div id="cdr3" class="section">
  <div class="section-header">
    <div class="section-number">2</div>
    <h2>CDR3 Diversity Analysis</h2>
  </div>
  <div class="section-intro">
    The <strong>CDR3</strong> (Complementarity-Determining Region 3) is the primary antigen-binding loop of VHH nanobodies.
    VHH CDR3s are typically <em>longer and more diverse</em> than VH CDR3s in conventional antibodies
    (VHH mean ~19 AA vs VH ~12 AA — Spinelli et al. 2022), enabling access to deep epitope cavities and enzyme active sites.
    CDR3 length diversity correlates with repertoire breadth; antigen-enriched samples may show length skewing toward specific structural motifs.
    Shorter CDR3s (&lt;12 AA) tend to form flat β-strand paratopes; longer CDR3s (&gt;16 AA) adopt protruding finger-like conformations.
    <br><cite>Spinelli et al. 2022. Front. Immunol. doi:10.3389/fimmu.2022.927966 · Muyldermans 2013. Annu. Rev. Biochem.</cite>
  </div>
  <div class="grid-2">
    {make_fig("fig5", fig5_traces, fig5_layout,
      "Total number of distinct CDR3 sequences per sample. Higher numbers indicate broader repertoire diversity.",
      "Unique CDR3 Sequences per Sample")}
    {make_fig("fig3", fig3_traces, fig3_layout,
      "Violin plot of CDR3 length distribution per sample. Box shows IQR; line shows median. VHH CDR3s typically peak at 15–23 AA.",
      "CDR3 Length Distribution (Violin)")}
  </div>
  {make_fig("fig4", fig4_traces, fig4_layout,
    "Area chart of CDR3 length frequency per sample. Peaks at specific lengths may indicate dominant structural motifs or antigen-driven convergence. "
    "Short CDR3s (&lt;12 AA): flat β-strand paratopes. Long CDR3s (&gt;16 AA): protruding loops for cavity/groove binding.",
    "CDR3 Length Frequency Profile")}
</div>

<!-- ══════════════════════════════════════════════════
     SECTION 3: AMINO ACID COMPOSITION
══════════════════════════════════════════════════ -->
<div id="aminoacids" class="section">
  <div class="section-header">
    <div class="section-number">3</div>
    <h2>CDR3 Amino Acid Composition</h2>
  </div>
  <div class="section-intro">
    The amino acid (AA) composition of CDR3 loops reflects antigen-driven selection.
    VHH CDR3s are characteristically enriched in <strong>Tyrosine (Y)</strong> and <strong>Serine (S)</strong> — enabling hydrogen bonding and Van der Waals contacts.
    <strong>Cysteine (C)</strong> pairs can form intra-CDR3 disulfide bonds with CDR1, creating a structural scaffold that extends paratope reach.
    Charged residues (R, D, E, K) mediate electrostatic complementarity. Convergent AA enrichment across samples suggests shared antigen pressure.
    <br><cite>Desmyter et al. 2002. J. Biol. Chem. doi:10.1074/jbc.D200025200 · Mitchell &amp; Colwell 2018. Proteins.</cite>
  </div>
  {make_fig("fig6", fig6_traces, fig6_layout,
    "Heatmap of amino acid frequency (% of all CDR3 residues) per sample. Dark = high frequency. Key VHH residues: Y, S, G, T. Cysteine (C) involvement suggests disulfide-mediated loop stabilization.",
    "Amino Acid Frequency Heatmap")}
  {make_fig("fig7", fig7_traces, fig7_layout,
    "Side-by-side comparison of amino acid usage per sample. Differences across samples may reflect divergent antigen-binding strategies or B-cell selection biases.",
    "Side-by-Side Amino Acid Usage")}
</div>

<!-- ══════════════════════════════════════════════════
     SECTION 4: CLONALITY & SHM
══════════════════════════════════════════════════ -->
<div id="clonality" class="section">
  <div class="section-header">
    <div class="section-number">4</div>
    <h2>Clonality &amp; Intra-clonal Homogeneity</h2>
  </div>
  <div class="section-intro">
    The mean identity of the members of a cluster to its representative is a <strong>descriptive measure of
    intra-clonal homogeneity</strong>: clusters close to 100% are internally uniform, clusters closer to the
    clustering threshold group sequences that differ more from their representative.
    <strong>This quantity is truncated by construction</strong>: CD-HIT admits a sequence into a cluster only if
    its identity to the representative is at least {identity_pct}, so no value below that threshold can be
    observed and the distribution is bounded on the left by the clustering parameter itself.
    It is therefore <em>not</em> an estimate of somatic hypermutation, and it is not equivalent to the V-gene
    mutation frequency reported by MiXCR or IMGT/V-QUEST, which is measured against an inferred germline
    reference rather than against an empirically chosen cluster member; the analogy is conceptual only.
    Sequence-level differences within a cluster may equally reflect somatic mutation, PCR and sequencing error,
    or repeated sampling of the same molecule, which this workflow does not deduplicate and therefore cannot tell apart.
    <br><cite>Li &amp; Godzik 2006 (CD-HIT). Bioinformatics. doi:10.1093/bioinformatics/btl158 · Bolotin et al. 2015 (MiXCR). Nature Methods. doi:10.1038/nmeth.3364</cite>
  </div>
  <div class="grid-2">
    {make_fig("fig9", fig9_traces, fig9_layout,
      f"Distribution of the mean identity of cluster members to their representative, one value per cluster "
      f"(n = {n_informative:,} multi-member clusters; singletons are excluded because the quantity is undefined for them). "
      f"Values are bounded below by the {identity_pct} clustering threshold, so the panel describes homogeneity within that window only.",
      "Mean Member Identity to Cluster Representative")}
    {make_fig("fig8", fig8_traces, fig8_layout,
      "Bubble chart: Shannon index of the clone abundances (x) against clonality, 1 − Pielou evenness (y); "
      "bubble size = unique CDR3 count. Samples move to the lower right as sequences are spread evenly over many "
      "clones, and to the upper left as they concentrate into a few dominant lineages. Both axes are computed on "
      "cluster sizes, so they are comparable across samples only alongside the depth-normalised richness "
      "(clusters per 1000 sequences) reported in the table above.",
      "Repertoire Diversity Landscape")}
  </div>
</div>

<!-- ══════════════════════════════════════════════════
     SECTION 5: METHODS
══════════════════════════════════════════════════ -->
<div id="methods" class="section">
  <div class="section-header">
    <div class="section-number">5</div>
    <h2>Methods Summary</h2>
  </div>
  <div class="methods-grid">

    <div class="methods-card">
      <h3>🔬 Pre-processing &amp; Translation</h3>
      <p>Paired-end FASTQ reads were adapter-trimmed with <strong>Cutadapt</strong> and merged with <strong>FLASH</strong> (max overlap 300 bp). 
      Merged reads were translated in-silico using conserved VHH framework start/end primer motifs (ATG…TVSS pattern; Deschaght 2017).</p>
      <p>CDR3 regions were extracted with <strong>nanocdr-x</strong>, a deep learning model trained on camelid VHH sequences.</p>
    </div>

    <div class="methods-card">
      <h3>🧩 CD-HIT Clustering</h3>
      <p>Translated amino acid sequences were clustered with <strong>CD-HIT</strong>:</p>
      <ul>
        <li>Identity threshold: {identity_threshold/100:.2f} ({identity_pct})</li>
        <li>Word size: {word_size}</li>
        <li>Global identity alignment (CD-HIT default <code>-G 1</code>); no minimum alignment coverage is imposed</li>
        <li>Four size thresholds analysed: 1, ≥5, ≥100, ≥1000</li>
      </ul>
      <p>This follows the validated approach of Deschaght et al. (2017) for large nanobody repertoires.</p>
    </div>

    <div class="methods-card">
      <h3>📐 Diversity Metrics</h3>
      <p><strong>Shannon H' (abundance):</strong> −Σ p<sub>i</sub> ln p<sub>i</sub> over the clone frequencies p<sub>i</sub> = cluster size / clustered sequences. <strong>Clonality</strong> is 1 − H'/ln(S), with S the number of clusters: 0 when every clone is equally abundant, approaching 1 as the library concentrates into one lineage (undefined, and reported as NA, when S = 1). <strong>Simpson</strong> is the Gini–Simpson index 1 − Σ p<sub>i</sub>², <strong>Gini</strong> the Gini coefficient of the cluster sizes, <strong>D50</strong> the smallest number of clones covering half of the sequences, and <strong>Top1/Top10</strong> the share of the library held by the largest one or ten clones.</p>
      <p><strong>Clusters per 1000 sequences:</strong> richness normalised by sequencing depth. Absolute cluster counts are not comparable between libraries sequenced at different depths; this ratio is.</p>
      <p><strong>CDR3 length evenness (H'):</strong> Shannon entropy of the CDR3 <em>length</em> distribution. It measures how evenly paratope lengths are spread and is reported separately from the abundance-based indices above, which it must not be confused with.</p>
      <p><strong>% Expanded Clonotypes:</strong> fraction of clusters with ≥5 members; reflects the proportion of somatically expanded B-cell lineages.</p>
      <p><strong>Intra-clonal homogeneity:</strong> distribution, across clusters, of the mean percent identity of the non-representative members to the cluster representative. Singleton clusters are excluded, the quantity being undefined for them. It is a descriptive summary bounded below by the {identity_pct} clustering threshold, not a germline-referenced mutation rate.</p>
    </div>

    <div class="methods-card">
      <h3>📚 Key References</h3>
      <ul class="ref-list">
        <li>Deschaght et al. 2017. <em>Front. Immunol.</em> doi:10.3389/fimmu.2017.00420</li>
        <li>Spinelli et al. 2022. <em>Front. Immunol.</em> doi:10.3389/fimmu.2022.927966</li>
        <li>Muyldermans S. 2013. <em>Annu. Rev. Biochem.</em> doi:10.1146/annurev-biochem-071812-105327</li>
        <li>Bolotin et al. 2015. <em>Nature Methods</em>. doi:10.1038/nmeth.3364</li>
        <li>Desmyter et al. 2002. <em>J. Biol. Chem.</em> doi:10.1074/jbc.D200025200</li>
        <li>Mitchell &amp; Colwell 2018. <em>Proteins</em>. doi:10.1002/prot.25497</li>
        <li>Li &amp; Godzik 2006. <em>Bioinformatics</em>. doi:10.1093/bioinformatics/btl158 (CD-HIT)</li>
      </ul>
    </div>

  </div>
</div>

</div><!-- /container -->

<div class="footer">
  Nanobody Repertoire Report &nbsp;·&nbsp; lescailab/nanorepertoire pipeline &nbsp;·&nbsp; {now}<br>
  <span style="opacity:.5">Based on Deschaght 2017 framework · Charts powered by Plotly.js 2.27</span>
</div>

</body>
</html>"""

print(f"Writing report to {output_html} ...")
with open(output_html, "w", encoding="utf-8") as f:
    f.write(html)
print(f"✅ Done! Open: {output_html}")
