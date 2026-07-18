from __future__ import annotations

import html
import json
import platform
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .data import PreparedData
from .exporters import export_all
from .search import SearchResult


def _json_safe(value):
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _feature_scores(data: PreparedData, result: SearchResult) -> pd.DataFrame:
    rows = []
    for feature in data.features.columns:
        containing = [item for item in result.interactions if feature in item.features]
        significant = [item for item in containing if item.significant]
        values = data.features[feature].to_numpy(dtype=float)
        if values.min() >= 0 and values.max() <= 2:
            frequency = float(values.mean() / 2.0)
            maf = min(frequency, 1.0 - frequency)
        else:
            maf = np.nan
        rows.append(
            {
                "feature": feature,
                "interaction_count": len(containing),
                "significant_interaction_count": len(significant),
                "fi_contribution": sum(item.fi_gain / item.order for item in containing),
                "max_fi_gain": max((item.fi_gain for item in containing), default=0.0),
                "min_p_value": min((item.p_value for item in containing), default=1.0),
                "min_q_value": min((item.q_value for item in containing), default=1.0),
                "minor_allele_frequency": maf,
                "gene": "",
                "pathway": "",
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["significant_interaction_count", "fi_contribution"], ascending=[False, False]
    )


def _summary(result: SearchResult, input_path: str | None) -> dict:
    return {
        "schema_version": "1.0",
        "software": {"name": "FIGHI", "version": __version__},
        "created_at": datetime.now(timezone.utc).isoformat(),
        "input": input_path,
        "trait": result.trait,
        "phenotype_mapping": result.phenotype_mapping,
        "quality_control": result.qc,
        "screened_atom_count": len(result.atoms),
        "evaluated_interactions": len(result.interactions),
        "significant_interactions": len(result.significant),
        "evaluated_by_order": result.evaluated_by_order,
        "runtime_seconds": result.runtime_seconds,
        "stopped_reason": result.stopped_reason,
        "warnings": result.warnings,
        "configuration": result.config.to_dict(),
        "interpretation": (
            "Significance is based on a one-degree-of-freedom efficient score test of the "
            "highest-order product, conditional on covariates and all lower-order terms "
            "within that candidate. Adjusted p-values are global across evaluated candidates."
        ),
    }


def _html_report(outdir: Path, result: SearchResult, summary: dict, plot_paths: dict[str, str]) -> Path:
    top = result.to_frame().head(50)
    rows = []
    for _, row in top.iterrows():
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(row['hyperedge']))}</td>"
            f"<td>{int(row['order'])}</td>"
            f"<td>{row['fi_gain']:.4g}</td>"
            f"<td>{row['p_value']:.3g}</td>"
            f"<td>{row['q_value']:.3g}</td>"
            f"<td>{'Yes' if row['significant'] else 'No'}</td>"
            "</tr>"
        )
    warnings = "".join(f"<li>{html.escape(item)}</li>" for item in result.warnings)
    images = "".join(
        f'<figure><img src="plots/{Path(path).name}" alt="{html.escape(name)}"><figcaption>{html.escape(name.replace("_", " ").title())}</figcaption></figure>'
        for name, path in plot_paths.items()
    )
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>FIGHI analysis report</title>
<style>
:root{{--ink:#0f172a;--muted:#475569;--line:#dbe3ec;--accent:#0f766e;--soft:#f1f5f9;}}
*{{box-sizing:border-box}}body{{margin:0;background:#f8fafc;color:var(--ink);font:16px/1.55 system-ui,sans-serif}}
main{{max-width:1120px;margin:auto;padding:42px 24px 72px}}header{{background:linear-gradient(135deg,#0f172a,#115e59);color:white;padding:36px;border-radius:18px}}
h1{{margin:0 0 8px;font-size:2.2rem}}h2{{margin-top:36px}}.meta{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:14px;margin-top:24px}}
.metric{{background:white;color:var(--ink);padding:14px 16px;border-radius:12px}}.metric b{{display:block;font-size:1.6rem;color:var(--accent)}}
.card{{background:white;border:1px solid var(--line);border-radius:14px;padding:22px;margin-top:18px;overflow:auto}}
table{{border-collapse:collapse;width:100%;font-size:.92rem}}th,td{{padding:10px;border-bottom:1px solid var(--line);text-align:left}}th{{background:var(--soft)}}
.plots{{display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:18px}}figure{{margin:0;background:white;border:1px solid var(--line);border-radius:14px;padding:12px}}
img{{width:100%;height:auto}}figcaption{{color:var(--muted);padding:8px}}code{{background:var(--soft);padding:2px 5px;border-radius:4px}}
</style></head><body><main>
<header><div>FIGHI v{__version__}</div><h1>Interaction discovery report</h1><p>Fisher-information-guided evidence with explicit statistical control.</p>
<div class="meta"><div class="metric"><b>{result.qc['samples_analyzed']}</b>samples</div><div class="metric"><b>{result.qc['features_analyzed']}</b>QC-passed features</div><div class="metric"><b>{len(result.interactions)}</b>interactions tested</div><div class="metric"><b>{len(result.significant)}</b>significant interactions</div></div></header>
<section class="card"><h2>Interpretation</h2><p>{html.escape(summary['interpretation'])}</p><p>Fisher-information gain is an evidence-ranking quantity; it is not, by itself, a significance decision. The <code>q_value</code> and configured alpha determine significance.</p></section>
<h2>Diagnostics</h2><div class="plots">{images or '<div class="card">No diagnostic plots were generated.</div>'}</div>
<section class="card"><h2>Top evaluated interactions</h2><table><thead><tr><th>Features</th><th>Order</th><th>FI gain</th><th>p</th><th>adjusted p</th><th>Significant</th></tr></thead><tbody>{''.join(rows) or '<tr><td colspan="6">No candidates evaluated.</td></tr>'}</tbody></table></section>
<section class="card"><h2>Quality and limitations</h2><ul>{warnings or '<li>No run warnings.</li>'}</ul><p>Screening, hierarchical expansion, model assumptions, population structure, LD, batch effects, and independent replication all affect interpretation. See the Methods and Limitations documentation before making biological claims.</p></section>
</main></body></html>"""
    path = outdir / "fighi_report.html"
    path.write_text(document, encoding="utf-8")
    return path


def write_outputs(
    outdir: str | Path,
    result: SearchResult,
    data: PreparedData,
    input_path: str | None = None,
    command: list[str] | None = None,
    plots: bool = True,
) -> dict[str, str]:
    target = Path(outdir)
    target.mkdir(parents=True, exist_ok=True)
    frame = result.to_frame().drop(columns=["features"], errors="ignore")
    results_path = target / "fighi_results.csv"
    frame.to_csv(results_path, index=False)
    significant_path = target / "fighi_significant_interactions.csv"
    frame.loc[frame["significant"]].to_csv(significant_path, index=False)
    feature_path = target / "fighi_feature_scores.csv"
    _feature_scores(data, result).to_csv(feature_path, index=False)

    summary = _summary(result, input_path)
    summary_path = target / "fighi_summary.json"
    summary_path.write_text(json.dumps(_json_safe(summary), indent=2), encoding="utf-8")
    model = {
        "schema_version": "1.0",
        "software_version": __version__,
        "trait": result.trait,
        "configuration": result.config.to_dict(),
        "screened_features": result.atoms,
        "significant_interactions": [item.to_dict(result.trait) for item in result.significant],
    }
    model_path = target / "fighi_model.json"
    model_path.write_text(json.dumps(_json_safe(model), indent=2), encoding="utf-8")
    provenance = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "command": shlex.join(command) if command else None,
        "python": sys.version,
        "platform": platform.platform(),
        "software_version": __version__,
    }
    provenance_path = target / "fighi_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    graph_paths = export_all(target, result.interactions, result.config.graph_top)
    if plots:
        from .plotting import make_plots

        plot_paths = make_plots(target, result)
    else:
        plot_paths = {}
    report_path = _html_report(target, result, summary, plot_paths)
    paths = {
        "results": str(results_path),
        "significant_interactions": str(significant_path),
        "feature_scores": str(feature_path),
        "summary": str(summary_path),
        "model": str(model_path),
        "provenance": str(provenance_path),
        "report": str(report_path),
        **graph_paths,
        **plot_paths,
    }
    return paths
