from __future__ import annotations

import html
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .errors import BenchmarkError
from .statistics import adjust_pvalues
from .utilities import sha256_file, write_json

_FORMATS = {"fighi", "plink_epistasis", "generic"}


def _load_manifest(path: str | Path) -> tuple[dict[str, Any], Path]:
    manifest_path = Path(path).expanduser().resolve()
    if not manifest_path.is_file():
        raise BenchmarkError(f"Benchmark manifest does not exist: {manifest_path}")
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BenchmarkError(f"Could not read benchmark manifest: {exc}") from exc
    if not isinstance(payload, dict):
        raise BenchmarkError("Benchmark manifest must contain a JSON object")
    return payload, manifest_path


def validate_benchmark_manifest(path: str | Path) -> dict[str, Any]:
    """Validate a benchmark manifest and return an auditable summary."""
    payload, manifest_path = _load_manifest(path)
    problems: list[str] = []
    analysis_id = payload.get("analysis_id")
    if not isinstance(analysis_id, str) or not analysis_id.strip():
        problems.append("analysis_id must be a non-empty string")
    correction = payload.get("correction", "fdr_bh")
    if correction not in {"fdr_bh", "bonferroni", "none"}:
        problems.append("correction must be fdr_bh, bonferroni, or none")
    alpha = payload.get("alpha", 0.05)
    if not isinstance(alpha, (int, float)) or not 0 < float(alpha) <= 1:
        problems.append("alpha must satisfy 0 < alpha <= 1")
    top_n = payload.get("top_n", 100)
    if not isinstance(top_n, int) or isinstance(top_n, bool) or top_n < 1:
        problems.append("top_n must be a positive integer")
    methods = payload.get("methods")
    if not isinstance(methods, list) or not methods:
        problems.append("methods must be a non-empty list")
        methods = []

    names: list[str] = []
    for index, method in enumerate(methods):
        label = f"methods[{index}]"
        if not isinstance(method, dict):
            problems.append(f"{label} must be an object")
            continue
        name = method.get("name")
        if not isinstance(name, str) or not name.strip():
            problems.append(f"{label}.name must be a non-empty string")
        elif not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
            problems.append(
                f"{label}.name may contain only letters, numbers, dot, dash, underscore"
            )
        else:
            names.append(name)
        command = method.get("command")
        enabled = method.get("enabled", True)
        if not isinstance(enabled, bool):
            problems.append(f"{label}.enabled must be true or false")
        if enabled and (
            not isinstance(command, list)
            or not command
            or not all(isinstance(x, str) for x in command)
        ):
            problems.append(f"{label}.command must be a non-empty string array")
        result = method.get("result")
        if enabled and not isinstance(result, dict):
            problems.append(f"{label}.result must be an object")
        elif isinstance(result, dict) and result.get("format", "generic") not in _FORMATS:
            problems.append(f"{label}.result.format must be one of {sorted(_FORMATS)}")
        elif enabled and isinstance(result, dict) and not isinstance(result.get("path"), str):
            problems.append(f"{label}.result.path must be a string")
        timeout = method.get("timeout_seconds")
        if timeout is not None and (
            not isinstance(timeout, (int, float))
            or isinstance(timeout, bool)
            or float(timeout) <= 0
        ):
            problems.append(f"{label}.timeout_seconds must be positive")
        if "env" in method and not isinstance(method["env"], dict):
            problems.append(f"{label}.env must be an object")
        if "cwd" in method and not isinstance(method["cwd"], str):
            problems.append(f"{label}.cwd must be a string")
    if len(names) != len(set(names)):
        problems.append("method names must be unique")

    shared = payload.get("shared_inputs", {})
    if shared is not None and not isinstance(shared, dict):
        problems.append("shared_inputs must be an object")
        shared = {}
    input_files: dict[str, dict[str, Any]] = {}
    for key, value in (shared or {}).items():
        if not isinstance(value, str) or not value:
            problems.append(f"shared_inputs.{key} must be a path string")
            continue
        resolved = (
            (manifest_path.parent / value).resolve()
            if not Path(value).is_absolute()
            else Path(value)
        )
        if not resolved.is_file():
            problems.append(f"shared input does not exist: {resolved}")
        else:
            input_files[key] = {"path": str(resolved), "sha256": sha256_file(resolved)}

    if problems:
        raise BenchmarkError("Invalid benchmark manifest:\n- " + "\n- ".join(problems))
    return {
        "valid": True,
        "manifest": str(manifest_path),
        "analysis_id": analysis_id,
        "enabled_methods": [method["name"] for method in methods if method.get("enabled", True)],
        "disabled_methods": [
            method["name"] for method in methods if not method.get("enabled", True)
        ],
        "shared_inputs": input_files,
    }


def _placeholders(
    payload: dict[str, Any], manifest_path: Path, outdir: Path, run_dir: Path
) -> dict[str, str]:
    values = {
        "manifest_dir": str(manifest_path.parent),
        "outdir": str(outdir),
        "run_dir": str(run_dir),
        "python_executable": sys.executable,
    }
    for key, value in payload.get("shared_inputs", {}).items():
        path = Path(value)
        values[key] = str(
            (manifest_path.parent / path).resolve() if not path.is_absolute() else path
        )
    return values


def _format_value(value: str, placeholders: dict[str, str]) -> str:
    try:
        return value.format_map(placeholders)
    except KeyError as exc:
        raise BenchmarkError(f"Unknown command placeholder: {exc.args[0]}") from exc


def _terminate_process(process: subprocess.Popen) -> None:
    """Terminate a timed-out comparator and its process group where supported."""
    if process.poll() is not None:
        return
    try:
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGTERM)
        else:  # pragma: no cover - Windows-specific
            process.terminate()
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGKILL)
            else:  # pragma: no cover - Windows-specific
                process.kill()
        except ProcessLookupError:
            return
        process.wait()


def _read_measurement(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def _executable_available(value: str, cwd: Path) -> bool:
    path = Path(value)
    if path.is_absolute():
        return path.is_file()
    if len(path.parts) > 1:
        return (cwd / path).is_file()
    return shutil.which(value) is not None


def _version(method: dict[str, Any], placeholders: dict[str, str], cwd: Path) -> str:
    command = method.get("version_command")
    if not isinstance(command, list) or not command:
        return "not recorded"
    argv = [_format_value(str(value), placeholders) for value in command]
    if not _executable_available(argv[0], cwd):
        return "unavailable"
    try:
        result = subprocess.run(
            argv, cwd=cwd, text=True, capture_output=True, timeout=30, check=False
        )
    except (OSError, subprocess.TimeoutExpired):
        return "unavailable"
    output = (result.stdout or result.stderr).strip().splitlines()
    return output[0][:500] if output else f"exit {result.returncode}"


def _run_method(
    method: dict[str, Any],
    payload: dict[str, Any],
    manifest_path: Path,
    outdir: Path,
) -> dict[str, Any]:
    name = method["name"]
    run_dir = outdir / "runs" / name
    run_dir.mkdir(parents=True, exist_ok=True)
    placeholders = _placeholders(payload, manifest_path, outdir, run_dir)
    argv = [_format_value(value, placeholders) for value in method["command"]]
    cwd_value = _format_value(str(method.get("cwd", "{manifest_dir}")), placeholders)
    cwd = Path(cwd_value).resolve()
    if not cwd.is_dir():
        raise BenchmarkError(f"Working directory for {name} does not exist: {cwd}")
    record: dict[str, Any] = {
        "method": name,
        "status": "running",
        "command": argv,
        "cwd": str(cwd),
        "version": _version(method, placeholders, cwd),
        "target": method.get("target", "not specified"),
        "notes": method.get("notes", ""),
        "return_code": None,
        "wall_seconds": None,
        "user_seconds": None,
        "system_seconds": None,
        "peak_rss_kb": None,
    }
    if not _executable_available(argv[0], cwd):
        record["status"] = "unavailable"
        record["error"] = f"Executable not found: {argv[0]}"
        return record

    environment = os.environ.copy()
    for key, value in method.get("env", {}).items():
        environment[str(key)] = _format_value(str(value), placeholders)
    stdout_path = run_dir / "stdout.log"
    stderr_path = run_dir / "stderr.log"
    measurement_path = run_dir / "resources.json"
    measured_argv = [
        sys.executable,
        str(Path(__file__).with_name("measurement.py")),
        str(measurement_path),
        *argv,
    ]
    start = time.perf_counter()
    try:
        with (
            stdout_path.open("w", encoding="utf-8") as stdout,
            stderr_path.open("w", encoding="utf-8") as stderr,
        ):
            process = subprocess.Popen(
                measured_argv,
                cwd=cwd,
                env=environment,
                stdout=stdout,
                stderr=stderr,
                start_new_session=os.name == "posix",
            )
            try:
                return_code = process.wait(timeout=method.get("timeout_seconds"))
            except subprocess.TimeoutExpired:
                _terminate_process(process)
                raise
        record["return_code"] = return_code
        record["status"] = "completed" if return_code == 0 else "failed"
    except subprocess.TimeoutExpired:
        record["status"] = "timeout"
        record["error"] = "Configured timeout exceeded"
    except OSError as exc:
        record["status"] = "failed"
        record["error"] = str(exc)
    measured = _read_measurement(measurement_path)
    for key in ("wall_seconds", "user_seconds", "system_seconds", "peak_rss_kb"):
        if key in measured:
            record[key] = measured[key]
    record["resource_source"] = measured.get("source", "parent_wall_clock_only")
    if record["wall_seconds"] is None:
        record["wall_seconds"] = time.perf_counter() - start
    record["stdout"] = str(stdout_path)
    record["stderr"] = str(stderr_path)
    record["resources"] = str(measurement_path) if measurement_path.exists() else None
    record["result_spec"] = method["result"]
    record["placeholders"] = placeholders
    return record


def _column(frame: pd.DataFrame, requested: str) -> str:
    if requested in frame:
        return requested
    lookup = {str(name).lower(): str(name) for name in frame.columns}
    if requested.lower() in lookup:
        return lookup[requested.lower()]
    raise BenchmarkError(f"Result column not found: {requested}")


def _normalize_hyperedge(values: list[Any]) -> str:
    features = sorted({str(value).strip() for value in values if str(value).strip()})
    if len(features) < 2:
        raise BenchmarkError("An interaction result must contain at least two distinct features")
    return "|".join(features)


def _parse_results(record: dict[str, Any]) -> pd.DataFrame:
    spec = record["result_spec"]
    placeholders = record["placeholders"]
    result_path = Path(_format_value(str(spec.get("path", "")), placeholders))
    if not result_path.is_absolute():
        result_path = Path(record["cwd"]) / result_path
    if not result_path.is_file():
        raise BenchmarkError(f"Expected result file was not produced: {result_path}")
    result_format = spec.get("format", "generic")
    if result_format == "plink_epistasis":
        frame = pd.read_csv(result_path, sep=r"\s+", comment="#")
        feature_columns = spec.get("feature_columns", ["SNP1", "SNP2"])
        p_column = spec.get("p_column", "P")
    elif result_format == "fighi":
        frame = pd.read_csv(result_path)
        feature_columns = []
        p_column = spec.get("p_column", "p_value")
    else:
        separator = spec.get("separator", ",")
        separator = "\t" if separator in {"tab", "\\t"} else separator
        frame = pd.read_csv(result_path, sep=separator, comment=spec.get("comment"))
        feature_columns = spec.get("feature_columns", [])
        p_column = spec.get("p_column", "p_value")
    if frame.empty:
        return pd.DataFrame(
            columns=[
                "method",
                "hyperedge",
                "order",
                "p_value",
                "rank_score",
                "rank_value",
                "effect",
            ]
        )

    if result_format == "fighi":
        edge_column = _column(frame, spec.get("hyperedge_column", "hyperedge"))
        hyperedges = (
            frame[edge_column]
            .astype(str)
            .map(lambda value: _normalize_hyperedge(re.split(r"\s*[|×]\s*", value)))
        )
    else:
        if not feature_columns:
            raise BenchmarkError("Generic parser requires result.feature_columns")
        resolved_features = [_column(frame, name) for name in feature_columns]
        hyperedges = frame[resolved_features].apply(
            lambda row: _normalize_hyperedge(row.tolist()), axis=1
        )
    if p_column:
        resolved_p = _column(frame, p_column)
        p_values = pd.to_numeric(frame[resolved_p], errors="coerce")
    else:
        p_values = pd.Series(np.nan, index=frame.index)
    score_column = spec.get("score_column")
    if score_column:
        rank_scores = pd.to_numeric(frame[_column(frame, score_column)], errors="coerce")
    else:
        rank_scores = pd.Series(np.nan, index=frame.index)
    score_ascending = bool(spec.get("score_ascending", False))
    effect_column = spec.get("effect_column")
    if effect_column:
        effects = pd.to_numeric(frame[_column(frame, effect_column)], errors="coerce")
    else:
        effects = pd.Series(np.nan, index=frame.index)
    parsed = pd.DataFrame(
        {
            "method": record["method"],
            "hyperedge": hyperedges,
            "order": hyperedges.str.count(r"\|") + 1,
            "p_value": p_values,
            "rank_score": rank_scores,
            "effect": effects,
            "source_file": str(result_path),
        }
    )
    valid_p = np.isfinite(parsed["p_value"]) & parsed["p_value"].between(0, 1)
    valid_score = np.isfinite(parsed["rank_score"])
    parsed = parsed.loc[valid_p | valid_score].copy()
    parsed["rank_value"] = np.where(
        valid_p.loc[parsed.index],
        parsed["p_value"],
        parsed["rank_score"] if score_ascending else -parsed["rank_score"],
    )
    return (
        parsed.sort_values("rank_value")
        .drop_duplicates("hyperedge", keep="first")
        .reset_index(drop=True)
    )


def _truth(payload: dict[str, Any], manifest_path: Path) -> set[str]:
    truth_value = payload.get("shared_inputs", {}).get("truth_file")
    if not truth_value:
        return set()
    truth_path = Path(truth_value)
    if not truth_path.is_absolute():
        truth_path = (manifest_path.parent / truth_path).resolve()
    truth_payload = json.loads(truth_path.read_text(encoding="utf-8"))
    values = truth_payload.get("truth_hyperedges")
    if values is None:
        values = truth_payload.get("truth_interactions", [])
    truth: set[str] = set()
    for value in values:
        if isinstance(value, list):
            truth.add(_normalize_hyperedge(value))
        else:
            truth.add(_normalize_hyperedge(re.split(r"\s*[|×]\s*", str(value))))
    return truth


def _average_precision(relevance: np.ndarray) -> float:
    positives = int(relevance.sum())
    if positives == 0:
        return float("nan")
    precision = np.cumsum(relevance) / np.arange(1, len(relevance) + 1)
    return float(np.sum(precision * relevance) / positives)


def _method_summaries(
    interactions: pd.DataFrame, truth: set[str], alpha: float, top_n: int
) -> pd.DataFrame:
    columns = [
        "method",
        "tests",
        "significant",
        "minimum_p_value",
        "truth_size",
        "true_positives",
        "precision",
        "recall",
        "top_n_recall",
        "average_precision",
    ]
    rows = []
    for method, frame in interactions.groupby("method", sort=True):
        ranked = frame.sort_values("rank_value")
        significant = ranked["q_harmonized"] <= alpha
        discovered = set(ranked.loc[significant, "hyperedge"])
        relevant = ranked["hyperedge"].isin(truth).to_numpy(dtype=int)
        true_positive = len(discovered & truth) if truth else np.nan
        rows.append(
            {
                "method": method,
                "tests": len(ranked),
                "significant": int(significant.sum()),
                "minimum_p_value": float(ranked["p_value"].min())
                if ranked["p_value"].notna().any()
                else np.nan,
                "truth_size": len(truth) if truth else np.nan,
                "true_positives": true_positive,
                "precision": float(true_positive / len(discovered))
                if truth and discovered
                else (0.0 if truth else np.nan),
                "recall": float(true_positive / len(truth)) if truth else np.nan,
                "top_n_recall": float(
                    len(set(ranked.head(top_n)["hyperedge"]) & truth) / len(truth)
                )
                if truth
                else np.nan,
                "average_precision": _average_precision(relevant) if truth else np.nan,
            }
        )
    return pd.DataFrame(rows, columns=columns)


def _overlap(interactions: pd.DataFrame, top_n: int) -> pd.DataFrame:
    ranked = {
        method: set(frame.nsmallest(top_n, "rank_value")["hyperedge"])
        for method, frame in interactions.groupby("method")
    }
    rows = []
    for left, right in combinations(sorted(ranked), 2):
        union = ranked[left] | ranked[right]
        intersection = ranked[left] & ranked[right]
        rows.append(
            {
                "method_a": left,
                "method_b": right,
                "top_n": top_n,
                "intersection": len(intersection),
                "union": len(union),
                "jaccard": len(intersection) / len(union) if union else 1.0,
            }
        )
    return pd.DataFrame(
        rows,
        columns=["method_a", "method_b", "top_n", "intersection", "union", "jaccard"],
    )


def _plots(outdir: Path, runs: pd.DataFrame, summary: pd.DataFrame) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_dir = outdir / "plots"
    plot_dir.mkdir(exist_ok=True)
    paths: list[str] = []
    completed = runs.loc[runs["status"] == "completed"].copy()
    if not completed.empty:
        figure, axes = plt.subplots(1, 2, figsize=(10, 4.2))
        wall = pd.to_numeric(completed["wall_seconds"], errors="coerce").fillna(0.0)
        memory = pd.to_numeric(completed["peak_rss_kb"], errors="coerce").fillna(0.0)
        axes[0].bar(completed["method"], wall, color="#0f766e")
        axes[0].set_ylabel("Wall time (s)")
        axes[0].tick_params(axis="x", rotation=30)
        axes[1].bar(completed["method"], memory / 1024.0, color="#334155")
        axes[1].set_ylabel("Peak RSS (MiB)")
        axes[1].tick_params(axis="x", rotation=30)
        figure.tight_layout()
        path = plot_dir / "resources.png"
        figure.savefig(path, dpi=160)
        plt.close(figure)
        paths.append(str(path))
    if not summary.empty:
        figure, axis = plt.subplots(figsize=(7, 4.2))
        axis.bar(summary["method"], summary["significant"], color="#1d4ed8")
        axis.set_ylabel("Harmonized significant interactions")
        axis.tick_params(axis="x", rotation=30)
        figure.tight_layout()
        path = plot_dir / "discoveries.png"
        figure.savefig(path, dpi=160)
        plt.close(figure)
        paths.append(str(path))
    return paths


def _report(
    outdir: Path,
    payload: dict[str, Any],
    runs: pd.DataFrame,
    summaries: pd.DataFrame,
    overlap: pd.DataFrame,
    plots: list[str],
) -> Path:
    def table(frame: pd.DataFrame) -> str:
        if frame.empty:
            return "<p>No rows.</p>"
        return frame.fillna("").to_html(
            index=False, border=0, escape=True, float_format=lambda value: f"{value:.5g}"
        )

    images = "".join(
        f'<img src="plots/{html.escape(Path(path).name)}" alt="Benchmark plot">' for path in plots
    )
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>FIGHI benchmark</title>
<style>body{{max-width:1180px;margin:auto;padding:36px 22px;font:15px/1.5 system-ui;color:#0f172a}}
h1,h2{{color:#115e59}}table{{border-collapse:collapse;width:100%;display:block;overflow:auto}}
th,td{{padding:8px 10px;border:1px solid #dbe3ec;text-align:left}}th{{background:#f1f5f9}}
img{{max-width:520px;width:100%;margin:10px}}.notice{{padding:16px;background:#fff7ed;border-left:4px solid #f97316}}</style>
</head><body><h1>FIGHI benchmark: {html.escape(str(payload["analysis_id"]))}</h1>
<p class="notice">Methods may test different statistical targets. Compare inputs, adjustment,
phenotype support and model assumptions alongside speed and overlap; agreement is not proof of truth.</p>
<h2>Execution</h2>{table(runs.drop(columns=["command"], errors="ignore"))}
<h2>Harmonized results</h2>{table(summaries)}<h2>Top-result overlap</h2>{table(overlap)}
<h2>Diagnostics</h2>{images or "<p>No plots available.</p>"}</body></html>"""
    path = outdir / "benchmark_report.html"
    path.write_text(document, encoding="utf-8")
    return path


def run_benchmark(
    manifest: str | Path,
    outdir: str | Path,
    *,
    strict: bool = False,
    plots: bool = True,
) -> dict[str, str]:
    """Run enabled methods, harmonize their outputs, and write comparison artifacts."""
    validation = validate_benchmark_manifest(manifest)
    payload, manifest_path = _load_manifest(manifest)
    target = Path(outdir).expanduser().resolve()
    target.mkdir(parents=True, exist_ok=True)
    write_json(target / "validated_inputs.json", validation)

    records: list[dict[str, Any]] = []
    interactions: list[pd.DataFrame] = []
    for method in payload["methods"]:
        if not method.get("enabled", True):
            records.append(
                {
                    "method": method["name"],
                    "status": "disabled",
                    "target": method.get("target", "not specified"),
                    "notes": method.get("notes", ""),
                }
            )
            continue
        record = _run_method(method, payload, manifest_path, target)
        if record["status"] == "completed":
            try:
                interactions.append(_parse_results(record))
            except BenchmarkError as exc:
                record["status"] = "parse_failed"
                record["error"] = str(exc)
        records.append(record)

    run_frame = pd.DataFrame(records)
    run_public = run_frame.drop(columns=["result_spec", "placeholders"], errors="ignore")
    if "command" in run_public:
        run_public["command"] = run_public["command"].map(
            lambda value: json.dumps(value) if isinstance(value, list) else value
        )
    run_path = target / "benchmark_runs.csv"
    run_public.to_csv(run_path, index=False)

    if interactions:
        combined = pd.concat(interactions, ignore_index=True)
        correction = payload.get("correction", "fdr_bh")
        combined["q_harmonized"] = np.nan
        for indices in combined.groupby("method").groups.values():
            valid_indices = [
                index for index in indices if np.isfinite(combined.at[index, "p_value"])
            ]
            if valid_indices:
                combined.loc[valid_indices, "q_harmonized"] = adjust_pvalues(
                    combined.loc[valid_indices, "p_value"].to_numpy(), correction
                )
        combined["significant_harmonized"] = combined["q_harmonized"] <= float(
            payload.get("alpha", 0.05)
        )
        combined = combined.sort_values(["method", "rank_value"]).reset_index(drop=True)
    else:
        combined = pd.DataFrame(
            columns=[
                "method",
                "hyperedge",
                "order",
                "p_value",
                "rank_score",
                "rank_value",
                "effect",
                "source_file",
                "q_harmonized",
                "significant_harmonized",
            ]
        )
    interaction_path = target / "benchmark_interactions.csv"
    combined.to_csv(interaction_path, index=False)

    truth = _truth(payload, manifest_path)
    top_n = int(payload.get("top_n", 100))
    summaries = _method_summaries(combined, truth, float(payload.get("alpha", 0.05)), top_n)
    summary_path = target / "benchmark_method_summary.csv"
    summaries.to_csv(summary_path, index=False)
    overlaps = _overlap(combined, top_n)
    overlap_path = target / "benchmark_overlap.csv"
    overlaps.to_csv(overlap_path, index=False)
    plot_paths = _plots(target, run_public, summaries) if plots else []
    report_path = _report(target, payload, run_public, summaries, overlaps, plot_paths)
    summary_json = write_json(
        target / "benchmark_summary.json",
        {
            "schema_version": "1.0",
            "analysis_id": payload["analysis_id"],
            "correction": payload.get("correction", "fdr_bh"),
            "alpha": payload.get("alpha", 0.05),
            "top_n": top_n,
            "truth_hyperedges": sorted(truth),
            "method_status": {str(row["method"]): str(row["status"]) for row in records},
            "caution": (
                "Comparator methods can differ in phenotype support, null model, covariate "
                "adjustment, search space and inferential target. Harmonized p-value adjustment "
                "does not make unlike tests equivalent."
            ),
        },
    )
    failures = run_public.loc[~run_public["status"].isin(["completed", "disabled"])]
    if strict and not failures.empty:
        failed = ", ".join(f"{row.method} ({row.status})" for row in failures.itertuples())
        raise BenchmarkError(f"Benchmark completed with unavailable or failed methods: {failed}")
    return {
        "runs": str(run_path),
        "interactions": str(interaction_path),
        "method_summary": str(summary_path),
        "overlap": str(overlap_path),
        "summary": str(summary_json),
        "report": str(report_path),
    }


def write_benchmark_template(
    path: str | Path,
    *,
    phenotype: str = "case",
    trait: str = "binary",
    covariates: str = "age,PC1",
    max_order: int = 2,
) -> Path:
    """Write a portable template with FIGHI enabled and comparator adapters documented."""
    payload = {
        "schema_version": "1.0",
        "analysis_id": "simulation_pairwise_seed17",
        "alpha": 0.05,
        "correction": "fdr_bh",
        "top_n": 100,
        "shared_inputs": {
            "data_file": "simulation.tsv.gz",
            "candidate_file": "candidates.txt",
            "sample_file": "samples.txt",
            "truth_file": "truth.json",
        },
        "methods": [
            {
                "name": "fighi",
                "enabled": True,
                "target": "highest-order product conditional on covariates and all within-candidate lower terms",
                "notes": (
                    "Template uses a fixed all-pairs universe and discovery-fraction 0 for a "
                    "same-sample search benchmark. Do not call this confirmatory after adaptive selection."
                ),
                "command": [
                    "{python_executable}",
                    "-m",
                    "fighi",
                    "run",
                    "{data_file}",
                    "--phenotype",
                    phenotype,
                    "--id-column",
                    "IID",
                    "--covariates",
                    covariates,
                    "--trait",
                    trait,
                    "--feature-file",
                    "{candidate_file}",
                    "--screen-method",
                    "all",
                    "--discovery-fraction",
                    "0",
                    "--max-order",
                    str(max_order),
                    "--outdir",
                    "{run_dir}/analysis",
                    "--no-plots",
                ],
                "version_command": ["{python_executable}", "-m", "fighi", "--version"],
                "result": {
                    "path": "{run_dir}/analysis/fighi_results.csv",
                    "format": "fighi",
                },
            },
            {
                "name": "plink_fast_epistasis",
                "enabled": False,
                "target": "PLINK fast pairwise epistasis statistic",
                "notes": "Enable after replacing PLINK_PREFIX and confirming phenotype/sample coding.",
                "command": [
                    "plink",
                    "--bfile",
                    "PLINK_PREFIX",
                    "--extract",
                    "{candidate_file}",
                    "--keep",
                    "{sample_file}",
                    "--fast-epistasis",
                    "--out",
                    "{run_dir}/plink",
                ],
                "version_command": ["plink", "--version"],
                "result": {
                    "path": "{run_dir}/plink.epi.cc",
                    "format": "plink_epistasis",
                },
            },
            {
                "name": "mdr_or_gmdr",
                "enabled": False,
                "target": "Cross-validated non-parametric interaction ranking",
                "notes": (
                    "Replace the complete command for the installed MDR/GMDR version; record "
                    "cross-validation, permutation, covariate handling and phenotype restrictions."
                ),
                "command": ["MDR_OR_GMDR_EXECUTABLE", "REPLACE_WITH_VERIFIED_ARGUMENTS"],
                "result": {
                    "path": "{run_dir}/mdr_results.tsv",
                    "format": "generic",
                    "separator": "tab",
                    "feature_columns": ["feature1", "feature2"],
                    "p_column": None,
                    "score_column": "cross_validation_score",
                    "score_ascending": False,
                },
            },
            {
                "name": "boost",
                "enabled": False,
                "target": "Fast pairwise case-control interaction screening",
                "notes": "Replace the command and map the installed BOOST implementation's output.",
                "command": ["BOOST_EXECUTABLE", "REPLACE_WITH_VERIFIED_ARGUMENTS"],
                "result": {
                    "path": "{run_dir}/boost_results.tsv",
                    "format": "generic",
                    "separator": "tab",
                    "feature_columns": ["SNP1", "SNP2"],
                    "p_column": "P",
                },
            },
            {
                "name": "beam",
                "enabled": False,
                "target": "Bayesian epistasis association mapping",
                "notes": "Declare priors, search/MCMC settings and whether the mapped value is a p-value.",
                "command": ["BEAM_EXECUTABLE", "REPLACE_WITH_VERIFIED_ARGUMENTS"],
                "result": {
                    "path": "{run_dir}/beam_results.tsv",
                    "format": "generic",
                    "separator": "tab",
                    "feature_columns": ["feature1", "feature2"],
                    "p_column": None,
                    "score_column": "posterior_score",
                    "score_ascending": False,
                },
            },
            {
                "name": "wtccc_exhaustive_or_pathway_aware",
                "enabled": False,
                "target": "Version-specific exhaustive or pathway-informed interaction search",
                "notes": "Record tool/version, pathway release, candidate construction and exact statistic.",
                "command": ["COMPARATOR_EXECUTABLE", "REPLACE_WITH_VERIFIED_ARGUMENTS"],
                "result": {
                    "path": "{run_dir}/comparator_results.tsv",
                    "format": "generic",
                    "separator": "tab",
                    "feature_columns": ["feature1", "feature2"],
                    "p_column": "p_value",
                },
            },
            {
                "name": "information_theoretic",
                "enabled": False,
                "target": "Entropy, mutual-information or information-gain interaction ranking",
                "notes": "Record discretization, estimator, bias correction and significance procedure.",
                "command": ["INFORMATION_TOOL", "REPLACE_WITH_VERIFIED_ARGUMENTS"],
                "result": {
                    "path": "{run_dir}/information_results.tsv",
                    "format": "generic",
                    "separator": "tab",
                    "feature_columns": ["feature1", "feature2"],
                    "p_column": None,
                    "score_column": "information_score",
                    "score_ascending": False,
                },
            },
        ],
    }
    return write_json(path, payload)
