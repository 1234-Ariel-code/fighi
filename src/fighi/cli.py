from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .annotation import annotate_results
from .benchmark import run_benchmark, validate_benchmark_manifest, write_benchmark_template
from .config import AnalysisConfig
from .conversion import convert_tped
from .data import prepare_file
from .demo import create_demo_dataset
from .errors import FIGHIError, InputValidationError
from .plink import prepare_plink_dataset
from .reporting import write_outputs
from .search import FIGHI
from .simulation import SimulationConfig, write_simulation


def _comma_list(value: str | None) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()] if value else []


def _require_output_directory(path: str | Path, overwrite: bool) -> Path:
    target = Path(path)
    if target.exists() and not target.is_dir():
        raise InputValidationError(f"Output path exists and is not a directory: {target}")
    if target.exists() and any(target.iterdir()) and not overwrite:
        raise InputValidationError(
            f"Output directory is not empty: {target}. Choose a new directory or pass --overwrite."
        )
    target.mkdir(parents=True, exist_ok=True)
    return target


def _add_analysis_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--trait", choices=["auto", "binary", "linear"], default="auto")
    parser.add_argument("--max-order", type=int, default=2, help="Highest interaction order (2-6)")
    parser.add_argument(
        "--top-m", type=int, default=50, help="Number of features retained by screening"
    )
    parser.add_argument(
        "--screen-method",
        choices=["hybrid", "marginal", "variance", "all"],
        default="hybrid",
        help="Pre-search feature screening strategy",
    )
    parser.add_argument("--min-maf", type=float, default=0.01)
    parser.add_argument("--max-missing", type=float, default=0.10)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--correction", choices=["fdr_bh", "bonferroni", "none"], default="fdr_bh")
    parser.add_argument(
        "--discovery-fraction",
        type=float,
        default=0.40,
        help="Fraction used only for screening/candidate discovery; 0 enables exploratory all-sample mode",
    )
    parser.add_argument("--beam-width", type=int, default=100)
    parser.add_argument("--max-candidates-per-order", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--stability-repeats", type=int, default=0)
    parser.add_argument("--stability-fraction", type=float, default=0.80)
    parser.add_argument("--stability-threshold", type=float, default=0.70)
    parser.add_argument("--graph-top", type=int, default=100)


def _config(args: argparse.Namespace) -> AnalysisConfig:
    return AnalysisConfig(
        trait=args.trait,
        max_order=args.max_order,
        top_m=args.top_m,
        screen_method=args.screen_method,
        min_maf=args.min_maf,
        max_missing=args.max_missing,
        alpha=args.alpha,
        correction=args.correction,
        discovery_fraction=args.discovery_fraction,
        beam_width=args.beam_width,
        max_candidates_per_order=args.max_candidates_per_order,
        seed=args.seed,
        stability_repeats=args.stability_repeats,
        stability_fraction=args.stability_fraction,
        stability_threshold=args.stability_threshold,
        graph_top=args.graph_top,
    ).validate()


def _prepare(args: argparse.Namespace, config: AnalysisConfig):
    return prepare_file(
        args.input,
        phenotype_column=args.phenotype,
        config=config,
        id_column=args.id_column,
        covariate_columns=_comma_list(args.covariates),
        feature_file=args.feature_file,
        delimiter=args.delimiter,
        allow_non_genotype=args.allow_non_genotype,
    )


def _run(args: argparse.Namespace) -> int:
    outdir = _require_output_directory(args.outdir, args.overwrite)
    config = _config(args)
    data = _prepare(args, config)
    print(
        f"[FIGHI] {data.qc['samples_analyzed']} samples; "
        f"{data.qc['features_analyzed']} QC-passed features; trait={data.trait}",
        flush=True,
    )
    result = FIGHI(config).run(data)
    paths = write_outputs(
        outdir,
        result,
        data,
        input_path=str(Path(args.input).resolve()),
        command=sys.argv,
        plots=not args.no_plots,
    )
    print(
        f"[FIGHI] Evaluated {len(result.interactions):,} interactions; "
        f"{len(result.significant):,} significant at alpha={config.alpha:g}."
    )
    print(f"[FIGHI] Report: {paths['report']}")
    return 0


def _validate(args: argparse.Namespace) -> int:
    config = _config(args)
    data = _prepare(args, config)
    payload = {
        "valid": True,
        "trait": data.trait,
        "phenotype_mapping": data.phenotype_mapping,
        "quality_control": data.qc,
        "warnings": data.warnings,
    }
    print(json.dumps(payload, indent=2))
    return 0


def _demo(args: argparse.Namespace) -> int:
    root = _require_output_directory(args.outdir, args.overwrite)
    data_path = create_demo_dataset(root / "fighi_demo.csv", samples=args.samples, seed=args.seed)
    analysis_path = root / "analysis"
    run_args = argparse.Namespace(
        input=str(data_path),
        phenotype="case",
        id_column="IID",
        covariates="age",
        feature_file=None,
        delimiter="auto",
        allow_non_genotype=False,
        outdir=str(analysis_path),
        overwrite=True,
        no_plots=False,
        trait="binary",
        max_order=2,
        top_m=20,
        screen_method="all",
        min_maf=0.01,
        max_missing=0.10,
        alpha=0.05,
        correction="fdr_bh",
        discovery_fraction=0.40,
        beam_width=100,
        max_candidates_per_order=100_000,
        seed=args.seed,
        stability_repeats=args.stability_repeats,
        stability_fraction=0.80,
        stability_threshold=0.70,
        graph_top=100,
    )
    return _run(run_args)


def _convert(args: argparse.Namespace) -> int:
    output = Path(args.output)
    if output.exists() and not args.overwrite:
        raise InputValidationError(
            f"Output already exists: {output}; pass --overwrite to replace it"
        )
    path = convert_tped(
        args.tped,
        args.tfam,
        output,
        phenotype_path=args.phenotype_file,
        phenotype_name=args.phenotype_name,
    )
    print(f"[FIGHI] Converted data: {path}")
    return 0


def _annotate(args: argparse.Namespace) -> int:
    outdir = _require_output_directory(args.outdir, args.overwrite)
    paths = annotate_results(
        args.feature_scores,
        args.snp_gene_map,
        outdir,
        gmt_files=args.gmt,
        alpha=args.alpha,
    )
    print(json.dumps(paths, indent=2))
    return 0


def _simulate(args: argparse.Namespace) -> int:
    outdir = _require_output_directory(args.outdir, args.overwrite)
    config = SimulationConfig(
        samples=args.samples,
        features=args.features,
        trait=args.trait,
        scenario=args.scenario,
        effect=args.effect,
        prevalence=args.prevalence,
        noise_sd=args.noise_sd,
        min_maf=args.min_maf,
        max_maf=args.max_maf,
        ld_rho=args.ld_rho,
        seed=args.seed,
    )
    paths = write_simulation(outdir, config)
    print(json.dumps(paths, indent=2))
    return 0


def _prepare_plink(args: argparse.Namespace) -> int:
    outdir = _require_output_directory(args.outdir, args.overwrite)
    paths = prepare_plink_dataset(
        args.input,
        args.candidates,
        outdir,
        input_type=args.input_type,
        sample_file=args.samples,
        phenotype_file=args.phenotype_file,
        phenotype_column=args.phenotype_column,
        iid_column=args.iid_column,
        covariate_columns=_comma_list(args.covariates),
        plink2=args.plink2,
        threads=args.threads,
        memory_mb=args.memory_mb,
        dry_run=args.dry_run,
    )
    print(json.dumps(paths, indent=2))
    return 0


def _benchmark_validate(args: argparse.Namespace) -> int:
    print(json.dumps(validate_benchmark_manifest(args.manifest), indent=2))
    return 0


def _benchmark_template(args: argparse.Namespace) -> int:
    output = Path(args.output)
    if output.exists() and not args.overwrite:
        raise InputValidationError(
            f"Output already exists: {output}; pass --overwrite to replace it"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    path = write_benchmark_template(
        output,
        phenotype=args.phenotype,
        trait=args.trait,
        covariates=args.covariates,
        max_order=args.max_order,
    )
    print(f"[FIGHI] Benchmark template: {path}")
    return 0


def _benchmark_run(args: argparse.Namespace) -> int:
    outdir = _require_output_directory(args.outdir, args.overwrite)
    paths = run_benchmark(
        args.manifest,
        outdir,
        strict=args.strict,
        plots=not args.no_plots,
    )
    print(json.dumps(paths, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fighi",
        description="Fisher-information-guided, statistically controlled interaction discovery",
    )
    parser.add_argument("--version", action="version", version=f"FIGHI {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run = subparsers.add_parser("run", help="Run an interaction analysis")
    run.add_argument("input", help="Sample-by-feature CSV or TSV")
    run.add_argument("--phenotype", required=True, help="Phenotype column")
    run.add_argument(
        "--id-column", default="auto", help="ID column; default auto-detects common names"
    )
    run.add_argument("--covariates", help="Comma-separated covariate columns")
    run.add_argument("--feature-file", help="Text file listing candidate feature columns")
    run.add_argument(
        "--delimiter", default="auto", help="auto, comma, tab (\\t), or another delimiter"
    )
    run.add_argument("--allow-non-genotype", action="store_true")
    run.add_argument("--outdir", required=True)
    run.add_argument("--overwrite", action="store_true")
    run.add_argument("--no-plots", action="store_true")
    _add_analysis_arguments(run)
    run.set_defaults(handler=_run)

    validate = subparsers.add_parser(
        "validate", help="Validate and summarize input without searching"
    )
    validate.add_argument("input")
    validate.add_argument("--phenotype", required=True)
    validate.add_argument("--id-column", default="auto")
    validate.add_argument("--covariates")
    validate.add_argument("--feature-file")
    validate.add_argument("--delimiter", default="auto")
    validate.add_argument("--allow-non-genotype", action="store_true")
    _add_analysis_arguments(validate)
    validate.set_defaults(handler=_validate)

    demo = subparsers.add_parser("demo", help="Generate data and run a complete demonstration")
    demo.add_argument("--outdir", default="fighi_demo")
    demo.add_argument("--samples", type=int, default=500)
    demo.add_argument("--seed", type=int, default=17)
    demo.add_argument("--stability-repeats", type=int, default=0)
    demo.add_argument("--overwrite", action="store_true")
    demo.set_defaults(handler=_demo)

    convert = subparsers.add_parser("convert-tped", help="Convert PLINK TPED/TFAM to FIGHI CSV")
    convert.add_argument("--tped", required=True)
    convert.add_argument("--tfam", required=True)
    convert.add_argument("--phenotype-file", help="Optional PLINK FID/IID/phenotype file")
    convert.add_argument("--phenotype-name", default="phenotype")
    convert.add_argument("--output", required=True)
    convert.add_argument("--overwrite", action="store_true")
    convert.set_defaults(handler=_convert)

    annotate = subparsers.add_parser(
        "annotate", help="Map features to genes and run GMT enrichment"
    )
    annotate.add_argument("--feature-scores", required=True)
    annotate.add_argument("--snp-gene-map", required=True)
    annotate.add_argument(
        "--gmt", action="append", default=[], help="GMT file; repeat for multiple files"
    )
    annotate.add_argument("--alpha", type=float, default=0.05)
    annotate.add_argument("--outdir", required=True)
    annotate.add_argument("--overwrite", action="store_true")
    annotate.set_defaults(handler=_annotate)

    simulate = subparsers.add_parser(
        "simulate", help="Generate a reproducible truth-known interaction benchmark"
    )
    simulate.add_argument("--outdir", required=True)
    simulate.add_argument("--samples", type=int, default=1000)
    simulate.add_argument("--features", type=int, default=100)
    simulate.add_argument("--trait", choices=["binary", "linear"], default="binary")
    simulate.add_argument(
        "--scenario",
        choices=["null", "main", "pairwise", "pairwise_main", "threeway", "structure"],
        default="pairwise",
    )
    simulate.add_argument("--effect", type=float, default=1.0)
    simulate.add_argument("--prevalence", type=float, default=0.5)
    simulate.add_argument("--noise-sd", type=float, default=1.0)
    simulate.add_argument("--min-maf", type=float, default=0.10)
    simulate.add_argument("--max-maf", type=float, default=0.45)
    simulate.add_argument("--ld-rho", type=float, default=0.0)
    simulate.add_argument("--seed", type=int, default=17)
    simulate.add_argument("--overwrite", action="store_true")
    simulate.set_defaults(handler=_simulate)

    prepare_plink = subparsers.add_parser(
        "prepare-plink", help="Export selected PLINK variants into a FIGHI-ready table"
    )
    prepare_plink.add_argument("--input", required=True, help="BED/PGEN prefix or VCF/BCF path")
    prepare_plink.add_argument(
        "--input-type", choices=["auto", "bed", "pgen", "vcf"], default="auto"
    )
    prepare_plink.add_argument("--candidates", required=True, help="One variant ID per line")
    prepare_plink.add_argument("--samples", help="Optional PLINK --keep file")
    prepare_plink.add_argument("--phenotype-file", help="IID-keyed phenotype/covariate table")
    prepare_plink.add_argument("--phenotype-column", default="phenotype")
    prepare_plink.add_argument("--iid-column", default="IID")
    prepare_plink.add_argument("--covariates", help="Comma-separated columns in --phenotype-file")
    prepare_plink.add_argument("--plink2", default="plink2")
    prepare_plink.add_argument("--threads", type=int, default=1)
    prepare_plink.add_argument("--memory-mb", type=int)
    prepare_plink.add_argument("--dry-run", action="store_true")
    prepare_plink.add_argument("--outdir", required=True)
    prepare_plink.add_argument("--overwrite", action="store_true")
    prepare_plink.set_defaults(handler=_prepare_plink)

    benchmark_template = subparsers.add_parser(
        "benchmark-template", help="Write an editable JSON comparator manifest"
    )
    benchmark_template.add_argument("--output", required=True)
    benchmark_template.add_argument("--phenotype", default="case")
    benchmark_template.add_argument("--trait", choices=["binary", "linear"], default="binary")
    benchmark_template.add_argument("--covariates", default="age,PC1")
    benchmark_template.add_argument("--max-order", type=int, choices=range(2, 7), default=2)
    benchmark_template.add_argument("--overwrite", action="store_true")
    benchmark_template.set_defaults(handler=_benchmark_template)

    benchmark_validate = subparsers.add_parser(
        "benchmark-validate", help="Validate a comparator manifest without executing it"
    )
    benchmark_validate.add_argument("manifest")
    benchmark_validate.set_defaults(handler=_benchmark_validate)

    benchmark_run = subparsers.add_parser(
        "benchmark-run", help="Execute and compare methods declared in a JSON manifest"
    )
    benchmark_run.add_argument("manifest")
    benchmark_run.add_argument("--outdir", required=True)
    benchmark_run.add_argument(
        "--strict", action="store_true", help="Fail if any enabled method fails"
    )
    benchmark_run.add_argument("--no-plots", action="store_true")
    benchmark_run.add_argument("--overwrite", action="store_true")
    benchmark_run.set_defaults(handler=_benchmark_run)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except FIGHIError as exc:
        print(f"FIGHI error: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("FIGHI interrupted by user", file=sys.stderr)
        return 130
