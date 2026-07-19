from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom

from .errors import InputValidationError
from .statistics import bh_adjust


def _find_column(frame: pd.DataFrame, candidates: list[str]) -> str | None:
    lookup = {str(name).lower(): str(name) for name in frame.columns}
    return next((lookup[name.lower()] for name in candidates if name.lower() in lookup), None)


def _load_mapping(path: str | Path) -> tuple[dict[str, set[str]], set[str]]:
    source = Path(path)
    separator = "\t" if source.suffix.lower() in {".tsv", ".tab"} else ","
    frame = pd.read_csv(source, sep=separator, low_memory=False)
    feature_column = _find_column(frame, ["feature", "snp", "rsid", "variant"])
    gene_column = _find_column(frame, ["gene", "symbol", "gene_symbol"])
    if not feature_column or not gene_column:
        raise InputValidationError(
            "SNP-gene mapping needs feature/SNP/rsid and gene/symbol columns"
        )
    mapping: dict[str, set[str]] = defaultdict(set)
    universe: set[str] = set()
    for feature, gene in zip(frame[feature_column], frame[gene_column], strict=False):
        if pd.isna(feature) or pd.isna(gene):
            continue
        for symbol in str(gene).replace("|", ";").split(";"):
            symbol = symbol.strip()
            if symbol:
                mapping[str(feature)].add(symbol)
                universe.add(symbol)
    return dict(mapping), universe


def _load_gmt(paths: list[str | Path]) -> dict[str, set[str]]:
    terms: dict[str, set[str]] = {}
    for path in paths:
        with Path(path).open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                name = fields[0].strip()
                if name:
                    terms[name] = {value.strip() for value in fields[2:] if value.strip()}
    return terms


def annotate_results(
    feature_scores: str | Path,
    snp_gene_map: str | Path,
    outdir: str | Path,
    gmt_files: list[str | Path] | None = None,
    alpha: float = 0.05,
) -> dict[str, str]:
    if not 0 < alpha < 1:
        raise InputValidationError("Annotation alpha must be between 0 and 1")
    source = Path(feature_scores)
    if not source.is_file():
        raise InputValidationError(f"Feature score file does not exist: {source}")
    frame = pd.read_csv(source)
    feature_column = _find_column(frame, ["feature", "snp", "rsid", "variant"])
    if not feature_column:
        raise InputValidationError("Feature score file needs a feature or SNP column")
    mapping, mapping_universe = _load_mapping(snp_gene_map)
    frame["gene"] = [
        ";".join(sorted(mapping.get(str(feature), set()))) for feature in frame[feature_column]
    ]
    significant_column = _find_column(frame, ["significant_interaction_count"])
    if significant_column:
        selected_features = frame.loc[frame[significant_column] > 0, feature_column].astype(str)
    else:
        selected_features = frame[feature_column].astype(str)
    selected_genes = (
        set().union(*(mapping.get(name, set()) for name in selected_features))
        if len(selected_features)
        else set()
    )

    enrichment_columns = [
        "term",
        "overlap",
        "term_size",
        "selected_size",
        "universe_size",
        "p_value",
        "q_value",
        "genes",
    ]
    enrichment = pd.DataFrame(columns=enrichment_columns)
    significant_terms: dict[str, set[str]] = {}
    if gmt_files:
        pathways = _load_gmt(gmt_files)
        pathway_universe = set().union(*pathways.values()) if pathways else set()
        universe = mapping_universe & pathway_universe
        selected = selected_genes & universe
        rows = []
        for term, members in pathways.items():
            members = members & universe
            overlap_genes = selected & members
            if not overlap_genes or not universe or not selected:
                continue
            p_value = float(
                hypergeom.sf(
                    len(overlap_genes) - 1,
                    len(universe),
                    len(members),
                    len(selected),
                )
            )
            rows.append(
                {
                    "term": term,
                    "overlap": len(overlap_genes),
                    "term_size": len(members),
                    "selected_size": len(selected),
                    "universe_size": len(universe),
                    "p_value": p_value,
                    "genes": ";".join(sorted(overlap_genes)),
                }
            )
        if rows:
            enrichment = pd.DataFrame(rows)
            enrichment["q_value"] = bh_adjust(enrichment["p_value"].to_numpy())
            enrichment = enrichment.loc[:, enrichment_columns].sort_values(["q_value", "p_value"])
            significant_terms = {
                row.term: set(str(row.genes).split(";"))
                for row in enrichment.itertuples()
                if row.q_value <= alpha
            }
    feature_pathways = []
    for gene_text in frame["gene"]:
        genes = set(str(gene_text).split(";")) if gene_text else set()
        terms = [term for term, members in significant_terms.items() if genes & members]
        feature_pathways.append(";".join(sorted(terms)))
    frame["pathway"] = feature_pathways

    target = Path(outdir)
    target.mkdir(parents=True, exist_ok=True)
    annotated_path = target / "fighi_feature_scores_annotated.csv"
    enrichment_path = target / "fighi_pathway_enrichment.csv"
    frame.to_csv(annotated_path, index=False)
    enrichment.to_csv(enrichment_path, index=False)
    return {"annotated_features": str(annotated_path), "enrichment": str(enrichment_path)}
