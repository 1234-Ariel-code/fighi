from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig
from .errors import InputValidationError


@dataclass(slots=True)
class PreparedData:
    features: pd.DataFrame
    standardized: pd.DataFrame
    phenotype: np.ndarray
    trait: str
    covariates: np.ndarray | None
    covariate_names: list[str]
    sample_ids: list[str] | None
    qc: dict
    warnings: list[str] = field(default_factory=list)
    phenotype_mapping: dict[str, float] | None = None


def _read_feature_names(path: str | Path) -> list[str]:
    names = [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines()]
    return [name for name in names if name and not name.startswith("#")]


def _separator(path: Path, delimiter: str | None) -> str:
    if delimiter and delimiter != "auto":
        if delimiter == "comma":
            return ","
        if delimiter == "tab":
            return "\t"
        return "\t" if delimiter == "\\t" else delimiter
    return "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","


def load_table(path: str | Path, delimiter: str | None = "auto") -> pd.DataFrame:
    input_path = Path(path)
    if not input_path.is_file():
        raise InputValidationError(f"Input file does not exist: {input_path}")
    try:
        frame = pd.read_csv(input_path, sep=_separator(input_path, delimiter), low_memory=False)
    except Exception as exc:
        raise InputValidationError(f"Could not read {input_path}: {exc}") from exc
    if frame.empty:
        raise InputValidationError("Input table contains no rows")
    if frame.columns.duplicated().any():
        duplicates = frame.columns[frame.columns.duplicated()].tolist()
        raise InputValidationError(f"Duplicate column names are not allowed: {duplicates[:10]}")
    return frame


def _prepare_phenotype(series: pd.Series, requested_trait: str) -> tuple[np.ndarray, str, dict | None]:
    non_missing = series.dropna()
    unique = list(pd.unique(non_missing))
    inferred_binary = len(unique) == 2
    trait = "binary" if requested_trait == "auto" and inferred_binary else requested_trait
    if trait == "auto":
        trait = "linear"

    if trait == "binary":
        if len(unique) != 2:
            raise InputValidationError(
                f"Binary phenotype must have exactly two observed levels; found {len(unique)}"
            )
        unique_strings = {str(value): value for value in unique}
        mapping: dict[object, float]
        if set(unique_strings) == {"0", "1"}:
            mapping = {unique_strings["0"]: 0.0, unique_strings["1"]: 1.0}
        elif set(unique_strings) == {"1", "2"}:
            mapping = {unique_strings["1"]: 0.0, unique_strings["2"]: 1.0}
        else:
            ordered = sorted(unique, key=lambda value: str(value))
            mapping = {ordered[0]: 0.0, ordered[1]: 1.0}
        encoded = series.map(mapping).to_numpy(dtype=float)
        public_mapping = {str(key): value for key, value in mapping.items()}
        return encoded, trait, public_mapping

    numeric = pd.to_numeric(series, errors="coerce")
    newly_missing = numeric.isna() & series.notna()
    if newly_missing.any():
        examples = series[newly_missing].astype(str).unique()[:5].tolist()
        raise InputValidationError(f"Linear phenotype contains non-numeric values: {examples}")
    values = numeric.to_numpy(dtype=float)
    if np.nanstd(values) <= 0:
        raise InputValidationError("Phenotype has zero variance")
    return values, "linear", None


def _prepare_covariates(frame: pd.DataFrame) -> tuple[np.ndarray | None, list[str], list[str]]:
    if frame.shape[1] == 0:
        return None, [], []
    warnings: list[str] = []
    work = frame.copy()
    numeric_columns = work.select_dtypes(include=[np.number, "bool"]).columns.tolist()
    categorical_columns = [name for name in work.columns if name not in numeric_columns]
    for name in numeric_columns:
        numeric = pd.to_numeric(work[name], errors="coerce")
        if numeric.isna().any():
            numeric = numeric.fillna(numeric.median())
            warnings.append(f"Covariate '{name}' had missing values imputed with its median")
        work[name] = numeric
    for name in categorical_columns:
        values = work[name].astype("string")
        if values.isna().any():
            mode = values.mode(dropna=True)
            fill = mode.iloc[0] if not mode.empty else "missing"
            values = values.fillna(fill)
            warnings.append(f"Covariate '{name}' had missing values imputed with its mode")
        work[name] = values
    encoded = pd.get_dummies(work, columns=categorical_columns, drop_first=True, dtype=float)
    varying = encoded.nunique(dropna=False) > 1
    dropped = encoded.columns[~varying].tolist()
    if dropped:
        warnings.append(f"Dropped constant covariates: {', '.join(dropped)}")
    encoded = encoded.loc[:, varying].astype(float)
    if encoded.empty:
        return None, [], warnings
    means = encoded.mean(axis=0)
    scales = encoded.std(axis=0, ddof=0).replace(0.0, 1.0)
    encoded = (encoded - means) / scales
    return encoded.to_numpy(dtype=float), encoded.columns.astype(str).tolist(), warnings


def prepare_dataframe(
    frame: pd.DataFrame,
    phenotype_column: str,
    config: AnalysisConfig,
    id_column: str | None = None,
    covariate_columns: list[str] | None = None,
    feature_names: list[str] | None = None,
    allow_non_genotype: bool = False,
) -> PreparedData:
    config.validate()
    covariate_columns = covariate_columns or []
    required = [phenotype_column] + covariate_columns + ([id_column] if id_column else [])
    missing_required = [name for name in required if name not in frame.columns]
    if missing_required:
        raise InputValidationError(f"Columns not found: {', '.join(missing_required)}")

    work = frame.copy()
    rows_initial = len(work)
    phenotype_missing = work[phenotype_column].isna()
    if phenotype_missing.any():
        work = work.loc[~phenotype_missing].reset_index(drop=True)
    if len(work) < 20:
        raise InputValidationError("At least 20 samples with observed phenotype are required")

    phenotype, trait, phenotype_mapping = _prepare_phenotype(work[phenotype_column], config.trait)
    if trait == "binary":
        class_counts = np.bincount(phenotype.astype(int), minlength=2)
        if class_counts.min() < 5:
            raise InputValidationError("Binary analysis requires at least 5 samples in each class")

    excluded = {phenotype_column, *covariate_columns}
    if id_column:
        excluded.add(id_column)
    if feature_names is None:
        selected = [name for name in work.columns if name not in excluded]
    else:
        missing_features = [name for name in feature_names if name not in work.columns]
        if missing_features:
            raise InputValidationError(f"Requested features not found: {missing_features[:10]}")
        selected = list(dict.fromkeys(feature_names))
    if len(selected) < 2:
        raise InputValidationError("At least two candidate feature columns are required")

    features = work[selected].apply(pd.to_numeric, errors="coerce")
    conversion_failures = []
    for name in selected:
        bad = features[name].isna() & work[name].notna()
        if bad.any():
            conversion_failures.append(name)
    if conversion_failures:
        raise InputValidationError(
            "Feature columns must be numeric genotype dosages. Non-numeric values found in: "
            + ", ".join(conversion_failures[:10])
        )

    warnings: list[str] = []
    missing_rates = features.isna().mean(axis=0)
    high_missing = missing_rates[missing_rates > config.max_missing].index.tolist()
    if high_missing:
        features = features.drop(columns=high_missing)
        warnings.append(f"Dropped {len(high_missing)} features above the missingness threshold")
    all_missing = features.columns[features.isna().all()].tolist()
    if all_missing:
        features = features.drop(columns=all_missing)
    features = features.fillna(features.mean(axis=0))

    if not allow_non_genotype and not features.empty:
        minimum = float(np.nanmin(features.to_numpy()))
        maximum = float(np.nanmax(features.to_numpy()))
        if minimum < -1e-8 or maximum > 2.0 + 1e-8:
            raise InputValidationError(
                "Genotype features must be dosages in [0, 2]. Use --allow-non-genotype "
                "only when intentionally analysing generic numeric features."
            )

    variances = features.var(axis=0, ddof=0)
    constant = variances[variances <= config.min_variance].index.tolist()
    if constant:
        features = features.drop(columns=constant)
        warnings.append(f"Dropped {len(constant)} constant or near-constant features")

    maf_dropped: list[str] = []
    if not allow_non_genotype:
        allele_frequency = features.mean(axis=0) / 2.0
        mafs = np.minimum(allele_frequency, 1.0 - allele_frequency)
        maf_dropped = mafs[mafs < config.min_maf].index.tolist()
        if maf_dropped:
            features = features.drop(columns=maf_dropped)
            warnings.append(f"Dropped {len(maf_dropped)} features below min_maf")
    if features.shape[1] < 2:
        raise InputValidationError("Fewer than two features remained after quality control")

    standardized = (features - features.mean(axis=0)) / features.std(axis=0, ddof=0)
    covariates, covariate_names, covariate_warnings = _prepare_covariates(
        work[covariate_columns]
    )
    warnings.extend(covariate_warnings)
    sample_ids = work[id_column].astype(str).tolist() if id_column else None
    qc = {
        "samples_input": rows_initial,
        "samples_analyzed": len(work),
        "samples_dropped_missing_phenotype": int(phenotype_missing.sum()),
        "features_requested": len(selected),
        "features_analyzed": features.shape[1],
        "features_dropped_missingness": len(high_missing),
        "features_dropped_constant": len(constant) + len(all_missing),
        "features_dropped_maf": len(maf_dropped),
        "phenotype_column": phenotype_column,
        "id_column": id_column,
        "covariates": covariate_names,
    }
    if trait == "binary":
        qc["class_counts"] = {
            "0": int((phenotype == 0).sum()),
            "1": int((phenotype == 1).sum()),
        }
    return PreparedData(
        features=features,
        standardized=standardized,
        phenotype=phenotype,
        trait=trait,
        covariates=covariates,
        covariate_names=covariate_names,
        sample_ids=sample_ids,
        qc=qc,
        warnings=warnings,
        phenotype_mapping=phenotype_mapping,
    )


def prepare_file(
    path: str | Path,
    phenotype_column: str,
    config: AnalysisConfig,
    id_column: str | None = None,
    covariate_columns: list[str] | None = None,
    feature_file: str | Path | None = None,
    delimiter: str | None = "auto",
    allow_non_genotype: bool = False,
) -> PreparedData:
    feature_names = _read_feature_names(feature_file) if feature_file else None
    frame = load_table(path, delimiter)
    if id_column == "auto":
        id_column = next(
            (name for name in ["IID", "iid", "sample_id", "sample", "ID", "id"] if name in frame.columns),
            None,
        )
    return prepare_dataframe(
        frame,
        phenotype_column=phenotype_column,
        config=config,
        id_column=id_column,
        covariate_columns=covariate_columns,
        feature_names=feature_names,
        allow_non_genotype=allow_non_genotype,
    )
