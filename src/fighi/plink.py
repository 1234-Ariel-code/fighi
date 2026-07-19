from __future__ import annotations

import shlex
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from .errors import ExternalToolError, InputValidationError
from .utilities import sha256_file, write_json

_PLINK_METADATA = {"FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE", "#FID", "#IID"}


def _candidate_ids(path: str | Path) -> list[str]:
    candidate_path = Path(path)
    if not candidate_path.is_file():
        raise InputValidationError(f"Candidate file does not exist: {candidate_path}")
    identifiers = []
    for line in candidate_path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            identifiers.append(stripped.split()[0])
    if not identifiers:
        raise InputValidationError("Candidate file contains no variant identifiers")
    if len(set(identifiers)) != len(identifiers):
        raise InputValidationError("Candidate file contains duplicate variant identifiers")
    return identifiers


def _detect_input(prefix: Path, requested: str) -> tuple[str, list[str]]:
    if requested == "auto":
        if prefix.suffix.lower() in {".vcf", ".bcf"} or prefix.name.lower().endswith(".vcf.gz"):
            requested = "vcf"
        elif prefix.with_suffix(".bed").is_file():
            requested = "bed"
        elif prefix.with_suffix(".pgen").is_file():
            requested = "pgen"
        else:
            raise InputValidationError(
                "Could not infer PLINK input. Provide a BED/PGEN prefix or VCF and set --input-type."
            )

    if requested == "bed":
        base = (
            prefix.with_suffix("") if prefix.suffix.lower() in {".bed", ".bim", ".fam"} else prefix
        )
        expected = [base.with_suffix(ext) for ext in (".bed", ".bim", ".fam")]
        missing = [str(path) for path in expected if not path.is_file()]
        if missing:
            raise InputValidationError(f"Missing BED input files: {', '.join(missing)}")
        return requested, ["--bfile", str(base)]
    if requested == "pgen":
        base = (
            prefix.with_suffix("")
            if prefix.suffix.lower() in {".pgen", ".pvar", ".psam"}
            else prefix
        )
        pvar = base.with_suffix(".pvar")
        pvar_zst = Path(str(pvar) + ".zst")
        expected = [base.with_suffix(".pgen"), base.with_suffix(".psam")]
        missing = [str(path) for path in expected if not path.is_file()]
        if not pvar.is_file() and not pvar_zst.is_file():
            missing.append(f"{pvar} (or {pvar_zst})")
        if missing:
            raise InputValidationError(f"Missing PGEN input files: {', '.join(missing)}")
        return requested, ["--pfile", str(base)]
    if requested == "vcf":
        if not prefix.is_file():
            raise InputValidationError(f"VCF/BCF input does not exist: {prefix}")
        flag = "--bcf" if prefix.suffix.lower() == ".bcf" else "--vcf"
        return requested, [flag, str(prefix)]
    raise InputValidationError("input-type must be auto, bed, pgen, or vcf")


def build_plink2_export_command(
    input_path: str | Path,
    candidate_file: str | Path,
    output_prefix: str | Path,
    *,
    input_type: str = "auto",
    sample_file: str | Path | None = None,
    plink2: str = "plink2",
    threads: int = 1,
    memory_mb: int | None = None,
) -> tuple[list[str], str]:
    """Build a PLINK 2 additive export command without executing a shell."""
    prefix = Path(input_path).expanduser().resolve()
    candidates = Path(candidate_file).expanduser().resolve()
    _candidate_ids(candidates)
    detected, input_args = _detect_input(prefix, input_type)
    command = [plink2, *input_args, "--extract", str(candidates)]
    if sample_file:
        samples = Path(sample_file).expanduser().resolve()
        if not samples.is_file():
            raise InputValidationError(f"Sample file does not exist: {samples}")
        command.extend(["--keep", str(samples)])
    command.extend(["--threads", str(max(1, threads))])
    if memory_mb is not None:
        if memory_mb < 64:
            raise InputValidationError("PLINK memory allocation must be at least 64 MB")
        command.extend(["--memory", str(memory_mb)])
    command.extend(["--export", "A", "include-alt", "--out", str(Path(output_prefix).resolve())])
    return command, detected


def _read_auxiliary(path: str | Path, iid_column: str) -> pd.DataFrame:
    table_path = Path(path)
    if not table_path.is_file():
        raise InputValidationError(f"Auxiliary table does not exist: {table_path}")
    try:
        frame = pd.read_csv(table_path, sep=None, engine="python", dtype={iid_column: str})
    except Exception as exc:
        raise InputValidationError(f"Could not read auxiliary table {table_path}: {exc}") from exc
    if iid_column not in frame:
        raise InputValidationError(f"Auxiliary table lacks IID column '{iid_column}'")
    frame[iid_column] = frame[iid_column].astype(str)
    if frame[iid_column].duplicated().any():
        raise InputValidationError(f"Auxiliary table contains duplicate {iid_column} values")
    return frame


def _rename_dosage_columns(columns: list[str], candidate_ids: list[str]) -> dict[str, str]:
    candidate_set = set(candidate_ids)
    longest_first = sorted(candidate_ids, key=len, reverse=True)
    mapping: dict[str, str] = {}
    used: set[str] = set()
    for column in columns:
        if column in candidate_set:
            variant = column
        else:
            variant = next(
                (identifier for identifier in longest_first if column.startswith(f"{identifier}_")),
                "",
            )
        if variant:
            if variant in used:
                raise InputValidationError(
                    f"PLINK export contains multiple dosage columns for candidate '{variant}'"
                )
            mapping[column] = variant
            used.add(variant)
    return mapping


def prepare_plink_dataset(
    input_path: str | Path,
    candidate_file: str | Path,
    outdir: str | Path,
    *,
    input_type: str = "auto",
    sample_file: str | Path | None = None,
    phenotype_file: str | Path | None = None,
    phenotype_column: str = "phenotype",
    iid_column: str = "IID",
    covariate_columns: list[str] | None = None,
    plink2: str = "plink2",
    threads: int = 1,
    memory_mb: int | None = None,
    dry_run: bool = False,
) -> dict[str, str]:
    """Export selected PLINK variants and create a sample-major FIGHI input table."""
    target = Path(outdir)
    target.mkdir(parents=True, exist_ok=True)
    output_prefix = target / "plink_export"
    command, detected = build_plink2_export_command(
        input_path,
        candidate_file,
        output_prefix,
        input_type=input_type,
        sample_file=sample_file,
        plink2=plink2,
        threads=threads,
        memory_mb=memory_mb,
    )
    command_path = target / "plink2_command.txt"
    command_path.write_text(shlex.join(command) + "\n", encoding="utf-8")
    candidates = _candidate_ids(candidate_file)
    initial_manifest = {
        "schema_version": "1.0",
        "status": "dry_run" if dry_run else "running",
        "input_type": detected,
        "input": str(Path(input_path).expanduser().resolve()),
        "candidate_file": str(Path(candidate_file).expanduser().resolve()),
        "candidate_file_sha256": sha256_file(candidate_file),
        "sample_file": str(Path(sample_file).expanduser().resolve()) if sample_file else None,
        "sample_file_sha256": sha256_file(sample_file) if sample_file else None,
        "command": command,
        "dosage_convention": "PLINK 2 --export A include-alt (alternate-allele counts)",
    }
    manifest_path = write_json(target / "prepare_plink_manifest.json", initial_manifest)
    if dry_run:
        return {"command": str(command_path), "manifest": str(manifest_path)}

    executable = shutil.which(command[0]) if not Path(command[0]).is_file() else command[0]
    if executable is None:
        raise ExternalToolError(
            f"PLINK 2 executable not found: {command[0]}. Install PLINK 2 or pass --plink2."
        )
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    (target / "plink2.stdout.log").write_text(completed.stdout, encoding="utf-8")
    (target / "plink2.stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        initial_manifest.update({"status": "failed", "return_code": completed.returncode})
        write_json(manifest_path, initial_manifest)
        raise ExternalToolError(
            f"PLINK 2 failed with exit code {completed.returncode}; see {target / 'plink2.stderr.log'}"
        )

    raw_path = output_prefix.with_suffix(".raw")
    if not raw_path.is_file():
        raise ExternalToolError(f"PLINK 2 completed but did not create expected export: {raw_path}")
    raw = pd.read_csv(raw_path, sep=r"\s+", dtype={"IID": str, "#IID": str})
    exported_iid = "IID" if "IID" in raw else "#IID" if "#IID" in raw else None
    if exported_iid is None:
        raise ExternalToolError("PLINK export does not contain an IID column")
    raw[exported_iid] = raw[exported_iid].astype(str)

    genotype_columns = [str(column) for column in raw.columns if column not in _PLINK_METADATA]
    rename = _rename_dosage_columns(genotype_columns, candidates)
    if not rename:
        raise InputValidationError("No requested candidates were present in the PLINK export")
    frame = raw[[exported_iid, *rename.keys()]].rename(columns={exported_iid: iid_column, **rename})

    covariates = covariate_columns or []
    if phenotype_file:
        auxiliary = _read_auxiliary(phenotype_file, iid_column)
        requested = [iid_column, phenotype_column, *covariates]
        missing = [name for name in requested if name not in auxiliary]
        if missing:
            raise InputValidationError(
                f"Auxiliary phenotype/covariate table lacks columns: {', '.join(missing)}"
            )
        frame = frame.merge(auxiliary[requested], on=iid_column, how="left", validate="one_to_one")
    elif "PHENOTYPE" in raw:
        phenotype = pd.to_numeric(raw["PHENOTYPE"], errors="coerce").replace(-9, pd.NA)
        frame.insert(1, phenotype_column, phenotype)
        if covariates:
            raise InputValidationError("Covariates require --phenotype-file")
    else:
        raise InputValidationError(
            "No phenotype was available. Supply --phenotype-file and --phenotype-column."
        )

    if frame[phenotype_column].isna().any():
        frame = frame.loc[frame[phenotype_column].notna()].reset_index(drop=True)
    if frame.empty:
        raise InputValidationError("No samples remain after matching non-missing phenotypes")
    ordered_candidates = [identifier for identifier in candidates if identifier in frame.columns]
    ordered = [iid_column, phenotype_column, *covariates, *ordered_candidates]
    frame = frame[ordered]

    data_path = target / "fighi_input.tsv.gz"
    frame.to_csv(data_path, sep="\t", index=False, compression="gzip")
    selected_path = target / "fighi_candidates.txt"
    selected_path.write_text("\n".join(ordered_candidates) + "\n", encoding="utf-8")
    samples_path = target / "fighi_samples.txt"
    samples_path.write_text(
        "\n".join(f"0\t{value}" for value in frame[iid_column].astype(str)) + "\n",
        encoding="utf-8",
    )
    initial_manifest.update(
        {
            "status": "completed",
            "return_code": 0,
            "samples": len(frame),
            "requested_candidates": len(candidates),
            "exported_candidates": len(ordered_candidates),
            "missing_candidates": [item for item in candidates if item not in ordered_candidates],
            "phenotype_column": phenotype_column,
            "covariate_columns": covariates,
            "output_sha256": sha256_file(data_path),
        }
    )
    write_json(manifest_path, initial_manifest)
    return {
        "data": str(data_path),
        "candidates": str(selected_path),
        "samples": str(samples_path),
        "manifest": str(manifest_path),
        "command": str(command_path),
    }
