from __future__ import annotations

import csv
import gzip
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from .errors import InputValidationError


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def _sample_ids(tfam: Path) -> list[str]:
    identifiers = []
    with _open_text(tfam) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                identifiers.append(fields[1])
    if not identifiers:
        raise InputValidationError("TFAM contains no sample identifiers")
    if len(set(identifiers)) != len(identifiers):
        raise InputValidationError("TFAM sample identifiers must be unique")
    return identifiers


def _scan_tped(tped: Path, sample_count: int) -> list[str]:
    identifiers = []
    with _open_text(tped) as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 4 + 2 * sample_count:
                raise InputValidationError(
                    f"TPED line {line_number} has {len(fields)} fields; expected {4 + 2 * sample_count}"
                )
            identifiers.append(fields[1])
    if len(set(identifiers)) != len(identifiers):
        raise InputValidationError("TPED variant identifiers must be unique")
    if len(identifiers) < 2:
        raise InputValidationError("TPED must contain at least two variants")
    return identifiers


def _phenotypes(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with _open_text(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 3 or fields[2].upper() in {"PHENO", "PHENOTYPE"}:
                continue
            if fields[2] != "-9":
                values[fields[1]] = fields[2]
    return values


def convert_tped(
    tped_path: str | Path,
    tfam_path: str | Path,
    output_path: str | Path,
    phenotype_path: str | Path | None = None,
    phenotype_name: str = "phenotype",
) -> Path:
    tped, tfam, output = Path(tped_path), Path(tfam_path), Path(output_path)
    if not tped.is_file() or not tfam.is_file():
        raise InputValidationError("Both TPED and TFAM files must exist")
    sample_ids = _sample_ids(tfam)
    variants = _scan_tped(tped, len(sample_ids))
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.NamedTemporaryFile(prefix="fighi-tped-", suffix=".dat", delete=False)
    temporary.close()
    matrix = np.memmap(
        temporary.name, dtype=np.float32, mode="w+", shape=(len(sample_ids), len(variants))
    )
    try:
        with _open_text(tped) as handle:
            for variant_index, line in enumerate(handle):
                fields = line.split()
                alleles = fields[4:]
                observed = [allele for allele in alleles if allele != "0"]
                reference = Counter(observed).most_common(1)[0][0] if observed else None
                for sample_index in range(len(sample_ids)):
                    first, second = alleles[2 * sample_index : 2 * sample_index + 2]
                    if reference is None or "0" in {first, second}:
                        matrix[sample_index, variant_index] = np.nan
                    else:
                        matrix[sample_index, variant_index] = float(
                            (first != reference) + (second != reference)
                        )
        phenotype_values = _phenotypes(Path(phenotype_path)) if phenotype_path else None
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            header = ["IID"]
            if phenotype_values is not None:
                header.append(phenotype_name)
            header.extend(variants)
            writer.writerow(header)
            for sample_index, identifier in enumerate(sample_ids):
                row: list[object] = [identifier]
                if phenotype_values is not None:
                    row.append(phenotype_values.get(identifier, ""))
                row.extend("" if np.isnan(value) else int(value) for value in matrix[sample_index])
                writer.writerow(row)
    finally:
        del matrix
        Path(temporary.name).unlink(missing_ok=True)
    return output

