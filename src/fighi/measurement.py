from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

try:
    import resource
except ImportError:  # pragma: no cover - unavailable on Windows
    resource = None  # type: ignore[assignment]


def _usage() -> tuple[float | None, float | None, float | None]:
    if resource is None:
        return None, None, None
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    peak_rss_kb = float(usage.ru_maxrss)
    if sys.platform == "darwin":  # macOS reports bytes; Linux reports KiB.
        peak_rss_kb /= 1024.0
    return float(usage.ru_utime), float(usage.ru_stime), peak_rss_kb


def main(argv: list[str] | None = None) -> int:
    """Execute one comparator in a fresh process and record isolated resource use."""
    arguments = list(sys.argv[1:] if argv is None else argv)
    if len(arguments) < 2:
        print("usage: python -m fighi.measurement METRICS_JSON COMMAND [ARG ...]", file=sys.stderr)
        return 2
    metrics_path = Path(arguments[0])
    command = arguments[1:]
    start = time.perf_counter()
    return_code = 127
    error = None
    try:
        completed = subprocess.run(command, check=False)
        return_code = completed.returncode
    except OSError as exc:
        error = str(exc)
    user_seconds, system_seconds, peak_rss_kb = _usage()
    payload = {
        "schema_version": "1.0",
        "source": "fighi_isolated_rusage" if resource else "fighi_wall_clock",
        "wall_seconds": time.perf_counter() - start,
        "user_seconds": user_seconds,
        "system_seconds": system_seconds,
        "peak_rss_kb": peak_rss_kb,
        "return_code": return_code,
        "error": error,
    }
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    metrics_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return return_code


if __name__ == "__main__":  # pragma: no cover - exercised through benchmark integration
    raise SystemExit(main())
