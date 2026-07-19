class FIGHIError(Exception):
    """Base class for expected, user-facing FIGHI errors."""


class InputValidationError(FIGHIError):
    """Raised when input data cannot support the requested analysis."""


class CandidateLimitError(FIGHIError):
    """Raised before a candidate search would exceed the configured safety limit."""


class ExternalToolError(FIGHIError):
    """Raised when an optional external program is missing or fails."""


class BenchmarkError(FIGHIError):
    """Raised when a benchmark manifest or comparator run is invalid."""
