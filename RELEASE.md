# FIGHI 1.0.0 release summary

Release date: 2026-07-17

This package is the production rewrite of the supplied FIGHI prototype. The detailed before/after review is in `docs/PROJECT_AUDIT.md`.

## Validated in this release

- 11 automated tests pass.
- The wheel builds and installs in a clean virtual environment.
- The installed `fighi` command runs without source-tree access.
- The seeded demo evaluates 190 pairs and recovers the planted `rs_demo_03 × rs_demo_07` interaction after independent holdout inference and global FDR correction.
- JSON/Cytoscape artifacts parse, GraphML is well formed, and the HTML report is created.
- The supplied legacy toy CSV validates and runs through the new CLI.
- Output-directory protection preserves unrelated user files unless overwrite is explicitly requested.

## Publication gates outside the software build

- Run the study-specific calibration and power suite described in `docs/VALIDATION.md`.
- Obtain an independent statistical-method review.
- Confirm the final author list, affiliations, repository namespace, package-index name, and citation/DOI metadata.
- Publish the same tested wheel artifact through the checklist in `docs/RELEASING.md`.

The implementation is release-ready; scientific claims for a manuscript must remain limited to validation actually performed on the intended cohorts and simulations.
