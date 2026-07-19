# Security policy

Version 1.x receives security fixes. Report vulnerabilities privately through GitHub's security-advisory feature; do not open a public issue before a fix is available.

FIGHI processes sensitive research data locally and does not transmit input data. The annotation command uses local mapping and GMT files. The project does not install dependencies at runtime and does not deserialize pickle files.

Users should:

- install only trusted releases in an isolated environment;
- keep genotype and phenotype files outside public repositories;
- review output permissions on shared systems;
- treat command lines and provenance files as potentially revealing file paths;
- pin and scan dependencies in regulated environments.

