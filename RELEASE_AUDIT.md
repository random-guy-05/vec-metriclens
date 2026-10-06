# Release audit

Release: **v1.0.0**  
Audit date: **2026-10-05**

## Release gate

GitHub Actions installs `anndata` and `veckit==0.1.2` on Python **3.10, 3.11, and 3.12** and runs:

- the complete pytest suite;
- a real synthetic T1 scorer/attribution integration test;
- Ruff;
- bytecode compilation.

No Challenge data is bundled.

## Coverage

Tests cover:

- published floor/ceiling skill transforms;
- signed zero-ideal metrics;
- task-weight aggregation;
- gene-level change attribution;
- spatial scale diagnostics;
- rejection of mismatched T2 board/setting combinations;
- end-to-end real-scorer report generation.

MetricLens diagnostics are explanatory and are not presented as additional official metrics.

## Current verification

The latest public GitHub Actions matrix is green. Published weights and validation anchors are documented in `docs/METRIC_CONTRACT.md` and tied to the current source snapshot.
