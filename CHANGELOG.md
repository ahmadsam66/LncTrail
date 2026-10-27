# Changelog

## v0.2.3

- Create `lnctrail_results.tsv` and `lnctrail_candidate_sites.tsv` before the prediction stage.
- Process target transcripts in configurable batches (`--batch-size`, default 250).
- Append a complete query–target audit row for each pair in a completed batch and flush the data to disk.
- Append candidate-site details per batch and update `run.log` plus `run_metadata.json` with progress.
- Finalize the progressive TSVs with global prioritization fields after all batches finish.
- Preserve the full result inventory, ranked gene/site shortlists, HTML report, optional IntaRNA backend, optional RIME integration, and scientific limitations.
- Add a mocked IntaRNA batch test confirming that completed rows are visible before the following batch starts.

## v0.2.2

- GitHub-focused project packaging, README, infographic, sample inputs, tests and CI workflow.
- GitHub-oriented README, infographic, example inputs, tests, and CI workflow.
