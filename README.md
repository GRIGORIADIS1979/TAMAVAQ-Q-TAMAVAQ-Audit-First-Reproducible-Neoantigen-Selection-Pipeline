# TAMAVAQ / Q-TAMAVAQ reproducible algorithms

A GitHub-style, reproducible, corrected algorithm package for the two TAMAVAQ manuscripts:

- `Tamavaq™ Journal Of Immunology-iTAMAVAQ-Manuscript6.docx`
- `SI Appendix I Tamavaq™ Q-TAMAVAQ6.docx`

This repository converts the manuscript concept into a **deterministic, auditable computational workflow**. It treats the quantum-walk/Grover language as a **selector and scoring abstraction** that can be simulated classically, logged, and tested. It does **not** claim clinical efficacy.

## Corrected scientific scope

TAMAVAQ is implemented here as a reproducible computational triage instrument:

- deterministic candidate indexing;
- NeoFox-like feature annotation and validation hooks;
- calibrated BIO probability field;
- unit-disciplined PHYS evidence field;
- Fubini--Study / kernel graph construction;
- continuous-time quantum-walk (CTQW) transport simulated with linear algebra;
- explicit predicate bits: `BIO`, `PHYS`, `MFG`, `PROV`, `VALID`, `MARK`;
- constraint-aware panel assembly;
- NeoDesign-style construct assembly;
- append-only batch-state and teleport-log JSON artifacts.

The output is a **version-bound shortlist and construct proposal**, not a medical recommendation.

## Repository layout

```text
src/tamavaq/
  pipeline.py           # end-to-end reproducible pipeline
schema/
  candidate.schema.json # input table contract
  batch_state.schema.json
  teleport_log.schema.json
examples/
  candidates.csv        # toy example; replace with real annotated candidates
  config.yaml
  run_example.sh
tests/
  test_pipeline.py
  test_reproducibility.py
docs/
  ALGORITHMS.md         # corrected algorithm specification
  REPRODUCIBILITY.md    # exact replay contract
  CITATIONS.md          # how cited methods are used
```

## Minimal install

```bash
python -m venv .venv
source .venv/bin/activate
pip install numpy pandas scipy pyyaml jsonschema pytest
```

Optional external tools for production feature annotation:

- NeoFox for feature annotation.
- NetMHCpan / MHCflurry / MixMHCpred / MixMHC2pred as available.
- neoDesign for protein construct optimization and lambda evaluation.
- NeoAgDT or an ILP solver for tumor-heterogeneity-aware exact optimization.

## Run the toy example

```bash
python -m tamavaq.pipeline \
  --candidates examples/candidates.csv \
  --config examples/config.yaml \
  --outdir runs/example
```

Expected outputs:

```text
runs/example/
  selected_panel.csv
  optimized_construct.txt
  batch_state.json
  teleport_log.jsonl
  scores.csv
```

## Input columns

Required columns:

```text
peptide, allele, gene, source_variant, calibrated_immunogenicity,
delta_g_kcal_mol, manufacturable, provenance_ok, expression_tpm,
vaf, clonality, anchor_type
```

Optional columns used if present:

```text
neofox_features_json, mhc_rank, cleavage_score, stability_score,
foreignness_score, tumor_cell_fraction
```

## Key design corrections

1. **No hidden ranking:** final selection is the result of logged predicates and constraints.
2. **No uncalibrated BIO gate:** immunogenicity must be calibrated or explicitly marked as uncalibrated in `batch_state.json`.
3. **No arbitrary docking-score physics:** PHYS uses `delta_g_kcal_mol` or a logged mapping from docking score to units.
4. **No opaque quantum claim:** CTQW and Grover-like amplification are simulated and versioned; hardware execution is optional and not required for reproducibility.
5. **No post-hoc feasibility edits:** DISTINCT/diversity/HLA quotas are enforced inside the panel assembly step.
6. **No silent construct edits:** NeoDesign-like ordering/linker choices are logged.

## Reproducibility contract

A run is replayable when the following are identical:

- input CSV bytes hash;
- config YAML bytes hash;
- package version or git commit hash;
- deterministic candidate index map;
- feature and calibration versions;
- CTQW parameters;
- predicate thresholds;
- panel constraints;
- construct assembly options.

