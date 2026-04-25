# Corrected TAMAVAQ / Q-TAMAVAQ algorithms

## Algorithm 1 — Candidate canonicalization and deterministic indexing

**Input:** candidate table `C`, config `θ`  
**Output:** indexed candidate table `C*`, index-map hash

1. Validate required columns.
2. Canonicalize each candidate id as:
   `sha256(peptide|allele|gene|source_variant|anchor_type)`.
3. Sort by `(canonical_id, peptide, allele)`.
4. Assign integer index `i = 0..n-1`.
5. Hash the index map and store in `batch_state.json`.

Correction: never use input row number as a stable identity.

## Algorithm 2 — NeoFox-compatible feature evidence object

**Input:** indexed candidates `C*`  
**Output:** evidence matrix `X`, feature manifest

For every candidate, build an evidence object:

```text
E_i = {
  peptide/allele/gene fields,
  antigen-presentation features,
  recognition/foreignness features,
  expression and VAF/clonality fields,
  calibrated immunogenicity,
  unit-disciplined ΔG,
  manufacturability/provenance flags,
  optional NeoFox feature JSON
}
```

Correction: keep feature blocks separately logged. Do not collapse biological evidence into a single undocumented score.

## Algorithm 3 — Calibration gate

**Input:** uncalibrated or calibrated immunogenicity score `p_i`  
**Output:** `BIO_i`

1. If model is calibrated, use `p_i` directly.
2. If uncalibrated, either fit a calibration model on validation data or mark run as exploratory.
3. Set `BIO_i = 1[p_i >= τ_bio]`.
4. Store calibration method, metrics, or `exploratory_uncalibrated=true`.

Correction: ranking performance is not a substitute for calibration.

## Algorithm 4 — Unit-disciplined physics gate

**Input:** `delta_g_kcal_mol_i` or docking score + logged mapping  
**Output:** `PHYS_i`

1. Prefer measured/estimated `ΔG` in kcal/mol.
2. If using docking score, apply a versioned regression mapping and uncertainty model.
3. Set `PHYS_i = 1[ΔG_i <= τ_dg]`.
4. Log physical unit, model version, and uncertainty.

Correction: never compare arbitrary docking scores as physical free energies without a logged mapping.

## Algorithm 5 — Fubini--Study / property-kernel graph

**Input:** normalized feature vectors `z_i`, kernel scale `σ`  
**Output:** kernel `K`, distance `D_FS`, affinity `W`, Laplacian `L`

1. Normalize `z_i` to unit norm.
2. Compute fidelity kernel `K_ij = |z_i · z_j|^2`.
3. Clip `K_ij` to `[0, 1]`.
4. Compute `D_FS_ij = arccos(sqrt(K_ij))`.
5. Compute heat weights `W_ij = exp(-D_FS_ij^2 / (2σ^2))` for `i != j`, else `0`.
6. Compute degree `D = diag(sum_j W_ij)`.
7. Compute Laplacian `L = D - W`.

Correction: when strict k-mer similarity collapses, use predicate-aligned property fingerprints.

## Algorithm 6 — CTQW transport scoring

**Input:** Laplacian `L`, prior distribution `π0`, time `t`  
**Output:** transported score `π_t`

1. Initialize amplitude `ψ0_i = sqrt(π0_i)`.
2. Compute `ψ_t = exp(-i L t) ψ0`.
3. Set `π_t_i = |ψ_t_i|^2`.
4. Normalize numerically.
5. Store eigenvalues, entropy, and Laplacian hash.

Correction: CTQW is used as a transparent transport smoother; it does not create biological evidence.

## Algorithm 7 — Predicate oracle in classical reproducible form

**Input:** candidate evidence, CTQW score, thresholds  
**Output:** predicate table

```text
BIO_i   = calibrated_immunogenicity_i >= τ_bio
PHYS_i  = delta_g_kcal_mol_i <= τ_dg
MFG_i   = manufacturable_i == true
PROV_i  = provenance_ok_i == true
VALID_i = no missing required fields and peptide alphabet valid
GEOM_i  = ctqw_score_i >= τ_ctqw
MARK_i  = BIO_i & PHYS_i & MFG_i & PROV_i & VALID_i & GEOM_i
```

Correction: make the oracle a pure function of logged predicate bits.

## Algorithm 8 — Grover-style amplification simulator

**Input:** `MARK`, iteration count `r`  
**Output:** amplified sampling weights

For reproducibility, use analytic amplitude amplification rather than unlogged random sampling:

1. Let `m = sum(MARK)` and `n = len(MARK)`.
2. If `m = 0`, return base transport scores.
3. Let `θ = asin(sqrt(m/n))` and target amplitude mass `sin((2r+1)θ)^2`.
4. Redistribute that mass over marked candidates proportional to transported scores.
5. Redistribute the remaining mass over unmarked candidates.

Correction: shallow Grover is an auditable weighting policy, not a claim of quantum hardware advantage.

## Algorithm 9 — NeoAgDT-style constrained panel optimization

**Input:** scores, predicate table, budget `k`, quotas  
**Output:** selected panel `S`

1. Restrict candidates to `MARK=1`.
2. Enforce DISTINCT by canonical id/peptide sequence.
3. Enforce HLA/category/anchor quotas.
4. If tumor-cell coverage matrix is available, solve an ILP maximizing expected covered tumor-cell response.
5. Otherwise solve deterministic greedy selection:
   - sort by amplified score, expression, clonality, then canonical id;
   - add candidates that preserve quotas until `k` is reached;
   - emit a feasibility warning if fewer than `k` candidates pass.

Correction: do not select the top `k` raw scores when constraints exist.

## Algorithm 10 — NeoDesign-style construct assembly

**Input:** selected peptides, linker library, risk scorer  
**Output:** ordered construct sequence, linker manifest

1. Build a graph where nodes are selected peptides.
2. Edge cost = junctional neoantigen risk + domain penalty + looseness/compactness penalty + linker penalty.
3. Find an ordering by dynamic programming for small `k`, deterministic greedy for large `k`.
4. At each junction, test direct concatenation.
5. If direct junction fails risk threshold, screen linker library and pick the minimal-risk linker.
6. Output final polyvalent protein sequence.
7. Store ordering, linker choices, and a lambda-evaluation placeholder for downstream mRNA design.

Correction: separate peptide selection from construct design; junctional risk is handled after panel feasibility.

