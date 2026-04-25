# Reproducibility and audit contract

## Immutable inputs

Each run stores SHA-256 hashes of:

- candidate table;
- configuration file;
- candidate index map;
- thresholds;
- graph/Laplacian arrays;
- selected output table.

## Batch state

`batch_state.json` records:

- batch id and timestamp;
- software versions;
- config hash;
- candidate table hash;
- deterministic candidate index map hash;
- calibration status and metrics;
- predicate thresholds;
- CTQW parameters;
- Laplacian hash;
- predicate vectors;
- selected candidate ids;
- construct assembly choices.

## Teleport log

`teleport_log.jsonl` is an append-only event ledger for module handoffs. In this corrected implementation, “teleportation” is treated as a transparent boundary-log abstraction. Each event records:

- module context;
- deterministic pseudo Bell outcomes generated from the event hash seed;
- Pauli-frame label for audit replay;
- predicate snapshot;
- checksum binding the event to `batch_state.json`.

## Determinism

The pipeline must use:

- stable sorting for candidate indexes;
- explicit random seeds;
- pinned config values;
- deterministic tie-breaking by canonical candidate id.

