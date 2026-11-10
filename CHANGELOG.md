# Changelog

## 0.5.0 — 2026-11-10

- All identification entry points use the current migrated sphingolipid engine,
  including the historical MS2PipelineConfig constructor. Retired backend
  selectors and coverage/intensity/count parameters now fail explicitly.
- Retire the notebook generators and experimental evidence classifier/scorer;
  preserve archived result data containers and the existing SphinGOlipID frontend.
- Preserve complete external evidence references when combining text results.
- Document the eight reference-only historical types without enabling new rules.
- Scientific engine, fragment libraries, gates, weights and confirmation unchanged.

## 0.4.0 — 2026-11-10

- Import the attributed independent scientific-engine snapshot, sphingolipid
  libraries, type/name crosswalk and production adapters; retain the original UI.
