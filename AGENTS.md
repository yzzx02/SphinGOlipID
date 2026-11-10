# SphinGOlipID development

Continue independent development on `codex/lipidgate-sphingolipid-migration`
or a subsequent `codex/` branch. Do not merge into main without the user's
explicit request. The imported reference baseline is
`codex/sphingolipid-hg2-b2-baseline-20261109`; preserve its upstream manifest.

The authoritative identification entry points are `engine.py`,
`reference_pipeline.py`, the GUI and the three CLI commands documented in
README. They use the isolated reference engine by default. `backend=legacy`
is an explicit compatibility route, not the scientific production default.
Keep and improve SphinGOlipID's existing frontend. Do not migrate LipidGate's
frontend, layouts or Qt application; adapt only the existing UI's backend calls
and necessary parameter/result fields.

Preserve admission, numeric scores, candidate order and full audit provenance
when editing storage, exports or interfaces. Retain the matched-fragment-count
Top1 split. Do not delete d/h or t alternatives by isomer preference.

Other never gates and counts only raw peaks strictly above 0.5% of the original
spectrum base peak. Chain confirmation is separate: d/t use the actual
composition-dependent LCB-H2O/LCB-2H2O, inclusive at 1%. Never use four fixed
masses or infer the hydroxyl location from a shared peak. t LCB-CH6O3 support
is removed; the distinct d-series channel remains.

Cer H uses LCB/Other 60/40 and a half-LCB gate. Strongest LCB scores there;
other independently matched LCB peaks >0.5% provide auxiliary Other support.
One centroid scores once. HexCer and SHexCer reuse auxiliary support with
HG/LCB/Other 60/20/20. SHexCer/+O use the strongest HG by default and no LCB
gate. LacCer/Hex2Cer have the two whole-disaccharide-plus-water HG channels,
either admitting, strongest HG, 60/20/20. Ordinary Cer Na is total-only with
mandatory water loss, 60/40. Cer1P H/Na require phosphate loss, 60/40; H uses
shared d/t chain confirmation, Na reports totals.

GM1/GD1/GT1/GD3 require two distinct Y/Z HG observations including Y0/Z0,
two distinct B observations, NeuAc and the existing extended-HexNAc signature
where required. HG/B/Other stays 50/25/25; HG and B use the strongest physical
peak. Existing negative and legacy linkage-isomer policies remain separate.

Current curated chain combinations are authoritative. The historical 85,710
rows are only a type/name crosswalk; never recreate a full grid to reach that
count. Positive explicit acyls, including attached acyls, have C>=14 and DB<=3;
only the exact fifteen marked reference source exceptions are allowed. LCB
DB<=2 except C18 may DB3; keep m18:3, with no LCB carbon-length filter. Never
infer acyls from totals. Source-only modifications without distinguishing
evidence remain explicitly labelled MS1 references, not confirmed MS2 types.

Libraries live in `src/sphingolipid_toolkit/data/ms2`. Curate a separate
destination, preserve unaffected blocks and order, apply chain-limit curation,
freshly normalize and verify indexes, then update the library/engine manifests.
Do not rebind an old rules fingerprint to different code. Snapshot provenance
is in `config/engine_revision.json`; changes to vendored scientific code require
a new checkpoint and scientific regression verification.

Run `python -m pytest -q` after functional changes. Keep algorithm and library
commits independently traceable. Use Git LFS for MSP, SQLite and compressed
reference-name assets; state their sizes before pushing. Never commit raw
samples, research outputs, caches or builds. Preserve existing unrelated edits.
