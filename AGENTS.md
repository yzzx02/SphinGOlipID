# SphinGOlipID development

Continue independent development on `codex/lipidgate-sphingolipid-migration`
or a subsequent `codex/` branch. Do not merge into main without the user's
explicit request. The imported reference baseline is
`codex/sphingolipid-hg2-b2-baseline-20261109`; preserve its upstream manifest.

The authoritative identification entry points are `engine.py`,
`reference_pipeline.py`, the GUI and the three CLI commands documented in
README. They all use the isolated current reference engine. The old notebook
backend, experimental scorer and legacy selector are retired; reject retired
score/intensity/count parameters rather than silently accepting them. The
compatibility MS2PipelineConfig constructor also uses the current engine.
Archived evidence schema objects are read-only historical vocabulary, never a
source of current fragments, gates or scores.
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
mandatory water loss, 60/40. Cer1P H/Na require phosphate loss. Protonated
Cer1P and SM use HG/LCB/Other 60/20/20 with composition-specific chain evidence.
Their sodium adducts use HG/Other 60/40 and report total composition.

For positive GM2 and GM3, combine the curated core Y/Z and B channels into one
HG gate. Require at least half of the channels for the candidate precursor,
using distinct observed peaks. Keep chain confirmation separate.

GM1/GD1/GT1/GD3 require two distinct Y/Z HG observations including Y0/Z0,
two distinct B observations, NeuAc and the existing extended-HexNAc signature
where required. Score the strongest diagnostic glycan peak in HG and use
HG/LCB/Other 60/20/20. Preserve the independent glycan admission requirements.

Lyso-SM requires phosphocholine m/z 184.0733. The optional 124.9998 ion
belongs to HG, base-specific dehydration ions to LCB, and the choline ions
104.1070 and 86.0964 and precursor dehydration to Other. Hex-IPC, Gb4 and
doubly ammoniated GD1a/GD1b/GT1b/GQ1b use HG/LCB/Other. Preserve the
acyl-specific FAH pool for positive SL. Use named CAEP and PG-Cer chains only
from the curated evidence list, with subclass gates and an independent LCB
observation. Preserve the sum-composition entries and negative-ion policies.

Current curated chain combinations are authoritative. The historical 85,710
rows are only a type/name crosswalk; never recreate a full grid to reach that
count. Positive explicit acyls, including attached acyls, have C>=14 and DB<=3;
only the exact fifteen marked reference source exceptions are allowed. LCB
DB<=3 at any carbon length for explicitly curated positive LCBs. Never
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


Local literature checkpoint: `sphingolipid-literature-20261110-r1`, documented
in `config/scientific_revision_literature_20261110.json`. Keep the attributed
migration manifest unchanged. The subsequent positive LCB triene revision is
documented in `config/scientific_revision_lcb_triene_20261115.json`.
Native DMS coverage is only N,N-dimethylSo(d18:1) [M+H]+. Require at least two
distinct HG observations among 58.0651/84.0808/110.0964, including 110.0964.
Use strongest HG, HG/Other 60/40, raw Other >0.5%, no Other/LCB gate or fallback;
retain existing support saturation and candidate ranking. Never open a DMS
homolog grid or reinterpret these masses as the ordinary d/t LCB interface.
LM1's reference formulas include the missing NeuAc residue and remain separate
from nLM1/GM1 linkage claims. TrimethylSo reference masses are intrinsic [M]+;
the single common 60.08 amine does not enable MS2 coverage. The other rare
types remain archived references, outside default identification candidates.
