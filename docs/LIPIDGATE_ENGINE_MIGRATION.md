# SphinGOlipID engine migration

The source remains on `codex/sphingolipid-identification`. The independent
destination is the SphinGOlipID repository, on
`codex/lipidgate-sphingolipid-migration`. Neither branch is merged into main.
The fixed HG2/B2 baseline and original positive/negative libraries are retained.

The destination carries an attributed, isolated snapshot of the reference MS2
engine, its shared chain confirmation, and bounded SQLite library reader. It
does not depend on a LipidGate installation or its desktop interface. Snapshot
manifests record original and namespaced module hashes; scientific regression
oracles are ported separately. Existing local destination edits are preserved.

## Library extension

The historical MassHunter MS1 catalog contains 85,710 rows. It is a type/name
crosswalk, not a required record count. The current curated library and its
eligible chain combinations remain authoritative. Historical names are
reference candidates, not observed or confirmed lipid identities.
Sugar-linkage and epimer names are retained as reference aliases. They do not
become separately confirmed MS2 structures.

`extend_positive_catalog_glycosphingolipid_library.py` appends 8,488 records
across 44 additional canonical families: 32 glycan-composition families
(6,528 records); CAEP (268), N-methylCAEP (169), PG-Cer (169); four acetyl-HexCer
families (1,148); HexSPB/Hex3SPB/SHexSPB (71); MIP2C (134); and the d18:0
3-ketosphinganine reference (one). All new spectra are protonated total/species
compositions. Existing approved Cer or free-SPB records must support every
backbone formula. No new LCB/acyl grid is inferred from a historical sum.
MIP2C uses current eligible phyto-Cer compositions rather than recreating the
old list's unsupported hydroxylated combinations. Acetyl-HexCer aliases retain
fixed class modifications, not a new short-fatty-acyl grid or FA(2:0) assignment.
Acyl C14 and chain double-bond policies remain in force. Subsequent chain-limit
curation removes zero records. Every pre-existing block and its relative order
is verified unchanged. The positive source has 1,807,359 raw MSP blocks.

| New spectra | Admission | Pools and weights | Reporting |
| --- | --- | --- | --- |
| Glycan-composition totals | Y0 and Z0, at least two B products; sialic families additionally need their NeuAc/NeuGc/KDN B signature | HG/B/Other 50/25/25 | Total composition |
| CAEP, N-methylCAEP | Phosphonate head ion and neutral head loss | HG/Other 60/40 | Total composition |
| PG-Cer | Glycerophosphate neutral loss | HG/Other 60/40 | Total composition |
| Ac1/Ac2/Ac5/Ac6HexCer (FMC aliases) | HOAc loss and glycan+HOAc loss, plus one glycan-loss channel (at least three HG) | HG/Other 60/40; strongest two HG | Total composition; no acetyl-site assignment |
| MIP2C | Y0 and Z0, plus an inositol-phosphate-loss channel (at least three HG) | HG/Other 60/40; strongest two HG | Total composition |
| HexSPB/Hex3SPB | Full glycan loss and glycan+water loss (both HG) | HG/Other 60/40; strongest HG | Species composition |
| SHexSPB | Sulfate loss and sulfohexose loss (both HG) | HG/Other 60/40; strongest HG | Species composition |
| KetoSPB d18:0 | Formaldehyde loss | HG/Other 60/40 | 3KDS candidate; no carbonyl-site proof |

Two-channel HG and B pools use the strongest matched physical peak; when
three HG are required, the existing strongest-two policy applies. Other requires
raw intensity strictly above 0.5% and never admits a candidate. New records
contain no LCB/acyl scoring pool and cannot upgrade to chain composition. Existing Cer, SM,
HexCer, ganglioside, sodium-adduct and negative-mode policies are unchanged.

Y/Z masses follow the retained protonated ceramide; B masses follow residue
elemental compositions. These are conservative theoretical templates, not new
experimentally validated spectra. Glycan composition does not establish linkage
or monosaccharide stereochemistry. Formulas or names without a supported rule
remain reference-only and are listed in the coverage audit, rather than receiving
a fabricated permissive gate. S1P's three-reference restriction is preserved.

The glycan-fragment interpretation follows the positive CID data in
[the intact glycosphingolipid study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11024888/)
and [the plasma glycosphingolipid study](https://www.mdpi.com/2218-1989/11/3/140).
The phosphonate head ions at nominal 126 and 140 are supported by
[the original CAEP/CMAEP study](https://doi.org/10.1016/0005-2760(90)90154-P),
with modern CAEP characterization described in
[Wang et al.](https://doi.org/10.1016/j.foodchem.2019.125425).
PG-Cer's 172.0137 Da neutral loss is retained from the destination's previous
explicit glycerophosphate-loss rule; its present total-only template has not
been validated against a PG-Cer standard in this migration.

FMC family definitions follow the primary structural characterization of
[acetylated GalCer](https://pmc.ncbi.nlm.nih.gov/articles/PMC3035502/); retained
destination neutral-loss templates now use exact HOAc (60.0211 Da), not 60 Da.
The published spectra are lithium adducts, so this does not independently
validate all present protonated predictions. MIP2C masses follow its elemental
head residue; its positive predictions have not been validated with a standard.
The [3KDS standard study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5748491/)
reports 300.2 to 270.3, consistent with CH2O loss, not an assumed CO loss.
Only d18:0 is enabled; generic ketoSo and oxo-Cer grids are not extrapolated.

All historical types remain in the separate reference catalog, including
N-methyl/dimethyl sphingoid bases, ketoSo, 3-keto Cer and ambiguous O-acetyl
names. These do not receive fabricated distinguishing fragments. TrimethylSo's
intrinsic charge is not silently treated as [M+H]+; its adduct needs verification.
LM1's historical name/formula conflict is retained explicitly. Myriocin is marked
as a non-sphingolipid. The type audit distinguishes reference retention, MS2 type
availability and individual composition coverage; these are not measured recall.

## Reproducibility

`config/library_revision_sphingolipid_catalog_migration_20261110.json` records
the independent source library, curation audit and code fingerprint. Library
indexes are rebuilt with fresh normalization and blockwise roundtrip checks.
The original baseline manifest is never rebound to the new library.

The upstream `extend_verified_catalog_index.py` may reuse an original-prefix normalization
only when its complete effective rules fingerprint matches. It verifies the
entire original raw prefix, freshly normalizes every appended record, checks
stable IDs and equality of all scientific fields, and rebuilds search plans
with blockwise roundtrip verification. Appended curation/parent tags use separate
MSP headers so they cannot be parsed as part of polarity. A different rules
fingerprint or changed prefix is rejected; index hashes are never merely rebound.

The historical name catalog is stored separately from fragment libraries.
Its coverage status distinguishes a reference name from an available MS2
composition. It is not an admission filter or a source of chain confirmation.
`build_sphingolipid_name_catalog.py` reproduces this catalog and the per-type
crosswalk using the reviewed alias configuration and a freshly normalized index.
An isobaric ordinary Cer does not count as coverage of an esterified subtype.

In SphinGOlipID, `scripts/build_sphingolipid_library.py` freshly normalizes a
curated MSP into a new, separate SQLite file and writes a verification audit.
It uses only the independent package, preserves stable source IDs and derived
records, filters to sphingolipids, and rejects an existing destination. After
review, install the new source/index together with an updated library manifest;
never replace only its rules hash.
