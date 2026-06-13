# Input Format

This page describes the v0.3 backend inputs consumed by `run_batch(config)`.

## Raw MS/MS text files

Default file names:

```text
hilic-msms-1.txt
hilic-msms-2.txt
...
hilic-msms-6.txt
```

The current parser supports the Agilent-style export used by the converted notebook. It expects spectrum blocks that contain:

- `index:...`
- `target m/z,<value>, m/z`
- `scan start time,<value>,...`
- two `binaryDataArray` blocks, first for fragment m/z and second for fragment intensity

`io_utils.parse_ms2_text_to_dataframe()` returns the normalized columns:

```text
file_name
scan_id
precursor_mz
rt
fragment_mz
fragment_intensity
```

## Precursor target Excel files

Default file names:

```text
hilic-msms-1.xlsx
hilic-msms-2.xlsx
...
hilic-msms-6.xlsx
```

Required columns:

```text
target
下限
上限
RT
Abund
```

The columns above must be numeric. Validation errors are raised before the legacy MS2 algorithm starts.

## MS1 theoretical library

Required columns:

```text
理论值
classy
name
structure
```

`理论值` must be numeric. The `classy`, `name`, and `structure` columns feed the legacy rule-based fragment generator.

