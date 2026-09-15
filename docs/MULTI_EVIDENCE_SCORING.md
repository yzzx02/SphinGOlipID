# Multi-evidence annotation/scoring

本模块基于 `542fa7f`，在现有 theoretical fragment generation 之上增加可审计的证据分类、门控和评分。
它目前通过显式调用的 adapter 返回独立结果，未接管 batch runner 的旧输出、Top-N 或历史表。
不得从最终 matched-only 表重建完整谱来计算新分数；需要候选完整理论片段和同一扫描的完整观测谱。

## 证据类别与来源

最终仅使用五类，`EvidenceType` 与 `EvidenceRole` 分开定义：

| fragment_type | evidence_role | 用途 |
|---|---|---|
| Precursor | precursor | 现有 MS1 exact-mass/adduct 匹配 prerequisite，排除在 MS/MS matching 和 weighted score 外 |
| HG | class_diagnostic | 现有头基、亚类诊断离子；如 SM 的 184.0733、NeuAc 的 292.10/274.09 |
| LCB | chain_specific | 现有 `LCB_fragment`/`LCB_mz` 字典中的确切标签；支持 LCB 碳数、不饱和度、羟基类型 |
| NL | structure_informative | 头基、糖基、主链顺序、分支及修饰的有来源中性丢失；包括作者指定的 Cer 例外 |
| common | supporting | 普通脱水、小型非特异离子和仅用于解释的未知标签；永不满足结构 gate |

`evidence.py` 中的集中规则按生成函数来源及其现有标签分类。生成函数加了只观察输出的 decorator；
只有 `capture_fragment_origins()` 上下文开启时才收集来源。LCB 加入时记录字典来源。
生成公式、质量、顺序和原始标签没有改动。旧表使用 exact-label fallback，不能仅凭一个 m/z 猜化学类别。
未知标签保持 `common`、`classification_status=UNSPECIFIED`、`scoring_eligible=False`。

同一候选完全相同的 theoretical m/z 合并成一个证据节点，保留所有 label。
跨类别解释采用 **HG > LCB > NL > common** 的固定优先级；不同但接近的 m/z 不按此优先级归并。
配对继续使用 `542fa7f` 的共享一对一 matcher：最多配对数 → 最小总绝对 ppm → 强度决胜。
分类不能让观测峰重复计分。MS1 Precursor 不作为一个 MS/MS fragment 参与配对。

## 作者确认的 Cer 例外

作者在本次实施中明确指定：普通 Cer 的 `M+H-H2O`、`M+H-2H2O` 可以作为 diagnostic NL。
代码识别其现有生成标签 `Cer(...)-1H2O` 和 `Cer(...)-2H2O`，以及上述标准表示，标注
`classification_source=author:Cer_dehydration`。这些片段的质量仍由原 `learn_fuc` 生成。
该决定覆盖本次最初“所有普通 H2O loss 都是 common”的要求，但只限 Cer 的一、二次脱水。
其他类别的普通脱水及 Cer 的更多次脱水仍为 common。

这是作者指定的 annotation policy，不是本次验证得到的普适亚类特异性结论。
一、二次脱水可能是相关过程；两条观测峰不等于两种独立化学特异性证据。其选择性仍需标准品/干扰物验证。
普通 Cer 因而使用 Precursor + diagnostic NL → Species level，再加 LCB → Molecular species level。

## Gate 与 score 独立

采用作者进一步确认的 **每个必需证据组分别至少达到 50%**。对于组的理论合格数量 `N_gate`：

```text
required = ceil(0.50 * N_gate)
PASS = precursor_pass AND N_gate > 0 AND qualified_matched >= required
```

HG 和 LCB 计 scoring-eligible 且具有对应诊断角色的理论片段；NL gate 只计其中 diagnostic NL 子集。
ordinary informative NL 可评分但不进入 diagnostic NL gate 分母或分子。未匹配理论片段留在分母。
三条 LCB 需至少两条，四条需至少两条，两条 diagnostic NL 需至少一条。缺少理论证据的组失败，不能除零后通过。
LCB 的标签身份还必须与候选名称中的 LCB 一致；其他 LCB 的匹配不能代替它。

合格观测证据要求 `normalized_intensity > 0` 且 `>= minimum_normalized_intensity`。
最低强度默认 **0.0**，表示不另设未校准的相对强度 cutoff，但零信号不会通过。
调用者仍应传入原 pipeline 的 raw `min_fragment_intensity`；示例使用 0，仅用于展示数学行为。
`minimum_group_fraction=0.50` 与最低强度都集中在 `EvidenceScoringConfig`。
敏感性表仅扫描归一化强度下限 0、0.01、0.05、0.10，固定 50% 门控；这不是阈值训练或新科学 cutoff。

| 类别规则 | Species level | Molecular species level |
|---|---|---|
| 有可靠 HG | Precursor + HG | 已满足 Species level + LCB |
| 无独立 HG，包括 Cer 作者例外和 neutral GSL | Precursor + diagnostic NL | 已满足 Species level + LCB |
| So | Precursor + appropriate LCB/backbone | 同一单链身份得到支持；不要求 FA |
| S1P、Lyso-SM、Lyso-sulfo | Precursor + HG | 已满足 Species level + LCB |
| Glu-So、Gb3-So | Precursor + diagnostic NL | 已满足 Species level + LCB |
| 未审阅/三链 OCer | UNSPECIFIED | UNSPECIFIED |

单链规则是独立声明。现有 `so_fuc()` 未追加 LCB 字典，`So()` 本身为 `pass`；因此当前生成结果不能靠
generic water 自动满足 backbone gate。本层不会为了通过 gate 新造片段。

稳定序列化值为 `candidate`、`species`、`molecular_species`；显示术语分别为
MS1 candidate、**Species level**、**Molecular species level**。RT 通过 `add_rt_validation()` 独立附加，
不修改注释层级、gate、权重或 MS/MS score。

## Group score 与 k

对同一完整 MS/MS spectrum，先用所有观测峰中的 base peak 归一化：

```text
I_ij = raw observed intensity / spectrum base-peak intensity
q_ij = I_ij / (I_ij + k_i)
S_i  = sum(q_ij over scoring-eligible theoretical fragments) / N_i
```

`I` 在 0–1 内；不直接把 raw intensity 代入饱和函数。未匹配的理论评分片段 `q=0`。
`N_i` 是该类**全部 scoring-eligible 理论 m/z 数量**，绝不是 matched count。
base peak 包括未匹配峰；整体信号乘 10 不改变新分数（前提是 raw matcher cutoff 未改变保留峰集合）。

| 配置 | 默认 k |
|---|---:|
| k_hg | 0.10 |
| k_lcb | 0.10 |
| k_nl_diagnostic | 0.10 |
| k_nl（ordinary informative） | 0.05 |
| k_common | 0.05 |

同一 I 下 k=0.10 比 k=0.05 得分低；高诊断价值的证据需要更高相对强度才能接近饱和。
权重总和为 100，但因为 q 的定义，即使 I=1，q 也小于 1；分数没有再人为拉伸到 100。
这些分数不是概率、FDR 或已经校准的鉴定置信度。

## Scoring-eligible policy

eligibility 在观察谱之前确定，不能因匹配好坏而删分母。

| 片段来源 | 默认评分集合 |
|---|---|
| 明确 HG | 精确列入 `HG_LABELS` 的离子 |
| LCB | 候选对应字典中现有离子，按完全相同 m/z 去重 |
| 固定头基/修饰 NL | 无附加水丢失的已有离散 NL；其 water satellite 保留解释但不评分 |
| GSL 主链/分支 | 当前生成集合中的首个主链丢失和完整糖链丢失两个无水拓扑锚点；HexCer 同质量只计一次 |
| 其他糖链路径 | 保留 NL 分类与匹配解释，不默认进入评分分母 |
| Cer 作者规则 | 一、二次脱水进入 diagnostic NL 评分集合 |
| common | 已知 generic water 和集中列出的 supporting ions |
| Unknown / OCer 待核验片段 | 不评分、不 gate |

GSL 两锚点是透明、**待作者验证可观测性的初始 policy**，不是论文已给出的必然可观测碎片。
它限制 branched-glycan 枚举路径造成的分母膨胀，但不保证不同结构的分数已具备统计可比性。
两个锚点若有多个解释路径，先按 m/z 去重。所有其他理论片段仍可匹配、输出解释。
旧 glycan 表没有生成来源时无法可靠还原锚点，默认只用于解释，可通过明确的
`eligibility_overrides={exact_label: bool}` 配置已审阅的代表集合。覆盖不会创造新片段。

SM 86.0964、104.1066 和 PG-Cer 93.0546 等小离子保守记 common，仅主要完整 HG ion 作 gate。
PI-Cer `[phosphate+H]+` 和 N-CAEP 的旧标签拼写保留原样，不通过本层改其质量或化学名称。

## 三套集中配置的模板

| 模板 | HG | LCB | NL | common |
|---|---:|---:|---:|---:|
| HG-dominant | 60 | 20 | 15 | 5 |
| NL-dominant glycosphingolipid | — | 35 | 55 | 10 |
| LCB-dominant | — | 60 | 30 | 10 |

模板权重严格沿用作者给定值，每套合计 100。类别映射见下表，也可查
`audit/multi_evidence_class_mapping.csv`。未知类别不自动套 GM3 权重。
若模板含有某组但此候选没有该组的 scoring-eligible 理论片段，则报告
`missing_theoretical_groups:...`，新总分为空；不以 0 补齐继续宣称模板有效，也不重分配权重。
So/S1P/Lyso/糖基 So 尚无作者指定的单链权重，模板明确为 UNSPECIFIED，保留独立 gate 结果。

`multi_evidence_score` 保留可计算的原始 weighted sum，供 gate-fail shadow 诊断；
`confidence_score` 仅在通过至少 Species level gate 后提供。没有通过 gate 的候选不获新排名。
新排序保留已有行序，另填 `multi_evidence_rank`；并列按 annotation 字典序确定。
旧 `匹配度分数`、`总分数`、`相对强度` 原样保留，另提供 `legacy_match_score`、`legacy_total_score`。

## FA assignment

对两链 Cer/SM/GSL，只有 LCB evidence gate 通过且 total composition 与候选 LCB/FA 的碳数、双键数互补一致时，
输出 `fa_assignment_mode=inferred_from_total_composition_minus_lcb`。
这表示 **precursor total composition + experimentally supported LCB → complementary fatty-acyl composition**。
它不是独立 FA 碎片确认；本实现始终不输出 `fa_directly_fragment_confirmed=True`。
缺少/不一致的总组成显式为 UNSPECIFIED；单链为 `not_applicable_single_chain`；三链不套两链相减。
不声称已解析 FA 双键位置、立体结构或羟基位置。

## 调用与产物

```python
from sphingolipid_toolkit.multi_evidence_adapter import (
    generate_evidence_fragments_for_candidates, evaluate_legacy_candidate,
)
from sphingolipid_toolkit.multi_evidence_scoring import PrecursorEvidence

# candidates 是已有 MS1 检索结果；本示例不读取任何历史数据。
theory = generate_evidence_fragments_for_candidates(candidates, observed_precursor, rt, abundance)
precursor = PrecursorEvidence.from_existing_match(
    existing_ms1_pass, observed_precursor, matched_library_mz, existing_adduct,
)
result = evaluate_legacy_candidate(
    "SM", "SM(d18:1/16:0)", theory.query("anno == 'SM(d18:1/16:0)'"),
    full_single_spectrum, precursor, legacy_row=legacy_result_row,
    ppm_tolerance=existing_fragment_ppm, min_fragment_intensity=existing_raw_cutoff,
    total_composition="SM(d34:1)",
)
# result.summary: 新旧分数、gate、层级和来源；result.fragments: 逐条 metadata/q。
```

MS1 pass 必须由现有匹配路径传入，adapter 不用 MS/MS 反向推断 MS1 pass，也不新增 MS1 tolerance。
`enabled=False` 直接返回旧结果行副本，不读谱、不增加字段。默认 batch runner 不调用本 adapter。
生成 adapter 使用新建临时目录，避免覆盖已有 `out_DB.csv`。

`python scripts/multi_evidence_shadow.py` 仅构造人工谱与实际规则生成的理论片段，并输出：

- `audit/multi_evidence_scoring_examples.csv`
- `audit/multi_evidence_gate_sensitivity.csv`
- `audit/multi_evidence_fragment_classification.csv`
- `audit/multi_evidence_class_mapping.csv`

这些表不能解释为历史 1,723 / 1,674 / 1,173 结果的影响统计。

## 完整 subclass → gate/template mapping

下表由集中 registry 导出。所有 HG/LCB/NL 条件均指该组分别达到 50%；所有 gate 还需 Precursor PASS。
<!-- CLASS_MAPPING -->

| lipid_class | has_HG | species_gate | molecular_species_gate | scoring_template | HG examples | LCB evidence | NL examples | common examples |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SM | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | [phosphocholine+H]+ | LCB_fragment / LCB_mz | C5H14NO4P; C3H9N | [C2H5NO+H]+; other precursor dehydration |
| PE-Cer | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | [phosphoethanolamine+H]+ | LCB_fragment / LCB_mz | phosphoethanolamine; C2H5N | [C2H5NO+H]+; other precursor dehydration |
| PI-Cer | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | [phosphate+H]+ | LCB_fragment / LCB_mz | phosphate; C6H10O5 | [C2H5NO+H]+; other precursor dehydration |
| PG-Cer | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | [phosphoglycerol+H]+ | LCB_fragment / LCB_mz | phosphoglycerol; C3H8O3 | [C2H5NO+H]+; other precursor dehydration |
| CAEP | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | [C2H8NO3P+H]+ | LCB_fragment / LCB_mz | C2H8NO3P | [C2H5NO+H]+; other precursor dehydration |
| N-CAEP | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | [C3H10NO3P+H]+] | LCB_fragment / LCB_mz | C3H10NO3P | [C2H5NO+H]+; other precursor dehydration |
| Cer1P | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | H4PO4 | LCB_fragment / LCB_mz | H3PO4 | [C2H5NO+H]+; other precursor dehydration |
| GM1 | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| GM1a | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| GM1b | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| GM2 | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| GM3 | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| GD1 | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| GD1a | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| GD1b | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| GD2 | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| GD3 | True | Precursor AND HG | Precursor AND HG AND LCB | HG-dominant | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| HexCer | False | Precursor AND NL | Precursor AND NL AND LCB | NL-dominant glycosphingolipid | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| GlcCer | False | Precursor AND NL | Precursor AND NL AND LCB | NL-dominant glycosphingolipid | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| GalCer | False | Precursor AND NL | Precursor AND NL AND LCB | NL-dominant glycosphingolipid | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| LacCer | False | Precursor AND NL | Precursor AND NL AND LCB | NL-dominant glycosphingolipid | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| Gb3 | False | Precursor AND NL | Precursor AND NL AND LCB | NL-dominant glycosphingolipid | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| Gb4 | False | Precursor AND NL | Precursor AND NL AND LCB | NL-dominant glycosphingolipid | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| GA1 | False | Precursor AND NL | Precursor AND NL AND LCB | NL-dominant glycosphingolipid | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| GA2 | False | Precursor AND NL | Precursor AND NL AND LCB | NL-dominant glycosphingolipid | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| type I B antigen | False | Precursor AND NL | Precursor AND NL AND LCB | NL-dominant glycosphingolipid | NeuAc+H; NeuAc+H-H2O | LCB_fragment / LCB_mz | first/complete glycan loss | [C2H5NO+H]+; other precursor dehydration |
| Cer | False | Precursor AND NL | Precursor AND NL AND LCB | LCB-dominant | — | LCB_fragment / LCB_mz | M+H-H2O; M+H-2H2O (author policy) | [C2H5NO+H]+; other precursor dehydration |
| FMC_1 | False | Precursor AND NL | Precursor AND NL AND LCB | NL-dominant glycosphingolipid | — | LCB_fragment / LCB_mz | Gal; HOAc | [C2H5NO+H]+; other precursor dehydration |
| FMC_3 | False | Precursor AND NL | Precursor AND NL AND LCB | NL-dominant glycosphingolipid | — | LCB_fragment / LCB_mz | Gal_OAc; HOAc | [C2H5NO+H]+; other precursor dehydration |
| FMC_5 | False | Precursor AND NL | Precursor AND NL AND LCB | NL-dominant glycosphingolipid | — | LCB_fragment / LCB_mz | Gal_4OAc; HOAc | [C2H5NO+H]+; other precursor dehydration |
| EO-Cer | False | Precursor AND NL | Precursor AND NL AND LCB | LCB-dominant | — | LCB_fragment / LCB_mz | E 18:2 | [C2H5NO+H]+; other precursor dehydration |
| EO-GlcCer | False | Precursor AND NL | Precursor AND NL AND LCB | NL-dominant glycosphingolipid | — | LCB_fragment / LCB_mz | E 18:2; Glc | [C2H5NO+H]+; other precursor dehydration |
| So | False | Precursor AND LCB | Precursor AND LCB | UNSPECIFIED | — | UNSPECIFIED: no LCB append in current single/three-chain path | — | [C2H5NO+H]+; other precursor dehydration |
| S1P | True | Precursor AND HG | Precursor AND HG AND LCB | UNSPECIFIED | H4PO4 | UNSPECIFIED: no LCB append in current single/three-chain path | H3PO4 | [C2H5NO+H]+; other precursor dehydration |
| Lyso-SM | True | Precursor AND HG | Precursor AND HG AND LCB | UNSPECIFIED | [phosphocholine+H]+ | UNSPECIFIED: no LCB append in current single/three-chain path | C5H14NO4P; C3H9N | [C2H5NO+H]+; other precursor dehydration |
| Lyso-sulfo | True | Precursor AND HG | Precursor AND HG AND LCB | UNSPECIFIED | [HSO4+H]+ | UNSPECIFIED: no LCB append in current single/three-chain path | C6H10O8S; SO3 | [C2H5NO+H]+; other precursor dehydration |
| Glu-So | False | Precursor AND NL | Precursor AND NL AND LCB | UNSPECIFIED | — | UNSPECIFIED: no LCB append in current single/three-chain path | C6H10O5; C6H10O5H2O | [C2H5NO+H]+; other precursor dehydration |
| Gb3-So | False | Precursor AND NL | Precursor AND NL AND LCB | UNSPECIFIED | — | UNSPECIFIED: no LCB append in current single/three-chain path | C6H10O5; 2*C6H10O5; 3*C6H10O5 | [C2H5NO+H]+; other precursor dehydration |
| 1-O-acetyl-Cer | False | UNSPECIFIED | UNSPECIFIED | UNSPECIFIED | — | UNSPECIFIED: no LCB append in current single/three-chain path | — | [C2H5NO+H]+; other precursor dehydration |
