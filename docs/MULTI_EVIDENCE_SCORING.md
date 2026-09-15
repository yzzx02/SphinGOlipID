# Multi-evidence scoring: fixed 60 / 20 / 20

本版基于 `bd428de`，按作者最终决定统一为三个 scoring pools。新层仍是显式调用的 adapter，
未接管默认 batch 输出、legacy 排名或历史表。五类 fragment evidence 被映射为三个 scoring pools：
HG、LCB 或 diagnostic NL 根据 subclass 充当 primary/secondary structural evidence，
ordinary NL 与 common 构成数量型 Support pool。

## 与 LipidGate 实现的对应

已读取 `yzzx02/LipidGate` commit `747dfacc4789c0d74e0f00fa2c7696ee8c423c08` 的
[`src/lipidgate/ms2/scoring.py`](https://github.com/yzzx02/LipidGate/blob/747dfacc4789c0d74e0f00fa2c7696ee8c423c08/src/lipidgate/ms2/scoring.py)。
核对 `_saturation_fragment_quality`、`_calculate_pool_scores`、`_non_precursor_base_intensity`
和 `_non_precursor_quality_overrides`。本版使用相同的结构质量函数、max aggregation、Support count saturation
和非母离子归一化窗口。参考记录：`audit/lipidgate_scoring_reference.json`。

**有意不复制** LipidGate 的 dynamic active-pool weight redistribution：本项目按作者要求固定 60/20/20。
本项目的 subclass registry、50% gates、scoring eligibility、一对一 matcher 和 annotation policy 保持独立。
这里的函数对齐不表示复制了 LipidGate 的全部类别规则或 scoring pipeline。

## 五类 metadata 与三个 scoring pools

| fragment_type | evidence_role | scoring_pool |
|---|---|---|
| Precursor | precursor | NA；现有 MS1 prerequisite，不进入 MS/MS score |
| HG | class_diagnostic | 根据 subclass 进入 primary/secondary，否则 NA |
| LCB | chain_specific | 根据 subclass 进入 primary/secondary，否则 NA |
| NL | diagnostic_nl | 根据 subclass 进入 primary/secondary；**永不进入 Support** |
| NL | supporting_nl | support；不参与结构 gate |
| common | supporting | support；不参与结构 gate |

`fragment_type != scoring_pool`。未映射为结构池的 diagnostic NL 保留化学类型和解释信息，不改成 supporting NL。
例如 GM3 的结构池为 HG/LCB，它的 diagnostic glycan NL 不重复计入 Support。
`scoring_eligible=False` 的片段仍可 matching/explanation，但不参与任何评分或 gate 分母。

`evidence.py` 保留生成来源 capture、exact legacy labels 和所有原始 id/mz 字段。不凭 m/z 数值猜化学类别。
未知来源保持 `UNSPECIFIED` 且不评分。旧 NL role `structure_informative` 可根据原 diagnostic_strength
迁移为 diagnostic_nl/supporting_nl，不改变质量。

## 集中 subclass registry

`STRUCTURAL_EVIDENCE_REGISTRY` 定义 primary_evidence、secondary_evidence、species_groups、molecular_groups。
原 `CLASS_RULES` 名称作为 registry alias 保留，但旧多模板权重及类别 k 已从当前评分代码移除。

| 类别组 | primary | secondary | Species level | Molecular species level |
|---|---|---|---|---|
| SM、GM/GD、PE/PI/PG-Cer、CAEP、N-CAEP、Cer1P | HG | LCB | Precursor + HG | Species + LCB |
| HexCer/GlcCer/GalCer、LacCer/Hex2Cer、Gb3/Gb4、neutral GA、type-I-B、FMC、EO-GlcCer | diagnostic NL | LCB | Precursor + diagnostic NL | Species + LCB |
| Cer；保留的 EO-Cer 映射 | LCB | diagnostic NL | Precursor + diagnostic NL | Species + LCB |
| So/S1P/Lyso/Glu-So/Gb3-So | UNSPECIFIED | UNSPECIFIED | 保留各自 gate | 保留各自 gate，不推 FA |
| 三链 OCer / 未注册类 | UNSPECIFIED | UNSPECIFIED | 规则未定者不升级 | 不套两链相减 |

单链类没有作者最终指定的 scoring policy，故 `scoring_policy=UNSPECIFIED`，结构池、权重及总分均为 NA。
其支持计数和已有 gate/level 仍可单独报告。实际 `So()` 为 pass，`so_fuc()` 不追加 LCB 字典，
所以不能凭现有 generic water/common 输出伪造 backbone support。未为单链新造权重或离子。

## Cer 作者指定 policy

Cer 的 `M-H2O`、`M-2H2O`（这里 M 表示已选 precursor ion），以及原标签
`Cer(...)-1H2O` / `Cer(...)-2H2O`、`M+H-H2O` / `M+H-2H2O`，属于
`NL + diagnostic_nl`，进入 Cer 的 secondary structural pool。
来源明确为 `author:Cer_dehydration`。

这是 **author-defined Cer structural evidence policy**。这些脱水是当前 annotation policy 中的结构支持证据，
**不能单独证明 Cer subclass identity，也不能排除 HexCer in-source fragmentation**，不声称它们是普适特异性
Cer diagnostic ions。相关 dehydration 不等于独立化学特异性。其他类别普通脱水及 Cer 更多次脱水仍是 common。

## Manuscript-compatible formulas

```text
Score = 60*S_primary + 20*S_secondary + 20*S_Supp

q_ij = (1+k_i)*I_ij / (I_ij+k_i)
S_primary   = max(q_ij of matched eligible primary evidence),   k=0.10
S_secondary = max(q_ij of matched eligible secondary evidence), k=0.05

support_target = min(N_Supp,3)
S_Supp = min(n_matched/support_target,1), if support_target > 0
```

结构池不平均、不除以 theoretical count。相同 I（0<I<1）时 primary 的 k=0.10 得分低于 secondary 的 k=0.05。
I=0 时 q=0，I=1 时 q=1，I/q 均限制在 0–1。理论池存在但没有匹配时，结构池为 0。
理论池根本不存在时为 NA，报告 missing_structural_theory，总分不可用，不假装完成评价。

| 集中配置 | 默认值 |
|---|---:|
| PRIMARY_POOL_SATURATION_HALF_INTENSITY | 0.10 |
| SECONDARY_POOL_SATURATION_HALF_INTENSITY | 0.05 |
| PRECURSOR_CLUSTER_EXCLUSION_DA | 4.1 |
| SUPPORT_SATURATION_COUNT | 3 |
| POOL_WEIGHTS primary / secondary / support | 60 / 20 / 20 |

GM3/SM：`60*S_HG + 20*S_LCB + 20*S_Supp`。
HexCer/LacCer：`60*S_NL + 20*S_LCB + 20*S_Supp`。
Cer：`60*S_LCB + 20*S_NL + 20*S_Supp`。
这些权重没有用实验结果训练或调整；分数不是概率/FDR。

## Intensity normalization

需要完整、单一扫描的 MS/MS spectrum，包括未匹配峰。
基本 `normalized_intensity = intensity / spectrum base peak`。
结构质量使用 `normalized_structural_intensity`：

1. 从观测谱中排除与已选 precursor 的 m/z 距离 **<=4.1 Da** 的峰，再找剩余峰的最大强度。
2. 对窗口外的已匹配结构峰，以 non-precursor base peak 归一化，限制到 0–1。
3. 与 LipidGate 一样，窗口内的非母离子匹配不使用该 override，保留 spectrum-relative intensity。
4. 不存在正值 non-precursor base peak 时，fallback 到 spectrum-relative intensity，并记录状态。
5. 缺少 actual precursor m/z 时明确报告 missing_precursor_mz_spectrum_relative_fallback，不猜 m/z。

`PrecursorEvidence.from_existing_match()` 保存 actual observed m/z；也可传 `precursor_mz=`；
legacy theory 的唯一 `target` 是最后的兼容来源。Precursor 类型不作为 MS/MS 节点，
恰好位于已选 precursor m/z 的 centroid 即使被错误标签称为结构片段，也不贡献 gate/score。

极强 precursor/isotope 峰不会压低窗口外结构峰的质量。Support 计数不使用这些强度质量值。
Gate 的最低归一化强度保留上一版 spectrum-relative 定义，不由结构评分的 override 偷改。
默认 minimum_normalized_intensity=0，但必须有正信号；调用方原 raw matcher cutoff 仍照常生效。

## Support 只按匹配数量

`N_support_total` 为 scoring-eligible common + supporting_nl 的唯一 theoretical m/z 数量。
`N_support_matched` 为一对一匹配得到的正信号、独立 support 节点数量。HG、LCB、diagnostic NL、Precursor 均排除。

| theory / matched | S_supp |
|---|---:|
| 8 / 1 | 1/3 |
| 8 / 2 | 2/3 |
| 8 / 3 或更多 | 1 |
| 2 / 1 | 0.5 |
| 2 / 2 | 1 |
| 1 / 1 | 1 |

Support 不计算 q，逐片段输出中的 support q 为 NA。只要匹配集合未变，改变支持峰强度不改变 S_supp。
改变强度可能影响一对一配对/原 cutoff 或 non-precursor base peak；这些是不同于 Support count 公式的上游影响。

无 support theory 时：`support_status=no_support_theory`，`S_supp/S_Supp=NA`，Support 贡献为 0，
20 分不转给结构池；结构两池都存在时总分最多 80。不会除零，也不会把缺失池伪装成评价为零的池。

## Gate 与 score 完全分开

每个 required HG/LCB/diagnostic NL 组分别满足：

```text
N_gate > 0
qualified_matched >= ceil(0.50*N_gate)
```

理论 gate 分母只含 scoring/gate-eligible evidence；ordinary NL/common 永不参与。
HG 两条需一条，LCB 三条需两条，diagnostic NL 四条需两条。LCB 标签身份需与候选一致。
Precursor PASS 是总前提，Molecular species level 必须先满足 Species level。

`1/4` primary 匹配即使 q=1 仍 gate fail；`2/4` LCB 匹配可 gate pass，池分数却是两者 q 的最大值，
不是 2/4、平均值或归一化覆盖率。Cer 的 species gate 是 secondary 的 diagnostic NL，不是机械要求 primary 先通过。

`multi_evidence_score` 保留 raw weighted sum 供 audit；`confidence_score` 仅在至少 Species gate 通过时提供。
不通过结构 gate 的 candidate 不获得新排名。新排序不覆盖 legacy 排序。
`add_rt_validation()` 仅附加正交 RT 状态，不改变 gate、score 或 annotation level。

## Eligibility 与相同质量节点

保留 bd428de 的生成 provenance 和 eligibility policy：明确 HG、LCB 字典、无附加水的固定诊断 loss，
以及 GSL 首个/完整糖链无水 loss 锚点进入各自合格集合。其他 glycan sequential/branch/topology losses
仍为 diagnostic_nl；未指定为合格的保留解释，可用显式 exact-label eligibility override 审阅后启用。
不因分数高低修改 eligibility，也不因为 max(q) 而把 explanation-only diagnostic NL 移入 Support。
糖链锚点可观测性仍待验证；本次未改变 parser、生成算法或片段质量。

同一 annotation 内仅在 **theoretical m/z 完全相同** 时，以当前 subclass 的
`primary > secondary > support > unmapped explanation` 选择唯一计分角色。
例如 Cer 同质量 LCB/diagnostic NL 优先 LCB；HexCer 同质量 diagnostic NL/LCB 优先 diagnostic NL。
保留原标签与 `evidence_annotations` JSON，避免丢失其他解释。
同一 centroid 不进入两个结构池或 structural+support；相近但不同 m/z 仍由原一对一 matcher 决定配对。
没有更改 `542fa7f` 的最大配对数、最小总绝对 ppm、强度决胜语义。

## 输出字段与兼容

新增/更新 primary_evidence_type、secondary_evidence_type、primary_gate_pass、secondary_gate_pass、
primary_matched_count、primary_theoretical_count、secondary_matched_count、secondary_theoretical_count、
S_primary、S_secondary、S_supp、固定权重、new_score/multi_evidence_score，以及 S_HG/S_LCB/S_NL/S_Supp aliases。
没有对应 pool 的 alias 为 NA。matched_count 是实际计分节点计数；更严格强度 minimum 下的 gate 分子/要求
另在 gate_summary 中完整记录。CSV 同时包含作者要求的 `S_supp` 和 `S_Supp`，读取时应保留大小写区分。

旧 `匹配度分数`、`总分数`、`相对强度` 及 legacy_match_score/legacy_total_score 保留。
序列化 level 保持 candidate/species/molecular_species；显示正式术语 Species level / Molecular species level。
FA 仍为 total composition 减去已支持 LCB 的 complementary inference，需碳数/双键数相符；
从不称为独立 FA fragment confirmation；单链不推 FA，三链保持 UNSPECIFIED。

```python
from sphingolipid_toolkit.multi_evidence_scoring import PrecursorEvidence
from sphingolipid_toolkit.multi_evidence_adapter import evaluate_legacy_candidate

precursor = PrecursorEvidence.from_existing_match(
    existing_ms1_pass, actual_precursor_mz, library_mz, actual_adduct,
)
result = evaluate_legacy_candidate(
    lipid_class, annotation, candidate_theory, full_single_spectrum, precursor,
    legacy_row=legacy_row, total_composition=total_composition,
    ppm_tolerance=existing_fragment_ppm, min_fragment_intensity=existing_raw_cutoff,
)
```

## Synthetic audit

运行 `python scripts/multi_evidence_shadow.py` 仅构造人工谱和原规则理论片段，固定覆盖 GM3、SM、HexCer、LacCer、Cer，
另保留 GM3 missing-HG 和 So UNSPECIFIED 示例。输出 scoring examples、gate sensitivity、fragment classification、class mapping。
对照脚本隔离加载 Git `bd428de` 的旧 evidence scoring，仅用于 audit，不是当前 production scoring 的可选模板。
同时保留 legacy score 和 bd428de_score/rank/level。旧版本示例表另存为 `multi_evidence_scoring_examples_bd428de.csv`。
敏感性表固定 50% 与全部权重/k，仅展示既定强度 minimum 网格；不据此挑选阈值。

**没有重跑、覆盖或验证历史 1,723 / 1,674 / 1,173 结果。**

## 完整 class mapping

<!-- CLASS_MAPPING -->

| lipid_class | primary_evidence | secondary_evidence | species_gate | molecular_species_gate | scoring_policy | scoring_weights |
| --- | --- | --- | --- | --- | --- | --- |
| SM | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| PE-Cer | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| PI-Cer | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| PG-Cer | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| CAEP | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| N-CAEP | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| Cer1P | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| GM1 | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| GM1a | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| GM1b | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| GM2 | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| GM3 | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| GD1 | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| GD1a | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| GD1b | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| GD2 | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| GD3 | HG | LCB | Precursor AND HG | Precursor AND HG AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| HexCer | diagnostic_nl | LCB | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| GlcCer | diagnostic_nl | LCB | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| GalCer | diagnostic_nl | LCB | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| LacCer | diagnostic_nl | LCB | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| Hex2Cer | diagnostic_nl | LCB | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| Gb3 | diagnostic_nl | LCB | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| Gb4 | diagnostic_nl | LCB | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| GA1 | diagnostic_nl | LCB | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| GA2 | diagnostic_nl | LCB | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| type I B antigen | diagnostic_nl | LCB | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| Cer | LCB | diagnostic_nl | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| FMC_1 | diagnostic_nl | LCB | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| FMC_3 | diagnostic_nl | LCB | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| FMC_5 | diagnostic_nl | LCB | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| EO-Cer | LCB | diagnostic_nl | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| EO-GlcCer | diagnostic_nl | LCB | Precursor AND NL | Precursor AND NL AND LCB | structural_60_20_20 | 60 / 20 / 20 |
| So | UNSPECIFIED | UNSPECIFIED | Precursor AND LCB | Precursor AND LCB | UNSPECIFIED | UNSPECIFIED |
| S1P | UNSPECIFIED | UNSPECIFIED | Precursor AND HG | Precursor AND HG AND LCB | UNSPECIFIED | UNSPECIFIED |
| Lyso-SM | UNSPECIFIED | UNSPECIFIED | Precursor AND HG | Precursor AND HG AND LCB | UNSPECIFIED | UNSPECIFIED |
| Lyso-sulfo | UNSPECIFIED | UNSPECIFIED | Precursor AND HG | Precursor AND HG AND LCB | UNSPECIFIED | UNSPECIFIED |
| Glu-So | UNSPECIFIED | UNSPECIFIED | Precursor AND NL | Precursor AND NL AND LCB | UNSPECIFIED | UNSPECIFIED |
| Gb3-So | UNSPECIFIED | UNSPECIFIED | Precursor AND NL | Precursor AND NL AND LCB | UNSPECIFIED | UNSPECIFIED |
| 1-O-acetyl-Cer | UNSPECIFIED | UNSPECIFIED | UNSPECIFIED | UNSPECIFIED | UNSPECIFIED | UNSPECIFIED |
