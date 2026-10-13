# Fixed three-pool scoring implementation report

基线 `bd428de`；本版按作者最终确认的 60/20/20 统一 multi-evidence score，保留显式 opt-in adapter。
没有修改 fragment chemistry/masses、glycan parser/generation、一对一 matcher、RT Linear-first 或历史输出。

## 移除/替换的旧逻辑

删除当前生产模块中的旧 HG/LCB/NL/common 多模板权重与类别 k 配置，移除 structural average(q)/theoretical-count denominator，
停止对 Support 使用强度质量函数。旧 NL role 细分为 diagnostic_nl/supporting_nl。
旧版本公式仅能在隔离的 `bd428de` Git audit 对照中执行，不再作为当前评分模板提供。

保留五类化学 metadata、来源、eligibility、50% gate、annotation level、FA complementary inference 与 legacy score coexistence。
相同 m/z 的唯一计分解释改按 subclass primary/secondary/support 顺序选择，底层 centroid matching 算法完全未动。

## 最终 mapping

| 类别组 | primary | secondary | weights |
|---|---|---|---|
| HG-dominant：SM、GM/GD、PE/PI/PG-Cer、CAEP、N-CAEP、Cer1P | HG | LCB | 60/20/20 |
| neutral GSL：HexCer、LacCer/Hex2Cer、Gb3/Gb4、GA、type-I-B 等 | diagnostic NL | LCB | 60/20/20 |
| Cer；沿用 EO-Cer 的结构 mapping | LCB | diagnostic NL | 60/20/20 |
| So、S1P、Lyso-SM、Lyso-sulfo、Glu-So、Gb3-So | UNSPECIFIED | UNSPECIFIED | UNSPECIFIED |
| OCer 三链与未注册类 | UNSPECIFIED | UNSPECIFIED | UNSPECIFIED |

完整 registry 表由脚本导出到 `audit/multi_evidence_class_mapping.csv` 和
[MULTI_EVIDENCE_SCORING.md](MULTI_EVIDENCE_SCORING.md)。没有为单链发明新权重。
Cer 一、二次脱水保留 author-defined structural evidence policy；不能单独证明 Cer subclass，也不能排除 HexCer in-source fragmentation。

## LipidGate 对齐与有意差异

实际读取参考 commit `747dfacc4789c0d74e0f00fa2c7696ee8c423c08` 的 `src/lipidgate/ms2/scoring.py`：

- q=(1+k)I/(I+k)，primary k=.10、secondary k=.05；q endpoints 和 clipping 一致。
- 结构池采用 matched max(q)，不平均、不除理论数量。
- Support 采用 min(matched/min(theory,3),1)，不使用强度质量。
- non-precursor base peak 排除 ±4.1 Da，窗口外匹配使用 override，无有效 base peak 则 spectrum-relative fallback。

本项目按作者要求不复制 LipidGate 的 active-pool dynamic weight redistribution；weights 固定。
无 Support theory 的值为 NA，状态 no_support_theory，贡献 0、不重新分配 20 分。
无结构 theory 的对应 pool/alias 为 NA，总分不可用。实际 precursor centroid 本身额外明确排除于独立 MS/MS 证据。
因此不能声称整个 LipidGate scoring pipeline 与本项目相同；只对齐上述数学函数与归一化机制。
参考版本、函数和配置记录在 `audit/lipidgate_scoring_reference.json`。

## Gates、重复计分和 annotation

Gates 仍按每个 required group 的 ceil(0.50*N_gate) 独立判断。最强单峰 q=1 不会绕过 1/4 的覆盖不足。
HG/LCB/diagnostic NL 可以 gate；supporting NL/common 不可以。
Cer species gate 依赖 secondary diagnostic NL，molecular 还需 LCB。
Gate 的 spectrum-relative minimum 未被结构质量的 non-precursor override 改写。

每个 exact theoretical m/z 只有一个计分 pool。primary > secondary > support，并保留多 label provenance。
同一个已匹配 centroid 不重复计入两个 structural pools 或 structural+support。
GM3 未映射到 HG/LCB 结构池的 diagnostic NL 保留解释，不算 Support。
回归测试覆盖了这些排他条件；现有最大配对数/最小 ppm 的 matcher 原文件未改动。

## 人工 shadow comparison

固定人工示例覆盖 GM3、SM、HexCer、LacCer、Cer；另测试 missing HG 与单链 So。
`multi_evidence_scoring_examples.csv` 保留各池类型、matched/theory、gate、三个 pool score、new_score、
legacy_match_score/legacy_total_score、bd428de_score/rank/level 及 annotation_level。
旧示例快照单独保留 `multi_evidence_scoring_examples_bd428de.csv`。

同一 SM precursor 的人工候选中，`SM(d18:1/16:0)` 的集中强 LCB 证据现在因 max(q) 排名第一，
`SM(d16:1/18:0)` 第二；bd428de 的平均质量评分顺序相反。legacy 顺序在该例中也为前者第一。
这是预先固定的人工输入呈现的公式影响，没有用真实结果调整 k/权重，也不是历史 rank 改变数量。
GM3 missing-HG 例仍不获得结构层级/新排名，即使 raw score 可计算。
极强 precursor 测试确认窗口外结构分数不受 precursor 强度压低；Support 计数不受强度质量函数影响。

## 未指定与未验证事项

单链 scoring policy 仍为 UNSPECIFIED；`so_fuc` 未追加的 LCB 不会被本层补造。
OCer 三链仍不套两链 FA 相减。GSL 首个/完整糖链锚点可观测性、Cer dehydration 选择性、
score calibration 和 threshold 外部验证仍未完成；本次不修改这些化学规则或用最终鉴定表训练它们。

## 验证与历史边界

完整 `python -m pytest -q`：**248 passed，30 个已有 pandas GroupBy.apply deprecation warnings**。
测试覆盖固定 mapping、q endpoints、max 与 50% 覆盖率分离、Support 三条饱和及无理论分支、
Support intensity independence、diagnostic NL 不重复计数、subclass-specific exact-mass priority、
non-precursor normalization/fallback/window boundary、precursor 本身排除、legacy 共存及全部既有 regression。
`git diff --check` 通过；对 `bd428de` 的逐文件 diff 确认 legacy core、glycan parser、fragment assignment、
legacy scoring、两个 production runners、rt_iup 和 rt_validation 均无修改。

没有重跑真实样品，没有覆盖历史 1,723 / 1,674 / 1,173 表，也没有声称这些历史鉴定已被新评分验证。
当前默认 batch pipeline 未启用新 evidence score 排序。
