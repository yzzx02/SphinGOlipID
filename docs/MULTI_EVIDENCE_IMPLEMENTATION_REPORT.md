# Multi-evidence implementation report

基线：`542fa7f`。本次增加证据 metadata、集中 class registry、独立 gates/score、兼容 adapter 和人工 shadow audit。
未修改任何 fragment chemistry/mass、`542fa7f` 一对一 assignment、legacy score 方程、Top-N 或 RT 模型选择。
`ms2_legacy_core.py` 仅增加 source-capture decorators 和 LCB 来源记录调用，默认无 capture 时执行原生成函数。

## 作者在实施中确认的变更

1. 每个必需 HG/LCB/diagnostic NL 组分别至少匹配 **50%** 的合格理论片段，向上取整。
2. 普通 Cer 的一、二次 precursor dehydration 作为 **Cer 专属 diagnostic NL**。

这两个决定已进入集中配置/分类、回归测试和文档。Cer 的来源标注为 `author:Cer_dehydration`，
不把它包装成已验证的普适诊断特异性，也不改动原脱水质量。

## Class-specific functions 的分类

完整逐片段清单：`audit/multi_evidence_fragment_classification.csv`。

| 现有函数/路径 | Evidence type | 判定与限制 |
|---|---|---|
| SM | HG / NL / common | 184.0733 主要 HG；已有头基 losses 为 NL；86.0964、104.1066 保守 supporting |
| PE_cer | HG / NL | 现有 phosphoethanolamine ion 与两个已有 losses；小胺 NL 可评分但不独立 gate |
| PI_cer | HG / NL | 保留原 `[phosphate+H]+` 标签及其既有质量，未修化学命名 |
| PG_cer | HG / NL / common | 主 phosphoglycerol HG；小型 glycerol ion 为 common |
| CAEP / N_CAEP | HG / NL | 保留各自既有离子及中性丢失，包括原标签拼写 |
| CerP / S1P | HG / NL | H4PO4 HG；H3PO4 loss NL；带水组合 NL 解释，不进入默认分母 |
| GSL_fragments | HG / NL | NeuAc 两个现有 diagnostic ions 为 HG；全部可解释糖链/分支 losses 为 NL |
| FMC_1 / FMC_3 / FMC_5 | NL | Gal/acetylated-Gal/HOAc losses；water satellites 不默认评分 |
| EO_cer / EO_Glccer | NL | 既有固定 E 18:2 / Glc loss；不称独立 FA-ion confirmation |
| Lyso_SM | HG / NL / common | phosphocholine HG，小衍生离子 supporting；不借用 Cer 的 FA 推断 |
| Lyso_sulfo | HG / NL | 既有 sulfate-related ion/losses |
| Glu_So / Gb3_So | NL | 既有糖基损失；无新增链证据 |
| So | 无专属输出 | 函数为 pass，不把 generic dehydration 偷换为 So backbone 证据 |
| O_FA_cer | 无专属输出 | 函数为 pass；随后旧 OCer 分支另行枚举 |
| OCer | NL / UNSPECIFIED | 固定三链 acyl-loss 标签，排除 gate/score，等待 author review |
| almost_fuc / LCB dictionary | LCB | 精确字典标签，保留 m/d/t 类型与碳数/双键数身份 |
| learn_fuc generic dehydration | common 或 Cer 专属 NL | 仅 Cer 一、二次脱水按作者决定改证据类型 |
| 60.044 `[C2H5NO+H]+` | common | 支持评分，不能通过任何结构 gate |

剩余无法可靠分类或判定强度的标签保持 common/UNSPECIFIED 且不评分；不会靠未知标签数量抬升 gate。
OCer 虽可识别为 loss，具体链赋值/代表性仍不明确，单列 TODO。旧 glycan 表缺少来源时可识别 NL，
但代表锚点无法安全还原时仅解释，不自行将所有路径纳入评分。

## 类别 gates 与 scoring templates

完整表位于 [MULTI_EVIDENCE_SCORING.md](MULTI_EVIDENCE_SCORING.md) 与
`audit/multi_evidence_class_mapping.csv`，逐类列出 has_HG、species/molecular gate、模板、HG/LCB/NL/common 示例。

- SM、PE/PI/PG-Cer、CAEP、N-CAEP、Cer1P，以及 registry 中有 NeuAc HG 的 GM/GD 类：HG gate，HG-dominant。
- neutral HexCer/GlcCer/GalCer/LacCer/Gb3/Gb4/GA/type I B、FMC、EO-GlcCer：无独立 HG，diagnostic NL gate，NL-dominant glycosphingolipid。
- Cer 与 EO-Cer：无独立 HG，diagnostic NL gate，LCB-dominant。
- 上述两链类获得 Molecular species level 还须 LCB gate；common 不能补足任何 gate。
- So 使用单链 LCB/backbone gate，S1P/Lyso-SM/Lyso-sulfo 使用 HG + LCB，Glu-So/Gb3-So 使用 NL + LCB。
  其权重未被作者定义且当前生成路径缺少独立 LCB append，模板保留 UNSPECIFIED，不发明第四套权重。
- 三链 OCer 与未注册类别：规则/模板明确 UNSPECIFIED，不套 GM3 默认。

## Synthetic rank 与敏感性

`multi_evidence_shadow.py` 只调用当前原有 generation/matching/scoring 函数，使用人工 precursor、RT 与观测强度。
原始 MS1 pass 是明确标注的 synthetic declaration；没有打开 mzML、历史库或最终结果表。
“旧分数”由保持 `542fa7f` 语义的 production legacy scorer 计算，非手写近似公式。

人工同一 SM precursor 对照中，两候选都匹配 2/3 LCB。`SM(d18:1/16:0)` 的 raw LCB 强度更集中，
legacy rank 为 1；`SM(d16:1/18:0)` 的支持分布更均衡，按既定饱和函数的新 rank 为 1。
这是固定输入展示的排序差异；没有为获得该排序调整权重。

GM3 缺失 HG 的人工例子仍保留 raw `multi_evidence_score` 供诊断，但 confidence/rank 被 gate 阻止。
Cer 的一、二次脱水与 LCB 可按作者规则得到 molecular_species；So 不能凭现有 common 输出升级。
归一化最低强度敏感性表固定 50% 组门控和全部权重，仅展示 0/0.01/0.05/0.10 下门控变化，
不会据此选择“最好看”的 cutoff。原始 weighted score 不由 gate 强度阈值改变。

## 仍需作者验证的内容

1. Cer dehydration 的区分类别能力、相关性及在干扰物中的特异性。
2. GSL 首个/完整糖链丢失两个无水锚点的可观测性；这只是初始结构化代表集合，尚无实验校准。
3. SM/PG 小型子离子目前仅 supporting；若有更强的 subclass-specific 证据，可显式审阅后改 registry。
4. 单链类 scoring weights 与已存在数据中的真正 backbone evidence；本次不新增理论碎片规则。
5. OCer 固定三链枚举及旧标签的解释；不能用两链 total-minus-LCB 冒充三链分配。
6. 归一化门控强度下限、50% 覆盖率和 confidence score 的外部验证；目前分数不是概率或 FDR。
7. 未注册 subclasses 和不同 glycan isomer 必须提供明确规则/结构，不能依赖名称前缀猜测。

## 对历史结果的影响范围

历史 **1,723 / 1,674 / 1,173** 均未重跑、读取为新评分输入或覆盖。
当前 batch runner 的输出结构、legacy score 和排序保持原样；新层为显式 API/独立 shadow 路径。
未来若启用新 confidence ranking 或 gate 过滤，可能改变 score/rank/annotation level/最终保留集合；
本次人工例子已证明“可以改变”，但没有真实历史变化数量。
真实验证若另行授权，只能写新的 audit/validation 路径，不能直接覆写历史表。

## Verification

最终完整运行 `python -m pytest -q`：**219 passed，30 个已有 pandas GroupBy.apply deprecation warnings**。
新增分类、门控、模板、评分和 shadow 回归测试，覆盖作者 Cer 例外、各组 50%、源 metadata、
exact-mass priority、一对一计数、未匹配理论分母、全谱强度归一化、k 行为、single-chain/FA 限制、
RT 正交验证、关闭 adapter 的 legacy 字段保真和人工排名变化。

另以 Python AST 对 `542fa7f` 和当前 `ms2_legacy_core.py` 作静态比较：去除新增的 metadata import、
decorator 与 `record_lcb_origin()` 调用后，两者 AST 完全相同。由此核对了全部质量表、化学函数主体、
legacy matching/scoring 均未改变。`git diff --check` 通过。

人工 shadow 脚本已成功生成四张独立 audit 表；新增文件没有实验样品来源。代码处于当前审计分支，
新评分尚未作为默认 batch ranking/历史结果生成流程启用。
