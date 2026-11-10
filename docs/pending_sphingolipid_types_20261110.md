> 此文件记录前一检查点的暂缓原因。后续文献复核、LM1公式修正和DMS有限补入见 [本轮处理](sphingolipid_literature_review_20261110.md)。

# 八个暂缓启用 MS2 的历史类型

本表是只读核对，不新增别名、库谱、碎片、门控或权重。八项在当前
`config/sphingolipid_type_coverage.csv` 中的 `canonical_class` 均为空；
下表的“相关现有类型”表示可用于后续研究的相近骨架，**不是已建立的直接对应**。
名称继续保存在独立 MS1 参考目录中，不作为已确认 MS2 类型。

| 原名称 | 相关现有类型（均未直接映射） | 暂时不加 MS2 规则的具体原因 |
| --- | --- | --- |
| 1-O-acetyl-Cer | Cer 中的酯化/O-acyl Cer 规则 | 字面 acetyl 是 C2 乙酰基，旧 `OCer` 代码却枚举可变长链 O-acyl。名称、组成和旧实现含义未统一；现有长链酯化模板不能直接替代 1-O-acetyl 专属证据。也不因此新增短脂酰链例外。 |
| 1-O-acetyl-ω-linoleoyloxy-Cer | Cer 中的 ω/O-acyl Cer 相关模板 | 同时涉及 1-O-acetyl 和 ω-linoleoyloxy 两处修饰，且有上述 acetyl/acyl 歧义。普通酯化 Cer 的碎片不能自动证明两种修饰及其位置；不是含糖的 AHexCer 别名。 |
| 3-keto Cer | Cer | 3 位醇氧化为羰基与增加一处 C=C 均少两个 H，可产生同分异构组成。普通 Cer 的 LCB 失水对不能单独确认 3-keto；本轮未验证其专属诊断离子及门控。 |
| LM1 | 含唾液酸的 Hex3HexNAcNeuAcCer 组成与 GM1 组成族相关；旧表分子式却对应中性的 Gb4/GA1 组成族 | 旧名称与分子式冲突，不能根据名称套 GM1，也不能根据旧分子式直接报 LM1。需先核准糖组成、质量和结构来源。 |
| N,N,N-trimethylSo | SPB | 字面结构涉及季铵固有电荷，需核实分子式电荷及离子形式；不能把它直接当普通中性 SPB 的 [M+H]+。本轮未建立并验证对应离子形式和特征碎片。 |
| N,N-dimethylSo | SPB | 例如 N,N-dimethylSo(d18:1) 与普通 SPB(d20:1) 同分子式。共有失水峰不足以确认 N 上的两个甲基；本轮没有验证 N-二甲基专属诊断片段。 |
| N-methylSo | SPB | 例如 N-methylSo(d18:1) 与普通 SPB(d19:1) 同分子式。普通 SPB 的失水证据不能确认 N-甲基化位置；本轮未验证专属诊断片段。 |
| ketoSo | KetoSPB 相关，但现有仅启用 ketoSa(d18:0) | 当前 KetoSPB 仅一条饱和 3-ketosphinganine 参考，不能外推整套不饱和 ketoSo。旧 ketoSo(d18:1) 的 C18H35NO2 又与普通 SPB(d18:2) 同分子式；需要独立诊断证据。 |

## 两个可复核的例子

旧目录 `LM1(d34:1)` 为 **C60H110N2O23**，对应中性 Hex3HexNAc-Cer
组成。含一个 NeuAc 的 lacto/neolacto 单唾液酸糖鞘脂应再包含 NeuAc
残基 C11H17NO8；按同一骨架推算为 **C71H127N3O31**。
这是对旧名称/分子式冲突的核对，未将旧条目改写为任何新类别。
原始结构研究所述 nLM1 包含 NeuAc-Gal-GlcNAc-Gal-Glc-Cer，
支持“不能缺失 NeuAc”的检查；并不由此证明 LM1 与 GM1 糖连接相同。
[原始结构研究](https://academic.oup.com/femspd/article/45/2/171/646673)。

已启用的 ketoSa(d18:0) 参考，其正离子 300.2→270.3 与失去 CH2O 一致。
该饱和标准的结果不能直接作为所有不饱和 ketoSo 的验证。
[3KDS 标准实验](https://pmc.ncbi.nlm.nih.gov/articles/PMC5748491/)。

其余暂缓项不是“文献中不存在碎片”，而是本次迁移尚未取得足以支撑
独立规则的核实证据。不以普通失水峰兜底确认修饰位置，也不为了旧清单
条数补造可检索类型。当前 104 个已启用类型的规则、本轮八项状态和原库均不改。
