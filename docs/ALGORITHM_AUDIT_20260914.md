# 投稿前算法实现审计与锁定（2026-09-14）

基线：GitHub `yzzx02/SphinGOlipID` 的 `main`，提交 `3e0b33f374c8f8e5abee7c983bbd74efd74045c8`。修正位于独立分支 `codex/algorithm-submission-lock`。本报告锁定已实现行为及已知限制，不宣称所有 manuscript claim 都已实现。

本次未运行任何真实样品分析，未重新计算、筛选或替换历史 1,723 条 annotation 或 1,173 条 RT 结果。理论 formula library、class-specific fragmentation tables、legacy fragment generation/matching/scoring/Top-N、历史结果文件均未修改；没有发布、License 或 Citation 工作。历史库只抽查糖链编码列作为测试模板，未生成或统计理论分子库。

## 1. 明确 bug 与影响矩阵

| 位置 | 修改前逻辑 | 修改后逻辑 | 是否可能影响已有结果 |
|---|---|---|---|
| `result_cleaning.deduplicate_by_rt_score_intensity` | 按 RT 遍历，只比较后续行；高分早峰和低分晚峰可能同时保留 | 每组按 score 降序、intensity 降序、RT 升序选择；与全部已保留 RT 比较，差值 ≤ tolerance 时不再保留；完全同分同强度同 RT 用行内容决胜，不依赖输入索引 | **会影响重新调用此 cleaning 工具得到的保留集合、条数及顺序**。没有对历史结果调用该工具 |
| `scoring._sequence_length` / `add_standard_score_columns` | list/tuple/set 用长度，字符串通常计为 1；ndarray 可能触发 ambiguous truth value | `parse_fragment_mz_values` 解析真实数值条目，再计数 | 修正标准化 `matched_fragment_count`；对已有数值 list 的正常计数不变。legacy 中文分数、强度与 Top-N 不变 |
| `scoring.match_fragments` | 返回全部 observed × theoretical 匹配，多个近邻 centroid 会重复贡献同一理论峰 | 按候选及 theoretical m/z 去重，并保留最大 observed intensity；使用 production 相同的 inclusive ppm bounds | **normalized API 的计数和强度可能降低**，ppm 浮点边界归属可能修正。两个 production runner 不调用此 API，legacy final result 不受此项修改影响 |
| `rt_iup` model summaries | 中文摘要，没有稳定英文 provenance 字段 | 追加 `fit_strategy`, `fit_type`, `r2`, `representative_point_count`, `inlier_count`, `x_min`, `x_max` | 只新增模型摘要字段，中文字段和数值决策不变；未来导出模型表会多这些列，不改历史文件 |
| config / CLI / 两条 MS2 workflow | 20/10 ppm 默认分别散落于各入口 | 在 `config.py` 集中定义两个历史默认；配置值一路传到实际边界，日志记录 effective fragment_ppm | 默认数值不变；新增日志。没有根据文稿猜测最终样品实际使用哪一个值 |

字符串解析支持 list、tuple、一维 ndarray、数值 scalar、Python-list/tuple 字符串、NumPy array display、逗号/分号/空白分隔数值。使用 `ast.literal_eval`，不执行字符串代码。None、NaN、空字符串、非有限条目不计数；无效的非数值内容明确抛出 `ValueError`，不再假定它是一个碎片。重复 m/z 仍按列表条目保留计数，因为 legacy 记录的是理论匹配条目。已存在的英文 alias 不被强制覆盖。

## 2. 去重修复前后的具体差异

- A：RT 10.00，score 0.90；B：RT 10.05，score 0.70；tolerance 0.1。旧逻辑可能同时保留 A、B；新逻辑保留 A。
- 若 B 的 score 更高，则保留 B；同分按总碎片强度，强度相同再按更早 RT。
- RT 相差大于 tolerance 的候选保留为不同 chromatographic feature，不把所有可能 isomer 合并。
- 10.00、10.08、10.16 形成接近链时，不做传递性整簇合并。若两端证据优于中间，允许保留相距 0.16 的两端；若中间最强，保留中间。最终同组任意两个已保留、已知 RT 的差都大于 tolerance。
- 缺失 RT 无法建立邻近关系，独立保留；缺失 group key 的处理沿用 pandas groupby 的原默认行为。负 tolerance 明确报错。空输入保留列结构。

注意：新旧集合并非必然是单向删减。旧逻辑可能删除 A，因为后续 B 更强；但 B 又被 C 删除。贪心逻辑可能保留最强 C 和与 C 距离足够远的 A。因此不能提前承诺历史 1,723 条保持相同条数，也不能给出未经验证的变动数字。

在本次基线源码中，`run_batch()` 和 `run_targeted_mzml_batch()` 都不自动调用此 cleaning 函数。统一 config 的 `enable_deduplication` 并不代表该步骤已经接入这两个 runner。已有独立工作流 `scripts/prepare_missing_fraction_workbooks.py::deduplicate` 是按 feature_id + lipid_name 合并重复谱图的另一套逻辑，本次不修改。

## 3. RT 实际算法与 provenance

普通 `rt_iup.fit_line()` 路径：

1. 提取同细类、同不饱和度的原始 carbon–RT 点，对同 carbon number 取 median RT。
2. 对 Linear / Quadratic 分别以 RANSAC 识别初始 inlier mask。Linear 的 RANSAC 最小抽样为 2，但有效拟合至少需要 3 个不同碳数；Quadratic 最小抽样为 3，有效拟合至少需要 4 个不同碳数。
3. 只在当前 inliers 上做 OLS，按绝对 residual ≤ `rt_window_min` 继续剔除；不重新加入已剔除点。内点集不变时停止，原有迭代上限为 5。
4. 在最后保留集合上创建新的 OLS estimator 并重新拟合，参数不直接使用 RANSAC estimator。迭代存在上限，因此“stable”不能解释为无限次迭代收敛的数学保证。
5. 在最终内点上计算 R²，要求 R² ≥ 0.99、RT 随 carbon number 的正向趋势；二次曲线另有下降尾部检查。
6. 若两个候选都有效，Quadratic 相对 Linear 的 R² 改善小于 0.002 时选 Linear；达到 0.002 时允许选 Quadratic。若 Linear 无有效候选而 Quadratic 有效，则选后者。
7. ECN 对每条曲线独立评价；跨不饱和度的平行性在 ECN 路径中仅记录审计，不自动删除单独有效的 ECN 曲线。随后可运行 IUP 排序、点级检查及原有 rescue。

默认值保持：`r2_threshold=0.99`、`rt_window_min=0.20`、`quadratic_min_distinct_x=4`、`quadratic_min_r2_gain=0.002`。没有新增“独立 iterative OLS 与 RANSAC 竞争”的模型家族。

新增摘要 `fit_strategy="RANSAC-seeded OLS refit"` 表示拟合阶段。IUP 后续可能调整截距，已有 `IUP原始参数`、位移记录等字段用于描述这一验证操作；不能把位移后的参数解释为又一次无约束 OLS。

**必须明确的现存差异：**`build_rescue_line()` 的模型排序优先考虑平行性、相对斜率差，再比较 R²、内点数和模型复杂度，不使用普通 `fit_line()` 的 0.002 改善门槛。已有 multi-start IUP 候选选择也不是普通两模型比较。新增测试锁定了 rescue 在改善小于 0.002 时仍可能选择 Quadratic 的现状。本次没有改变这条影响历史 RT 结果的路径；manuscript 不应将普通拟合门槛泛化到全部 rescue 决策。若要求统一，需单独评估最终保留集合的影响。

另一个 workflow gap：targeted batch 的 `apply_rt_fit_filter()` 实际调用较早的 `rt_validation.fit_ransac_models()`，不是 `rt_iup.py` 的完整 ECN/IUP 步骤；legacy text runner 本身也不会自动执行完整 RT/IUP。文章必须说明最终数据实际调用了哪条路径。

## 4. Glycan sequence 与 branch 审计

### production 调用链

- legacy text：MS1 Excel 的 `structure`、`classy` → `core.isomer_search()` 选出的候选 → `run_one_file()` → `learn_fuc(classy, name, structure, ...)` → `GSL_fragments(structure, classy, M, name)`。
- targeted mzML：`_read_legacy_or_normalized_ms1_library()` → `search_ms1_candidates()` → `generate_legacy_fragments_for_candidates()` → 同一个 `learn_fuc()` / `GSL_fragments()`。
- `classy` 以字母开头时选择已有 subclass function；非字母类别且有 structure 时进入 glycan path。
- `GSL_fragments` 对 `structure.split(" ")` 得到的 residue-loss tokens 累加质量，输出 sequential precursor-loss fragment 及脱水伴随峰；遇到 NeuAc 另加 292.10、274.09 诊断峰。所有 residue masses 原样保留，包括 NeuAc 291.09 和 HexNAc 203.08 的 legacy 精度。

### 已验证的真实编码

历史 `MS1 DB_new  3.0.xlsx`（2025-07，Sheet1）中：

| 名称模板 | structure（引号内前导空格必须保留） | classy |
|---|---|---|
| LacCer(d18:0) | `" -Gal -Glc"` | 空 |
| GM3(d18:0) | `" -NeuAc -Gal -Glc"` | 空 |
| GM1(d18:0) | `" -Gal -GalNAc -NeuAc -Gal -Glc"` | `"0,3"` |

原文件 SHA256：`7719a35322067ceac1bca20447d1db2b3c092374c363c4975cab53c670e9c88e`。仅将这些既有编码写成独立测试模板，使用人工 precursor m/z、RT 和小型候选；不读取该库运行鉴定。

`MS2 Match (1).ipynb`（2026-06）的 GSL 函数与 production 相同；notebook SHA256：`7ef7bc565b5cfa98a4a3f7a4942865da8f6f9848968fd75a8be809d3e1beb26e`。这些历史源文件不属于 main 的已跟踪内容，报告只记录来源和哈希，不上传整个库或 notebook。

逗号分隔 `classy` 在 production 中解释为 split 后 token 的整数位置。GM1 的位置 3 选择 NeuAc，产生 branch-first losses。位置 0 对应前导空 token，原有 `[:y-1]` 会重新生成部分 sequential fragments；随后两个 production workflow 都在同 annotation 内按理论 m/z 去重。已对原始重复输出及生产路径实际可见的分支 m/z 写回归测试。

例如人工 precursor 1500.0，GM1 branch-first 额外得到 1208.91（失 NeuAc）及 1046.8572（再失 Gal），并各有脱水峰。简单 GlcCer/LacCer、GM3、GM1 及 d18:1 LCB 诊断峰均有数值 oracle 测试。

### claim 的准确边界

**已有有限的 branch-position-aware 处理，不能说 production 只支持线性糖链；但 `#` branch-aware glycan parsing is currently missing。**`#` token 会触发未知 residue/解析错误，并未被编码为糖树。当前实现也不是任意嵌套糖树、连接位点或立体化学解析器。

targeted MassHunter CSV fallback 通过名称 `_derive_legacy_structure()` 推导一条序列，`_derive_legacy_classy()` 对这些糖脂返回空字符串，因此没有复制历史 GM1 的位置元数据。该 fallback 的 GM1 等类别不能宣称运行了上述分支路径。新增测试明确暴露该 gap；本次没有按脂质名称猜测分支或悄悄改动类别映射。

本机另一历史 notebook `贪婪算法的数据库生成.ipynb` 包含括号形式的实验性 parser 示例，但不属于当前 production 调用链，且不等同于 `#` 语法。没有将这些不同实现直接拼接到 legacy。由于生产已有经过原始数据佐证的整数位置编码，本次选择锁定它，而非发明或替换语法。

## 5. Fragment matching 重复计数审计

`core.query()` 会返回窗口内全部理论行；`process_sheet()` 不执行一对一分配。`ms2_pipeline._finalize_results()` 和 targeted `finalize_spectrum_matches()` 都只在每个 annotation 内按理论 m/z 选择最大强度 centroid。

| 人工输入 | legacy / targeted 最终行为 |
|---|---|
| observed 100.0000，强度 50；theory 100.0001、100.0002 | count=2，matched intensity=100；两条匹配记录可以引用同一 observed m/z |
| observed 100.0000/100.0005，强度 50/80；theory 100.0001 | count=1，matched intensity=80；一个理论峰不因多个 centroid 重复加分 |
| 两条理论行 m/z 完全相同 | 同一 annotation 中只计一次 |

因此存在同一 observed peak 重复贡献给不同近邻理论峰的风险：单一观测峰可能满足 `min_matched_fragments=2`，并影响 coverage、强度加权分数及排名。其 intensity ratio 也不能无条件解释为独立观测峰总强度的比例。本次遵照要求**未将 legacy 改为一对一分配**；回归测试保留了这些数值。

修正后的 normalized API 与上述最终去重方式一致：只消除重复理论 m/z / 多 observed centroid，保留 distinct-theory one-to-many；多候选时使用 `anno`/`lipid_name`/`candidate_id` 或显式 `candidate_col` 区分身份。比较语义以相同输入理论表为前提：legacy text 建表阶段会四舍五入理论 m/z 到 4 位，targeted 保留生成值的精度，这个已有前处理差异未被擅自统一。

**是否参与 manuscript 最终结果：**这类 count/intensity 进入两个实际 production 的候选分数，具备影响最终候选的直接路径。main 的六组分整理脚本也读取这些已保存的分数和 fragment counts；但 main 不包含可核验的完整“最终 1,723 条”生成命令与参数链。本次没有逐谱重算，因此无法确认最终 1,723 中哪些行实际发生 one-to-many，不能编造受影响条数或声称风险不存在。

## 6. ppm 配置来源及不变项

`config.py` 统一持有 `LEGACY_FRAGMENT_PPM=20.0`、`TARGETED_FRAGMENT_PPM=10.0`。统一 config、旧 `MS2PipelineConfig`、两个 CLI 均引用对应常量，显式传入值优先。text 建表用 config 的 fraction 生成上下限；targeted 匹配同样使用 config 的 fraction。日志记录 `Effective fragment_ppm=...`，单 spectrum 调用也有 debug provenance。

新增 tests 从 config/CLI 一直验证到实际 fragment generation、matching 及 final count：15 ppm 偏差人工峰在 config=10 时排除、config=20 时命中，并测试窗口端点和强度 cutoff。不仅检查 config 属性。未改变 MS1 ppm、RT/IUP tolerance、legacy min score、Top-N 默认或最终论文 tolerance 的取值。

## 7. Regression tests 与验证结果

命令：`python -m pytest -q`（独立工作树，pytest 的 pythonpath 指向本树 `src`）。

- 修改前基线：**32 passed**。
- 修改后完整 suite：**113 passed, 22 warnings**，无失败、无跳过。
- warnings 均为现有 legacy pandas `GroupBy.apply` 的未来兼容性提示；本次不通过改变分组行为消除它们。
- 环境：Python 3.13.2，NumPy 2.2.3，pandas 2.3.3，scikit-learn 1.6.1。
- 新增六个测试模块：`test_result_cleaning.py`、`test_scoring.py`、`test_glycan_fragments.py`、`test_rt_model_regression.py`、`test_fragment_matching_regression.py`、`test_fragment_ppm_config.py`。
- 测试覆盖去重四个必需 case、行顺序不变性、强度/RT/full-precision metadata 决胜、列表及 Excel 字符串计数、LCB/GSL 数值、实际 branch metadata 传递、两个 production matcher 的 one-to-many 和 centroid 去重、legacy score/Top-N 排序、10/20 ppm、真实 RANSAC 去 outlier 后创建新的最终 OLS、median representatives、Linear/Quadratic 门槛和正趋势、R² 拒绝、summary aliases。

所有新增测试都是人工数据或历史编码模板下的人工质量，不访问原始 `.d` / mzML 样品或历史结果文件。已有 RT 形状 fixture 也只是测试代码内的小型表，不运行外部样品。测试只在 pytest 临时目录写合成文件；未渲染 Excel/Word。`ms2_legacy_core.py`、`formula_library.py` 和 `supplementary_data/` 对基线无 diff。

## 8. 是否建议重跑最终 QC 1,723 annotations

**建议安排一次独立、可回滚的 QC 复核；本次未执行，也未自动替换最终表。**先锁定最终使用的 source commit、理论结构元数据版本、实际 fragment_ppm、排序与后处理脚本、RT 路径，再决定复核范围：

1. 若最终整理调用本次修复的 RT near-duplicate cleaning，必须对比修复前后的保留集合；变化来源明确是 cleaning，而非 MS2 fragment rules。
2. 若最终处理中读入字符串形式的 `实际mz` 再创建标准化 count，需要核对修正 count 是否被后续过滤/representative selection 使用。若 count 仅作展示，则不应因此改写 legacy ranking。
3. 对 one-to-many 峰碰撞、GM1 等分支元数据缺失及 RT rescue 的模型选择范围做定向复核；这些是已有风险，不是本次悄悄修复后的新鉴定结果。
4. 在没有完整最终 provenance 的情况下，不承诺精确复现 1,723 / 1,173，也不把测试通过等同于真实数据无变化。若后续统一 branch metadata、匹配分配或 rescue 规则，应单独形成可审阅的结果差异清单。

本分支锁定的是“明确 bug 修复 + legacy 行为回归 + 已知 gap”，不是对上述未实现 claim 或真实样品结果不变的无条件认证。
