# 候选结果字段

主结果为 UTF-8 BOM 编码的 results.csv；数值保存 17 位有效数字，复核时使用 round-trip 浮点读取。机器可读定义见 [results.schema.json](results.schema.json)。

| 字段 | 含义 |
|---|---|
| candidate_id / rank | 候选编号与综合排序位置 |
| track / target | 赛道与人源 GPR84 靶点 |
| SMILES / full_inchikey | 标准化母体结构与完整结构身份 |
| ranking_score | 固定模型的未校准排序分数 |
| fixed_threshold_class | score ≥0.046076174825429916 时为 1 |
| applicability_domain | 训练相似度与描述符范围联合规则 |
| composite_score | 五分量固定权重综合分数 |
| BBB_probability / ADMET_component | 计算模型预测与 ADMET 综合分量 |
| DILI / hERG / AMES | 预测风险指标，较高值代表较高预测风险 |
| Vina_best_kcal_mol / Vina_seed_sd | Vina 最佳得分及三种子离散度，kcal/mol |
| structure_file / receptor_file | 仓库内相对路径，结构坐标单位 Å |
| supplier / product_id / product_url | 查询时的商品身份与来源 |
| model_version / run_version | 冻结模型版本与公开运行版本 |
| candidate_status / notes | 计算候选阶段与证据说明 |
| independent_acceptance | 冻结模型的独立评价标记，与候选细胞实验状态分别记录 |

数据中的空白表示缺失，缺失证据不用于补造药理或下游标签。模型分数、BBB/ADMET 和 Vina 指标分别按其计算模型解释；细胞功能与药效按对应实验体系和实测证据评价。
