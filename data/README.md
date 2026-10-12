# 数据来源与复现输入

本仓库提供正式模型训练成员、结构特征、发表分组与时间划分、开发 OOF、来源测量审计及项目计算筛选结果。各记录的原始来源、版本、用途及适用许可见 sources.csv 和 [资源与许可说明](../docs/THIRD_PARTY_NOTICES.md)。

## 训练数据

training/training_structures.csv、training/training_members.csv 和 training/development.npz 对应相同的 400 个结构及顺序。测量端点、单位、作用类型、物种、来源 URL/DOI 和排除原因在 training/audit；发表组与滚动时间窗口在 training/splits。

正式成员由 training_members.csv 定义。测量审计保留 2,674 条重建记录，包括因结构、测量或标签原因排除的记录，便于复核标签形成过程。数据截止 2025-12-31。

## 十万结构筛选输入

原始库来自 ZINC SmallWorld In-Stock-25Q2 的确定性十万结构选择，记录日期 2026-09-22。官方获取入口为 https://sw.docking.org 和 https://zinc20.docking.org/tranches/home/。

ZINC 允许分享筛选结果；大规模原始数据库快照的再分发须取得书面许可。本仓库公开的是项目生成的结构评分结果，原始 CSV 快照留在本地，按 [官方条款](https://wiki.docking.org/index.php/Terms_And_Conditions) 使用。

完整筛选默认只读取 results/official_100k_scores.parquet 的结构与编号列，重新计算十万条分数。结构内容、顺序和数量与冻结输入逐项校验，保存的分数不作为新的推理输入。该方式安装依赖后可离线复算。

如需将公开结果中的结构单独整理为本地输入：

```powershell
python prepare_library.py --from-results
python predict.py --output results/reproduction
```

如已通过官方渠道独立取得与冻结清单完全一致的原始 CSV 快照：

```powershell
python prepare_library.py --input "local_data/public_purchasable_100k.csv"
```

导入程序先核对 SHA-256、十万成员及顺序，随后保存在本地。导入文件和运行生成的结构输入已列入 .gitignore。官方服务可能需要交互式访问；原库下载不是默认演示与完整复算的前置步骤。冻结清单哈希、原始来源和访问日期见 library/library_manifest.json。

## 下游证据与结构

evidence 的公开计算结果覆盖 300 条 ADMET/BBB、60 条对接及对应身份查询。按完整 InChIKey 核对结构，保留全部 60 条证据及 Top1000 覆盖情况。输入、输出、权重签名、参数和受体准备信息随相应 CSV/JSON/结构文件提供。

原始第三方网络响应、商品页面和中间缓存保留在本地；公开表提供来源链接、状态及本项目计算/核对结果。供应信息按记录日期理解，后续使用可通过表中链接重新核对。

## 数据处理与用途

原始结构经去盐、电荷处理和保留立体化学的标准化。训练数据按发表与结构连接键分组；未标注商业结构用于候选筛选。来源表与测量表保留去重、质量检查、标签排除及历史用途。排除标记只提供筛选所需状态，仓库不分发受限身份映射或盲法标签。

第三方记录保留原许可和署名；根目录 MIT 许可适用于本项目自有代码及模型。复核时同时查阅来源记录和具体数据许可。
