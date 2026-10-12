# 智靶药研 GPR84 小分子 AI 虚拟筛选

面向缺血性脑卒中相关 GPR84 靶向研究，本项目将药理数据整理、机器学习排序、理化与结构警示过滤、ADMET/BBB 预测和分子对接证据整合为可复核的计算筛选流程，帮助确定后续实验的候选优先级。

本仓库开放正式模型、训练与推理代码、400 个训练结构、数据划分与测量来源、十万结构的计算筛选结果、候选三维结构及结果复核程序。公开版本 **v1.0.0**；正式模型 **GPR84-ML-v3-official-1**。

## 项目方法与成果

- **有来源的药理数据**：按人源靶点、作用类型、测量单位和结构身份整理标签；关联论文与结构连接键分组，保存去重、冲突处理和划分记录。
- **可检查的模型输入**：XGBoost 使用 MACCS 167 位指纹与 10 项理化描述符，共 177 维；训练成员、插补规则、参数、随机种子和权重均可校验。
- **分层筛选与证据整合**：十万结构评分后应用理化及 PAINS/BRENK 过滤，再结合 ADMET/BBB、三种子对接、新颖性、商品身份和结构多样性选择计算候选。
- **可复核的结果**：提供全库评分、Top1000、60 条对接证据的完整排序、5 个候选及结构文件，并支持逐项数值比对。

```mermaid
flowchart LR
    A[药理数据与来源核对] --> B[结构标准化与分组评价]
    B --> C[固定参数训练与模型冻结]
    C --> D[十万结构重新评分]
    D --> E[理化与结构警示过滤]
    E --> F[ADMET 对接与身份证据整合]
    F --> G[多样性选择与候选复核]
```

| 内容 | 结果与口径 |
|---|---|
| 训练结构 | 400 个，其中阳性 326、阴性 74 |
| 发表连通组交叉验证 ROC-AUC | 0.937635 |
| 同一评价 MCC | 0.754376 |
| 同一评价 Balanced Accuracy | 0.821257 |
| 滚动时间外推平均 ROC-AUC | 0.739567 |
| 全库计算评分 | 100,000 条结构 |
| 计算筛选候选 | 5 个，来自合格 Top1000 与完整下游证据的交集 |

指标为开发阶段的模型选择评价。完整评价方法、置信区间、适用范围和下游证据覆盖见 [Model Card](docs/MODEL_CARD.md)；公平比较的原始记录见 [历史开发记录](logs/history/README.md)。

## 环境安装

使用 **Windows、Python 3.12、CPU**。建议至少 8 GB 可用内存与 2 GB 工作空间；完整筛选默认 4 个进程，可用 `--workers 1` 降低并发。核心流程在依赖安装完成后可离线运行。

在仓库根目录打开 PowerShell：

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe train.py --check-only
```

核心依赖锁定为 NumPy 2.5.3、pandas 2.3.3、SciPy 1.18.1、scikit-learn 1.8.0、RDKit 2026.3.6、XGBoost 3.4.1、joblib 1.6.0；完整列表见 [requirements.txt](requirements.txt)。

## 快速演示

```powershell
.\.venv\Scripts\python.exe demo.py --output results/demo
.\.venv\Scripts\python.exe screen.py --input data/examples/smiles.csv --output results/example_scores.csv
```

演示重新计算 60 条真实证据结构的模型分数，结合公开全库结果中的 Top1000 成员执行证据整合和多样性选择。它展示与完整流程相同的决策规则，通常约一分钟；真实测量耗时见 [复现记录](docs/REPRODUCTION_RECORD.md)。

## 完整筛选与结果复核

```powershell
.\.venv\Scripts\python.exe predict.py --output results/reproduction --workers 4
.\.venv\Scripts\python.exe verify_submission.py --results-dir results/reproduction --output results/reproduction_verification --recompute-development
```

默认从本项目公开筛选结果的结构与编号列读取十万结构，**重新计算全部模型分数**，再执行过滤、证据融合和候选选择。输入结构与顺序须通过冻结身份校验。导入原始快照和准备独立输入的方式见 [数据说明](data/README.md)。

主结果为 `results/reproduction/results.csv`，包含候选编号、SMILES、模型与运行版本、预测指标、综合分数及结构路径。完整分数保存为 Parquet；字段、单位与缺失规则见 [结果字段说明](docs/RESULT_FIELDS.md)。

## 固定参数训练复现

```powershell
.\.venv\Scripts\python.exe train.py --reproduce --output results/training_reproduction
```

使用同一 400 个成员、冻结特征、训练中位数插补、500 轮及种子 42 生成独立复现模型与训练分数比较，正式权重保持原样。训练、分组评价和结果核对步骤见 [复现指南](docs/REPRODUCIBILITY.md)。

## 仓库导航

| 位置 | 内容 |
|---|---|
| `src/` | 结构评分、证据整合和结果生成 |
| `models/` | 正式权重、冻结化学处理代码与模型锁 |
| `data/training/` | 400 个成员、特征、来源测量、划分和 OOF |
| `data/evidence/` | 300 条 ADMET、60 条对接、来源链接与候选结构 |
| `results/` | 十万评分、Top1000、候选清单和证据覆盖 |
| `logs/` | 历史开发与本版本实际复现记录 |
| `docs/` | 模型、数据、字段、复现和共享说明 |

## 共享与引用

自有代码和正式模型采用 [MIT 许可](LICENSE)。第三方数据、资源与软件的来源和适用许可见 [资源与许可说明](docs/THIRD_PARTY_NOTICES.md)；ZINC 原始快照通过官方渠道独立获取，本仓库公开项目计算筛选结果。

引用方式见 [CITATION.cff](CITATION.cff)，版本变更见 [CHANGELOG.md](CHANGELOG.md)，使用建议与贡献方式见 [CONTRIBUTING.md](CONTRIBUTING.md)。欢迎通过 GitHub Issues 报告复现问题或提交方法改进。
