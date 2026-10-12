# 复现流程与结果核对

## 运行顺序

在仓库根目录使用 Python 3.12 配置 requirements.txt。随后按需要执行：

```powershell
python train.py --check-only
python screen.py --input data/examples/smiles.csv --output results/example_scores.csv
python demo.py --output results/demo
python train.py --reproduce --output results/training_reproduction
python predict.py --output results/reproduction --workers 4
python verify_submission.py --results-dir results/reproduction --output results/reproduction_verification --recompute-development
```

完整筛选读取十万结构并重新生成分数；演示使用 60 条真实证据结构和冻结 Top1000 成员展示同一证据整合过程。固定参数训练产生单独的复现模型和真实回执，不覆盖正式权重。

## 一致性标准

- 正式权重、关键化学代码和 400 个成员通过冻结哈希校验。
- 输入库的十万结构、编号和行顺序通过独立语义哈希校验；保存的筛选分数不供新推理使用。
- 全库评分、训练相似度、适用域、固定阈值分类和候选资格的列哈希与参考记录逐项一致。
- Top5 身份、顺序、模型分数精确一致，综合分数绝对差 ≤1e-12。
- 原始 OOF 按 float_precision=round_trip 读取，开发指标及两种各 2,000 次 bootstrap 可重算比对。

## 输入和输出规则

screen.py 支持 `SMILES`、`smiles`、`isomeric_smiles`，依次选择第一个存在的字段。空值、无效结构和重复行保留输入顺序；无效结构标记 invalid_SMILES，其余行继续处理。盐形式标准化后保持结构身份，可辨认立体化学予以保留。

分数等于固定阈值归阳性。输出目录为仓库内 results/、logs/validation/ 或仓库外独立目录；代码、配置、模型、输入 CSV 和不可变清单受到覆盖保护。版本或关键哈希不一致时终止并提供具体原因。

## 记录与维护

每次运行保存环境、参数、耗时、输入与输出哈希和证据覆盖。源文件与结果的交付快照为 checksums.sha256；运行时只读输入的校验清单为 IMMUTABLE_SHA256.json。结果和新日志会随复现生成，属于可变输出。

正式模型版本、阈值和排名权重保持冻结。文档或公开代码维护后可运行 `python tools/release_manifest.py` 更新发布清单；科学模型的关键锁仍独立校验。实际验收环境与各步骤测量见 [复现记录](REPRODUCTION_RECORD.md)。
