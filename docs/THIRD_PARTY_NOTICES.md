# 第三方资源与许可

自有代码、文档及正式模型采用根目录 MIT 许可。数据库记录、供应信息、第三方软件和原始研究材料分别遵守其来源条款；以下许可说明标明各自适用范围。查阅日期 2026-10-12。

| 资源 | 用途与版本 | 适用许可或获取规则 | 官方来源 |
|---|---|---|---|
| ChEMBL | ChEMBL 37，药理测量与交叉核对 | CC BY-SA 3.0；保留出处、版本和修改说明 | https://www.ebi.ac.uk/chembl/ |
| BindingDB | Q9NQS5 历史快照，测量核对 | BindingDB 自行整理记录为 CC BY 3.0；导入 ChEMBL 的记录按 CC BY-SA 3.0 | https://www.bindingdb.org/rwd/bind/info.jsp |
| IUPHAR/BPS GtoPdb | 2026.3，作用类型核对 | 数据库 ODbL；内容 CC BY-SA 4.0 | https://www.guidetopharmacology.org/download.jsp |
| PubChem / PubMed | 结构身份、文献与元数据核对 | 按 NCBI 使用政策及具体来源署名 | https://www.ncbi.nlm.nih.gov/home/about/policies/ |
| ACS Figshare 补充材料 | 测量来源核对与自行整理的事实记录 | 按具体材料页面许可；4c00359.s001 为 CC BY-NC 4.0；仓库不包含原始 SI 文件 | https://acs.figshare.com/articles/dataset/26134126 |
| ZINC / SmallWorld | In-Stock-25Q2 十万结构筛选 | 公开项目筛选结果；原始大规模快照不随仓库再分发 | https://wiki.docking.org/index.php/Terms_And_Conditions |
| RCSB PDB 9WIJ | 受体结构、原配体与对接准备 | PDB archive 数据 CC0；保留结构来源与引用 | https://www.rcsb.org/structure/9WIJ |
| 商品查询来源 | 候选身份及查询日期核对 | 提供链接与项目核对结果；原始网络响应及商品页面留在本地 | 见 data/evidence/05_supplier_evidence.csv |

各来源事实记录作为带来源的集合提供，原始数据库及研究材料的许可继续适用。不同来源的数据库权益与内容许可不由 MIT 覆盖；涉及 ACS 非商业材料的使用按其页面条件执行。需要原件时通过原来源独立获取。

## 软件依赖

| 软件 | 版本与用途 | 主要许可 | 来源 |
|---|---|---|---|
| RDKit | 2026.3.6，结构标准化与特征 | BSD-3-Clause | https://github.com/rdkit/rdkit |
| XGBoost | 3.4.1，正式分类模型 | Apache-2.0 | https://github.com/dmlc/xgboost |
| scikit-learn | 1.8.0，插补与评价 | BSD-3-Clause | https://github.com/scikit-learn/scikit-learn |
| NumPy / pandas / SciPy | 数值计算与数据处理，见 requirements.txt | 各项目 BSD 许可及所附第三方许可 | https://numpy.org/ 、https://pandas.pydata.org/ 、https://scipy.org/ |
| joblib / PyArrow | 模型与表格读写 | joblib BSD；Arrow Apache-2.0 | https://github.com/joblib/joblib 、https://arrow.apache.org/ |
| ADMET-AI | 历史预测 2.0.0，Chemprop 2.3.1、PyTorch 2.14.0 | 按各软件与预训练权重原许可 | https://github.com/swansonk14/admet_ai |
| AutoDock Vina / Meeko / molscrub | 历史对接 1.2.7 / 0.8.0 / 0.3.0 | 按各软件原许可 | https://github.com/ccsb-scripps/AutoDock-Vina |

历史 ADMET 与对接计算的输入、输出、版本、参数和签名保存在 data/evidence；默认流程复用这些已保存证据。仓库不打包第三方软件环境、ADMET 预训练权重、论文全文、原始 SI 或完整商品响应。

项目自有代码与模型的公开授权由项目负责人确认，MIT 许可与源文件同时交付。第三方来源引文及开发数据的用途信息保留在 data/sources.csv 和 data/training/audit。
