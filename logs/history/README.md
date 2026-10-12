# 历史开发记录

本目录保留已有模型开发、配置搜索、公平比较、时间外推及来源核对的原始记录。XGBoost_F4_result.json、XGBoost_F4_trials.csv 和 continuation_round_01_search.log 对应正式配置的开发选择；原始 OOF 位于 data/training。

comparison.csv、seed_stability.csv 及 historical_evaluation / historical_validation_status 为相应历史模型与数据口径的记录。来源修订、样本支持及数据接触状态写在原记录中；这些记录不作为当前正式模型新的独立验证成绩。

原始 OOF、旧导出和阈值边界精度审计保留，读取策略统一为 round-trip 浮点。现存开发日志与本版本新生成的训练复现日志按时间分别存放；新增记录见 ../reproduction 和 ../validation。
