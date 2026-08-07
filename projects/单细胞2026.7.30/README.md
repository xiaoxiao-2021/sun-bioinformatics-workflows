# 单细胞项目（2026-07-30）

## 数据来源

待补充。请记录数据集编号、下载地址、物种、组织、实验平台、样本分组和数据使用限制。

## 分析目标

待补充。请明确主要科学问题、拟分析的细胞类型、比较组和预期输出。

## 当前状态

项目目录和 README 已建立，具体数据来源、分析设计、代码入口及研究进度尚待整理。

## 使用方法

待补充。后续加入分析代码时，请在此记录运行顺序、输入文件、软件环境、关键参数和结果目录。项目特异代码应保存在本项目内；可复用的单细胞流程应整理到 `workflows/single_cell/`。

## 目录约定

- `data/raw/` 或 `datasets/<dataset>/raw/`：原始单细胞数据。
- `data/processed/` 或 `datasets/<dataset>/processed/`：可继续用于下游分析的整理后输入数据。
- `metadata/`：样本、分组、临床和批次信息。
- `scripts/`：项目特异分析代码。
- `results/statistics/`、`results/tables/` 和 `results/summary/`：数值、统计和汇总结果。
- `figures/`：PDF、PNG、SVG 等图形。
- `objects/`：Seurat RDS、RData、h5ad 等中间对象。
- `docs/`：分析记录、方法和结果解释。

仅在实际产生对应内容时创建目录，避免空目录堆积。
