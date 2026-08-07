# GSE182148 TAOK3 RNA m5C analysis

## 分析目标

基于GSE182148公开的RNA亚硫酸氢盐测序处理后数据，
检索TAOK3相关候选m5C位点，并比较Control与NSUN2 knockdown
条件下的位点甲基化水平。

## Research questions

1. TAOK3是否存在于公开的候选m5C位点表中？
2. TAOK3包含几个候选m5C位点？
3. NSUN2敲低后，相关位点的m5C水平是否下降？
4. 这些变化是否在多个胰腺癌细胞系中重复出现？

## 数据来源

- GEO accession: GSE182148
- Cell lines:
  - BxPC-3
  - MiaPaCa-2
  - PANC-1
- Groups:
  - Control
  - NSUN2 knockdown

## 当前状态

进行中。项目已建立分步骤 R 脚本和 Quarto 分析报告，现阶段代码覆盖数据检查、TAOK3 位点整理、汇总、可视化和结果表导出。最终生物学结论仍需结合数据质控和重复性验证后更新。

## 使用方法

在本项目根目录按顺序运行：

```bash
Rscript scripts/01_data_check.R
Rscript scripts/02_reshape_TAOK3.R
Rscript scripts/03_summary_TAOK3.R
Rscript scripts/04_visualize_TAOK3.R
Rscript scripts/05_export_result_table.R
```

也可使用 Quarto 渲染完整报告：

```bash
quarto render GSE182148_TAOK3_m5C.qmd
```

运行前请确认所需输入数据已按下述目录约定放置，并根据脚本中的参数说明检查文件名和字段。

## Directory structure

- `data/raw/`: GEO 下载的原始补充文件
- `data/processed/`: 清洗和重塑后、继续作为下游输入的 TAOK3 数据
- `metadata/`: 样本、分组及数据来源说明（需要时建立）
- `scripts/`: 分步骤分析代码
- `results/statistics/`: 位点和分组的统计比较结果
- `results/tables/`: 最终汇总结果表
- `figures/`: PDF、PNG 和 SVG 图形输出
- `objects/`: RDS、RData 或 h5ad 等中间对象（需要时建立）
- `docs/`: 数据理解、方法和结果解释记录
- `GSE182148_TAOK3_m5C.qmd`: 完整分析报告

其中 `TAOK3_sites.csv` 和 `TAOK3_sites_long.csv` 属于可继续分析的 processed data；`TAOK3_group_comparison.csv` 包含组间汇总与差值，属于 statistics result。

## Analysis workflow

1. 检查处理后数据结构
2. 确认样本与字段含义
3. 提取TAOK3相关位点
4. 计算或整理m5C水平
5. 比较Control与NSUN2 knockdown
6. 评估覆盖度、重复一致性及跨细胞系重复性
7. 输出表格和图片
