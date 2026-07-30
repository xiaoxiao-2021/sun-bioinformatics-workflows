# GSE182148 TAOK3 RNA m5C analysis

## Project purpose

基于GSE182148公开的RNA亚硫酸氢盐测序处理后数据，
检索TAOK3相关候选m5C位点，并比较Control与NSUN2 knockdown
条件下的位点甲基化水平。

## Research questions

1. TAOK3是否存在于公开的候选m5C位点表中？
2. TAOK3包含几个候选m5C位点？
3. NSUN2敲低后，相关位点的m5C水平是否下降？
4. 这些变化是否在多个胰腺癌细胞系中重复出现？

## Dataset

- GEO accession: GSE182148
- Cell lines:
  - BxPC-3
  - MiaPaCa-2
  - PANC-1
- Groups:
  - Control
  - NSUN2 knockdown

## Directory structure

- `data/raw/`: GEO下载的原始补充文件
- `data/processed/`: 整理后的TAOK3相关数据
- `scripts/`: 分步骤分析代码
- `results/tables/`: 输出表格
- `results/figures/`: 输出图片
- `docs/`: 数据理解和分析记录
- `GSE182148_TAOK3_m5C.qmd`: 完整分析报告

## Analysis workflow

1. 检查处理后数据结构
2. 确认样本与字段含义
3. 提取TAOK3相关位点
4. 计算或整理m5C水平
5. 比较Control与NSUN2 knockdown
6. 评估覆盖度、重复一致性及跨细胞系重复性
7. 输出表格和图片