# Projects

本目录用于保存具体科研课题。

每一个子目录代表一个相对独立的科研项目、课题或科学问题。

## 一个项目建议包含

project_name/
├── datasets/
├── metadata/
├── scripts/
├── results/
├── objects/
├── docs/
└── README.md

## 各目录含义

datasets/
保存该项目使用的数据。

推荐：
- raw/：原始数据，不修改
- processed/：清洗或整理后的分析输入

metadata/
保存样本信息、分组信息、临床信息、数据集说明等。

scripts/
保存该项目专属的分析代码。

results/
保存统计分析结果和最终输出：
- figures/
- tables/
- logs/

objects/
保存 RDS、RData、h5ad 等分析对象。

docs/
保存分析记录、方法说明和结果解释。

README.md
记录项目目标、数据来源、分析流程、当前状态和运行方法。

## 不应该放在这里

- Conda 实际环境
- 通用 workflow
- 与项目无关的工具
- GitHub 不适合保存的大型原始数据