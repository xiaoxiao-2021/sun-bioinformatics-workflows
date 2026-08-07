# sun-bioinformatics-workflows

这是一个面向博士阶段及长期科研积累的个人生物信息学科研平台，用于统一管理具体科研课题、可复用分析流程、通用脚本、软件环境记录和方法文档。

## 目录说明

| 目录 | 用途 |
|---|---|
| `projects/` | 存放具体科研课题。每个项目独立记录数据来源、分析目标、当前状态和使用方法。 |
| `workflows/` | 存放可跨项目复用、可配置和可重复运行的分析流程。 |
| `scripts/` | 存放 R、Python 和 Shell 通用工具函数及辅助脚本。 |
| `environments/` | 记录 Conda/Micromamba 等软件环境配置，不在仓库中保存已安装环境。 |
| `reference/` | 记录或组织基因组、注释和数据库等参考资源；大型资源文件不提交到 Git。 |
| `docs/` | 存放学习笔记、分析方法、标准操作流程和科研记录。 |
| `tools/` | 存放独立软件的使用说明、安装记录或本地工具入口；软件本体按许可和体积决定是否纳入版本管理。 |

## 使用原则

- 项目特异的数据和结果放在对应的 `projects/` 子目录中。
- 可复用逻辑从具体项目中抽离到 `workflows/` 或 `scripts/`。
- 软件依赖通过 `environments/` 中的配置文件记录，避免只依赖交互式环境。
- 原始测序数据、大型中间文件和本地数据库不提交到 Git。
- 每个项目和 workflow 的 README 应随研究进展持续更新。

## 科研项目数据与输出规范

项目内文件按语义归类，不按扩展名机械移动：

| 路径 | 内容 |
|---|---|
| `data/raw/` 或 `datasets/<dataset>/raw/` | 下载或接收的原始数据。 |
| `data/processed/` 或 `datasets/<dataset>/processed/` | 清洗、整理或格式转换后，仍会作为后续分析输入的数据。 |
| `metadata/` 或 `datasets/<dataset>/metadata/` | 样本信息、分组、临床信息和数据集登记信息。 |
| `scripts/` | 项目分析代码。 |
| `results/tables/` | 最终或面向汇报的结果表。 |
| `results/statistics/` | P 值、效应量、相关系数、差异分析等统计结果。 |
| `results/summary/` | 描述性统计和分析汇总。 |
| `figures/` | PDF、PNG、SVG 等全部图形输出。 |
| `objects/` | RDS、RData、h5ad 等中间分析对象。 |
| `docs/` | 分析记录、方法说明和结果解释。 |

`processed data` 不等于 `results`：即使某个统计结果会被下一步脚本读取，只要它表达的是检验、比较、排序或汇总结论，就仍应放在 `results/`。脚本应显式创建输出目录并写入项目相对路径，禁止依赖当前工作目录生成 `Rplots.pdf`，也禁止写入个人服务器的绝对路径。

本次目录重构只建立平台框架，不下载数据、不安装软件，也不创建或修改 Conda 环境。
