library(readr)
library(dplyr)
library(tidyr)

# 读取第一步提取的TAOK3数据
taok3_dat <- read_csv(
  "data/processed/TAOK3_sites.csv",
  show_col_types = FALSE
)

# 宽格式转成长格式
taok3_long <- taok3_dat |>
  pivot_longer(
  cols = matches(
    "^(BxPC3|MiaPaCa2|PANC1)_(CT|KD)\\d+_(coverage|methR)$"
  ),
  names_to = c(
    "cell_line",
    "group",
    "replicate",
    ".value"
  ),
  names_pattern =
    "^(BxPC3|MiaPaCa2|PANC1)_(CT|KD)(\\d+)_(coverage|methR)$"
) |>
  mutate(
    group = recode(
      group,
      CT = "Control",
      KD = "NSUN2_KD"
    ),
    cell_line = recode(
      cell_line,
      BxPC3 = "BxPC-3",
      MiaPaCa2 = "MiaPaCa-2",
      PANC1 = "PANC-1"
    ),
    replicate = as.integer(replicate),
    site = paste(chromosome, position, sep = ":")
  ) |>
  select(
    gene,
    site_id,
    site,
    chromosome,
    position,
    strand,
    cell_line,
    group,
    replicate,
    coverage,
    methR
  ) |>
  arrange(site_id, cell_line, group, replicate)

# 基础检查
dim(taok3_long)
head(taok3_long)

# 检查各细胞系和分组的样本数量
taok3_long |>
  count(cell_line, group)

# 查看各位点的样本数据
taok3_long |>
  select(
    site_id,
    cell_line,
    group,
    replicate,
    coverage,
    methR
  ) |>
  print(n = 45)

# 保存可继续用于后续分析的长格式数据
processed_dir <- "data/processed"
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(
  taok3_long,
  file.path(processed_dir, "TAOK3_sites_long.csv")
)
