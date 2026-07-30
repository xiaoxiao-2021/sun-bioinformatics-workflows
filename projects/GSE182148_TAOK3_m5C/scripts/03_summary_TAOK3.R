library(readr)
library(dplyr)
library(tidyr)

# 读取长格式数据
taok3_long <- read_csv(
  "data/processed/TAOK3_sites_long.csv",
  show_col_types = FALSE
) |>
  mutate(
    chromosome = as.character(chromosome),
    chromosome = if_else(
      grepl("^chr", chromosome),
      chromosome,
      paste0("chr", chromosome)
    ),
    site = paste0(chromosome, ":", position)
  )

# 防止某一组全部为NA时得到NaN
mean_or_na <- function(x) {
  if (all(is.na(x))) {
    NA_real_
  } else {
    mean(x, na.rm = TRUE)
  }
}

# 按位点、细胞系和分组汇总
group_summary <- taok3_long |>
  group_by(
    gene,
    site_id,
    site,
    cell_line,
    group
  ) |>
  summarise(
    total_replicates = n(),
    valid_replicates = sum(!is.na(methR)),
    mean_methR = mean_or_na(methR),
    mean_coverage = mean_or_na(coverage),
    .groups = "drop"
  )

# 将Control和KD放到同一行，便于比较
site_comparison <- group_summary |>
  pivot_wider(
    names_from = group,
    values_from = c(
      total_replicates,
      valid_replicates,
      mean_methR,
      mean_coverage
    )
  ) |>
  mutate(
    delta_methR = mean_methR_Control -
      mean_methR_NSUN2_KD
  ) |>
  arrange(site_id, cell_line)

# 查看结果
print(group_summary, n = Inf)
print(site_comparison, n = Inf)

# 保存结果
write_csv(
  site_comparison,
  "data/processed/TAOK3_group_comparison.csv"
)
