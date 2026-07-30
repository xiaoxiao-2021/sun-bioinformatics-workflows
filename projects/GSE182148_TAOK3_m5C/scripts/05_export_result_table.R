library(readr)
library(dplyr)

# 读取分组比较结果
comparison <- read_csv(
  "data/processed/TAOK3_group_comparison.csv",
  show_col_types = FALSE
)

# 重新补充准确的位点坐标
site_map <- read_csv(
  "data/processed/TAOK3_sites.csv",
  show_col_types = FALSE
) |>
  transmute(
    site_id,
    site = paste0("chr", chromosome, ":", position)
  ) |>
  distinct()

# 整理最终结果表
result_table <- comparison |>
  select(-any_of("site")) |>
  left_join(site_map, by = "site_id") |>
  transmute(
    gene,
    site_id,
    site,
    cell_line,

    control_valid_n = valid_replicates_Control,
    kd_valid_n = valid_replicates_NSUN2_KD,

    control_mean_methR = round(mean_methR_Control, 3),
    kd_mean_methR = round(mean_methR_NSUN2_KD, 3),
    delta_methR = round(delta_methR, 3),

    control_mean_coverage = round(mean_coverage_Control, 1),
    kd_mean_coverage = round(mean_coverage_NSUN2_KD, 1),

    evaluation = case_when(
      is.na(control_mean_methR) |
        is.na(kd_mean_methR) ~ "Not evaluable",
      TRUE ~ "Evaluable"
    ),

    change_direction = case_when(
      evaluation == "Not evaluable" ~ NA_character_,
      delta_methR > 0 ~ "Decreased after NSUN2 KD",
      delta_methR < 0 ~ "Increased after NSUN2 KD",
      delta_methR == 0 ~ "No change"
    )
  ) |>
  arrange(site_id, cell_line)

# 查看结果
print(result_table, n = Inf)

# 保存最终结果表
write_csv(
  result_table,
  "results/tables/TAOK3_m5C_result_table.csv"
)
