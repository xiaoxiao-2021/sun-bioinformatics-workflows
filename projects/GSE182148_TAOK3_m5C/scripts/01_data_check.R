library(readxl)
library(dplyr)
library(stringr)
library(readr)

# 读取Excel
file_path <- "data/raw/GSE182148_GEO_process_data.xlsx"

raw_dat <- read_excel(
  file_path,
  sheet = "NSUN2 m5C",
  col_names = FALSE,
  na = "ND"
)

# 15个样本
sample_names <- c(
  "BxPC3_CT1", "BxPC3_CT2", "BxPC3_CT3",
  "BxPC3_KD1", "BxPC3_KD2",
  "MiaPaCa2_CT1", "MiaPaCa2_CT2", "MiaPaCa2_CT3",
  "MiaPaCa2_KD1", "MiaPaCa2_KD2",
  "PANC1_CT1", "PANC1_CT2", "PANC1_CT3",
  "PANC1_KD1", "PANC1_KD2"
)

# 为每个样本生成coverage和methR列名
sample_columns <- as.vector(
  rbind(
    paste0(sample_names, "_coverage"),
    paste0(sample_names, "_methR")
  )
)

column_names <- c(
  "chromosome",
  "position",
  "strand",
  "site_id",
  "gene",
  sample_columns,
  "CTave",
  "KDave",
  "delta_methR"
)

# 删除前4行说明和表头
dat <- raw_dat[-c(1:4), ]

names(dat) <- column_names

# 整理数据类型
dat <- dat |>
  mutate(
    gene = str_trim(gene),
    across(
      -c(chromosome, strand, site_id, gene),
      as.numeric
    )
  )

# 查看数据
dim(dat)
names(dat)
head(dat)

# 提取TAOK3
taok3_dat <- dat |>
  filter(gene == "TAOK3")

taok3_dat |>
  select(
    chromosome,
    position,
    strand,
    site_id,
    gene,
    CTave,
    KDave,
    delta_methR
  )

# 保存提取结果
write_csv(
  taok3_dat,
  "data/processed/TAOK3_sites.csv"
)