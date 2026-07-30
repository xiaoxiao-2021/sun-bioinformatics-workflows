library(readr)
library(dplyr)
library(ggplot2)

# 读取样本级长格式数据
plot_dat <- read_csv(
  "data/processed/TAOK3_sites_long.csv",
  show_col_types = FALSE
) |>
  filter(!is.na(methR))

# 绘制每个样本的甲基化率
p <- ggplot(
  plot_dat,
  aes(
    x = group,
    y = methR
  )
) +
  geom_jitter(
    width = 0.08,
    height = 0,
    size = 2
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 4
  ) +
  stat_summary(
    fun = mean,
    geom = "line",
    aes(group = 1)
  ) +
  facet_grid(
    site_id ~ cell_line
  ) +
  labs(
    x = NULL,
    y = "m5C methylation rate",
    title = "TAOK3 m5C levels after NSUN2 knockdown"
  ) +
  theme_classic()

print(p)

# 保存图片
ggsave(
  "results/figures/TAOK3_m5C_by_cell_line.pdf",
  plot = p,
  width = 9,
  height = 7
)

ggsave(
  "results/figures/TAOK3_m5C_by_cell_line.png",
  plot = p,
  width = 9,
  height = 7,
  dpi = 300
)
