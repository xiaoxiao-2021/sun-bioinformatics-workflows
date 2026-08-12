library(readxl)

input_file <- "datasets/raw/基因表达水平检测结果-U251.xlsx"

print(excel_sheets(input_file))
#读取sheet1数据
dat <- read_excel(input_file, sheet = 1)

cat("数据维度：\n")
print(dim(dat))

cat("\n列名：\n")
print(colnames(dat))

cat("\n前6行：\n")
print(head(dat))

cat("\n数据结构：\n")
str(dat)

sample_cols <- c(
  "OE1", "OE2", "OE3",
  "NC1", "NC2", "NC3"
)

stopifnot(all(sample_cols %in% colnames(dat)))

expr <- as.matrix(dat[, sample_cols])
storage.mode(expr) <- "numeric"

cat("\nNA数量：\n")
print(sum(is.na(expr)))

cat("\ngene_id重复数量：\n")
print(sum(duplicated(dat$gene_id)))

cat("\n整个矩阵中0的比例：\n")
print(mean(expr == 0, na.rm = TRUE))

cat("\n整个矩阵中0的比例：\n")
print(mean(expr == 0, na.rm = TRUE))

detected_n <- rowSums(expr > 0, na.rm = TRUE)

cat("\n每个基因在多少个样本中检测到：\n")
print(table(detected_n))


cat("\n整数值比例：\n")
print(mean(abs(expr - round(expr)) < 1e-8, na.rm = TRUE))

cat("\n表达值范围：\n")
print(range(expr, na.rm = TRUE))

cat("\n非零表达值分位数：\n")

print(
  quantile(
    expr[expr > 0],
    probs = c(
      0,
      0.25,
      0.50,
      0.75,
      0.90,
      0.95,
      0.99,
      1
    ),
    na.rm = TRUE
  )
)

filter_summary <- c(
  detected_gt0_in3 = sum(rowSums(expr > 0) >= 3),
  expr_ge_0.1_in3 = sum(rowSums(expr >= 0.1) >= 3),
  expr_ge_0.5_in3 = sum(rowSums(expr >= 0.5) >= 3),
  expr_ge_1_in3   = sum(rowSums(expr >= 1) >= 3),
  expr_ge_1_in2   = sum(rowSums(expr >= 1) >= 2)
)

print(filter_summary)


