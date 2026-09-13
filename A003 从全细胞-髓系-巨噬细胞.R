############################################################
# A003. 从全细胞组成引入髓系细胞和巨噬细胞
#
# 分析内容：
# 1. 六大类细胞每个样本组成
# 2. 六大类细胞D与DD比较
# 3. D减DD组成效应
# 4. 总髓系细胞比例
# 5. 巨噬细胞和中性粒细胞单独展示
############################################################


############################################################
# A003.1 加载对象和包
############################################################

setwd(
  "C:/Users/phill/Documents/Lab/scRNAseq/GHRL25040179/0730-3 redo"
)

library(Seurat)
library(ggplot2)
library(patchwork)

if (!exists("allcell_original")) {
  
  allcell_original <- readRDS(
    "A002_allcell_original_publication_ready.rds"
  )
}

if (
  !"A002_celltype" %in%
  colnames(allcell_original@meta.data)
) {
  
  allcell_original$A002_celltype <- factor(
    as.character(
      allcell_original$Annotated_by_hand
    ),
    levels = c(
      "Malignant_Cells",
      "Macrophages",
      "Neutrophils",
      "CTLs/NK_Cells",
      "CAFs",
      "Endothelial_Cells"
    )
  )
}


############################################################
# A003.2 设置统一顺序和低饱和度配色
############################################################

A003_celltype_levels <- c(
  "Malignant_Cells",
  "Macrophages",
  "Neutrophils",
  "CTLs/NK_Cells",
  "CAFs",
  "Endothelial_Cells"
)

A003_celltype_colors <- c(
  "Malignant_Cells" = "#C98276",
  "Macrophages" = "#7D9AB5",
  "Neutrophils" = "#87A889",
  "CTLs/NK_Cells" = "#72A6A1",
  "CAFs" = "#D1A56E",
  "Endothelial_Cells" = "#A98EB5"
)

A003_group_colors <- c(
  "D" = "#B65C50",
  "DD" = "#739B8C"
)

A003_direction_colors <- c(
  "Higher_in_D" = "#B65C50",
  "Higher_in_DD" = "#739B8C"
)


############################################################
# A003.3 计算六大类细胞样本数量和百分比
############################################################

A003_count_table <- table(
  celltype = allcell_original$A002_celltype,
  sample = allcell_original$orig.ident
)

A003_percentage_table <- prop.table(
  A003_count_table,
  margin = 2
) * 100

A003_composition_data <- as.data.frame(
  A003_percentage_table
)

colnames(A003_composition_data) <- c(
  "celltype",
  "orig.ident",
  "percentage"
)

A003_composition_data$group <- ifelse(
  grepl(
    "^DD",
    A003_composition_data$orig.ident
  ),
  "DD",
  "D"
)

A003_composition_data$group <- factor(
  A003_composition_data$group,
  levels = c(
    "D",
    "DD"
  )
)

A003_composition_data$celltype <- factor(
  A003_composition_data$celltype,
  levels = A003_celltype_levels
)

write.csv(
  as.data.frame.matrix(
    A003_count_table
  ),
  "A003_allcell_major_sample_cell_count.csv",
  row.names = TRUE
)

write.csv(
  as.data.frame.matrix(
    A003_percentage_table
  ),
  "A003_allcell_major_sample_percentage.csv",
  row.names = TRUE
)

write.csv(
  A003_composition_data,
  "A003_allcell_major_sample_percentage_long.csv",
  row.names = FALSE
)


############################################################
# A003.4 绘制每个样本的六大类组成堆叠图
############################################################

############################################################
# A003修正
# 修复DD3堆叠图最上层被删除的问题
# 不需要重新计算任何统计
############################################################


############################################################
# A003.4修正版：六大类组成堆叠图
############################################################

p_A003_stacked <- ggplot(
  A003_composition_data,
  aes(
    x = orig.ident,
    y = percentage,
    fill = celltype
  )
) +
  geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = A003_celltype_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = c(
      0,
      25,
      50,
      75,
      100
    ),
    expand = expansion(
      mult = c(
        0,
        0.01
      )
    )
  ) +
  coord_cartesian(
    ylim = c(
      0,
      100.5
    ),
    clip = "off"
  ) +
  labs(
    title = "Cellular composition of the tumor microenvironment",
    x = NULL,
    y = "Percentage of all cells (%)",
    fill = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    legend.position = "right",
    legend.text = element_text(
      size = 9
    )
  )

ggsave(
  "A003_allcell_major_composition_stacked_bar_v2.pdf",
  p_A003_stacked,
  width = 8,
  height = 5.7
)


############################################################
# A003.12修正版：重新生成组合图
############################################################

p_A003_summary <- (
  p_A003_stacked |
    p_A003_effect
) /
  (
    p_A003_major_myeloid |
      p_A003_myeloid
  ) +
  plot_annotation(
    tag_levels = "A"
  )

ggsave(
  "A003_allcell_to_myeloid_summary_v2.pdf",
  p_A003_summary,
  width = 13,
  height = 10
)


############################################################
# 检查六个样本百分比总和
############################################################

A003_sample_percentage_check <- colSums(
  A003_percentage_table
)

cat(
  "\n六个样本的百分比总和，均应为100：\n"
)

print(
  A003_sample_percentage_check
)

cat(
  "\nA003堆叠图和组合图修正完成。\n"
)


############################################################
# A003.5 计算六大类D与DD组成差异
############################################################

A003_statistics_list <- list()

for (
  A003_current_celltype in
  A003_celltype_levels
) {
  
  A003_current_data <- A003_composition_data[
    A003_composition_data$celltype ==
      A003_current_celltype,
    ,
    drop = FALSE
  ]
  
  A003_D_values <- A003_current_data$percentage[
    A003_current_data$group == "D"
  ]
  
  A003_DD_values <- A003_current_data$percentage[
    A003_current_data$group == "DD"
  ]
  
  A003_t_test <- tryCatch(
    t.test(
      A003_D_values,
      A003_DD_values
    ),
    error = function(e) {
      NULL
    }
  )
  
  A003_wilcox_test <- tryCatch(
    wilcox.test(
      A003_D_values,
      A003_DD_values,
      exact = TRUE
    ),
    error = function(e) {
      NULL
    }
  )
  
  A003_statistics_list[[A003_current_celltype]] <- data.frame(
    celltype = A003_current_celltype,
    mean_D = mean(
      A003_D_values
    ),
    mean_DD = mean(
      A003_DD_values
    ),
    delta_D_minus_DD =
      mean(A003_D_values) -
      mean(A003_DD_values),
    log2_ratio_DD_vs_D = log2(
      (
        mean(A003_DD_values) + 0.01
      ) /
        (
          mean(A003_D_values) + 0.01
        )
    ),
    t_test_p = ifelse(
      is.null(A003_t_test),
      NA,
      A003_t_test$p.value
    ),
    wilcox_p = ifelse(
      is.null(A003_wilcox_test),
      NA,
      A003_wilcox_test$p.value
    ),
    stringsAsFactors = FALSE
  )
}

A003_statistics <- do.call(
  rbind,
  A003_statistics_list
)

rownames(A003_statistics) <- NULL

A003_statistics$t_test_FDR <- p.adjust(
  A003_statistics$t_test_p,
  method = "BH"
)

A003_statistics$direction <- ifelse(
  A003_statistics$delta_D_minus_DD > 0,
  "Higher_in_D",
  "Higher_in_DD"
)

write.csv(
  A003_statistics,
  "A003_allcell_major_DD_vs_D_statistics.csv",
  row.names = FALSE
)


############################################################
# A003.6 六大类细胞D与DD百分比箱线图
############################################################

p_A003_boxplot <- ggplot(
  A003_composition_data,
  aes(
    x = group,
    y = percentage,
    fill = group
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.78
  ) +
  geom_point(
    aes(
      color = group
    ),
    size = 2.3,
    position = position_jitter(
      width = 0.07
    )
  ) +
  facet_wrap(
    ~ celltype,
    ncol = 3,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = A003_group_colors
  ) +
  scale_color_manual(
    values = A003_group_colors
  ) +
  labs(
    title = "Major cell-type composition in D and DD tumors",
    x = NULL,
    y = "Percentage of all cells (%)",
    fill = NULL,
    color = NULL
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5
    ),
    strip.background = element_rect(
      fill = "#F2F0EA",
      color = "black",
      linewidth = 0.4
    ),
    strip.text = element_text(
      size = 9,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 9,
      face = "bold"
    ),
    legend.position = "none"
  )

ggsave(
  "A003_allcell_major_DD_vs_D_percentage_boxplot.pdf",
  p_A003_boxplot,
  width = 8.7,
  height = 6.5
)


############################################################
# A003.7 绘制D减DD的组成效应图
############################################################

A003_effect_data <- A003_statistics

A003_effect_data$celltype <- factor(
  A003_effect_data$celltype,
  levels = rev(
    A003_effect_data$celltype[
      order(
        A003_effect_data$delta_D_minus_DD
      )
    ]
  )
)

p_A003_effect <- ggplot(
  A003_effect_data,
  aes(
    x = delta_D_minus_DD,
    y = celltype,
    fill = direction
  )
) +
  geom_col(
    width = 0.68
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.45,
    color = "black"
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%+.2f",
        delta_D_minus_DD
      )
    ),
    hjust = ifelse(
      A003_effect_data$delta_D_minus_DD >= 0,
      -0.15,
      1.15
    ),
    size = 3.4
  ) +
  scale_fill_manual(
    values = A003_direction_colors
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.18,
        0.18
      )
    )
  ) +
  labs(
    title = "D-associated changes in major cell composition",
    subtitle = "Positive values indicate higher abundance in D",
    x = "Difference in percentage points (D - DD)",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 9.5,
      hjust = 0.5
    ),
    axis.text.y = element_text(
      size = 10
    ),
    legend.position = "bottom"
  )

ggsave(
  "A003_allcell_major_D_minus_DD_effect.pdf",
  p_A003_effect,
  width = 7.3,
  height = 5.2
)


############################################################
# A003.8 计算每个样本总髓系细胞比例
#
# 总髓系细胞：
# Macrophages + Neutrophils
############################################################

A003_sample_levels <- c(
  "D1",
  "D2",
  "D3",
  "DD1",
  "DD2",
  "DD3"
)

A003_myeloid_summary <- data.frame(
  orig.ident = A003_sample_levels,
  group = c(
    "D",
    "D",
    "D",
    "DD",
    "DD",
    "DD"
  ),
  Macrophages = NA_real_,
  Neutrophils = NA_real_,
  Total_myeloid = NA_real_,
  stringsAsFactors = FALSE
)

for (
  A003_current_sample in
  A003_sample_levels
) {
  
  A003_sample_column <- which(
    colnames(A003_percentage_table) ==
      A003_current_sample
  )
  
  A003_sample_row <- which(
    A003_myeloid_summary$orig.ident ==
      A003_current_sample
  )
  
  A003_mac_percentage <-
    A003_percentage_table[
      "Macrophages",
      A003_sample_column
    ]
  
  A003_neutrophil_percentage <-
    A003_percentage_table[
      "Neutrophils",
      A003_sample_column
    ]
  
  A003_myeloid_summary$Macrophages[
    A003_sample_row
  ] <- A003_mac_percentage
  
  A003_myeloid_summary$Neutrophils[
    A003_sample_row
  ] <- A003_neutrophil_percentage
  
  A003_myeloid_summary$Total_myeloid[
    A003_sample_row
  ] <- A003_mac_percentage +
    A003_neutrophil_percentage
}

A003_myeloid_summary$group <- factor(
  A003_myeloid_summary$group,
  levels = c(
    "D",
    "DD"
  )
)

write.csv(
  A003_myeloid_summary,
  "A003_allcell_myeloid_sample_summary.csv",
  row.names = FALSE
)


############################################################
# A003.9 总髓系比例统计
############################################################

A003_myeloid_D <-
  A003_myeloid_summary$Total_myeloid[
    A003_myeloid_summary$group == "D"
  ]

A003_myeloid_DD <-
  A003_myeloid_summary$Total_myeloid[
    A003_myeloid_summary$group == "DD"
  ]

A003_myeloid_t_test <- t.test(
  A003_myeloid_D,
  A003_myeloid_DD
)

A003_myeloid_wilcox <- wilcox.test(
  A003_myeloid_D,
  A003_myeloid_DD,
  exact = TRUE
)

A003_myeloid_statistics <- data.frame(
  cell_population = "Macrophages_and_Neutrophils",
  mean_D = mean(
    A003_myeloid_D
  ),
  mean_DD = mean(
    A003_myeloid_DD
  ),
  delta_D_minus_DD =
    mean(A003_myeloid_D) -
    mean(A003_myeloid_DD),
  t_test_p =
    A003_myeloid_t_test$p.value,
  wilcox_p =
    A003_myeloid_wilcox$p.value,
  stringsAsFactors = FALSE
)

write.csv(
  A003_myeloid_statistics,
  "A003_allcell_myeloid_DD_vs_D_statistics.csv",
  row.names = FALSE
)


############################################################
# A003.10 绘制总髓系细胞比例
############################################################

A003_myeloid_label <- paste0(
  "D - DD = ",
  round(
    A003_myeloid_statistics$delta_D_minus_DD,
    2
  ),
  " percentage points\n",
  "P = ",
  format(
    A003_myeloid_statistics$t_test_p,
    digits = 2
  )
)

p_A003_myeloid <- ggplot(
  A003_myeloid_summary,
  aes(
    x = group,
    y = Total_myeloid,
    fill = group
  )
) +
  geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_point(
    size = 3,
    position = position_jitter(
      width = 0.07
    )
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(
      A003_myeloid_summary$Total_myeloid
    ) * 1.14,
    label = A003_myeloid_label,
    size = 3.5
  ) +
  scale_fill_manual(
    values = A003_group_colors
  ) +
  scale_y_continuous(
    limits = c(
      0,
      max(
        A003_myeloid_summary$Total_myeloid
      ) * 1.30
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  labs(
    title = "Myeloid-cell abundance",
    subtitle = "Macrophages + Neutrophils",
    x = NULL,
    y = "Percentage of all cells (%)"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 10,
      hjust = 0.5
    ),
    axis.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    legend.position = "none"
  )

ggsave(
  "A003_allcell_myeloid_DD_vs_D_boxplot.pdf",
  p_A003_myeloid,
  width = 4.5,
  height = 5.5
)


############################################################
# A003.11 巨噬和中性粒细胞单独过渡图
############################################################

A003_major_myeloid_data <- A003_composition_data[
  A003_composition_data$celltype %in%
    c(
      "Macrophages",
      "Neutrophils"
    ),
  ,
  drop = FALSE
]

A003_major_myeloid_data$celltype <- factor(
  A003_major_myeloid_data$celltype,
  levels = c(
    "Macrophages",
    "Neutrophils"
  )
)

p_A003_major_myeloid <- ggplot(
  A003_major_myeloid_data,
  aes(
    x = group,
    y = percentage,
    fill = group
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_point(
    size = 2.7,
    position = position_jitter(
      width = 0.07
    )
  ) +
  facet_wrap(
    ~ celltype,
    nrow = 1
  ) +
  scale_fill_manual(
    values = A003_group_colors
  ) +
  labs(
    title = "Major myeloid populations",
    x = NULL,
    y = "Percentage of all cells (%)"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5
    ),
    strip.background = element_rect(
      fill = "#F2F0EA",
      color = "black",
      linewidth = 0.4
    ),
    strip.text = element_text(
      size = 10,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    legend.position = "none"
  )

ggsave(
  "A003_allcell_macrophage_neutrophil_DD_vs_D.pdf",
  p_A003_major_myeloid,
  width = 7.2,
  height = 4.8
)


############################################################
# A003.12 组合成全细胞到髓系的总结图
############################################################

p_A003_summary <- (
  p_A003_stacked |
    p_A003_effect
) /
  (
    p_A003_major_myeloid |
      p_A003_myeloid
  ) +
  plot_annotation(
    tag_levels = "A"
  )

ggsave(
  "A003_allcell_to_myeloid_summary.pdf",
  p_A003_summary,
  width = 13,
  height = 10
)


############################################################
# A003.13 保存结果
############################################################

A003_output <- list(
  count_table =
    A003_count_table,
  percentage_table =
    A003_percentage_table,
  composition_data =
    A003_composition_data,
  major_cell_statistics =
    A003_statistics,
  myeloid_sample_summary =
    A003_myeloid_summary,
  myeloid_statistics =
    A003_myeloid_statistics
)

saveRDS(
  A003_output,
  "A003_allcell_to_myeloid_results.rds"
)


############################################################
# A003.14 最终统一打印
############################################################

cat(
  "\n",
  "============================================================\n",
  "A003最终需要反馈的控制台结果\n",
  "============================================================\n"
)

cat(
  "\n1. 六大类细胞在六个样本中的数量：\n"
)

print(
  A003_count_table
)

cat(
  "\n2. 六大类细胞在六个样本中的百分比：\n"
)

print(
  round(
    A003_percentage_table,
    2
  )
)

cat(
  "\n3. 六大类细胞D与DD组成统计：\n"
)

print(
  A003_statistics,
  row.names = FALSE
)

cat(
  "\n4. 每个样本的巨噬、中性粒和总髓系比例：\n"
)

print(
  A003_myeloid_summary,
  row.names = FALSE
)

cat(
  "\n5. 总髓系细胞D与DD统计：\n"
)

print(
  A003_myeloid_statistics,
  row.names = FALSE
)

cat(
  "\n",
  "============================================================\n",
  "A003控制台结果打印结束\n",
  "============================================================\n"
)