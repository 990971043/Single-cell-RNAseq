############################################################
# A004. 整体免疫微环境细胞级功能程序分析
#
# 分析对象：
# 1. 巨噬细胞免疫抑制程序
# 2. 巨噬细胞修复/重塑程序
# 3. 中性粒细胞PMN-MDSC样程序
# 4. T/NK细胞杀伤程序
# 5. T/NK细胞耗竭程序
# 6. 全部免疫细胞IFN/dsDNA程序
# 7. 全部免疫细胞抗原呈递程序
#
# 统计单位：6个生物学样本
############################################################


############################################################
# A004.1 加载对象和包
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

DefaultAssay(allcell_original) <- "RNA"

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

allcell_original$group <- factor(
  as.character(
    allcell_original$group
  ),
  levels = c(
    "D",
    "DD"
  )
)


############################################################
# A004.2 设置统一低饱和度红绿配色
############################################################

A004_group_colors <- c(
  "D" = "#B65C50",
  "DD" = "#739B8C"
)

A004_direction_colors <- c(
  "Higher_in_D" = "#B65C50",
  "Higher_in_DD" = "#739B8C"
)

A004_expression_colors <- c(
  "#FFFFFF",
  "#DCE4E1",
  "#94B4AE",
  "#D2A08B",
  "#A94F45"
)


############################################################
# A004.3 定义免疫功能程序
############################################################

A004_programs <- list(
  Mac_Suppressive = c(
    "Mrc1",
    "Cd163",
    "Mertk",
    "Tgfb1",
    "Il10",
    "Cd274",
    "Hmox1",
    "Msr1",
    "Trem2",
    "Lgals9",
    "Arg1",
    "Retnla"
  ),
  
  Mac_Repair_Remodeling = c(
    "Spp1",
    "Vegfa",
    "Mmp9",
    "Mmp12",
    "Mmp14",
    "Fn1",
    "Areg",
    "Gas6",
    "Ccl8",
    "Cfh"
  ),
  
  Neutrophil_MDSC_like = c(
    "S100a8",
    "S100a9",
    "Lcn2",
    "Retnlg",
    "Mmp9",
    "Arg2",
    "Il1b",
    "Cxcr2",
    "Wfdc17",
    "Olr1",
    "Ptgs2"
  ),
  
  TNK_Cytotoxic = c(
    "Nkg7",
    "Gzmb",
    "Prf1",
    "Ccl5",
    "Ifng",
    "Fasl"
  ),
  
  TNK_Exhaustion = c(
    "Pdcd1",
    "Lag3",
    "Havcr2",
    "Tigit",
    "Ctla4",
    "Tox"
  ),
  
  Immune_IFN_dsDNA = c(
    "Isg15",
    "Ifit1",
    "Ifit2",
    "Ifit3",
    "Rsad2",
    "Mx1",
    "Oas1a",
    "Oas3",
    "Zbp1",
    "Cxcl10",
    "Stat1",
    "Irf7"
  ),
  
  Immune_APC = c(
    "H2-Aa",
    "H2-Ab1",
    "H2-Eb1",
    "Cd74",
    "B2m",
    "Tap1",
    "Tap2",
    "Ciita"
  )
)


############################################################
# A004.4 指定每个程序在哪类细胞中计算
############################################################

A004_program_lineage <- c(
  "Mac_Suppressive" =
    "Macrophages",
  
  "Mac_Repair_Remodeling" =
    "Macrophages",
  
  "Neutrophil_MDSC_like" =
    "Neutrophils",
  
  "TNK_Cytotoxic" =
    "CTLs/NK_Cells",
  
  "TNK_Exhaustion" =
    "CTLs/NK_Cells",
  
  "Immune_IFN_dsDNA" =
    "All_Immune",
  
  "Immune_APC" =
    "All_Immune"
)

A004_immune_celltypes <- c(
  "Macrophages",
  "Neutrophils",
  "CTLs/NK_Cells"
)


############################################################
# A004.5 检查每个程序实际存在的基因
############################################################

A004_gene_availability_list <- list()

for (
  A004_current_program in
  names(A004_programs)
) {
  
  A004_requested_genes <-
    A004_programs[[A004_current_program]]
  
  A004_present_genes <- intersect(
    A004_requested_genes,
    rownames(allcell_original)
  )
  
  A004_missing_genes <- setdiff(
    A004_requested_genes,
    rownames(allcell_original)
  )
  
  A004_gene_availability_list[[A004_current_program]] <- data.frame(
    program = A004_current_program,
    gene = A004_requested_genes,
    gene_status = ifelse(
      A004_requested_genes %in%
        A004_present_genes,
      "Present",
      "Missing"
    ),
    stringsAsFactors = FALSE
  )
  
  cat(
    "\n",
    A004_current_program,
    "\n存在基因：",
    paste(
      A004_present_genes,
      collapse = ", "
    ),
    "\n"
  )
  
  if (length(A004_missing_genes) > 0) {
    
    cat(
      "缺失基因：",
      paste(
        A004_missing_genes,
        collapse = ", "
      ),
      "\n"
    )
  }
}

A004_gene_availability <- do.call(
  rbind,
  A004_gene_availability_list
)

rownames(A004_gene_availability) <- NULL

write.csv(
  A004_gene_availability,
  "A004_immune_program_gene_availability.csv",
  row.names = FALSE
)


############################################################
# A004.6 逐个功能程序进行AddModuleScore
############################################################

set.seed(730)

A004_cell_score_list <- list()
A004_score_index <- 1

for (
  A004_current_program in
  names(A004_programs)
) {
  
  A004_current_lineage <-
    A004_program_lineage[
      A004_current_program
    ]
  
  if (
    A004_current_lineage ==
    "All_Immune"
  ) {
    
    A004_target_cells <- colnames(
      allcell_original
    )[
      as.character(
        allcell_original$A002_celltype
      ) %in% A004_immune_celltypes
    ]
    
  } else {
    
    A004_target_cells <- colnames(
      allcell_original
    )[
      as.character(
        allcell_original$A002_celltype
      ) == A004_current_lineage
    ]
  }
  
  A004_present_genes <- intersect(
    A004_programs[[A004_current_program]],
    rownames(allcell_original)
  )
  
  if (
    length(A004_present_genes) < 3
  ) {
    
    cat(
      "\n跳过程序：",
      A004_current_program,
      "，因为可用基因少于3个。\n"
    )
    
    next
  }
  
  A004_current_object <- subset(
    allcell_original,
    cells = A004_target_cells
  )
  
  A004_current_object <- AddModuleScore(
    A004_current_object,
    features = list(
      A004_present_genes
    ),
    name = "A004_score",
    assay = "RNA",
    ctrl = 50,
    seed = 730,
    search = FALSE
  )
  
  A004_current_meta <-
    A004_current_object@meta.data
  
  A004_cell_score_list[[A004_score_index]] <- data.frame(
    cell = rownames(
      A004_current_meta
    ),
    orig.ident =
      as.character(
        A004_current_meta$orig.ident
      ),
    group =
      as.character(
        A004_current_meta$group
      ),
    celltype =
      as.character(
        A004_current_meta$A002_celltype
      ),
    lineage =
      as.character(
        A004_current_lineage
      ),
    program =
      A004_current_program,
    score =
      A004_current_meta$A004_score1,
    stringsAsFactors = FALSE
  )
  
  A004_score_index <-
    A004_score_index + 1
}

A004_cell_scores <- do.call(
  rbind,
  A004_cell_score_list
)

rownames(A004_cell_scores) <- NULL

write.csv(
  A004_cell_scores,
  "A004_immune_program_cell_scores.csv",
  row.names = FALSE
)


############################################################
# A004.7 按样本汇总功能程序评分
############################################################

A004_sample_scores <- aggregate(
  score ~ orig.ident +
    group +
    lineage +
    program,
  data = A004_cell_scores,
  FUN = mean
)

A004_sample_scores$group <- factor(
  A004_sample_scores$group,
  levels = c(
    "D",
    "DD"
  )
)

A004_program_order <- c(
  "Mac_Suppressive",
  "Mac_Repair_Remodeling",
  "Neutrophil_MDSC_like",
  "Immune_IFN_dsDNA",
  "Immune_APC",
  "TNK_Cytotoxic",
  "TNK_Exhaustion"
)

A004_sample_scores$program <- factor(
  A004_sample_scores$program,
  levels = A004_program_order
)

write.csv(
  A004_sample_scores,
  "A004_immune_program_scores_by_sample.csv",
  row.names = FALSE
)


############################################################
# A004.8 计算D与DD样本级统计
############################################################

A004_statistics_list <- list()

for (
  A004_current_program in
  A004_program_order
) {
  
  A004_current_data <- A004_sample_scores[
    as.character(
      A004_sample_scores$program
    ) == A004_current_program,
    ,
    drop = FALSE
  ]
  
  if (nrow(A004_current_data) == 0) {
    next
  }
  
  A004_D_scores <- A004_current_data$score[
    A004_current_data$group == "D"
  ]
  
  A004_DD_scores <- A004_current_data$score[
    A004_current_data$group == "DD"
  ]
  
  A004_t_test <- tryCatch(
    t.test(
      A004_D_scores,
      A004_DD_scores
    ),
    error = function(e) {
      NULL
    }
  )
  
  A004_wilcox_test <- tryCatch(
    wilcox.test(
      A004_D_scores,
      A004_DD_scores,
      exact = TRUE
    ),
    error = function(e) {
      NULL
    }
  )
  
  A004_all_scores <- c(
    A004_D_scores,
    A004_DD_scores
  )
  
  A004_score_sd <- sd(
    A004_all_scores
  )
  
  A004_standardized_effect <- ifelse(
    is.na(A004_score_sd) |
      A004_score_sd == 0,
    0,
    (
      mean(A004_D_scores) -
        mean(A004_DD_scores)
    ) / A004_score_sd
  )
  
  A004_statistics_list[[A004_current_program]] <- data.frame(
    program = A004_current_program,
    lineage =
      unique(
        A004_current_data$lineage
      )[1],
    mean_D =
      mean(A004_D_scores),
    mean_DD =
      mean(A004_DD_scores),
    delta_D_minus_DD =
      mean(A004_D_scores) -
      mean(A004_DD_scores),
    standardized_effect_D_minus_DD =
      A004_standardized_effect,
    t_test_p = ifelse(
      is.null(A004_t_test),
      NA,
      A004_t_test$p.value
    ),
    wilcox_p = ifelse(
      is.null(A004_wilcox_test),
      NA,
      A004_wilcox_test$p.value
    ),
    stringsAsFactors = FALSE
  )
}

A004_statistics <- do.call(
  rbind,
  A004_statistics_list
)

rownames(A004_statistics) <- NULL

A004_statistics$t_test_FDR <- p.adjust(
  A004_statistics$t_test_p,
  method = "BH"
)

A004_statistics$direction <- ifelse(
  A004_statistics$delta_D_minus_DD > 0,
  "Higher_in_D",
  "Higher_in_DD"
)

write.csv(
  A004_statistics,
  "A004_immune_program_DD_vs_D_statistics.csv",
  row.names = FALSE
)


############################################################
# A004.9 样本级功能程序箱线图
############################################################

p_A004_program_boxplot <- ggplot(
  A004_sample_scores,
  aes(
    x = group,
    y = score,
    fill = group
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_point(
    size = 2.5,
    position = position_jitter(
      width = 0.07
    )
  ) +
  facet_wrap(
    ~ program,
    ncol = 3,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = A004_group_colors
  ) +
  labs(
    title = "Immune-state programs in D and DD tumors",
    subtitle = "Each point represents one biological sample",
    x = NULL,
    y = "Module score"
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
    plot.subtitle = element_text(
      size = 9.5,
      hjust = 0.5
    ),
    strip.background = element_rect(
      fill = "#F2F0EA",
      color = "black",
      linewidth = 0.4
    ),
    strip.text = element_text(
      size = 8.5,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 9,
      face = "bold"
    ),
    legend.position = "none"
  )

ggsave(
  "A004_immune_program_scores_DD_vs_D.pdf",
  p_A004_program_boxplot,
  width = 10,
  height = 8
)


############################################################
# A004.10 功能程序标准化效应图
############################################################

A004_effect_data <- A004_statistics

A004_effect_data$program <- factor(
  A004_effect_data$program,
  levels = rev(
    A004_effect_data$program[
      order(
        A004_effect_data$
          standardized_effect_D_minus_DD
      )
    ]
  )
)

p_A004_effect <- ggplot(
  A004_effect_data,
  aes(
    x = standardized_effect_D_minus_DD,
    y = program,
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
        standardized_effect_D_minus_DD
      )
    ),
    hjust = ifelse(
      A004_effect_data$
        standardized_effect_D_minus_DD >= 0,
      -0.15,
      1.15
    ),
    size = 3.3
  ) +
  scale_fill_manual(
    values = A004_direction_colors
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.20,
        0.20
      )
    )
  ) +
  labs(
    title = "D-associated changes in immune-state programs",
    subtitle = "Positive values indicate higher activity in D",
    x = "Standardized effect (D - DD)",
    y = NULL,
    fill = NULL
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
    plot.subtitle = element_text(
      size = 9.5,
      hjust = 0.5
    ),
    axis.text.y = element_text(
      size = 9
    ),
    legend.position = "bottom"
  )

ggsave(
  "A004_immune_program_D_minus_DD_effect.pdf",
  p_A004_effect,
  width = 7.7,
  height = 5.8
)


############################################################
# A004.11 制作样本级程序热图数据
############################################################

A004_heatmap_data <- A004_sample_scores

A004_heatmap_data$z_score <- NA_real_

for (
  A004_current_program in
  A004_program_order
) {
  
  A004_current_rows <- which(
    as.character(
      A004_heatmap_data$program
    ) == A004_current_program
  )
  
  if (length(A004_current_rows) == 0) {
    next
  }
  
  A004_current_values <-
    A004_heatmap_data$score[
      A004_current_rows
    ]
  
  A004_current_sd <- sd(
    A004_current_values
  )
  
  if (
    is.na(A004_current_sd) |
    A004_current_sd == 0
  ) {
    
    A004_heatmap_data$z_score[
      A004_current_rows
    ] <- 0
    
  } else {
    
    A004_heatmap_data$z_score[
      A004_current_rows
    ] <- as.numeric(
      scale(
        A004_current_values
      )
    )
  }
}

A004_heatmap_data$orig.ident <- factor(
  A004_heatmap_data$orig.ident,
  levels = c(
    "D1",
    "D2",
    "D3",
    "DD1",
    "DD2",
    "DD3"
  )
)

A004_heatmap_data$program <- factor(
  A004_heatmap_data$program,
  levels = rev(
    A004_program_order
  )
)

p_A004_heatmap <- ggplot(
  A004_heatmap_data,
  aes(
    x = orig.ident,
    y = program,
    fill = z_score
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.65
  ) +
  scale_fill_gradient2(
    low = "#739B8C",
    mid = "#F5F2EB",
    high = "#B65C50",
    midpoint = 0,
    limits = c(
      -2,
      2
    ),
    oob = scales::squish,
    name = "Z-score"
  ) +
  labs(
    title = "Immune-state program activity by sample",
    x = NULL,
    y = NULL
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
    axis.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 9
    ),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

ggsave(
  "A004_immune_program_sample_heatmap.pdf",
  p_A004_heatmap,
  width = 7.2,
  height = 5.3
)


############################################################
# A004.12 免疫功能marker D/DD dotplot
############################################################

A004_dot_marker_information <- data.frame(
  gene = c(
    "Mrc1",
    "Cd163",
    "Mertk",
    "Tgfb1",
    "Cd274",
    "Hmox1",
    
    "Spp1",
    "Vegfa",
    "Mmp9",
    "Mmp12",
    "Areg",
    "Gas6",
    
    "S100a8",
    "S100a9",
    "Lcn2",
    "Retnlg",
    "Arg2",
    "Cxcr2",
    
    "Isg15",
    "Ifit1",
    "Ifit3",
    "Rsad2",
    "Mx1",
    "Zbp1",
    
    "Nkg7",
    "Gzmb",
    "Prf1",
    "Ccl5",
    "Ifng",
    "Fasl"
  ),
  marker_class = rep(
    c(
      "Suppressive",
      "Repair/remodeling",
      "MDSC-like",
      "IFN/dsDNA",
      "Cytotoxic"
    ),
    each = 6
  ),
  stringsAsFactors = FALSE
)

A004_dot_marker_information <-
  A004_dot_marker_information[
    A004_dot_marker_information$gene %in%
      rownames(allcell_original),
    ,
    drop = FALSE
  ]

A004_immune_object <- subset(
  allcell_original,
  subset = A002_celltype %in%
    c(
      "Macrophages",
      "Neutrophils",
      "CTLs/NK_Cells"
    )
)

A004_immune_object$A004_celltype_group <- paste0(
  as.character(
    A004_immune_object$A002_celltype
  ),
  "_",
  as.character(
    A004_immune_object$group
  )
)

A004_dotplot_object <- DotPlot(
  A004_immune_object,
  features =
    A004_dot_marker_information$gene,
  group.by = "A004_celltype_group",
  assay = "RNA",
  scale = FALSE
)

A004_dot_data <- A004_dotplot_object$data

A004_dot_data$gene <- as.character(
  A004_dot_data$features.plot
)

A004_dot_data$celltype_group <-
  as.character(
    A004_dot_data$id
  )

A004_dot_data$average_expression <-
  A004_dot_data$avg.exp

A004_dot_data$positive_percentage <-
  A004_dot_data$pct.exp

A004_dot_data$marker_class <-
  A004_dot_marker_information$marker_class[
    match(
      A004_dot_data$gene,
      A004_dot_marker_information$gene
    )
  ]

A004_dot_data$relative_expression <- NA_real_

for (
  A004_current_gene in
  unique(A004_dot_data$gene)
) {
  
  A004_current_rows <- which(
    A004_dot_data$gene ==
      A004_current_gene
  )
  
  A004_current_values <-
    A004_dot_data$average_expression[
      A004_current_rows
    ]
  
  A004_current_min <- min(
    A004_current_values,
    na.rm = TRUE
  )
  
  A004_current_max <- max(
    A004_current_values,
    na.rm = TRUE
  )
  
  if (
    A004_current_max >
    A004_current_min
  ) {
    
    A004_dot_data$relative_expression[
      A004_current_rows
    ] <- (
      A004_current_values -
        A004_current_min
    ) / (
      A004_current_max -
        A004_current_min
    )
    
  } else {
    
    A004_dot_data$relative_expression[
      A004_current_rows
    ] <- 0
  }
}

A004_dot_data$celltype_group <- factor(
  A004_dot_data$celltype_group,
  levels = rev(
    c(
      "Macrophages_D",
      "Macrophages_DD",
      "Neutrophils_D",
      "Neutrophils_DD",
      "CTLs/NK_Cells_D",
      "CTLs/NK_Cells_DD"
    )
  )
)

A004_dot_data$gene <- factor(
  A004_dot_data$gene,
  levels =
    A004_dot_marker_information$gene
)

A004_dot_data$marker_class <- factor(
  A004_dot_data$marker_class,
  levels = unique(
    A004_dot_marker_information$marker_class
  )
)

write.csv(
  A004_dot_data,
  "A004_immune_state_marker_dotplot_data.csv",
  row.names = FALSE
)

p_A004_dotplot <- ggplot(
  A004_dot_data,
  aes(
    x = gene,
    y = celltype_group
  )
) +
  geom_point(
    aes(
      size = positive_percentage,
      color = relative_expression
    )
  ) +
  facet_grid(
    . ~ marker_class,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_size_area(
    name = "Expressing\ncells (%)",
    max_size = 8,
    limits = c(
      0,
      100
    )
  ) +
  scale_color_gradientn(
    name = "Relative mean\nexpression",
    colors = A004_expression_colors,
    limits = c(
      0,
      1
    ),
    oob = scales::squish
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    strip.background = element_rect(
      fill = "#F2F0EA",
      color = "black",
      linewidth = 0.4
    ),
    strip.text = element_text(
      size = 8.5,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8,
      face = "italic"
    ),
    axis.text.y = element_text(
      size = 8.5
    ),
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.4
    ),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  "A004_immune_state_marker_dotplot.pdf",
  p_A004_dotplot,
  width = 15,
  height = 5
)


############################################################
# A004.13 保存完整结果
############################################################

A004_output <- list(
  programs = A004_programs,
  gene_availability =
    A004_gene_availability,
  cell_scores =
    A004_cell_scores,
  sample_scores =
    A004_sample_scores,
  statistics =
    A004_statistics,
  heatmap_data =
    A004_heatmap_data,
  dotplot_data =
    A004_dot_data
)

saveRDS(
  A004_output,
  "A004_immune_state_program_results.rds"
)


############################################################
# A004.14 最终统一打印
############################################################

cat(
  "\n",
  "============================================================\n",
  "A004最终需要反馈的控制台结果\n",
  "============================================================\n"
)

cat(
  "\n1. 各程序实际使用的基因数：\n"
)

A004_used_gene_number <- aggregate(
  gene ~ program,
  data = A004_gene_availability[
    A004_gene_availability$gene_status ==
      "Present",
    ,
    drop = FALSE
  ],
  FUN = length
)

colnames(A004_used_gene_number)[2] <-
  "used_gene_number"

print(
  A004_used_gene_number,
  row.names = FALSE
)

cat(
  "\n2. 每个样本的免疫程序评分：\n"
)

print(
  A004_sample_scores,
  row.names = FALSE
)

cat(
  "\n3. D与DD免疫程序统计：\n"
)

print(
  A004_statistics,
  row.names = FALSE
)

cat(
  "\n",
  "============================================================\n",
  "A004控制台结果打印结束\n",
  "============================================================\n"
)