############################################################
# A005. 样本级伪批量GSVA免疫通路活性分析
#
# 分析层级：
# 1. 全部免疫细胞
# 2. 巨噬细胞
# 3. 中性粒细胞
# 4. T/NK细胞
#
# 每个细胞谱系 × 每个样本生成一个伪批量表达谱
# GSVA计算通路活性
# limma进行D与DD比较
############################################################


############################################################
# A005.1 加载对象和包
############################################################

setwd(
  "C:/Users/phill/Documents/Lab/scRNAseq/GHRL25040179/0730-3 redo"
)

library(Seurat)
library(Matrix)
library(edgeR)
library(limma)
library(GSVA)
library(GSEABase)
library(ggplot2)

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


############################################################
# A005.2 设置统一样本和配色
############################################################

A005_sample_levels <- c(
  "D1",
  "D2",
  "D3",
  "DD1",
  "DD2",
  "DD3"
)

A005_sample_group <- factor(
  c(
    "D",
    "D",
    "D",
    "DD",
    "DD",
    "DD"
  ),
  levels = c(
    "D",
    "DD"
  )
)

names(A005_sample_group) <-
  A005_sample_levels

A005_group_colors <- c(
  "D" = "#B65C50",
  "DD" = "#739B8C"
)

A005_direction_colors <- c(
  "Higher_in_D" = "#B65C50",
  "Higher_in_DD" = "#739B8C"
)


############################################################
# A005.3 定义功能基因集
############################################################

A005_gene_sets <- list(
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
  )
)

for (
  A005_current_set in
  names(A005_gene_sets)
) {
  
  A005_gene_sets[[A005_current_set]] <- intersect(
    A005_gene_sets[[A005_current_set]],
    rownames(allcell_original)
  )
}


############################################################
# A005.4 定义各细胞谱系使用的程序
############################################################

A005_lineage_programs <- list(
  All_Immune = c(
    "Mac_Suppressive",
    "Mac_Repair_Remodeling",
    "Neutrophil_MDSC_like",
    "Immune_IFN_dsDNA",
    "Immune_APC",
    "TNK_Cytotoxic",
    "TNK_Exhaustion"
  ),
  
  Macrophages = c(
    "Mac_Suppressive",
    "Mac_Repair_Remodeling",
    "Immune_IFN_dsDNA",
    "Immune_APC"
  ),
  
  Neutrophils = c(
    "Neutrophil_MDSC_like",
    "Immune_IFN_dsDNA",
    "Immune_APC"
  ),
  
  CTLs_NK_Cells = c(
    "TNK_Cytotoxic",
    "TNK_Exhaustion",
    "Immune_IFN_dsDNA",
    "Immune_APC"
  )
)


############################################################
# A005.5 读取RNA counts
############################################################

A005_counts <- tryCatch(
  GetAssayData(
    allcell_original,
    assay = "RNA",
    layer = "counts"
  ),
  error = function(e) {
    
    allcell_original <<- JoinLayers(
      allcell_original,
      assay = "RNA"
    )
    
    GetAssayData(
      allcell_original,
      assay = "RNA",
      layer = "counts"
    )
  }
)

cat(
  "\nA005 counts矩阵：",
  nrow(A005_counts),
  "genes ×",
  ncol(A005_counts),
  "cells\n"
)


############################################################
# A005.6 创建GSVA兼容函数
############################################################

A005_run_GSVA <- function(
    expression_matrix,
    gene_sets
) {
  
  A005_available_sets <- lapply(
    gene_sets,
    function(x) {
      intersect(
        x,
        rownames(expression_matrix)
      )
    }
  )
  
  A005_available_sets <- A005_available_sets[
    lengths(A005_available_sets) >= 3
  ]
  
  if (
    length(A005_available_sets) == 0
  ) {
    
    stop(
      "当前谱系没有至少包含3个基因的可用基因集。"
    )
  }
  
  if (
    "gsvaParam" %in%
    getNamespaceExports("GSVA")
  ) {
    
    A005_GSVA_parameter <- GSVA::gsvaParam(
      exprData = as.matrix(
        expression_matrix
      ),
      geneSets = A005_available_sets,
      kcdf = "Gaussian",
      minSize = 3,
      maxSize = nrow(
        expression_matrix
      )
    )
    
    A005_GSVA_score <- GSVA::gsva(
      A005_GSVA_parameter,
      verbose = FALSE
    )
    
  } else {
    
    A005_GSVA_score <- GSVA::gsva(
      expr = as.matrix(
        expression_matrix
      ),
      gset.idx.list =
        A005_available_sets,
      method = "gsva",
      kcdf = "Gaussian",
      min.sz = 3,
      max.sz = nrow(
        expression_matrix
      ),
      verbose = FALSE
    )
  }
  
  return(
    as.matrix(
      A005_GSVA_score
    )
  )
}


############################################################
# A005.7 逐个谱系建立样本伪批量表达谱
############################################################

A005_lineages <- names(
  A005_lineage_programs
)

A005_GSVA_score_list <- list()
A005_logCPM_list <- list()
A005_run_status_list <- list()

for (
  A005_current_lineage in
  A005_lineages
) {
  
  cat(
    "\n",
    "============================================================\n",
    "正在建立伪批量：",
    A005_current_lineage,
    "\n",
    "============================================================\n"
  )
  
  if (
    A005_current_lineage ==
    "All_Immune"
  ) {
    
    A005_current_cells <- colnames(
      allcell_original
    )[
      as.character(
        allcell_original$A002_celltype
      ) %in% c(
        "Macrophages",
        "Neutrophils",
        "CTLs/NK_Cells"
      )
    ]
    
  } else if (
    A005_current_lineage ==
    "CTLs_NK_Cells"
  ) {
    
    A005_current_cells <- colnames(
      allcell_original
    )[
      as.character(
        allcell_original$A002_celltype
      ) == "CTLs/NK_Cells"
    ]
    
  } else {
    
    A005_current_cells <- colnames(
      allcell_original
    )[
      as.character(
        allcell_original$A002_celltype
      ) == A005_current_lineage
    ]
  }
  
  A005_current_meta <-
    allcell_original@meta.data[
      A005_current_cells,
      ,
      drop = FALSE
    ]
  
  A005_current_sample_factor <- factor(
    A005_current_meta$orig.ident,
    levels = A005_sample_levels
  )
  
  A005_sample_design <- Matrix::sparse.model.matrix(
    ~ 0 + A005_current_sample_factor
  )
  
  colnames(A005_sample_design) <-
    A005_sample_levels
  
  A005_current_counts <- A005_counts[
    ,
    A005_current_cells,
    drop = FALSE
  ]
  
  A005_pseudobulk_counts <-
    A005_current_counts %*%
    A005_sample_design
  
  colnames(A005_pseudobulk_counts) <-
    A005_sample_levels
  
  A005_DGE <- edgeR::DGEList(
    counts = A005_pseudobulk_counts,
    group = A005_sample_group
  )
  
  A005_keep_gene <- rowSums(
    edgeR::cpm(
      A005_DGE
    ) > 1
  ) >= 2
  
  A005_DGE <- A005_DGE[
    A005_keep_gene,
    ,
    keep.lib.sizes = FALSE
  ]
  
  A005_DGE <- edgeR::calcNormFactors(
    A005_DGE,
    method = "TMM"
  )
  
  A005_logCPM <- edgeR::cpm(
    A005_DGE,
    log = TRUE,
    prior.count = 2
  )
  
  A005_current_program_names <-
    A005_lineage_programs[[A005_current_lineage]]
  
  A005_current_gene_sets <-
    A005_gene_sets[
      A005_current_program_names
    ]
  
  A005_GSVA_scores <- A005_run_GSVA(
    expression_matrix = A005_logCPM,
    gene_sets = A005_current_gene_sets
  )
  
  A005_GSVA_score_list[[A005_current_lineage]] <-
    A005_GSVA_scores
  
  A005_logCPM_list[[A005_current_lineage]] <-
    A005_logCPM
  
  A005_run_status_list[[A005_current_lineage]] <- data.frame(
    lineage = A005_current_lineage,
    cell_number =
      length(A005_current_cells),
    retained_gene_number =
      nrow(A005_logCPM),
    GSVA_program_number =
      nrow(A005_GSVA_scores),
    run_status = "Success",
    stringsAsFactors = FALSE
  )
}

A005_run_status <- do.call(
  rbind,
  A005_run_status_list
)

rownames(A005_run_status) <- NULL

write.csv(
  A005_run_status,
  "A005_pseudobulk_GSVA_run_status.csv",
  row.names = FALSE
)


############################################################
# A005.8 整理GSVA样本评分并运行limma
############################################################

A005_score_long_list <- list()
A005_limma_list <- list()
A005_result_index <- 1

for (
  A005_current_lineage in
  A005_lineages
) {
  
  A005_current_scores <-
    A005_GSVA_score_list[[A005_current_lineage]]
  
  A005_design <- model.matrix(
    ~ A005_sample_group
  )
  
  rownames(A005_design) <-
    A005_sample_levels
  
  A005_fit <- limma::lmFit(
    A005_current_scores,
    A005_design
  )
  
  A005_fit <- limma::eBayes(
    A005_fit,
    trend = TRUE,
    robust = TRUE
  )
  
  A005_limma_result <- limma::topTable(
    A005_fit,
    coef = 2,
    number = Inf,
    sort.by = "none"
  )
  
  A005_limma_result$program <-
    rownames(A005_limma_result)
  
  rownames(A005_limma_result) <- NULL
  
  A005_limma_result$lineage <-
    A005_current_lineage
  
  A005_limma_result$effect_D_minus_DD <-
    -A005_limma_result$logFC
  
  A005_limma_result$direction <- ifelse(
    A005_limma_result$effect_D_minus_DD > 0,
    "Higher_in_D",
    "Higher_in_DD"
  )
  
  A005_limma_list[[A005_current_lineage]] <-
    A005_limma_result
  
  for (
    A005_current_program in
    rownames(A005_current_scores)
  ) {
    
    A005_score_long_list[[A005_result_index]] <- data.frame(
      lineage = A005_current_lineage,
      program = A005_current_program,
      orig.ident =
        A005_sample_levels,
      group =
        as.character(
          A005_sample_group
        ),
      GSVA_score =
        as.numeric(
          A005_current_scores[
            A005_current_program,
            A005_sample_levels
          ]
        ),
      stringsAsFactors = FALSE
    )
    
    A005_result_index <-
      A005_result_index + 1
  }
}

A005_sample_scores <- do.call(
  rbind,
  A005_score_long_list
)

rownames(A005_sample_scores) <- NULL

A005_limma_results <- do.call(
  rbind,
  A005_limma_list
)

rownames(A005_limma_results) <- NULL

A005_limma_results$global_FDR <- p.adjust(
  A005_limma_results$P.Value,
  method = "BH"
)

A005_limma_results$lineage_program <- paste0(
  A005_limma_results$lineage,
  " | ",
  A005_limma_results$program
)

write.csv(
  A005_sample_scores,
  "A005_pseudobulk_GSVA_scores_by_sample.csv",
  row.names = FALSE
)

write.csv(
  A005_limma_results,
  "A005_pseudobulk_GSVA_DD_vs_D_limma.csv",
  row.names = FALSE
)


############################################################
# A005.9 GSVA评分样本箱线图
############################################################

A005_sample_scores$group <- factor(
  A005_sample_scores$group,
  levels = c(
    "D",
    "DD"
  )
)

A005_sample_scores$lineage_program <- paste0(
  A005_sample_scores$lineage,
  "\n",
  A005_sample_scores$program
)

p_A005_boxplot <- ggplot(
  A005_sample_scores,
  aes(
    x = group,
    y = GSVA_score,
    fill = group
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.80
  ) +
  geom_point(
    size = 2.2,
    position = position_jitter(
      width = 0.07
    )
  ) +
  facet_wrap(
    ~ lineage_program,
    ncol = 4,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = A005_group_colors
  ) +
  labs(
    title = "Pseudobulk GSVA immune-program activity",
    subtitle = "Each point represents one biological sample",
    x = NULL,
    y = "GSVA score"
  ) +
  theme_classic(
    base_size = 9
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
      linewidth = 0.35
    ),
    strip.text = element_text(
      size = 7.4,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 8,
      face = "bold"
    ),
    legend.position = "none"
  )

ggsave(
  "A005_pseudobulk_GSVA_scores_DD_vs_D.pdf",
  p_A005_boxplot,
  width = 13,
  height = 11
)


############################################################
# A005.10 绘制limma效应汇总图
############################################################

A005_effect_data <- A005_limma_results

A005_effect_data$lineage_program <- factor(
  A005_effect_data$lineage_program,
  levels = rev(
    A005_effect_data$lineage_program[
      order(
        A005_effect_data$effect_D_minus_DD
      )
    ]
  )
)

p_A005_effect <- ggplot(
  A005_effect_data,
  aes(
    x = effect_D_minus_DD,
    y = lineage_program,
    fill = direction
  )
) +
  geom_col(
    width = 0.68
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    color = "black"
  ) +
  scale_fill_manual(
    values = A005_direction_colors
  ) +
  labs(
    title = "Pseudobulk GSVA pathway changes",
    subtitle = "Positive values indicate higher activity in D",
    x = "GSVA effect (D - DD)",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(
    base_size = 9
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
      size = 7.8
    ),
    legend.position = "bottom"
  )

ggsave(
  "A005_pseudobulk_GSVA_D_minus_DD_effect.pdf",
  p_A005_effect,
  width = 9,
  height = 8.5
)


############################################################
# A005.11 样本级GSVA热图
############################################################

A005_heatmap_data <- A005_sample_scores

A005_heatmap_data$z_score <- NA_real_

for (
  A005_current_row in
  unique(A005_heatmap_data$lineage_program)
) {
  
  A005_current_rows <- which(
    A005_heatmap_data$lineage_program ==
      A005_current_row
  )
  
  A005_current_values <-
    A005_heatmap_data$GSVA_score[
      A005_current_rows
    ]
  
  A005_current_sd <- sd(
    A005_current_values
  )
  
  if (
    is.na(A005_current_sd) |
    A005_current_sd == 0
  ) {
    
    A005_heatmap_data$z_score[
      A005_current_rows
    ] <- 0
    
  } else {
    
    A005_heatmap_data$z_score[
      A005_current_rows
    ] <- as.numeric(
      scale(
        A005_current_values
      )
    )
  }
}

A005_heatmap_data$orig.ident <- factor(
  A005_heatmap_data$orig.ident,
  levels = A005_sample_levels
)

p_A005_heatmap <- ggplot(
  A005_heatmap_data,
  aes(
    x = orig.ident,
    y = lineage_program,
    fill = z_score
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
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
    title = "Pseudobulk GSVA activity by sample",
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      size = 9,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 7.5
    ),
    axis.ticks = element_blank()
  )

ggsave(
  "A005_pseudobulk_GSVA_sample_heatmap.pdf",
  p_A005_heatmap,
  width = 8,
  height = 8.5
)


############################################################
# A005.12 保存完整结果
############################################################

A005_output <- list(
  gene_sets =
    A005_gene_sets,
  lineage_programs =
    A005_lineage_programs,
  run_status =
    A005_run_status,
  logCPM =
    A005_logCPM_list,
  GSVA_score_matrices =
    A005_GSVA_score_list,
  sample_scores =
    A005_sample_scores,
  limma_results =
    A005_limma_results
)

saveRDS(
  A005_output,
  "A005_pseudobulk_GSVA_results.rds"
)


############################################################
# A005.13 最终统一打印
############################################################

cat(
  "\n",
  "============================================================\n",
  "A005最终需要反馈的控制台结果\n",
  "============================================================\n"
)

cat(
  "\n1. 各谱系伪批量GSVA运行状态：\n"
)

print(
  A005_run_status,
  row.names = FALSE
)

cat(
  "\n2. GSVA的D与DD limma结果：\n"
)

print(
  A005_limma_results[
    ,
    c(
      "lineage",
      "program",
      "effect_D_minus_DD",
      "P.Value",
      "adj.P.Val",
      "global_FDR",
      "direction"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat(
  "\n",
  "============================================================\n",
  "A005控制台结果打印结束\n",
  "============================================================\n"
)