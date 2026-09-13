############################################################
# A002. 原始全细胞大类发表版UMAP和marker dotplot
############################################################


############################################################
# A002.1 加载包并读取对象
############################################################

setwd(
  "C:/Users/phill/Documents/Lab/scRNAseq/GHRL25040179/0730-3 redo"
)

library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)

if (!exists("allcell_original")) {
  
  allcell_original <- readRDS(
    "0730_004_combined_harmony_mt_red_singlet_annotation_0.2.rds"
  )
}

DefaultAssay(allcell_original) <- "RNA"


############################################################
# A002.2 固定大类顺序和莫奈风格配色
############################################################

A002_celltype_levels <- c(
  "Malignant_Cells",
  "Macrophages",
  "Neutrophils",
  "CTLs/NK_Cells",
  "CAFs",
  "Endothelial_Cells"
)

allcell_original$A002_celltype <- factor(
  as.character(
    allcell_original$Annotated_by_hand
  ),
  levels = A002_celltype_levels
)

A002_celltype_colors <- c(
  "Malignant_Cells" = "#C98276",
  "Macrophages" = "#7895B2",
  "Neutrophils" = "#87A889",
  "CTLs/NK_Cells" = "#7AA6A8",
  "CAFs" = "#D1A56E",
  "Endothelial_Cells" = "#A98EB5"
)


############################################################
# A002.3 全细胞大类UMAP
############################################################

p_A002_UMAP <- DimPlot(
  allcell_original,
  reduction = "umap",
  group.by = "A002_celltype",
  cols = A002_celltype_colors,
  pt.size = 0.045,
  label = FALSE,
  shuffle = TRUE,
  seed = 730,
  raster = FALSE
) +
  coord_fixed() +
  labs(
    title = "Tumor microenvironment",
    color = NULL
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "right",
    legend.text = element_text(
      size = 10
    ),
    legend.key.height = grid::unit(
      0.58,
      "cm"
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )

ggsave(
  "A002_allcell_major_UMAP_publication.pdf",
  p_A002_UMAP,
  width = 7.6,
  height = 6.2
)


############################################################
# A002.4 D与DD分面全细胞UMAP
############################################################

A002_allcell_D <- subset(
  allcell_original,
  subset = group == "D"
)

A002_allcell_DD <- subset(
  allcell_original,
  subset = group == "DD"
)

p_A002_D <- DimPlot(
  A002_allcell_D,
  reduction = "umap",
  group.by = "A002_celltype",
  cols = A002_celltype_colors,
  pt.size = 0.04,
  label = FALSE,
  shuffle = TRUE,
  seed = 730,
  raster = FALSE
) +
  ggtitle("D") +
  coord_fixed() +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "none"
  )

p_A002_DD <- DimPlot(
  A002_allcell_DD,
  reduction = "umap",
  group.by = "A002_celltype",
  cols = A002_celltype_colors,
  pt.size = 0.04,
  label = FALSE,
  shuffle = TRUE,
  seed = 730,
  raster = FALSE
) +
  ggtitle("DD") +
  coord_fixed() +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "none"
  )

p_A002_D_DD <- (
  p_A002_D |
    p_A002_DD
) +
  plot_annotation(
    title = "Tumor microenvironment in D and DD tumors",
    theme = theme(
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      )
    )
  )

ggsave(
  "A002_allcell_major_UMAP_D_DD_publication.pdf",
  p_A002_D_DD,
  width = 10.2,
  height = 5.8
)


############################################################
# A002.5 设置六大类代表性marker
############################################################

A002_marker_information <- data.frame(
  gene = c(
    "Krt8",
    "Krt18",
    "Hmga2",
    "Mmp10",
    
    "Csf1r",
    "C1qa",
    "C1qc",
    "Mrc1",
    
    "S100a8",
    "S100a9",
    "Retnlg",
    "Lcn2",
    
    "Cd3d",
    "Cd3e",
    "Nkg7",
    "Gzmb",
    
    "Col1a1",
    "Col1a2",
    "Dcn",
    "Pdgfra",
    
    "Pecam1",
    "Cdh5",
    "Emcn",
    "Kdr"
  ),
  marker_class = rep(
    c(
      "Malignant",
      "Macrophage",
      "Neutrophil",
      "T/NK",
      "CAF",
      "Endothelial"
    ),
    each = 4
  ),
  stringsAsFactors = FALSE
)

A002_available_marker_information <-
  A002_marker_information[
    A002_marker_information$gene %in%
      rownames(allcell_original),
    ,
    drop = FALSE
  ]

cat(
  "\n用于A002 dotplot的marker：\n"
)

print(
  A002_available_marker_information,
  row.names = FALSE
)


############################################################
# A002.6 读取RNA counts和data
############################################################

A002_counts <- tryCatch(
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

A002_normalized_data <- GetAssayData(
  allcell_original,
  assay = "RNA",
  layer = "data"
)


############################################################
# A002.7 计算平均表达和至少2 UMI阳性率
############################################################

A002_dot_result_list <- list()
A002_dot_index <- 1

for (
  A002_current_celltype in
  A002_celltype_levels
) {
  
  A002_current_cells <- colnames(
    allcell_original
  )[
    as.character(
      allcell_original$A002_celltype
    ) == A002_current_celltype
  ]
  
  for (
    A002_current_gene in
    A002_available_marker_information$gene
  ) {
    
    A002_current_average <- mean(
      A002_normalized_data[
        A002_current_gene,
        A002_current_cells,
        drop = TRUE
      ]
    )
    
    A002_current_percentage <- mean(
      A002_counts[
        A002_current_gene,
        A002_current_cells,
        drop = TRUE
      ] >= 2
    ) * 100
    
    A002_current_marker_class <-
      A002_available_marker_information$marker_class[
        A002_available_marker_information$gene ==
          A002_current_gene
      ][1]
    
    A002_dot_result_list[[A002_dot_index]] <- data.frame(
      celltype = A002_current_celltype,
      gene = A002_current_gene,
      marker_class =
        A002_current_marker_class,
      average_expression =
        A002_current_average,
      percent_cells_2UMI =
        A002_current_percentage,
      stringsAsFactors = FALSE
    )
    
    A002_dot_index <- A002_dot_index + 1
  }
}

A002_dot_data <- do.call(
  rbind,
  A002_dot_result_list
)

rownames(A002_dot_data) <- NULL


############################################################
# A002.8 计算每个基因的相对平均表达
# 0为该基因在六类细胞中的最低平均表达
# 1为该基因在六类细胞中的最高平均表达
############################################################

A002_dot_data$relative_expression <- NA_real_

for (
  A002_current_gene in
  unique(A002_dot_data$gene)
) {
  
  A002_gene_rows <- which(
    A002_dot_data$gene ==
      A002_current_gene
  )
  
  A002_gene_values <-
    A002_dot_data$average_expression[
      A002_gene_rows
    ]
  
  A002_gene_min <- min(
    A002_gene_values,
    na.rm = TRUE
  )
  
  A002_gene_max <- max(
    A002_gene_values,
    na.rm = TRUE
  )
  
  if (
    A002_gene_max >
    A002_gene_min
  ) {
    
    A002_dot_data$relative_expression[
      A002_gene_rows
    ] <- (
      A002_gene_values -
        A002_gene_min
    ) / (
      A002_gene_max -
        A002_gene_min
    )
    
  } else {
    
    A002_dot_data$relative_expression[
      A002_gene_rows
    ] <- 0
  }
}

A002_dot_data$celltype <- factor(
  A002_dot_data$celltype,
  levels = rev(
    A002_celltype_levels
  )
)

A002_dot_data$gene <- factor(
  A002_dot_data$gene,
  levels =
    A002_available_marker_information$gene
)

A002_dot_data$marker_class <- factor(
  A002_dot_data$marker_class,
  levels = c(
    "Malignant",
    "Macrophage",
    "Neutrophil",
    "T/NK",
    "CAF",
    "Endothelial"
  )
)

write.csv(
  A002_dot_data,
  "A002_allcell_major_marker_dotplot_data.csv",
  row.names = FALSE
)


############################################################
# A002.9 绘制发表版marker dotplot
############################################################

p_A002_dotplot <- ggplot(
  A002_dot_data,
  aes(
    x = gene,
    y = celltype
  )
) +
  geom_point(
    aes(
      size = percent_cells_2UMI,
      color = relative_expression
    )
  ) +
  facet_grid(
    . ~ marker_class,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_size_area(
    name = "Cells with\n>=2 UMI (%)",
    max_size = 9,
    limits = c(
      0,
      100
    ),
    breaks = c(
      0,
      25,
      50,
      75,
      100
    )
  ) +
  scale_color_gradientn(
    name = "Relative mean\nexpression",
    colors = c(
      "#FFFFFF",
      "#D9E2E3",
      "#96B6B1",
      "#D5A28D",
      "#A94F45"
    ),
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
    base_size = 11
  ) +
  theme(
    strip.background = element_rect(
      fill = "#F3F0EA",
      color = "black",
      linewidth = 0.45
    ),
    strip.text = element_text(
      size = 10,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10,
      face = "italic"
    ),
    axis.text.y = element_text(
      size = 10
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.45
    ),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    )
  )

ggsave(
  "A002_allcell_major_marker_dotplot_publication.pdf",
  p_A002_dotplot,
  width = 13.5,
  height = 4.8
)


############################################################
# A002.10 保存带有A002字段的对象
############################################################

saveRDS(
  allcell_original,
  "A002_allcell_original_publication_ready.rds"
)


############################################################
# A002.11 最终统一打印
############################################################

cat(
  "\n",
  "============================================================\n",
  "A002最终需要反馈的控制台结果\n",
  "============================================================\n"
)

cat(
  "\n1. 六个大类细胞数：\n"
)

print(
  table(
    allcell_original$A002_celltype,
    useNA = "ifany"
  )
)

cat(
  "\n2. 每个marker平均表达最高的细胞类型：\n"
)

A002_marker_peak_table <- do.call(
  rbind,
  lapply(
    split(
      A002_dot_data,
      A002_dot_data$gene
    ),
    function(x) {
      x[
        which.max(
          x$average_expression
        ),
        c(
          "gene",
          "marker_class",
          "celltype",
          "average_expression",
          "percent_cells_2UMI"
        ),
        drop = FALSE
      ]
    }
  )
)

rownames(A002_marker_peak_table) <- NULL

print(
  A002_marker_peak_table,
  row.names = FALSE
)

write.csv(
  A002_marker_peak_table,
  "A002_allcell_marker_peak_celltype.csv",
  row.names = FALSE
)

cat(
  "\n",
  "============================================================\n",
  "A002控制台结果打印结束\n",
  "============================================================\n"
)