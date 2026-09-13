############################################################
# 0730-3 redo
# A001-B001. 读取原始全细胞和巨噬细胞对象
# 核实原始注释、分群及TAM2状态
############################################################


############################################################
# A001.1 设置工作目录并加载包
############################################################

setwd(
  "C:/Users/phill/Documents/Lab/scRNAseq/GHRL25040179/0730-3 redo"
)

library(Seurat)
library(ggplot2)
library(patchwork)

cat(
  "\n当前工作目录：",
  getwd(),
  "\n"
)


############################################################
# A001.2 读取原始全细胞对象
############################################################

allcell_original <- readRDS(
  "0730_004_combined_harmony_mt_red_singlet_annotation_0.2.rds"
)

cat(
  "\n原始全细胞对象读取完成：\n"
)

print(
  allcell_original
)


############################################################
# A001.3 核实全细胞对象内容
############################################################

cat(
  "\n全细胞对象assay：\n"
)

print(
  Assays(allcell_original)
)

cat(
  "\n全细胞对象reductions：\n"
)

print(
  Reductions(allcell_original)
)

cat(
  "\n全细胞对象metadata字段：\n"
)

print(
  colnames(allcell_original@meta.data)
)

cat(
  "\n全细胞原始大类注释：\n"
)

print(
  table(
    allcell_original$Annotated_by_hand,
    useNA = "ifany"
  )
)

cat(
  "\n全细胞六个样本分布：\n"
)

print(
  table(
    allcell_original$orig.ident,
    allcell_original$group
  )
)

cat(
  "\n全细胞resolution 0.2分群：\n"
)

print(
  table(
    allcell_original$RNA_snn_res.0.2,
    useNA = "ifany"
  )
)


############################################################
# A001.4 保存全细胞核实表
############################################################

A001_allcell_celltype_number <- as.data.frame(
  table(
    celltype =
      allcell_original$Annotated_by_hand
  )
)

A001_allcell_sample_celltype <- as.data.frame.matrix(
  table(
    allcell_original$Annotated_by_hand,
    allcell_original$orig.ident
  )
)

write.csv(
  A001_allcell_celltype_number,
  "A001_allcell_original_celltype_number.csv",
  row.names = FALSE
)

write.csv(
  A001_allcell_sample_celltype,
  "A001_allcell_original_sample_celltype_count.csv",
  row.names = TRUE
)


############################################################
# B001.1 读取原始巨噬细胞对象
############################################################

mac_original <- readRDS(
  "0730_007_mac_0801_0.2.rds"
)

DefaultAssay(mac_original) <- "RNA"

cat(
  "\n原始巨噬细胞对象读取完成：\n"
)

print(
  mac_original
)


############################################################
# B001.2 核实巨噬细胞对象内容
############################################################

cat(
  "\n巨噬细胞对象assay：\n"
)

print(
  Assays(mac_original)
)

cat(
  "\n巨噬细胞对象reductions：\n"
)

print(
  Reductions(mac_original)
)

cat(
  "\n巨噬细胞对象metadata字段：\n"
)

print(
  colnames(mac_original@meta.data)
)

cat(
  "\n巨噬细胞六个样本分布：\n"
)

print(
  table(
    mac_original$orig.ident,
    mac_original$group
  )
)

cat(
  "\n原始resolution 0.2分群：\n"
)

print(
  table(
    mac_original$RNA_snn_res.0.2,
    useNA = "ifany"
  )
)


############################################################
# B001.3 严格恢复原始TAM1-TAM8命名
############################################################

B001_cluster_to_TAM <- c(
  "0" = "TAM1",
  "1" = "TAM2",
  "2" = "TAM3",
  "3" = "TAM4",
  "4" = "TAM5",
  "5" = "TAM6",
  "6" = "TAM7",
  "7" = "TAM8"
)

mac_original$B001_TAM <- unname(
  B001_cluster_to_TAM[
    as.character(
      mac_original$RNA_snn_res.0.2
    )
  ]
)

mac_original$B001_TAM <- factor(
  mac_original$B001_TAM,
  levels = c(
    "TAM1",
    "TAM2",
    "TAM3",
    "TAM4",
    "TAM5",
    "TAM6",
    "TAM7",
    "TAM8"
  )
)

Idents(mac_original) <- "B001_TAM"

cat(
  "\n恢复后的TAM1-TAM8细胞数：\n"
)

print(
  table(
    mac_original$B001_TAM,
    useNA = "ifany"
  )
)


############################################################
# B001.4 核实原注释是否已经存在
############################################################

if (
  "Annotated_by_hand" %in%
  colnames(mac_original@meta.data)
) {
  
  cat(
    "\n对象中原有的Annotated_by_hand：\n"
  )
  
  print(
    table(
      mac_original$Annotated_by_hand,
      useNA = "ifany"
    )
  )
  
  cat(
    "\n原注释与本次恢复命名的对应关系：\n"
  )
  
  print(
    table(
      original_annotation =
        mac_original$Annotated_by_hand,
      restored_annotation =
        mac_original$B001_TAM,
      useNA = "ifany"
    )
  )
}


############################################################
# B001.5 计算每个样本的TAM1-TAM8数量和百分比
############################################################

B001_TAM_sample_count <- table(
  TAM = mac_original$B001_TAM,
  sample = mac_original$orig.ident
)

B001_TAM_sample_percentage <- prop.table(
  B001_TAM_sample_count,
  margin = 2
) * 100

cat(
  "\n原始8个TAM在六个样本中的细胞数：\n"
)

print(
  B001_TAM_sample_count
)

cat(
  "\n每个样本内部8个TAM的百分比：\n"
)

print(
  round(
    B001_TAM_sample_percentage,
    2
  )
)

write.csv(
  as.data.frame.matrix(
    B001_TAM_sample_count
  ),
  "B001_original_8TAM_sample_cell_count.csv",
  row.names = TRUE
)

write.csv(
  round(
    as.data.frame.matrix(
      B001_TAM_sample_percentage
    ),
    4
  ),
  "B001_original_8TAM_sample_percentage.csv",
  row.names = TRUE
)


############################################################
# B001.6 单独打印原始TAM2的样本分布
############################################################

B001_TAM2_count <- B001_TAM_sample_count[
  "TAM2",
  ,
  drop = TRUE
]

B001_TAM2_percentage <-
  B001_TAM_sample_percentage[
    "TAM2",
    ,
    drop = TRUE
  ]

B001_TAM2_summary <- data.frame(
  orig.ident = names(
    B001_TAM2_count
  ),
  group = ifelse(
    grepl(
      "^DD",
      names(B001_TAM2_count)
    ),
    "DD",
    "D"
  ),
  TAM2_cell_number = as.numeric(
    B001_TAM2_count
  ),
  TAM2_percentage = as.numeric(
    B001_TAM2_percentage
  ),
  stringsAsFactors = FALSE
)

cat(
  "\n原始TAM2在六个样本中的数量和百分比：\n"
)

print(
  B001_TAM2_summary,
  row.names = FALSE
)

write.csv(
  B001_TAM2_summary,
  "B001_original_TAM2_sample_summary.csv",
  row.names = FALSE
)


############################################################
# B001.7 检查Creb5在原始对象中是否存在
############################################################

B001_Creb5_present <-
  "Creb5" %in% rownames(mac_original)

cat(
  "\nCreb5是否存在于原始巨噬对象：",
  B001_Creb5_present,
  "\n"
)

if (B001_Creb5_present) {
  
  B001_Creb5_expression <- FetchData(
    mac_original,
    vars = c(
      "Creb5",
      "B001_TAM",
      "orig.ident",
      "group"
    ),
    layer = "data"
  )
  
  B001_Creb5_summary <- aggregate(
    B001_Creb5_expression$Creb5,
    by = list(
      TAM =
        B001_Creb5_expression$B001_TAM
    ),
    FUN = function(x) {
      c(
        average_expression = mean(x),
        positive_percentage =
          mean(x > 0) * 100
      )
    }
  )
  
  cat(
    "\nCreb5在原始8个TAM中的表达概况：\n"
  )
  
  print(
    B001_Creb5_summary,
    row.names = FALSE
  )
}


############################################################
# B001.8 保存原始复现对象
############################################################

saveRDS(
  mac_original,
  "B001_original_mac_8TAM_restored.rds"
)


############################################################
# B001.9 最终统一打印需要反馈的结果
############################################################

cat(
  "\n",
  "============================================================\n",
  "A001-B001最终需要反馈的控制台结果\n",
  "============================================================\n"
)

cat(
  "\n1. 原始全细胞大类细胞数：\n"
)

print(
  table(
    allcell_original$Annotated_by_hand,
    useNA = "ifany"
  )
)

cat(
  "\n2. 原始巨噬细胞resolution 0.2分群：\n"
)

print(
  table(
    mac_original$RNA_snn_res.0.2,
    useNA = "ifany"
  )
)

cat(
  "\n3. 恢复后的TAM1-TAM8细胞数：\n"
)

print(
  table(
    mac_original$B001_TAM,
    useNA = "ifany"
  )
)

cat(
  "\n4. TAM1-TAM8在六个样本中的细胞数：\n"
)

print(
  B001_TAM_sample_count
)

cat(
  "\n5. 每个样本内部TAM1-TAM8百分比：\n"
)

print(
  round(
    B001_TAM_sample_percentage,
    2
  )
)

cat(
  "\n6. TAM2单独汇总：\n"
)

print(
  B001_TAM2_summary,
  row.names = FALSE
)

if (B001_Creb5_present) {
  
  cat(
    "\n7. Creb5在八个TAM中的表达概况：\n"
  )
  
  print(
    B001_Creb5_summary,
    row.names = FALSE
  )
}

cat(
  "\n",
  "============================================================\n",
  "A001-B001控制台结果打印结束\n",
  "============================================================\n"
)