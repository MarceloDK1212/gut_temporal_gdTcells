library(Seurat)
library(tidyverse)
library(dplyr)
library(readxl)
library(colorspace)

# download gut_temporal_gdTcells_processed.rds from here: https://drive.google.com/drive/folders/1EWbUZrY9am-fp6OK6hSwn5LkSeahgUv2?usp=drive_link
# download alternative_volcanoes.xlsx from here: https://drive.google.com/drive/folders/1EWbUZrY9am-fp6OK6hSwn5LkSeahgUv2?usp=drive_link
# figure aesthetics were further modified using Adobe Illustrator

gut_gdT_obj <- readRDS("gut_temporal_gdTcells_processed.rds")

# FIGURE 2A
DimPlot(gut_gdT_obj, group.by = "RNA_snn_res.0.6")

# FIGURE 2B
FeaturePlot(gut_gdT_obj, features = "Trgv6", order = TRUE)
FeaturePlot(gut_gdT_obj, features = "Trgv4", order = TRUE)
FeaturePlot(gut_gdT_obj, features = "Trdv2-2", order = TRUE)
FeaturePlot(gut_gdT_obj, features = "Trdv5", order = TRUE)

# FIGURE 2C
FeaturePlot(gut_gdT_obj, features = "Rorc", order = TRUE)
FeaturePlot(gut_gdT_obj, features = "Tbx21", order = TRUE)
FeaturePlot(gut_gdT_obj, features = "Il23r", order = TRUE)
FeaturePlot(gut_gdT_obj, features = "Cd27", order = TRUE)

# FIGURE 2D
DimPlot(gut_gdT_obj, split.by = "timepoint")

# Figure 2E
markers <- FindAllMarkers(gut_gdT_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
significant_markers <- dplyr::filter(markers, p_val_adj < 0.05)
clusters_to_plot <- c(0, 4, 5, 10, 11)
markers_subset <- significant_markers %>% dplyr::filter(cluster %in% clusters_to_plot)
top10_subset <- markers_subset %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pseudo_gut_gdT_obj <- AggregateExpression(gut_gdT_obj, assays = "RNA", return.seurat = T, group.by = "RNA_snn_res.0.6")
clusters_to_plot <- c("g0", "g5", "g4", "g10", "g11")
cells_to_plot <- WhichCells(pseudo_gut_gdT_obj, idents = clusters_to_plot)
mycols <- rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
DoHeatmap(pseudo_gut_gdT_obj, features = top10_subset$gene, cells = cells_to_plot, draw.lines = FALSE, group.bar.height = 0.025) + scale_fill_gradientn(colours = mycols) 

# Figure 2F
data <- read_excel("C:/Users/magre/Downloads/alternative_volcanoes.xlsx", sheet = "Vg6") #excel sheet with pre-selected differentially expressed genes from Table S2 and their log2FC
data_clean <- data %>% mutate(across(-1, ~ as.numeric(str_replace(., ",", "."))))
data_long <- data_clean %>% pivot_longer(cols = -1, names_to = "Comparison", values_to = "FoldChange")
desired_order <- c("d7vsd0", "d21vsd7", "d42vsd21", "d42vsd7")
data_long$Comparison <- factor(data_long$Comparison, levels = desired_order)
pal <- qualitative_hcl(43, palette = "Dark 3")
ggplot(data_long, aes(x = Comparison, y = FoldChange, color = ...1)) +
  geom_point(size = 4, alpha = 1, position = position_jitter(width = 0.15), na.rm = TRUE) +
  geom_text(aes(label = ...1), size = 3, hjust = 0.5, show.legend = F) +
  scale_color_manual(values = pal) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "none")
# note: Jitter added to avoid overplotting (position_jitter, width = 0.15). No random seed was set, so point positions may vary slightly between runs

# Figure 2G
FeaturePlot(gut_gdT_obj, features = "Il22", order = TRUE, split.by = "timepoint", cols = c("grey", "red"))
FeaturePlot(gut_gdT_obj, features = "Bhlhe40", order = TRUE, split.by = "timepoint", cols = c("grey", "red"))

# Figure 2H
data <- read_excel("C:/Users/magre/Downloads/alternative_volcanoes.xlsx", sheet = "Vg4") #excel sheet with pre-selected differentially expressed genes from Table S2 and their log2FC
data_clean <- data %>% mutate(across(-1, ~ as.numeric(str_replace(., ",", "."))))
data_long <- data_clean %>% pivot_longer(cols = -1, names_to = "Comparison", values_to = "FoldChange")
desired_order <- c("d7vsd0", "d21vsd7", "d42vsd21", "d42vsd7")
data_long$Comparison <- factor(data_long$Comparison, levels = desired_order)
pal <- qualitative_hcl(57, palette = "Dark 3")
ggplot(data_long, aes(x = Comparison, y = FoldChange, color = ...1)) +
  geom_point(size = 4, alpha = 1, position = position_jitter(width = 0.05), na.rm = TRUE) +
  geom_text(aes(label = ...1), size = 3, hjust = 0.5, show.legend = F) +
  scale_color_manual(values = pal) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "none")
# note: Jitter added to avoid overplotting (position_jitter, width = 0.05). No random seed was set, so point positions may vary slightly between runs

# Figure 2I
FeaturePlot(gut_gdT_obj, features = "Gzma", order = TRUE, split.by = "timepoint", cols = c("grey", "red"))
FeaturePlot(gut_gdT_obj, features = "Itga1", order = TRUE, split.by = "timepoint", cols = c("grey", "red"))
