library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(Matrix)
library(zellkonverter)
library(SingleCellExperiment)
library(scales)

setwd("~/Downloads")

h5ad_file <- "SG_publication_allmice_DVC.h5ad"
sce <- readH5AD(h5ad_file)

# convert to Seurat
allmice.data <- as.Seurat(sce, counts = "X", data = NULL)
saveRDS(allmice.data, file = "allmice_from_h5ad.rds")

print(allmice.data)
dim(allmice.data)
head(colnames(allmice.data@meta.data))
head(rownames(allmice.data))

# qc checks (mitochondrial + glp1r/gfral)
allmice.data[["percent.mt"]] <- PercentageFeatureSet(allmice.data, pattern = "^(mt-|MT-)")
summary(allmice.data$percent.mt)

VlnPlot(allmice.data, features = "percent.mt", ncol = 1)

gene_names <- c("Glp1r", "Gfral")
present_genes <- gene_names[gene_names %in% rownames(allmice.data)]
print(present_genes)

# subset to Glp1r+, Gfral+
gene_expression <- FetchData(allmice.data, vars = c("Glp1r", "Gfral"))

cells_to_keep <- rownames(gene_expression)[
  gene_expression$Glp1r > 0 | gene_expression$Gfral > 0
]

filtered_Glp1r_Gfral <- allmice.data[, cells_to_keep]
saveRDS(filtered_Glp1r_Gfral, file = "filtered_Glp1r_Gfral_rawsubset.rds")

# pre-processing
filtered_Glp1r_Gfral <- NormalizeData(filtered_Glp1r_Gfral, verbose = FALSE)
filtered_Glp1r_Gfral <- FindVariableFeatures(filtered_Glp1r_Gfral, verbose = FALSE)
filtered_Glp1r_Gfral <- ScaleData(
  filtered_Glp1r_Gfral,
  features = rownames(filtered_Glp1r_Gfral),
  verbose = FALSE
)
filtered_Glp1r_Gfral <- RunPCA(filtered_Glp1r_Gfral, verbose = FALSE)

ElbowPlot(filtered_Glp1r_Gfral)

# adriana's filters
filtered_Glp1r_Gfral <- FindNeighbors(filtered_Glp1r_Gfral, dims = 1:6, verbose = FALSE)
filtered_Glp1r_Gfral <- FindClusters(filtered_Glp1r_Gfral, resolution = 0.15, verbose = FALSE)
filtered_Glp1r_Gfral <- RunUMAP(filtered_Glp1r_Gfral, dims = 1:6, verbose = FALSE)

DimPlot(filtered_Glp1r_Gfral, reduction = "umap", label = TRUE)
saveRDS(filtered_Glp1r_Gfral, file = "filtered_Glp1r_Gfral_processed.rds")

# metadata labels
## receptor group labels
glp1r_cells <- WhichCells(filtered_Glp1r_Gfral, expression = Glp1r > 0)
gfral_cells <- WhichCells(filtered_Glp1r_Gfral, expression = Gfral > 0)

filtered_Glp1r_Gfral$cell_type <- ifelse(
  Cells(filtered_Glp1r_Gfral) %in% glp1r_cells & Cells(filtered_Glp1r_Gfral) %in% gfral_cells,
  "Both",
  ifelse(
    Cells(filtered_Glp1r_Gfral) %in% glp1r_cells,
    "Glp1r",
    "Gfral"
  )
)

filtered_Glp1r_Gfral$cell_type <- factor(
  filtered_Glp1r_Gfral$cell_type,
  levels = c("Glp1r", "Gfral", "Both")
)

table(filtered_Glp1r_Gfral$cell_type)

## feeding condition (adlib / fast / refed)
table(filtered_Glp1r_Gfral$orig.ident)

filtered_Glp1r_Gfral$diet_category <- case_when(
  filtered_Glp1r_Gfral$orig.ident == "adlib" ~ "Ad libitum",
  filtered_Glp1r_Gfral$orig.ident == "fast"  ~ "Fasted",
  filtered_Glp1r_Gfral$orig.ident == "refed" ~ "Refed",
  TRUE ~ "Other"
)

filtered_Glp1r_Gfral$diet_category <- factor(
  filtered_Glp1r_Gfral$diet_category,
  levels = c("Ad libitum", "Fasted", "Refed", "Other")
)

table(filtered_Glp1r_Gfral$diet_category)

## annotate with atlas identity columns
colnames(filtered_Glp1r_Gfral@meta.data)

table(filtered_Glp1r_Gfral$identity_layer1)
table(filtered_Glp1r_Gfral$identity_layer2)
table(filtered_Glp1r_Gfral$identity_layer3)

cluster_id2 <- filtered_Glp1r_Gfral@meta.data %>%
  select(seurat_clusters, identity_layer2) %>%
  group_by(seurat_clusters, identity_layer2) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(proportion = count / sum(count)) %>%
  arrange(seurat_clusters, desc(proportion))

print(cluster_id2)

cluster_id3 <- filtered_Glp1r_Gfral@meta.data %>%
  select(seurat_clusters, identity_layer3) %>%
  group_by(seurat_clusters, identity_layer3) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(proportion = count / sum(count)) %>%
  arrange(seurat_clusters, desc(proportion))

print(cluster_id3)

if (length(levels(filtered_Glp1r_Gfral)) == 6) {
  new.cluster.ids <- c(
    "Mixed Neurons",
    "Sall3 Neurons",
    "Monoamine Neurons",
    "Excitatory Neurons",
    "Astrocyte",
    "Oligodendrocyte"
  )
  names(new.cluster.ids) <- levels(filtered_Glp1r_Gfral)
  filtered_Glp1r_Gfral <- RenameIdents(filtered_Glp1r_Gfral, new.cluster.ids)
}

filtered_Glp1r_Gfral$annotated_cluster <- Idents(filtered_Glp1r_Gfral)

DimPlot(filtered_Glp1r_Gfral, reduction = "umap", label = TRUE) +
  ggtitle("Annotated Glp1r/Gfral receptor-positive cells")

table(filtered_Glp1r_Gfral$annotated_cluster)

# neuron-only object for downstream analyses
neurons <- subset(filtered_Glp1r_Gfral, subset = identity_layer1 == "neuron")
saveRDS(neurons, file = "filtered_Glp1r_Gfral_neurons.rds")

table(neurons$cell_type)
table(neurons$diet_category)
table(Idents(neurons))

# composition across feeding states
## receptor-group composition across diet
comp_receptor <- filtered_Glp1r_Gfral@meta.data %>%
  group_by(diet_category, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(diet_category) %>%
  mutate(proportion = count / sum(count))

print(comp_receptor)

p_comp_receptor_prop <- ggplot(comp_receptor, aes(x = diet_category, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(
    title = "Composition of receptor-positive groups across feeding states",
    x = "Feeding state",
    y = "Proportion"
  ) +
  theme_bw()

p_comp_receptor_count <- ggplot(comp_receptor, aes(x = diet_category, y = count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Counts of receptor-positive groups across feeding states",
    x = "Feeding state",
    y = "Cell count"
  ) +
  theme_bw()

p_comp_receptor_prop
p_comp_receptor_count

ggsave("composition_receptor_proportion.png", p_comp_receptor_prop, width = 7, height = 5, dpi = 300)
ggsave("composition_receptor_counts.png", p_comp_receptor_count, width = 7, height = 5, dpi = 300)


## CLEANED GRAPH (remove "Other")
comp_receptor_clean <- comp_receptor %>%
  filter(diet_category %in% c("Ad libitum", "Fasted", "Refed")) %>%
  mutate(
    diet_category = factor(diet_category, levels = c("Ad libitum", "Fasted", "Refed")),
    cell_type = factor(cell_type, levels = c("Glp1r", "Gfral", "Both"))
  )

p1 <- ggplot(comp_receptor_clean, aes(x = diet_category, y = proportion, fill = cell_type)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.6) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "Glp1r" = "#E76F51",
      "Gfral" = "#2A9D8F",
      "Both"  = "#457B9D"
    )
  ) +
  labs(
    title = "Receptor-positive cell composition across feeding states",
    subtitle = "Relative proportion of Glp1r+, Gfral+, and double-positive cells",
    x = NULL,
    y = "Proportion of cells",
    fill = "Cell group"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 12, color = "grey30"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

p1
ggsave("composition_receptor_proportion_clean.png", p1, width = 8, height = 5, dpi = 300)


comp_receptor_counts_clean <- comp_receptor %>%
  filter(diet_category %in% c("Ad libitum", "Fasted", "Refed")) %>%
  mutate(
    diet_category = factor(diet_category, levels = c("Ad libitum", "Fasted", "Refed")),
    cell_type = factor(cell_type, levels = c("Glp1r", "Gfral", "Both"))
  )

p2 <- ggplot(comp_receptor_counts_clean, aes(x = diet_category, y = count, fill = cell_type)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.6) +
  scale_fill_manual(
    values = c(
      "Glp1r" = "#E76F51",
      "Gfral" = "#2A9D8F",
      "Both"  = "#457B9D"
    )
  ) +
  labs(
    title = "Counts of receptor-positive cells across feeding states",
    subtitle = "Absolute number of cells in each receptor-defined group",
    x = NULL,
    y = "Cell count",
    fill = "Cell group"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 12, color = "grey30"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

p2
ggsave("composition_receptor_counts_pretty.png", p2, width = 8, height = 5, dpi = 300)

p1 + p2

# chi-square test
receptor_table <- table(filtered_Glp1r_Gfral$diet_category, filtered_Glp1r_Gfral$cell_type)
print(receptor_table)

chisq_receptor <- chisq.test(receptor_table)
print(chisq_receptor)

chisq.test(
  table(neurons$diet_category, neurons$cell_type)[c("Ad libitum","Fasted","Refed"),]
)


# feeding condition analysis
comp_subtype_clean <- comp_subtype %>%
  filter(diet_category %in% c("Ad libitum", "Fasted", "Refed")) %>%
  filter(annotated_cluster %in% c(
    "Mixed Neurons",
    "Sall3 Neurons",
    "Monoamine Neurons",
    "Excitatory Neurons"
  )) %>%
  mutate(
    diet_category = factor(diet_category, levels = c("Ad libitum", "Fasted", "Refed")),
    cell_type = factor(cell_type, levels = c("Glp1r", "Gfral", "Both")),
    annotated_cluster = factor(
      annotated_cluster,
      levels = c("Mixed Neurons", "Sall3 Neurons", "Monoamine Neurons", "Excitatory Neurons")
    )
  )

p_subtype_clean <- ggplot(
  comp_subtype_clean,
  aes(x = diet_category, y = prop, fill = cell_type)
) +
  geom_col(width = 0.8, color = "white", linewidth = 0.6) +
  facet_wrap(~ annotated_cluster, ncol = 2) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "Glp1r" = "#E76F51",
      "Gfral" = "#2A9D8F",
      "Both"  = "#457B9D"
    )
  ) +
  labs(
    title = "Receptor composition within neuronal subtypes",
    subtitle = "Glp1r+, Gfral+, and double-positive cells across feeding states",
    x = NULL,
    y = "Proportion of cells",
    fill = "Cell group"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 12, color = "grey30"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

p_subtype_clean

mixed <- subset(neurons, annotated_cluster == "Mixed Neurons")
chisq.test(table(mixed$diet_category, mixed$cell_type)[c("Ad libitum","Fasted","Refed"), ])

exc <- subset(neurons, annotated_cluster == "Excitatory Neurons")
chisq.test(table(exc$diet_category, exc$cell_type)[c("Ad libitum","Fasted","Refed"), ])

mono <- subset(neurons, annotated_cluster == "Monoamine Neurons")
chisq.test(table(mono$diet_category, mono$cell_type)[c("Ad libitum","Fasted","Refed"), ])

neurons_no_oligo <- neurons %>%
  subset(annotated_cluster != "Oligodendrocyte")

p_glp1r <- VlnPlot(
  neurons_no_oligo,
  features = "Glp1r",
  group.by = "annotated_cluster",
  pt.size = 0
) +
  labs(
    title = "Glp1r expression across cell types",
    x = "Cell type",
    y = "Expression level"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

p_gfral <- VlnPlot(
  neurons_no_oligo,
  features = "Gfral",
  group.by = "annotated_cluster",
  pt.size = 0
) +
  labs(
    title = "Gfral expression across cell types",
    x = "Cell type",
    y = "Expression level"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

p_glp1r | p_gfral

avg_expr <- AverageExpression(
  neurons_no_oligo,
  features = c("Glp1r", "Gfral"),
  group.by = "annotated_cluster",
  assays = DefaultAssay(neurons_no_oligo)
)

avg_expr
names(avg_expr)



# pseudotime attempt
library(SingleCellExperiment)
library(slingshot)

sce <- as.SingleCellExperiment(neurons)

sce <- slingshot(
  sce,
  clusterLabels = sce$annotated_cluster,
  reducedDim = "UMAP"
)

library(SingleCellExperiment)
library(slingshot)

sce <- as.SingleCellExperiment(neurons)

# Check reduced dimensions
dim(reducedDims(sce)$UMAP)
sum(is.na(reducedDims(sce)$UMAP))
head(reducedDims(sce)$UMAP)

# Check cluster labels
table(is.na(sce$annotated_cluster))
table(sce$annotated_cluster, useNA = "ifany")

umap_mat <- Embeddings(neurons, "umap")

keep_cells <- complete.cases(umap_mat) & !is.na(neurons$annotated_cluster)

neurons_clean <- subset(neurons, cells = colnames(neurons)[keep_cells])

# Rebuild SCE
sce <- as.SingleCellExperiment(neurons_clean)

# Explicitly set UMAP and cluster labels
reducedDims(sce)$UMAP <- Embeddings(neurons_clean, "umap")
colData(sce)$cluster <- neurons_clean$annotated_cluster

# Sanity checks
sum(is.na(reducedDims(sce)$UMAP))
table(is.na(colData(sce)$cluster))
table(colData(sce)$cluster)

receptor bias=corr(gene,Glp1r)−corr(gene,Gfral)

library(dplyr)
library(ggplot2)

plot_df <- corr_df %>%
  filter(!gene %in% c("Glp1r","Gfral")) %>%
  slice_max(abs(receptor_bias), n = 20)

ggplot(plot_df,
       aes(x = reorder(gene, receptor_bias),
           y = receptor_bias,
           fill = receptor_bias > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c("#2C7BB6","#D7191C"),
                    labels=c("Gfral-associated","Glp1r-associated"),
                    name="Program") +
  theme_classic(base_size = 14) +
  labs(
    title = "Genes associated with receptor-linked transcriptional programs",
    x = "Gene",
    y = "Correlation difference\n(Glp1r correlation − Gfral correlation)"
  ) +
  theme(legend.position="top")








# pathway/module
inflammatory_genes <- c("Nfkbia","Jun","Fos","Stat3","Irf1","Socs3","Cxcl10","B2m")
oxphos_genes <- c("Atp5f1a","Cox4i1","Ndufa1","Ndufb8","Uqcr10","Sdhd","Atp5me")
kras_genes <- c("Dusp1","Fos","Jun","Egr1","Spry2","Rasd1","Pik3r3")

present <- function(g) g[g %in% rownames(filtered_Glp1r_Gfral)]

present(inflammatory_genes)
present(oxphos_genes)
present(kras_genes)

filtered_Glp1r_Gfral <- AddModuleScore(
  filtered_Glp1r_Gfral,
  features = list(present(inflammatory_genes)),
  name = "InflammatoryScore"
)

filtered_Glp1r_Gfral <- AddModuleScore(
  filtered_Glp1r_Gfral,
  features = list(present(oxphos_genes)),
  name = "OxphosScore"
)

filtered_Glp1r_Gfral <- AddModuleScore(
  filtered_Glp1r_Gfral,
  features = list(present(kras_genes)),
  name = "KrasScore"
)

colnames(filtered_Glp1r_Gfral@meta.data)

VlnPlot(filtered_Glp1r_Gfral,
        features = "OxphosScore1",
        group.by = "diet_category",
        pt.size = 0)

VlnPlot(filtered_Glp1r_Gfral,
        features = "KrasScore1",
        group.by = "diet_category",
        pt.size = 0)




## 9b) Cluster composition across diet
comp_cluster <- filtered_Glp1r_Gfral@meta.data %>%
  group_by(diet_category, annotated_cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(diet_category) %>%
  mutate(proportion = count / sum(count))

print(comp_cluster)

p_comp_cluster_prop <- ggplot(comp_cluster, aes(x = diet_category, y = proportion, fill = annotated_cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(
    title = "Annotated cluster composition across feeding states",
    x = "Feeding state",
    y = "Proportion"
  ) +
  theme_bw()

p_comp_cluster_count <- ggplot(comp_cluster, aes(x = diet_category, y = count, fill = annotated_cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Annotated cluster counts across feeding states",
    x = "Feeding state",
    y = "Cell count"
  ) +
  theme_bw()

p_comp_cluster_prop
p_comp_cluster_count

ggsave("composition_cluster_proportion.png", p_comp_cluster_prop, width = 8, height = 5, dpi = 300)
ggsave("composition_cluster_counts.png", p_comp_cluster_count, width = 8, height = 5, dpi = 300)

# Stats
cluster_table <- table(filtered_Glp1r_Gfral$diet_category, filtered_Glp1r_Gfral$annotated_cluster)
print(cluster_table)

chisq_cluster <- chisq.test(cluster_table)
print(chisq_cluster)

## 9c) Neuron-only cluster composition across diet
neurons$annotated_cluster <- Idents(neurons)

comp_neurons <- neurons@meta.data %>%
  group_by(diet_category, annotated_cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(diet_category) %>%
  mutate(proportion = count / sum(count))

p_comp_neurons <- ggplot(comp_neurons, aes(x = diet_category, y = proportion, fill = annotated_cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(
    title = "Neuron subtype composition across feeding states",
    x = "Feeding state",
    y = "Proportion"
  ) +
  theme_bw()

p_comp_neurons
ggsave("composition_neurons_proportion.png", p_comp_neurons, width = 8, height = 5, dpi = 300)

############################
# 10) NEW ANALYSIS 2: Module / pathway scoring
############################
# Since your teammate slides suggest pathway-level differences are more informative
# than sparse DEGs, score a few pathway-like gene modules per cell.

# IMPORTANT:
# These are compact gene sets for a quick class-project analysis.
# If some genes are missing, we keep only those present.

present_only <- function(gene_vec, obj) {
  unique(gene_vec[gene_vec %in% rownames(obj)])
}

# Small curated gene lists
inflammatory_genes <- c("Nfkbia", "Jun", "Fos", "Stat3", "Irf1", "Socs3", "Cxcl10", "B2m")
oxphos_genes      <- c("Atp5f1a", "Cox4i1", "Ndufa1", "Ndufb8", "Uqcr10", "Sdhd", "Atp5me")
kras_genes        <- c("Dusp1", "Fos", "Jun", "Egr1", "Spry2", "Rasd1", "Pik3r3")
ifng_genes        <- c("Stat1", "Irf1", "Ifit3", "Isg15", "Cxcl10", "Tap1", "B2m")
xeno_genes        <- c("Gstm1", "Gpx1", "Aldh1a1", "Nqo1", "Hmox1")

inflammatory_present <- present_only(inflammatory_genes, filtered_Glp1r_Gfral)
oxphos_present       <- present_only(oxphos_genes, filtered_Glp1r_Gfral)
kras_present         <- present_only(kras_genes, filtered_Glp1r_Gfral)
ifng_present         <- present_only(ifng_genes, filtered_Glp1r_Gfral)
xeno_present         <- present_only(xeno_genes, filtered_Glp1r_Gfral)

print(inflammatory_present)
print(oxphos_present)
print(kras_present)
print(ifng_present)
print(xeno_present)

# Add module scores
if (length(inflammatory_present) >= 2) {
  filtered_Glp1r_Gfral <- AddModuleScore(
    filtered_Glp1r_Gfral,
    features = list(inflammatory_present),
    name = "InflammatoryScore"
  )
}

if (length(oxphos_present) >= 2) {
  filtered_Glp1r_Gfral <- AddModuleScore(
    filtered_Glp1r_Gfral,
    features = list(oxphos_present),
    name = "OxphosScore"
  )
}

if (length(kras_present) >= 2) {
  filtered_Glp1r_Gfral <- AddModuleScore(
    filtered_Glp1r_Gfral,
    features = list(kras_present),
    name = "KrasScore"
  )
}

if (length(ifng_present) >= 2) {
  filtered_Glp1r_Gfral <- AddModuleScore(
    filtered_Glp1r_Gfral,
    features = list(ifng_present),
    name = "IFNgScore"
  )
}

if (length(xeno_present) >= 2) {
  filtered_Glp1r_Gfral <- AddModuleScore(
    filtered_Glp1r_Gfral,
    features = list(xeno_present),
    name = "XenoScore"
  )
}

# See created metadata columns
grep("Score1", colnames(filtered_Glp1r_Gfral@meta.data), value = TRUE)

############################
# 11) Plot module scores in all receptor-positive cells
############################

# By feeding state
if ("InflammatoryScore1" %in% colnames(filtered_Glp1r_Gfral@meta.data)) {
  p_inflam_diet <- VlnPlot(filtered_Glp1r_Gfral, features = "InflammatoryScore1", group.by = "diet_category", pt.size = 0)
  print(p_inflam_diet)
  ggsave("InflammatoryScore_by_diet.png", p_inflam_diet, width = 7, height = 5, dpi = 300)
}

if ("OxphosScore1" %in% colnames(filtered_Glp1r_Gfral@meta.data)) {
  p_oxphos_diet <- VlnPlot(filtered_Glp1r_Gfral, features = "OxphosScore1", group.by = "diet_category", pt.size = 0)
  print(p_oxphos_diet)
  ggsave("OxphosScore_by_diet.png", p_oxphos_diet, width = 7, height = 5, dpi = 300)
}

if ("KrasScore1" %in% colnames(filtered_Glp1r_Gfral@meta.data)) {
  p_kras_diet <- VlnPlot(filtered_Glp1r_Gfral, features = "KrasScore1", group.by = "diet_category", pt.size = 0)
  print(p_kras_diet)
  ggsave("KrasScore_by_diet.png", p_kras_diet, width = 7, height = 5, dpi = 300)
}

# By receptor group
if ("InflammatoryScore1" %in% colnames(filtered_Glp1r_Gfral@meta.data)) {
  p_inflam_celltype <- VlnPlot(filtered_Glp1r_Gfral, features = "InflammatoryScore1", group.by = "cell_type", pt.size = 0)
  print(p_inflam_celltype)
  ggsave("InflammatoryScore_by_celltype.png", p_inflam_celltype, width = 7, height = 5, dpi = 300)
}

if ("OxphosScore1" %in% colnames(filtered_Glp1r_Gfral@meta.data)) {
  p_oxphos_celltype <- VlnPlot(filtered_Glp1r_Gfral, features = "OxphosScore1", group.by = "cell_type", pt.size = 0)
  print(p_oxphos_celltype)
  ggsave("OxphosScore_by_celltype.png", p_oxphos_celltype, width = 7, height = 5, dpi = 300)
}

if ("KrasScore1" %in% colnames(filtered_Glp1r_Gfral@meta.data)) {
  p_kras_celltype <- VlnPlot(filtered_Glp1r_Gfral, features = "KrasScore1", group.by = "cell_type", pt.size = 0)
  print(p_kras_celltype)
  ggsave("KrasScore_by_celltype.png", p_kras_celltype, width = 7, height = 5, dpi = 300)
}

############################
# 12) Stats on module scores
############################
meta_all <- filtered_Glp1r_Gfral@meta.data

if ("InflammatoryScore1" %in% colnames(meta_all)) {
  print(kruskal.test(InflammatoryScore1 ~ diet_category, data = meta_all))
  print(kruskal.test(InflammatoryScore1 ~ cell_type, data = meta_all))
}

if ("OxphosScore1" %in% colnames(meta_all)) {
  print(kruskal.test(OxphosScore1 ~ diet_category, data = meta_all))
  print(kruskal.test(OxphosScore1 ~ cell_type, data = meta_all))
}

if ("KrasScore1" %in% colnames(meta_all)) {
  print(kruskal.test(KrasScore1 ~ diet_category, data = meta_all))
  print(kruskal.test(KrasScore1 ~ cell_type, data = meta_all))
}

############################
# 13) Repeat module scoring in neurons only
############################
# Keep neuron object aligned with parent metadata
neurons <- subset(filtered_Glp1r_Gfral, subset = identity_layer1 == "neuron")

if (length(inflammatory_present) >= 2) {
  neurons <- AddModuleScore(neurons, features = list(inflammatory_present), name = "InflammatoryScore")
}
if (length(oxphos_present) >= 2) {
  neurons <- AddModuleScore(neurons, features = list(oxphos_present), name = "OxphosScore")
}
if (length(kras_present) >= 2) {
  neurons <- AddModuleScore(neurons, features = list(kras_present), name = "KrasScore")
}

meta_neurons <- neurons@meta.data

if ("InflammatoryScore1" %in% colnames(meta_neurons)) {
  p_neuron_inflam_diet <- VlnPlot(neurons, features = "InflammatoryScore1", group.by = "diet_category", pt.size = 0)
  print(p_neuron_inflam_diet)
  ggsave("Neuron_InflammatoryScore_by_diet.png", p_neuron_inflam_diet, width = 7, height = 5, dpi = 300)
}

if ("OxphosScore1" %in% colnames(meta_neurons)) {
  p_neuron_oxphos_diet <- VlnPlot(neurons, features = "OxphosScore1", group.by = "diet_category", pt.size = 0)
  print(p_neuron_oxphos_diet)
  ggsave("Neuron_OxphosScore_by_diet.png", p_neuron_oxphos_diet, width = 7, height = 5, dpi = 300)
}

if ("KrasScore1" %in% colnames(meta_neurons)) {
  p_neuron_kras_celltype <- VlnPlot(neurons, features = "KrasScore1", group.by = "cell_type", pt.size = 0)
  print(p_neuron_kras_celltype)
  ggsave("Neuron_KrasScore_by_celltype.png", p_neuron_kras_celltype, width = 7, height = 5, dpi = 300)
}

if ("InflammatoryScore1" %in% colnames(meta_neurons)) {
  print(kruskal.test(InflammatoryScore1 ~ diet_category, data = meta_neurons))
}
if ("OxphosScore1" %in% colnames(meta_neurons)) {
  print(kruskal.test(OxphosScore1 ~ diet_category, data = meta_neurons))
}
if ("KrasScore1" %in% colnames(meta_neurons)) {
  print(kruskal.test(KrasScore1 ~ cell_type, data = meta_neurons))
}

############################
# 14) Save final objects
############################
saveRDS(filtered_Glp1r_Gfral, file = "filtered_Glp1r_Gfral_final_with_new_analyses.rds")
saveRDS(neurons, file = "neurons_final_with_module_scores.rds")

cat("Done.\n")