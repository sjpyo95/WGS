# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(data.table)
library(stats)
library(pheatmap)

# ✅ 1. Load data
df_wide <- read.csv("77ea-matrix.csv", header=TRUE, row.names=1, check.names=FALSE)

# ✅ 2. Remove variants starting with "MT:"
df_wide <- df_wide[!grepl("^MT:", rownames(df_wide)), ]

head(df_wide)
# ✅ 3. Convert data from Wide to Long format
df_long <- df_wide %>%
  tibble::rownames_to_column(var = "Variant") %>%
  pivot_longer(
    cols = -Variant,
    names_to = "SampleID",
    values_to = "Genotype"
  )

# ✅ 4. Define sample groups (Control vs Target)
df_long <- df_long %>%
  mutate(
    Group = ifelse(grepl("^OUT", SampleID), "OUT", "IN"))

head(df_long)
dim(df_long)

# ✅ 5. Compute genotype consistency within each group
consistency_metrics <- df_long %>%
  group_by(Variant, Group) %>%
  summarise(
    total_samples = n(),
    major_genotype = as.numeric(names(which.max(table(Genotype)))),
    major_freq = max(table(Genotype)) / n(),
    .groups = 'drop'
  )

# ✅ 6. Compare major genotypes between Control and Target
consistency_comparison <- consistency_metrics %>%
  pivot_wider(
    id_cols = Variant,
    names_from = Group,
    values_from = c(major_genotype, major_freq),
    names_sep = '_',
    values_fn = list(
      major_genotype = function(x) x[1],
      major_freq = mean
    )
  )

# ✅ 7. Filter stable variants (Concordance ≥ 65% in each group)
threshold <- 0.65
consistent_variants_union <- consistency_comparison %>%
  filter((!is.na(major_freq_control) & major_freq_control >= threshold) &
           (!is.na(major_freq_target) & major_freq_target >= threshold))

# ✅ 8. Identify variants with different major genotypes between Control and Target
dif_major <- consistent_variants_union %>%
  filter(major_genotype_control != major_genotype_target)

dif_major_var <- dif_major$Variant

df_dif <- df_long %>%
  filter(Variant %in% dif_major_var)

# ✅ 9. Statistical analysis: logFC computation + T-test + Wilcoxon rank-sum test
results_diff <- df_dif %>%
  group_by(Variant) %>%
  do({
    dat <- .
    meanControl <- mean(dat$Genotype[dat$Group == 'control'], na.rm = TRUE)
    meanTarget <- mean(dat$Genotype[dat$Group == 'target'], na.rm = TRUE)
    
    eps <- 1e-6  # Small value to prevent division by zero
    
    # ✅ New logFC calculation: Always use the larger value as numerator
    if (meanTarget >= meanControl) {
      logFC <- log2((meanTarget + eps) / (meanControl + eps))
    } else {
      logFC <- -log2((meanControl + eps) / (meanTarget + eps))
    }
    
    # Perform T-test
    t_test_res <- t.test(Genotype ~ Group, data = dat)
    
    # Perform Wilcoxon rank-sum test (non-parametric test)
    wilcox_res <- wilcox.test(Genotype ~ Group, data = dat, alternative = "two.sided")
    
    data.frame(
      mean_control = meanControl,
      mean_target = meanTarget,
      logFC = logFC,
      t_p_value = t_test_res$p.value,
      wilcox_p_value = wilcox_res$p.value
    )
  }) %>%
  ungroup()

# ✅ 10. Apply multiple testing correction (Benjamini-Hochberg FDR)
results_diff <- results_diff %>%
  mutate(
    t_p_value_adj = p.adjust(t_p_value, method = "BH"),
    wilcox_p_value_adj = p.adjust(wilcox_p_value, method = "BH")
  )

# ✅ 11. Select significant variants (logFC ≥ 1 & adjusted p-value < 0.05)
selected_res <- results_diff %>%
  filter((abs(logFC) >= 1) &
           (wilcox_p_value_adj < 0.05) & (t_p_value_adj < 0.05))

selected <- selected_res$Variant

# ✅ 12. Prepare data for heatmap
df_selected <- df_wide[rownames(df_wide) %in% selected,]

# ✅ 13. Order rows by logFC (Descending order)
selected_res <- selected_res %>%
  arrange(desc(logFC))  # Sorting in descending order

ordered_variants <- selected_res$Variant
df_selected <- df_selected[ordered_variants, , drop = FALSE]  # Reorder rows

mat_selected <- as.matrix(df_selected)

# ✅ 14. Separate Control and Target for Clustering
control_cols <- colnames(mat_selected)[grepl("^control", colnames(mat_selected))]
target_cols <- colnames(mat_selected)[!grepl("^control", colnames(mat_selected))]

# Compute clustering order separately for Control and Target groups
control_dist <- dist(t(mat_selected[, control_cols]))
control_clust <- hclust(control_dist)

target_dist <- dist(t(mat_selected[, target_cols]))
target_clust <- hclust(target_dist)

# Arrange column order based on clustering results
ordered_cols <- c(control_cols[control_clust$order], target_cols[target_clust$order])
mat_selected <- mat_selected[, ordered_cols]

# ✅ 15. Add logFC values as row annotations
logFC_values <- selected_res %>%
  select(Variant, logFC) %>%
  column_to_rownames(var = "Variant")

# ✅ 16. Define annotation for heatmap
anno_col <- data.frame(Group = ifelse(grepl('^control', colnames(mat_selected)), 'control', 'target'))
rownames(anno_col) <- colnames(mat_selected)

anno_row <- logFC_values  # logFC values as row annotations

# Define heatmap colors
my_colors <- c("#929292", "#F1C40F", "#E74C3C")
my_breaks <- c(-0.5, 0.5, 1.5, 2.5)

# ✅ 17. Save heatmap as PDF
#pdf("heatmap_selected_variants_updated.pdf", width = 15, height = 5)  # Start PDF
pheatmap(mat_selected,
         annotation_col = anno_col,
         annotation_row = anno_row,  # Show logFC values next to heatmap
         cluster_rows = FALSE,  # No clustering for rows, ordered by logFC
         cluster_cols = FALSE,  # Column clustering disabled (since we pre-sorted)
         color = my_colors,
         breaks = my_breaks,
         scale = 'none',
         main = 'Heatmap of Selected Variants (Genotype)')
dev.off()  # Finish PDF

# ✅ 18. Save final results
write.csv(results_diff, "genotype_comparison_with_wilcoxon_updated.csv", row.names = FALSE)
write.csv(selected_res, "selected_variants_updated.csv", row.names = FALSE)

print("✅ Analysis completed! Results saved -> genotype_comparison_with_wilcoxon_updated.csv")
print("✅ Heatmap PDF saved -> heatmap_selected_variants_updated.pdf")
