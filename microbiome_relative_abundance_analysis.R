library(tidyverse)

control_data <- read.csv('control_91.csv', header = TRUE, stringsAsFactors = FALSE)
target_data <- read.csv('target_31.csv', header = TRUE, stringsAsFactors = FALSE)
sample_group <- read.csv("sample_group.csv", header = TRUE, stringsAsFactors = FALSE)

head(control_data)
head(target_data)

control_data[control_data == ''] <- 'unknown'
target_data[target_data == ''] <- 'unknown'

#control_data <- control_data %>% mutate(across(everything(), as.numeric, .names = 'cleaned_{col))
#target_data <- target_data %>% mutate(across(everthing(), as.numeric, .names = 'cleaned_{col}'))

common_taxonomy <- intersect(control_data$Species, target_data$Species)

control_data <- control_data %>% filter(Species %in% common_taxonomy) %>% arrange(Species)
target_data <- target_data %>% filter(Species %in% common_taxonomy) %>% arrange(Species)

row.names(control_data) <- control_data$Species
row.names(target_data) <- target_data$Species

control_aligned <- control_data[,-c(1:7)] %>% mutate_all(as.numeric)
target_aligned <- target_data[,-c(1:7)] %>% mutate_all(as.numeric)

head(control_aligned)
head(target_aligned)
dim(control_aligned)
dim(target_aligned)

all(row.names(control_aligned) == row.names(target_aligned))

combined_data <- cbind(control_aligned, target_aligned)

head(combined_data)
dim(combined_data)

head(sample_group)

library(reshape2)

combined_data$species <- row.names(combined_data)

combined_melt <- melt(combined_data)
combined_melt <- combined_melt %>%
  left_join(sample_group, by = c('variable' = 'Samples'))

colnames(combined_melt) <- c('species', 'samples', 'ratio', 'group')
head(combined_melt)
dim(combined_melt)
combined_melt$ratio <- combined_melt$ratio*100

statistical_results <- combined_melt %>%
  group_by(species) %>%
  summarise(
    mean_control = mean(ratio[group == "Control"], na.rm = TRUE),
    mean_target = mean(ratio[group == "Target"], na.rm = TRUE),
    p_value_t_test = t.test(
      ratio[group == "Control"],
      ratio[group == "Target"]
    )$p.value,
    p_value_wilcox = wilcox.test(
      ratio[group == "Control"],
      ratio[group == "Target"]
    )$p.value
  ) %>%
  mutate(
    logFC = log2((mean_target + 1e-8) / (mean_control + 1e-8))
  )

significant_species <- statistical_results %>%
  filter((p_value_t_test < 0.05 | p_value_wilcox < 0.05) & abs(logFC) > 1)

significant_species <- significant_species[order(significant_species$logFC, decreasing = TRUE),]

cat("Number of significant species after all cutoffs:", nrow(significant_species))

write.csv(significant_species, 't_test_wilcoxon_test.csv', row.names = FALSE)

library(pheatmap)

# 필터링된 데이터로 Heatmap 데이터 생성
heatmap_data <- combined_data[rownames(combined_data) %in% significant_species$species, ]
heatmap_data <- heatmap_data[match(significant_species$species, rownames(heatmap_data)), ]

head(heatmap_data)
# Ensure heatmap_data is numeric
heatmap_data$species <- NULL

# 데이터 정규화
heatmap_data <- t(scale(t(heatmap_data)))

annotation_col <- sample_group %>%
  filter(Samples %in% colnames(heatmap_data)) %>%
  column_to_rownames("Samples")

annotation_row <- significant_species %>%
  select(species, logFC) %>%
  column_to_rownames('species')

head(annotation_col)
head(annotation_row)
heatmap_data <- heatmap_data[, colnames(heatmap_data) %in% rownames(annotation_col)]


# Heatmap 생성
pheatmap(
  heatmap_data,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Heatmap of Significant Species",
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = list(logFC = colorRampPalette(c("orange", "gray", "green"))(50)),
  fontsize_row = 8,
  fontsize_col = 8,
)

sig_species <- significant_species$species
head(sig_species)

filtered_data <- combined_melt %>%
  filter(species %in% sig_species)

head(filtered_data)

cv_results <- filtered_data %>%
  group_by(species, group) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE),
    cv_ratio = ifelse(mean_ratio > 0, sd_ratio / mean_ratio, NA)
  ) %>%
  pivot_wider(names_from = group, values_from = c(mean_ratio, sd_ratio, cv_ratio))
head(cv_results)
write.csv(cv_results, 'cv_results.csv', row.names = FALSE)

consistent_species <- cv_results %>%
  filter(cv_ratio_Control < 1.5 & cv_ratio_Target < 1.5)

cat("Number of consistent species:", nrow(consistent_species), "\n")

consistent_species

ggplot(filtered_data, aes(x = group, y = ratio, fill = group)) +
  geom_boxplot() +
  facet_wrap(~ species, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Group-wise Variability for Significant Species",
    x = "Group",
    y = "Relative Abundance (%)"
  )

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("KEGGREST")
library(KEGGREST)

# 각 미생물에 대한 KEGG 정보를 검색
kegg_results <- lapply(significant_species, function(species) {
  keggFind("genome", species)
})

# 결과 확인
kegg_results

devtools::install_github("dami82/easyPubMed")
library(easyPubMed)
species_query <- "Bifidobacterium adolescentis"
pubmed_results <- get_pubmed_ids(species_query)

# Fetch the articles
article_data <- fetch_pubmed_data(pubmed_results)