rm(list = ls())

library(ggplot2)
library(reshape2)
library(dplyr)
library(ggdendro)
library(cowplot)
library(tidyverse)
library(tidyr)
library(gridExtra)


mlst <- read.csv("mlst.csv")
virulence <- read.csv("virulence.csv")

# select cols from mlst df needed for plot 
mlst_selected <- mlst %>%
  select(SampleID, SequenceType, aroE, gdh, gki, recP, spi, xpt, ddl) %>%
  mutate_at(vars(aroE, gdh, gki, recP, spi, xpt, ddl), ~ as.numeric(as.character(.)))

# select cols from viurlence df needed for plot 
virulence_binary <- virulence %>%
  select(SampleID, GENE) %>%
  mutate(presence = 1) %>%
  distinct(SampleID, GENE, .keep_all = TRUE) %>%
  spread(key = GENE, value = presence, fill = 0)

# merge mlst and virulence dfs to calculate distances for dendrogram
merged_data <- merge(mlst_selected, virulence_binary, by = "SampleID")

# find distances and create dendrogram
dist_matrix <- dist(merged_data[, -c(1, 2)])  # exclude sampleID and seqtype
hc <- hclust(dist_matrix)
dendro_data <- dendro_data(hc)

# plot dendrogram 
dendrogram_plot <- ggplot() +
  geom_segment(data = segment(dendro_data), aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank())

# prep virulence data for heatmap
virulence_long <- melt(virulence_binary, id.vars = "SampleID")

# plot heatmap
heatmap_plot <- ggplot(virulence_long, aes(x = SampleID, y = variable, fill = as.factor(value))) +
  geom_tile() +
  scale_fill_manual(values = c("0" = "yellow", "1" = "red"), labels = c("0" = "Absent", "1" = "Present"),name=NULL) +
  theme_minimal() +
  labs(x = "SampleID", y = "Virulence Gene") +
  theme(axis.text.x = element_text(hjust = 1), 
        panel.grid = element_blank())

# add in table for mlst data below plot
mlst_table <- mlst_selected %>%
  select(-SequenceType) %>%
  arrange(SampleID)

table_grob <- tableGrob(mlst_table, rows = NULL)

# combine heatmap and dendrogram 
combined_plot <- plot_grid(dendrogram_plot, heatmap_plot, ncol = 1, align = "v", rel_heights = c(0.5, 3))

# add in table to combined plot
final_plot <- plot_grid(combined_plot, table_grob, ncol = 1, rel_heights = c(3.5, 1.5))
print(final_plot)

# save and export 
ggsave("virulence_mlst_plot.png",width=350,height=350,units="mm", bg = "white")


