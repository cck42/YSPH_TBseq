rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggtree)
library(ape)
library(phytools)
library(patchwork)
library(forcats)

ampliconfile   <- "data/amplicon_positions.csv"
assemblyfile   <- "data/assembly_lengths.csv"
divergencefile <- "data/fastani.out"

outpng        <- "specificity_plot.png"
outpdf        <- "specificity_plot.pdf"

# plot amplicon mappings --------------------------------------------------

maxlen =  6.5e6 #4.45e6

assembly_lengths <- read.csv(assemblyfile, header = FALSE, col.names = c("strain", "start", "end")) %>%
  #filter(!strain %in% c("M.fortuitum", "M.intracellulare")) %>%
  arrange(end) %>%
  mutate(strain = case_when(
    strain == "OXC141" ~ "S. pneumoniae",
    strain == "H37Rv" ~ "M. tuberculosis", 
    strain == "PA01" ~ "P. aeruginosa",
    TRUE ~ strain
  )) %>%
  mutate(strain = factor(strain, levels = unique(strain), ordered = TRUE)) %>%
  mutate(y = as.numeric(strain) - 0.5)

strain_order <- levels(assembly_lengths$strain)

ampmap <- read.csv(ampliconfile, col.names = c("strain", "start", "end", "panel")) %>%
  #filter(!strain %in% c("M.fortuitum", "M.intracellulare")) %>%
  mutate(strain = case_when(
    strain == "OXC141" ~ "S. pneumoniae",
    strain == "H37Rv" ~ "M. tuberculosis", 
    strain == "PA01" ~ "P. aeruginosa",
    TRUE ~ strain
  )) %>%
  mutate(strain = factor(strain, levels = strain_order, ordered = TRUE)) %>%
  mutate(panel = sapply(panel, function(primer_type) {
    if (primer_type == "fwd") {
      return(1)
    } else if (primer_type == "rev") {
      return(2)
    } else {
      return(NA)
    }
  })) %>%  
  mutate(y = as.numeric(strain) - 1 + (panel - 1) * 0.5)

ampplot <- ggplot(ampmap, aes(xmin = start, xmax = end, 
                              ymin = y + 0.05, ymax = y + 0.45)) + 
  geom_rect() +
  geom_segment(data = assembly_lengths, aes(x = start, xend = end, y = y, yend = y)) + 
  scale_y_continuous(
    breaks = unique(assembly_lengths$y),
    labels = levels(assembly_lengths$strain),
    limits = c(0, length(levels(assembly_lengths$strain))),
    expand = c(0, 0)
  ) +
  scale_x_continuous(limits = c(0, maxlen), breaks = seq(0, maxlen, 5e5), labels = \(x) x / 1e6, expand = c(0, 0)) + 
  xlab("Genome Position (Mbp)") +
  theme(
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    #panel.background = element_blank(), panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Predicted Amplicon Coverage")
ampplot

ggsave(outpng, width=500,height=350,dpi=400,units="mm")
ggsave(outpdf, width=500,height=350,dpi=400,units="mm")
