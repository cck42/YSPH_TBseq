library(tidyverse)
library(ggplot2)
library(ggtree)
library(ape)
library(phytools)
library(patchwork)

SP_treefile       <- "./data/SP_ml_tree.newick"
SP_ampliconfile   <- "./data/SP_amplicon_positions.csv"
SP_assemblyfile   <- "./data/SP_assembly_lengths.csv"
SP_divergencefile <- "./data/SP_divergence_fastani.out"
SP_pangenomefile  <- "./data/SP_gene_presence_absence.csv"
TB_treefile       <- "./data/TB_ml_tree.newick"
TB_ampliconfile   <- "./data/TB_amplicon_positions.csv"
TB_assemblyfile   <- "./data/TB_assembly_lengths.csv"
TB_divergencefile <- "./data/TB_divergence_fastani.out"
TB_pangenomefile  <- "./data/TB_gene_presence_absence.csv"

SP_treefile       <- "./data/supp/SP_ml_tree.newick"
SP_ampliconfile   <- "./data/supp/SP_amplicon_positions.csv"
SP_assemblyfile   <- "./data/supp/SP_assembly_lengths.csv"
SP_divergencefile <- "./data/supp/SP_divergence_fastani.out"
SP_pangenomefile  <- "./data/supp/SP_gene_presence_absence.csv"
TB_treefile       <- "./data/supp/TB_ml_tree.newick"
TB_ampliconfile   <- "./data/supp/TB_amplicon_positions.csv"
TB_assemblyfile   <- "./data/supp/TB_assembly_lengths.csv"
TB_divergencefile <- "./data/supp/TB_divergence_fastani.out"
TB_pangenomefile  <- "./data/supp/TB_gene_presence_absence.csv"

outpng        <- "./Fig2_combined_pangenome_plot.png"
outpdf        <- "./Fig2_combined_pangenome_plot.pdf"



# tree plot ---------------------------------------------------------------

get_tree_plot <- function(treefile) {
  
  # load in tree
  tree <-  read.tree(treefile)
  tree$tip.label <- gsub("_reference", "", tree$tip.label)
  
  # root with midpoint
  tree <- midpoint_root(tree)
  
  treesize = max(tree$edge.length)
  
  treeplot <- ggtree(tree,ladderize=F) +
    #geom_nodelab(aes(label=""), size=1, nudge_x=-0.01, nudge_y=0.25) +
    geom_tiplab(align=T, offset = treesize*0.01, linetype="dotted", 
                linesize = 0.4,size=3) + 
    scale_x_continuous(expand=c(0,treesize*0.35))
  treeplot
}

get_tree_order <- function(treefile) {
  # load in tree
  tree <-  read.tree(treefile)
  tree$tip.label <- gsub("_reference", "", tree$tip.label)
  tree <- midpoint_root(tree)
  #get tree order 
  istip <- tree$edge[,2]<=length(tree$tip.label)
  treeorder <- tree$tip.label[tree$edge[istip,2]]
  treeorder
}


# pangenome representation ------------------------------------------------


get_pangenome_plot <- function(pangenomefile, divergencefile, treefile) {
  mindivpc=90
  treeorder <- get_tree_order(treefile)
  # load in matrix data
  pangenome <- read.csv(pangenomefile, header = T)
  
  pangenome <- pangenome %>% 
    select(!c(2,3,6,7:14)) %>% 
    rename_with(function(x) {gsub("_reference","",x)}) %>% 
    rename_with(function(x) {gsub("\\.\\.","_",x)}) %>% 
    rename_with(tolower,c(1:3)) %>% 
    mutate(across(!c(1:3), function(x) {x != ""})) %>%
    mutate(gindex = seq.int(nrow(pangenome))) %>% 
    relocate(gindex) %>% 
    pivot_longer(!c(1:4),names_to = "strain",values_to="present") %>% 
    mutate(strain = factor(strain,ordered=T,levels=treeorder))
    
  
  # read in divergence / fastani output
  divtable <- read.table(divergencefile, sep="\t",header = F,
                         col.names = c("strain","ref","identity","fragments","matches")) %>% 
                    mutate(across(c(strain,ref),function(x) {gsub(".fasta","",gsub(".*\\/","",x,perl=T))})) %>%
                    mutate(across(c(strain,ref),function(x) {gsub(".fa","",gsub(".*\\/","",x,perl=T))})) %>%
                    mutate(across(c(strain,ref),function(x) {gsub("_reference","",x)}))
  
  
  # plot pangenome / divergence heatmap 
  pangenome <- merge(pangenome,divtable)
  panplot <- ggplot(subset(pangenome,present),aes(x=gindex,y=strain,fill=identity)) + 
    geom_tile() + 
    scale_fill_distiller(limits=c(mindivpc,100),direction=1,palette="YlOrRd",name="seq identity") + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_discrete(expand=c(0,0)) + 
    theme(legend.position="bottom",
          plot.background=element_blank(),
          panel.background=element_blank(),
          panel.grid = element_blank(),
          axis.title=element_blank(),axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          axis.text.y=element_blank(),
          plot.title = element_text(hjust = 0.5))+
    ggtitle("pangenome representation")
  panplot
}

get_amplicon_plot <- function(assemblyfile, ampliconfile, treefile) {
  # plot amplicon mappings --------------------------------------------------
  
  maxlen = 4.45e6
  
  treeorder <- get_tree_order(treefile)
  
  assembly_lengths <- read.csv(assemblyfile, header = FALSE, col.names = c("strain", "start", "end")) %>%
    mutate(strain = factor(strain, ordered = TRUE, levels = treeorder)) %>%
    mutate(y = as.numeric(strain) - 0.5)
  
  ampmap <- read.csv(ampliconfile, col.names = c("strain", "start", "end", "panel")) %>%
    mutate(strain = factor(strain, ordered = TRUE, levels = treeorder)) %>%
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
  
  ampplot <- ggplot(ampmap, aes(xmin = start, xmax = end, ymin = y + 0.05, ymax = y + 0.45)) + 
    geom_rect() + 
    geom_segment(data = assembly_lengths, aes(x = start, xend = end, y = y, yend = y)) + 
    scale_y_continuous(breaks = seq(0.5, length(treeorder), 1), labels = treeorder,
                       limits = c(0, length(treeorder)), expand = c(0, 0)) + 
    scale_x_continuous(limits = c(0, maxlen), breaks = seq(0, maxlen, 5e5), labels = \(x) x / 1e6, expand = c(0, 0)) + 
    xlab("genome position (mb)") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(), panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle("amplicon coverage")
  
  return(ampplot)
}

# assemble combined plots -------------------------------------------------


SP_treeplot <- get_tree_plot(SP_treefile)
SP_panplot <- get_pangenome_plot(SP_pangenomefile, SP_divergencefile, SP_treefile)
SP_ampplot <- get_amplicon_plot(SP_assemblyfile, SP_ampliconfile, SP_treefile)

TB_treeplot <- get_tree_plot(TB_treefile)
TB_panplot <- get_pangenome_plot(TB_pangenomefile, TB_divergencefile, TB_treefile)
TB_ampplot <- get_amplicon_plot(TB_assemblyfile, TB_ampliconfile, TB_treefile)

SP_treeplot + SP_panplot + SP_ampplot + 
TB_treeplot + TB_panplot + TB_ampplot + 
  plot_layout(widths=c(2.5,3,3),ncol=3)

ggsave(outpng, width=500,height=500,dpi=400,units="mm")
ggsave(outpdf, width=500,height=500,dpi=400,units="mm")


