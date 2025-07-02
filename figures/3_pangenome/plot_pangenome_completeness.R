library(tidyverse)
library(ggplot2)
library(ggtree)
library(ape)
library(phytools)
library(patchwork)

# treefile       <- "./data/SP_ml_tree.newick"
# ampliconfile   <- "./data/SP_amplicon_positions.csv"
# assemblyfile   <- "./data/SP_assembly_lengths.csv"
# divergencefile <- "./data/SP_divergence_fastani.out"
# pangenomefile  <- "./data/SP_gene_presence_absence.csv"
# outfile        <- "./data/SP_pangenome_plot.png"

treefile       <- "./data/TB_ml_tree.newick"
ampliconfile   <- "./data/TB_amplicon_positions.csv"
assemblyfile   <- "./data/TB_assembly_lengths.csv"
divergencefile <- "./data/TB_divergence_fastani.out"
pangenomefile  <- "./data/TB_gene_presence_absence.csv"
outfile        <- "./data/TB_pangenome_plot.png"



#outgroup="JYGP01"


# tree plot ---------------------------------------------------------------

# load in tree
tree <-  read.tree(treefile)
tree$tip.label <- gsub("_reference", "", tree$tip.label)

# root with midpoint
#tree <- ladderize(root(tree,outgroup))
tree <- midpoint_root(tree)

treesize = max(tree$edge.length)

treeplot <- ggtree(tree,ladderize=F) +
  #geom_nodelab(aes(label=""), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_tiplab(align=T, offset = treesize*0.01, linetype="dotted", 
              linesize = 0.4,size=3) + 
  scale_x_continuous(expand=c(0,treesize*0.35))
treeplot

#get tree order 
istip <- tree$edge[,2]<=length(tree$tip.label)
treeorder <- tree$tip.label[tree$edge[istip,2]]

# pangenome representation ------------------------------------------------

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
  scale_fill_gradient(high="darkred",low="white",limits=c(80,100),name="seq identity") + 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) + 
  theme(legend.position="bottom",
        panel.background=element_blank(),
        panel.grid = element_blank(),
        axis.title=element_blank(),axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("pangenome representation")
panplot


# plot amplicon mappings --------------------------------------------------
assembly_lengths <- read.csv(assemblyfile,header=F,col.names=c("strain","start","end")) %>%
  mutate(strain = factor(strain,ordered=T,levels=treeorder)) %>%
  mutate(y = as.numeric(strain)-0.5)

ampmap <- read.csv(ampliconfile,col.names = c("strain","start","end","panel")) %>%
              mutate(strain = factor(strain,ordered=T,levels=treeorder)) %>%
              mutate(y = as.numeric(strain)-1 + (panel-1)*0.5)

ampplot <- ggplot(ampmap,aes(xmin=start,xmax=end,ymin=y+0.05,ymax=y+0.45)) + geom_rect() + 
  geom_segment(data=assembly_lengths,aes(x=start,xend=end,y=y,yend=y)) + 
  scale_y_continuous(breaks=seq(0.5,length(treeorder),1),labels=treeorder,
                     limits=c(0,length(treeorder)),expand=c(0,0)) + 
  xlab("mb")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("amplicon coverage")


# assemble combined plots -------------------------------------------------


treeplot + panplot + ampplot + plot_layout(widths=c(2,3,3))
ggsave(outfile, width=350,height=150,dpi=400,units="mm")
