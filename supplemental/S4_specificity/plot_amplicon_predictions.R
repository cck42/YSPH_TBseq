library(tidyverse )
library(patchwork)

filelist <- list.files("./",pattern="TBseq_.*_p.")
#use do.call and rbind to read in all files and combine them into one data frame
amplicons <- do.call(rbind, lapply(filelist, function(x) {read.csv(x, header=TRUE, sep="\t") %>% mutate(filename=x)})) %>% 
  separate(filename, into=c("prefix", "genus","spp", "pool"), sep="_", remove=FALSE) %>%
  mutate(strain = paste(genus, spp, sep=" ")) %>%
  select(-prefix, -genus, -spp, -filename) %>%
  mutate(mismatches = pmax(p1mismatches,p2mismatches),
         summismatches = p1mismatches+p2mismatches)

amplicons


spporder = c("P aeruginosa","K pneumoniae","M canettii","M tuberculosis","S aureus","S odontolytica","S pneumoniae","H influenzae")


assembly_lengths <- read.table("reflengths.txt", header = FALSE, col.names = c("strain","length", "refid")) %>%
  mutate(strain = gsub("_"," ",strain)) %>%
  mutate(strain = factor(strain, ordered = TRUE, levels = spporder)) %>%
  mutate(y = as.numeric(strain) - 0.5)

ampmap <- merge(amplicons,assembly_lengths,by="strain") %>% 
  mutate(y=ifelse(pool=="p2",y-0.5,y))

maxlen=6.5e6
  
ampplot <- ggplot(subset(ampmap,mismatches <= 3 & summismatches <=5), aes(xmin = amplicon_start, xmax = amplicon_end, ymin = y + 0.05, ymax = y + 0.45)) + 
  geom_rect(color="grey",linewidth = 0.2) + 
  geom_segment(data = assembly_lengths, aes(x = 0, xend = length, y = y, yend = y),inherit.aes=F) + 
  scale_y_continuous(breaks = seq(0.5, length(spporder), 1), labels = spporder,
                     limits = c(0, length(spporder)), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, maxlen), breaks = seq(0, maxlen, 1e6), labels = \(x) x / 1e6, expand = c(0, 0)) + 
  xlab("genome position (mb)") +
  theme(axis.title.y = element_blank(),
        panel.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

ampplot
ggsave("amplicon_coverage.png", ampplot, width=10, height=6)





# read in depth files -----------------------------------------------------
aligndepth <- rbind(
                read.table("S_odont_align/TB166_depth.txt",col.names = c("ref","pos","depth")) %>% mutate(sample="TB seq"),
                read.table("S_odont_align/TB205_depth.txt",col.names = c("ref","pos","depth")) %>% mutate(sample="unamplified")) %>%
                mutate(sample=factor(sample, levels=c("unamplified","TB seq")))
#average across 1kb blocks
bsize=5000
aligndepth1k <- aligndepth %>% 
  mutate(pos = (floor(pos/bsize)*bsize)+(bsize/2)) %>%
  group_by(sample,ref,pos) %>% 
  summarize(depth = round(mean(depth)))
  
#make depth plots
depplot <- ggplot(aligndepth1k,aes(x=pos,y=depth)) + geom_bar(stat="identity") + 
  scale_x_continuous(limits = c(0, 2.5e6), breaks = seq(0, 2.5e6, 5e5), labels = \(x) x / 1e6, expand = c(0, 0)) + 
  facet_grid(sample ~ .,scale="free_y") + theme(axis.text.x=element_blank(),
                                                axis.ticks.x=element_blank(),
                                                axis.title.x=element_blank())
depplot

# just S. odontolytica -----------------------------------------------------
SOtab <- subset(ampmap,strain=="S odontolytica") %>% mutate(y=y-5)
SOlen <- subset(assembly_lengths,strain=="S odontolytica") %>% mutate(y=y-5)


SOplot <- ggplot(subset(SOtab,mismatches<=3 & summismatches<=5), aes(xmin = amplicon_start, xmax = amplicon_end, ymin = 0.05, ymax = 0.45)) + 
  geom_rect(alpha=0.3) + 
  geom_segment(data = subset(SOlen,strain=="S odontolytica"), aes(x = 0, xend = length, y = 0, yend = 0),inherit.aes=F) + 
  scale_x_continuous(limits = c(0, 2.5e6), breaks = seq(0, 2.5e6, 5e5), labels = \(x) x / 1e6, expand = c(0, 0)) + 
  scale_y_continuous(breaks=c(0.2),labels="predicted\namplicons")+
  xlab("genome position (mb)") +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

depplot + SOplot + plot_layout(heights=c(19,1))

ggsave("S_odont_depth_amplicons.png", width=10, height=6)


dust <- read.table("S_odontolytica.dust",col.names=c("chr","start","end"))

dplot <- ggplot(dust, aes(xmin = start, xmax = end, ymin = 0.05, ymax = 0.45)) + 
  geom_rect(alpha=1,fill="black",color="black") + 
  scale_x_continuous(limits = c(0, 2.5e6), breaks = seq(0, 2.5e6, 5e5), labels = \(x) x / 1e6, expand = c(0, 0)) + 
  scale_y_continuous(breaks=c(0.2),labels="low\ncomplexity")+
  xlab("genome position (mb)") +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))


depplot + SOplot + dplot + plot_layout(heights=c(9,1,1))
ggsave("S_odont_depth_amplicons_dust.png", width=10, height=6)

# 
# lims=c(0.22e6,0.24e6)
# (depplot + xlim(lims[1],lims[2])) +  (SOplot + xlim(lims[1],lims[2]))  + plot_layout(heights=c(9,1))
# lims=c(0.365e6,0.39e6)
# (depplot + xlim(lims[1],lims[2])) +  (SOplot + xlim(lims[1],lims[2]))  + plot_layout(heights=c(9,1))
# lims=c(1.42e6,1.46e6)
# (depplot + xlim(lims[1],lims[2])) +  (SOplot + xlim(lims[1],lims[2]))  + plot_layout(heights=c(9,1))
# lims=c(2.23e6,2.24e6)
# (depplot + xlim(lims[1],lims[2])) +  (SOplot + xlim(lims[1],lims[2]))  + plot_layout(heights=c(9,1))
