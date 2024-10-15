rm(list = ls())

library("tidyverse")
library(RColorBrewer)
library(patchwork)
library(ggtext)


# data import  ------------------------------------------------------------

metadataSP <- read.csv("data/S.pneumoniae.csv") %>% 
  select(c("Seq_ID", "Original_ID", "Sample_Type", 
           "Sample_Dilution", "Primer_Conc", 
           "NGS_Prep_Method", "Ct1")) %>% 
  rename("CT"="Ct1") %>%
  rename_with(tolower) %>%
  mutate("species"="S.pneumoniae")

metadataTB <- read.csv("data/M.tuberculosis.csv") %>% 
  rename("Sample_Dilution"="Template_dilution") %>% 
  select(c("Seq_ID", "Original_ID", "Sample_Type", 
           "Sample_Dilution",  "Primer_Conc",  
           "NGS_Prep_Method", "CT")) %>% 
  rename_with(tolower) %>%
  mutate("species"="M.tuberculosis")


heattabSP <- read.table("data/SP_CZID_heatmap_240812.csv",sep=",",header=1) %>% 
    mutate(ID=gsub("_S.+","",sample_name,perl=T)) %>%
    rename(seq_id = sample_name)

heattabTB <- read.table("data/TB_CZID_heatmap_241001.csv",sep=",",header=1) %>% 
  mutate(ID=gsub("_S.+","",sample_name,perl=T)) %>%
  rename(seq_id = sample_name)



covtabTB1 <- read.table("data/TB.coverage.tsv",sep="\t",header=T) %>% 
  rename_with(tolower)


covtabTB2 <- read.table("data/TB005-6_combinedcoverage.txt",sep="\t",header=T) %>% 
  rename_with(tolower) %>%
  rename(seq_id = sample)

covtabTB <- rbind(covtabTB1,covtabTB2)


covtabSP <- read.table("data/SP.coverage.tsv",sep="\t",header=T) %>% 
  rename_with(tolower) %>% 
  mutate(subsample = as.numeric(subsample)) %>% 
  filter(subsample==1.0)

#View(covtabSP)

# sample type comparison (SP) ---------------------------------------------

#matched samples
sampletype_samples <- c("Yale-SP00051", "Yale-SP00052", "Yale-SP00053", "Yale-SP00054",
                        "Yale-SP00055", "Yale-SP00056", "Yale-SP00057", "Yale-SP00058", "Yale-SP00059",
                        "Yale-CS00171", "Yale-CS00172", "Yale-CS00173", "Yale-CS00174", "Yale-CS00175",
                        "Yale-CS00176", "Yale-CS00162", "Yale-CS00163", "Yale-CS00164")


SPsampletab <- merge(heattabSP, metadataSP,by="seq_id", all.x = T) %>% 
                filter(seq_id %in% sampletype_samples)



cols = c("Streptococcus pneumoniae"="#006D2C","Streptococcus mitis"="#BDD7E7",
         "Streptococcus oralis"="#6BAED6","Streptococcus sp. oral taxon 061"="#3182BD",
         "Streptococcus sp. oral taxon 431"="#08519C")

# identify taxa that are not in the cols vector
non_streptococcus_taxa <- setdiff(unique(heattabSP$taxon_name), names(cols))
# create a vector of grey scale colors
grey_scale <- colorRampPalette(c("grey80", "grey20"))(length(non_streptococcus_taxa))
# assign grey scale colors to these taxa
grey_cols <- setNames(grey_scale, non_streptococcus_taxa)
# combine the grey scale colors with the existing cols vector
all_cols <- c(cols, grey_cols)


SPsampletab$taxon <- factor(SPsampletab$taxon_name, levels=rev(names(all_cols)), ordered=T)
SPsampletab$ngs_prep_method <- fct_recode(SPsampletab$ngs_prep_method,
                                       "no amplification"="Nextera XT (GLab)","amplicon"="COVIDseq")
SPsampletab$ngs_prep_method <- factor(SPsampletab$ngs_prep_method,levels=c("no amplification","amplicon"),ordered=T)




# plot
sampletypeplot <- ggplot(SPsampletab, aes(x = original_id, y = NT_rpm, fill = taxon)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        axis.title.x = element_blank()) + 
  ylab("aligned reads (NT rpm)") +
  facet_grid(~sample_type + ngs_prep_method, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = all_cols, 
                    breaks = c("Streptococcus pneumoniae", "Streptococcus mitis",
                               "Streptococcus oralis", "Streptococcus sp. oral taxon 061",
                               "Streptococcus sp. oral taxon 431")) +  
  theme(axis.text.x = ggtext::element_markdown())

sampletypeplot



# amplicon coverage colors ------------------------------------------------

ampcol <- scale_color_manual(values=c("no amplification"="blue","amplicon"="dark red"))
ampfil <- scale_fill_manual(values=c("no amplification"="blue","amplicon"="dark red"))


# SP amp vs unamp ---------------------------------------------------------

amptab = rbind(metadataSP,metadataTB) %>% 
  mutate(seq_id=gsub("Yale-","",seq_id)) %>% 
  #filter(sample_type=="NP swab") %>% 
  #filter(sample_type=="Saliva") %>% 
  filter(is.na(sample_dilution)) %>% 
  merge(covtabSP,by="seq_id") %>% 
  mutate(group=paste(original_id,ngs_prep_method))

amptab$ngs_prep_method <- fct_recode(amptab$ngs_prep_method,
                                     "no amplification"="Hybrid CovidSeq without amplification",
                                     "no amplification"="Nextera XT (GLab)",
                                     "amplicon"="Mpox/Hybrid CovidSeq",
                                     "amplicon"="COVIDseq",
                                     "amplicon"="Hybrid COVIDseq")

#View(amptab)

crosssamples <-  amptab %>% 
  filter(is.na(sample_dilution)) %>% 
  filter(sample_type %in% c("Saliva","NP swab","Culture Isolate")) %>% 
  filter(original_id %in% c("A889-NP", "B042-NP", "C677-NP",
                            "A889-S","B042-S","C677-S",
                            "W1527-CI","W3317-CI","W4034-CI")) %>% 
  pull(original_id)


amptabcf <- subset(amptab,original_id %in% crosssamples)  %>% 
#                 filter(!grepl("CS",seq_id)) #remove earlier CS run to avoid duplicates
                  filter(!seq_id %in% c("CS00171","CS00162","CS00172","CS00163","CS00173","CS00164"))


samporder <- amptabcf %>% 
  filter(ngs_prep_method=="no amplification") %>% 
  arrange(coverage) %>% 
  pull(original_id)

amptabcf$original_id <- factor(amptabcf$original_id,levels=samporder,ordered=T)
amptabcf$ngs_prep_method <- factor(amptabcf$ngs_prep_method,levels=c("no amplification","amplicon"),ordered=T)


spcovbox  <- ggplot(amptabcf,aes(x=ngs_prep_method,y=coverage,color=ngs_prep_method,group=ngs_prep_method)) + 
  geom_boxplot() + 
  geom_point() + 
  geom_line(aes(group=original_id),color="grey") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45,hjust=1),
        legend.position="none") + ampcol +
  ylim(0,100) +
  facet_grid(. ~ sample_type)

spcovbox


# dilution curves ---------------------------------------------------------

# spreps <- c("CS00030", "CS00047", "CS00048", "CS00049",
#           "CS00050", "CS00051", "SP00057", "SP00058", "SP00059")
# tbreps <- c("TB015","TB016","TB017", "TB018", "TB019", "TB020", "TB021", "TB022", "TB023", 
#           "TB090", "TB091", "TB092",
#           "TB093", "TB094", "TB095")

covtab <- rbind(covtabSP,covtabTB) %>% 
  mutate(seq_id=gsub("Yale-","",seq_id),
         subsample = as.numeric(subsample)) %>%
  filter(subsample==1)

metatab <- rbind(metadataSP,metadataTB) %>% 
  mutate(seq_id=gsub("Yale-","",seq_id)) %>%
  filter(primer_conc %in% c("100uM","200uM","0")) %>% 
  filter(!is.na(ct)) %>% 
  filter(!is.na(sample_dilution)) %>%
  filter(!original_id %in% c("111-2565-18"))

dilutiontab <- merge(metatab,covtab,by="seq_id") %>% mutate(group=paste(original_id,ngs_prep_method))

dilutiontab$ngs_prep_method <- fct_recode(dilutiontab$ngs_prep_method,
                                          "no amplification"="Hybrid CovidSeq without amplification",
                                          "amplicon"="Hybrid COVIDseq",
                                          "amplicon"="Mpox/Hybrid CovidSeq",)



dilutionplot <- ggplot(dilutiontab,aes(x=ct,y=coverage,group=group,color=ngs_prep_method)) + 
  geom_smooth(se = F,linewidth=0.3) + 
  ylab("genome coverage") + xlab("CT ( -> genome copies / ul)") +
  geom_point() + facet_grid(species ~ .) + ampcol

dilutionplot



dilutionplotTB <- ggplot(subset(dilutiontab,species=="M.tuberculosis"),aes(x=ct,y=coverage,group=group,color=ngs_prep_method)) + 
  geom_smooth(se = F,linewidth=0.3) + 
  ylab("genome coverage") + xlab("CT ( -> genome copies / ul)") +
  geom_point() + facet_grid(species ~ .) + ampcol

dilutionplotTB



dilutionplotSP <- ggplot(subset(dilutiontab,species=="S.pneumoniae"),aes(x=ct,y=coverage,group=group,color=ngs_prep_method)) + 
  geom_smooth(se = F,linewidth=0.3) + 
  ylab("genome coverage") + xlab("CT ( -> genome copies / ul)") +
  geom_point() + facet_grid(species ~ .) + ampcol

dilutionplotSP



# get mean cts for samples ------------------------------------------------

# meancttab <- rbind(metadataSP,metadataTB) %>% 
#   filter(!is.na(ct)) %>% 
#   filter(is.na(sample_dilution)) %>% 
#   filter(!sample_type %in% c("Culture isolate")) %>% 
#   group_by(species,sample_type) %>% 
#   summarize(mean = mean(ct,na.rm = T),
#             lq = quantile(ct,c(0.25),na.rm = T),
#             med = quantile(ct,c(0.5),na.rm = T),
#             uq = quantile(ct,c(0.75),na.rm = T),
#   )
# meancttab
# 
# dilutionplot2 <- ggplot(dilutiontab,aes(x=ct,y=coverage,group=group,color=ngs_prep_method)) + 
#   geom_rect(data=meancttab,aes(xmin=lq,xmax=uq,ymin=0,ymax=100),alpha=0.3, inherit.aes=F) +
#   geom_smooth(se = F,linewidth=0.3) + 
#   geom_point() + 
#   ylab("genome coverage") + xlab("CT ( -> genome copies / ul)") +
#   facet_grid(species ~ .)
# 
# dilutionplot2
# 


# TB amp vs unamp ---------------------------------------------------------

amptab = rbind(metadataSP,metadataTB) %>% 
         mutate(seq_id=gsub("Yale-","",seq_id)) %>% 
        filter(is.na(sample_dilution)) %>% 
        merge(covtab,by="seq_id") %>% 
        mutate(group=paste(original_id,ngs_prep_method))

amptab$ngs_prep_method <- fct_recode(amptab$ngs_prep_method,
                                       "no amplification"="Hybrid CovidSeq without amplification",
                                       "amplicon"="Mpox/Hybrid CovidSeq",)


# negsamples <-  amptab %>% 
#   filter(is.na(sample_dilution)) %>% 
#   filter(ngs_prep_method=="no amplification") %>% 
#                 pull(original_id)
# possamples = amptab %>% 
#   filter(is.na(sample_dilution)) %>% 
#   filter(ngs_prep_method=="amplicon") %>% 
#                 pull(original_id)
# crosssamples <- intersect(negsamples,possamples)

crosssamples <- c("Peru-41", "Peru-42", "Peru-43", "Peru-44", "Peru-45", "Peru-46", "Peru-47", "Peru-48", "Peru-49", "Peru-50")

amptabcf <- subset(amptab,original_id %in% crosssamples)


samporder <- amptabcf %>% 
                filter(ngs_prep_method=="no amplification") %>% 
                arrange(coverage) %>% 
                pull(original_id)

amptabcf$original_id <- factor(amptabcf$original_id,levels=samporder,ordered=T)


tbcovplot <- ggplot(amptabcf,aes(x=original_id,y=coverage,fill=ngs_prep_method)) + 
      geom_bar(stat="identity",position="dodge") + 
    theme(axis.text.x=element_text(angle=90),
          axis.title.x=element_blank(),
          legend.position="none") + ampfil

tbdepplot <- ggplot(amptabcf,aes(x=original_id,y=meandepth,fill=ngs_prep_method)) + 
  geom_bar(stat="identity",position="dodge") + 
  theme(axis.text.x=element_text(angle=90),
        axis.title.x=element_blank()) + ampfil

tbcovbox  <- ggplot(amptabcf,aes(x=ngs_prep_method,y=coverage,color=ngs_prep_method,group=ngs_prep_method)) + 
  geom_boxplot() + 
  geom_point() + 
  geom_line(aes(group=original_id),color="grey") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45,hjust=1),
        legend.position="none") + ampcol +
  ylim(0,100) +
  facet_grid(. ~ sample_type)




# sputum metagenomics (TB) ---------------------------------------------


#matched samples
sputum_samples <- c("Peru-41", "Peru-42", "Peru-43", "Peru-44", "Peru-45", "Peru-46", "Peru-47", "Peru-48", "Peru-49", "Peru-50")


TBsampletab <- merge(heattabTB, metadataTB,by="seq_id", all.x = T) %>% 
  filter(original_id %in% sputum_samples)

TBsampletab$taxon = paste(TBsampletab$genus_name,"spp.")
TBsampletab$taxon[TBsampletab$taxon_name=="Mycobacterium tuberculosis"] = "Mycobacterium tuberculosis"
TBsampletab$taxon[is.na(TBsampletab$taxon)] <- "NA"

#View(TBsampletab %>% group_by(taxon) %>% summarize(reads = sum(NR_r)))

cols = c("Mycobacterium tuberculosis"="#bd0026",
         "Mycobacterium spp."="#fdb462",
         "Schaalia spp."="#bebada",
         "Actinomyces spp."="#bc80bd",
         "Pseudomonas spp."="#ccebc5",
         "Neisseria spp."="#ffffb3",
         "Streptococcus spp."="#80b1d3",
         "Rothia spp."="#8dd3c7"
)
# identify taxa that are not in the cols vector
misctaxa <- setdiff(unique(TBsampletab$taxon), names(cols))
# create a vector of grey scale colors
grey_scale <- colorRampPalette(c("grey80", "grey20"))(length(misctaxa))
# assign grey scale colors to these taxa
grey_cols <- setNames(grey_scale, misctaxa)
# combine the grey scale colors with the existing cols vector
all_cols <- c(cols,grey_cols)


TBsampletab$taxon <- factor(TBsampletab$taxon, levels=rev(names(all_cols)), ordered=T)


TBsampletab$ngs_prep_method <- fct_recode(TBsampletab$ngs_prep_method,
                                          "no amplification"="Hybrid CovidSeq without amplification",
                                          "amplicon"="Mpox/Hybrid CovidSeq")
TBsampletab$ngs_prep_method <- factor(TBsampletab$ngs_prep_method,levels=c("no amplification","amplicon"),ordered=T)


TBsampletab$original_id <- factor(TBsampletab$original_id,levels=samporder,ordered=T)

# plot
sputumplot <- ggplot(TBsampletab, aes(x = original_id, y = NT_rpm, fill = taxon)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        axis.title.x = element_blank()) + 
  ylab("aligned reads (NT rpm)") +
  facet_grid(~sample_type + ngs_prep_method, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = all_cols, 
                    breaks = c(names(cols))) +  
  theme(axis.text.x=element_text(angle=45,hjust=1))

sputumplot





# merge plots -------------------------------------------------------------

playout <- "AADB
CCCC"

sampletypeplot + tbcovbox + dilutionplot + guide_area() + plot_layout(design = playout,guides="collect")
ggsave("2_combined_coverage_plot.png",dpi=400,units="mm",width=300,height=200)



sputumplot + tbcovbox + dilutionplotTB + guide_area() + plot_layout(design = playout,guides="collect")
ggsave("2_combined_coverage_plot_TB.png",dpi=400,units="mm",width=300,height=200)


playout <- "AADBB
CCCCC"

sampletypeplot + spcovbox + dilutionplotSP + guide_area() + plot_layout(design = playout,guides="collect")
ggsave("2_combined_coverage_plot_SP.png",dpi=400,units="mm",width=300,height=200)


spcomb <- sampletypeplot + spcovbox + guide_area() + plot_layout(design = "AACBB",guides="collect")
tbcomb <- sputumplot + tbcovbox + guide_area() + plot_layout(design = "AACB",guides="collect")
spcomb / tbcomb
ggsave("2_combined_coverage_plot_preprint.png",dpi=400,units="mm",width=300,height=200)

ggsave("2_combined_coverage_plot_preprint.svg",dpi=400,units="mm",width=300,height=200)


