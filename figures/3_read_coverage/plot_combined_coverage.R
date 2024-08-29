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


heattab <- read.table("data/SP_CZID_heatmap_240812.csv",sep=",",header=1) %>% 
    mutate(ID=gsub("_S.+","",sample_name,perl=T)) %>%
    rename(seq_id = sample_name)


covtabTB <- read.table("data/TB.coverage.tsv",sep="\t",header=T) %>% 
  rename_with(tolower)

covtabSP <- read.table("data/SP.coverage.tsv",sep="\t",header=T) %>% 
  rename_with(tolower)

# sample type comparison (SP) ---------------------------------------------


#matched samples
sampletype_samples <- c("Yale-SP00051", "Yale-SP00052", "Yale-SP00053", "Yale-SP00054",
                        "Yale-SP00055", "Yale-SP00056", "Yale-SP00057", "Yale-SP00058", "Yale-SP00059",
                        "Yale-CS00171", "Yale-CS00172", "Yale-CS00173", "Yale-CS00174", "Yale-CS00175",
                        "Yale-CS00176", "Yale-CS00162", "Yale-CS00163", "Yale-CS00164")


SPsampletab <- merge(heattab, metadataSP,by="seq_id", all.x = T) %>% 
                filter(seq_id %in% sampletype_samples)



cols = c("Streptococcus pneumoniae"="#006D2C","Streptococcus mitis"="#BDD7E7",
         "Streptococcus oralis"="#6BAED6","Streptococcus sp. oral taxon 061"="#3182BD",
         "Streptococcus sp. oral taxon 431"="#08519C")

# identify taxa that are not in the cols vector
non_streptococcus_taxa <- setdiff(unique(heattab$taxon_name), names(cols))
# create a vector of grey scale colors
grey_scale <- colorRampPalette(c("grey80", "grey20"))(length(non_streptococcus_taxa))
# assign grey scale colors to these taxa
grey_cols <- setNames(grey_scale, non_streptococcus_taxa)
# combine the grey scale colors with the existing cols vector
all_cols <- c(cols, grey_cols)


SPsampletab$taxon_name <- factor(SPsampletab$taxon_name, levels=rev(names(all_cols)), ordered=T)
SPsampletab$ngs_prep_method <- fct_recode(SPsampletab$ngs_prep_method,
                                       "no amplification"="Nextera XT (GLab)","amplicon"="COVIDseq")
SPsampletab$ngs_prep_method <- factor(SPsampletab$ngs_prep_method,levels=c("no amplification","amplicon"),ordered=T)




# plot
sampletypeplot <- ggplot(SPsampletab, aes(x = original_id, y = NT_rpm, fill = taxon_name)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        axis.title.x = element_blank()) + 
  ylab("aligned reads (NT)") +
  facet_grid(~sample_type + ngs_prep_method, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = all_cols, 
                    breaks = c("Streptococcus pneumoniae", "Streptococcus mitis",
                               "Streptococcus oralis", "Streptococcus sp. oral taxon 061",
                               "Streptococcus sp. oral taxon 431")) +  
  theme(axis.text.x = ggtext::element_markdown())

sampletypeplot



# sputum extraction (TB) --------------------------------------------------

sputumplot <- plot_spacer()



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
  geom_point() + facet_grid(species ~ .)

dilutionplot





# merge plots -------------------------------------------------------------

playout <- "AABD
CCCC"

sampletypeplot + sputumplot + dilutionplot + guide_area() + plot_layout(design = playout,guides="collect")

ggsave("2_combined_coverage_plot.png",dpi=400,units="mm",width=300,height=200)
