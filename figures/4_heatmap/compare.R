library(tidyverse)
library(readr)
library(cowplot)
library(rlang)
library(forcats)


tb003_mykrobe <- read.csv('/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/figures/4_heatmap/data/TB003_mykrobe.csv')
tb004_mykrobe <- read.csv('/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/figures/4_heatmap/data/TB004_mykrobe.csv')
tb005_mykrobe <- read.csv('/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/figures/4_heatmap/data/TB005_mykrobe.csv')
sample_metadata <-read.csv('/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/metadata/M.tuberculosis.csv',
                           na.strings=c("","NaN","NA"))


all_mykrobe <- tb003_mykrobe %>%
  bind_rows(tb004_mykrobe,tb005_mykrobe) 
mykrobe_combined <- all_mykrobe %>%
  left_join(sample_metadata, by=c(sample='Seq_ID')) %>%
  filter(Sample_source %in% c("Cayetano","Moldova"))


heatmap_res_DNA <- mykrobe_combined %>%
  filter(NGS_Run_ID=="TB003") %>%
  filter(Sample_Type=="DNA") %>%
  filter(Primer_Conc!='0')%>%
  #complete(drug,sample) %>%
  filter(!is.na(drug)) %>%
  mutate('Template_dilution'=as.character(format(Template_dilution,scientific=TRUE))) %>%
  ggplot() +
  geom_tile(aes(x=Template_dilution,y=drug,fill=susceptibility)) + 
  scale_fill_manual(
    values=c("#FFC561","#B20026","#006D2C"),
    name="Susceptibility",
    labels=c("Partial resistance","Resistant","Susceptible"),
    na.value='gray'
  ) +
  labs(x='Sample',
       y='Drug') +
  theme_half_open()+
  theme(axis.text.x=element_text(angle=-90,vjust=0.5,hjust=0.1),
        panel.background = element_rect(fill='gray'),
        legend.position='none') +
  facet_wrap("Original_ID",strip.position="bottom",nrow=1)
plot(heatmap_res_DNA)


heatmap_res_sputum <- mykrobe_combined %>%
  filter(Sample_Type=="Sputum extract") %>%
  complete(drug,sample) %>%
  filter(!is.na(drug)) %>%
  ggplot()+
  geom_tile(aes(x=sample,y=drug,fill=susceptibility))+ 
  scale_fill_manual(
    values=c("#FFC561","#B20026","#006D2C"),
    name="Susceptibility",
    labels=c("Partial resistance","Resistant","Susceptible"),
    na.value='gray'
  ) +
  labs(x='Sample',
       y='Drug') +
  theme_half_open()+
  theme(axis.text.x=element_text(angle=-90,vjust=0.5,hjust=0.1)) 
plot(heatmap_res_sputum)

reslegend <- get_legend(
  heatmap_res_sputum+
    theme(legend.box.margin = margin(0,0,0,0),
                           legend_position="bottom")
)

heatmap_res <- plot_grid(heatmap_res_DNA,
                         heatmap_res_sputum+theme(legend.position='none'), 
                         #labels=c("DNA","Sputum extracts"),
                         #label_y=1.05,
                         labels='AUTO',
                         ncol=1,align='h')
heatmap_res
heatmap_res_legend <- plot_grid(heatmap_res,reslegend,rel_widths = c(10,1.5))
heatmap_res_legend

save_plot("/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/figures/4_heatmap/Fig4_heatmap.png",
          heatmap_res_legend,
          base_height=7.75,base_width=12.15)


#Now I'm going to reformat the metadata a bit so it works nicely in nextstrain
TB004_coverage <- read.delim('/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/figures/4_heatmap/data/TB004_coverage.tsv')
TB005_coverage <-read.delim('/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/figures/4_heatmap/data/TB005_coverage.tsv')
TB003_coverage <- read.delim('/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/figures/3_read_coverage/data/TB.coverage.tsv')
all_coverage <- TB003_coverage %>%
  rename(sample=Seq_ID) %>%
  bind_rows(TB004_coverage,TB005_coverage) %>%
  filter(subsample==1.0)


excluded_runs<-c('TB001','TB002')
olddrugs <- c('Isoniazid','Streptomycin','Rifampicin','Ethambutol','Ethionamide','Pyrazinamide','Ofloxacin','Kanamycin')
newdrugs <- c('INH','STR','RIF','EMB','ETH','PZA','OFL','KAN')
metadata_ns <- sample_metadata %>%
  left_join(all_coverage,by=c(Seq_ID="sample")) %>% #combine with coverage data
  full_join(all_mykrobe,by=c(Seq_ID="sample" )) %>% #combine with mykrobe data
  filter(!(NGS_Run_ID %in% excluded_runs)) %>% #Remove some sequencing runs
  filter(Sample_Type!='Water') %>% #remove NTCs
  select(Seq_ID,drug,Original_ID,Sample_Type,Sample_source,collection_date,CT,starting_quant,Primer_Conc,
         NGS_Run_ID,numreads,covbases,coverage,meandepth,susceptibility,lineage)%>%
  complete(drug,Seq_ID) %>% #Explicitly write out NA for anything where drug susceptibility was not possible
  pivot_wider(names_from=drug,values_from=susceptibility) %>%
  filter(!is.na(Sample_source)) %>%
  filter(coverage>=75 & meandepth>10) %>%
  rename_with(~newdrugs[which(olddrugs==.x)],.cols=olddrugs) %>% #rename headers to abbreviations to match up
  mutate(across(all_of(newdrugs),.fns=~fct_collapse(.,'0' = c("S"),'1' = c("r","R")))) %>% #recode to 0 (no res) or 1 (any res) 
  rename(mykrobe_lineage=lineage) %>%
  mutate(continent=ifelse(Sample_source=="Moldova",'Europe','South_America'))%>%
  mutate(region=ifelse(Sample_source=="Moldova",'Eastern_Europe','South_America'))%>%
  mutate(country=ifelse(Sample_source=="Moldova",'Moldova','Peru')) %>%
  mutate(collection_date=ifelse(!is.na(collection_date),collection_date,'20XX-XX-XX')) %>%
  mutate(study="Glab")

#COmbined this with the reference nextstrain metadata
context_metadata <- read.delim('/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/pipelines/nextstrain/data/input_metadata.tsv') %>%
  mutate(across(all_of(newdrugs),as.factor))
merged_metadata_ns <- context_metadata %>%
  full_join(metadata_ns,by=c(strain='Seq_ID',date='collection_date','EMB','ETH','INH','KAN','OFL','PZA','RIF','STR',
                             'continent','region','country','study'))

#check these
glab_merged <- merged_metadata_ns %>%
  filter(country=="Peru")

write.csv(merged_metadata_ns,file='/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/pipelines/nextstrain/data/metadata_merged.csv')


