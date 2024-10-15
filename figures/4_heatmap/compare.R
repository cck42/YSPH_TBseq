library(tidyverse)
library(readr)
library(cowplot)
library(rlang)
library(forcats)
library(readxl)
library(stringr)


tb003_mykrobe <- read.csv('data/TB003_mykrobe.csv')
tb004_mykrobe <- read.csv('data/TB004_mykrobe.csv')
tb005_mykrobe <- read.csv('data/TB005_mykrobe.csv')
sample_metadata <-read.csv('data/M.tuberculosis.csv',
                           na.strings=c("","NaN","NA"))
sputum_extracts_summary <- read_excel('data/sputum_extraction_summary.xlsx',sheet=1) %>%
  mutate("Heat_inactivation"=as_factor(Heat_inactivation))%>%
  mutate("Chemical_inactivation"=as_factor(Chemical_inactivation))%>%
  mutate("liquefaction"=as_factor(liquefaction))%>%
  mutate("fastprep"=as_factor(fastprep))%>%
  mutate("pt_sample_num"=as_factor(pt_sample_num))

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
    values=c("#FFC561","#B20026","gray"),
    name="Susceptibility",
    labels=c("Partial resistance","Resistant","Susceptible"),
    na.value='#555555'
  ) +
  labs(x='Sample',
       y='Drug') +
  theme_half_open()+
  theme(axis.text.x=element_text(angle=-90,vjust=0.5,hjust=0.1),
        panel.background = element_rect(fill='#555555')) +
  facet_wrap("Original_ID",strip.position="bottom",nrow=1)
plot(heatmap_res_DNA)

save_plot("supp_dilutions-heatmap.svg",heatmap_res_DNA,base_width=12.5,base_height=4)

mykrobe_sputum_combined <- mykrobe_combined %>%
  right_join(sputum_extracts_summary, by=c('sample'="Yale-ID"))

heatmap_res_sputum_small <- mykrobe_sputum_combined %>%
  #filter(Sample_Type=="Sputum extract") %>%
  #complete(drug,sample,pt_sample_num) %>%
  filter(!is.na(drug)) %>%
  filter(liquefaction=="NALC 0.5% - NaOH"&fastprep=="40s") %>%
  ggplot()+
  geom_tile(aes(x=ID,y=drug,fill=susceptibility,height=0.95,width=0.95))+ 
  scale_fill_manual(
    values=c("#B20026","gray"),
    name="Susceptibility",
    labels=c("Resistant","Susceptible"),
    na.value='gray'
  ) +
  labs(x='Sample',
       y='Drug') +
  theme_half_open()+
  theme(axis.title=element_blank(),
        axis.text.x=element_text(angle=45,vjust=1,hjust=1),
        panel.background = element_rect(fill='white')) 
plot(heatmap_res_sputum_small)

ggsave("Fig4_heatmap.png",width=220,height=150,units="mm")

save_plot("Fig4_heatmap.svg",heatmap_res_sputum_small)


heatmap_res_sputum <- mykrobe_sputum_combined %>%
  #filter(Sample_Type=="Sputum extract") %>%
  #complete(expand(.,drug,sample)) %>%
  filter(!is.na(drug)) %>%
  #filter(!is.na(sample)) %>%
  ggplot()+
  geom_tile(aes(x=pt_sample_num,y=drug,fill=susceptibility))+ 
  scale_fill_manual(
    values=c("#FFC561","#B20026","grey"),
    name="Susceptibility",
    labels=c("Partial resistance","Resistant","Susceptible"),
    na.value='#555555'
  ) +
  labs(x='Sample',
       y='Drug') +
  theme_half_open()+
  theme(axis.text.x=element_text(angle=-90,vjust=0.5,hjust=0.1),
        panel.background = element_rect(fill='#555555')) +
  #facet_wrap(c('liquefaction','fastprep'),dir="v",shrink=FALSE,nrow=2,drop=TRUE,scales="free_x")
  facet_grid(extraction ~ batch, scales="free",switch="y",labeller=label_context) +
  theme(strip.text.y=element_blank(),
        panel.spacing.y=unit(2,"lines"))
plot(heatmap_res_sputum)
save_plot("heatmap_res_sputum_all.pdf",heatmap_res_sputum,
          base_width=12,base_height=8)

reslegend <- get_legend(
  heatmap_res_sputum+
    theme(legend.box.margin = margin(0,0,0,50))
)



#calculate mean and SD (SQ and coverage) for each extraction method
extr_stats <- sputum_extracts_summary %>%
  group_by(group) %>%
  summarize(meanquant=mean(starting_quant),std_dev_quant=sd(starting_quant),std_err_quant=(sd(starting_quant)/sqrt(n())),
            meancov=mean(coverage),std_dev_cov=sd(coverage),std_err_cov=(sd(coverage)/sqrt(n()))) %>%
  rename("method"=group)
write_tsv(extr_stats,"extraction_eff_by_method.tsv")

#SQ by extraction type
sq_median <- sputum_extracts_summary %>%
  ggplot() +
  geom_boxplot(aes(x=group,y=starting_quant,fill=group))+
  theme_half_open() +
  scale_fill_brewer(palette="Dark2",
    name="Extraction \nmethod"
  ) +
  scale_y_log10()+
  labs(x="Method",
       y="Starting quantity (GE/uL)")
plot(sq_median)

sq_avg <- extr_stats %>%
  ggplot() +
  geom_pointrange(aes(x=method,y=meanquant,ymin=meanquant-std_err_quant,ymax=meanquant+std_err_quant,color=method))+
  theme_half_open() +
  scale_color_brewer(palette="Dark2",
                    name="Extraction method"
  ) +
  scale_y_log10(limits= c(1,10000000))+
  labs(x="Method",
       y="Starting quantity (GE/uL)")
plot(sq_avg)

#Coverage by extraction type
cov_median <- sputum_extracts_summary %>%
  ggplot() +
  geom_boxplot(aes(x=group,y=coverage,fill=group))+
  theme_half_open() +
  scale_fill_brewer(palette="Dark2",
                    name="Extraction \nmethod"
  ) +
  labs(x="Method",
       y="Coverage (%)")
plot(cov_median)

cov_avg <- extr_stats %>%
  ggplot() +
  geom_pointrange(aes(x=method,y=meancov,ymin=meancov-std_err_cov,ymax=meancov+std_err_cov,color=method))+
  theme_half_open() +
  scale_color_brewer(palette="Dark2",
                     name="Extraction method"
  ) +
  ylim(0,100)+
  labs(x="Method",
       y="Coverage (%)")
plot(cov_avg)

extrlegend <- get_legend(
  sq_median+
    theme(legend.box.margin = margin(0,0,300,50))
)

extr_summ <- plot_grid(sq_avg+theme(legend.position='none'),
                     sq_median+theme(legend.position='none'),
                     cov_avg+theme(legend.position='none'),
                     cov_median+theme(legend.position='none'),
                     #extrlegend,
                     nrow=2,
                     labels=c("Mean (+/- SE)","Median (+/- 1.5 IQR)",
                              "Mean (+/- SE)","Median (+/- 1.5 IQR)"),
                     rel_widths=c(1,1,.25),
                     label_x=.05)
extr_summ
save_plot("extr_summary.jpg",extr_summ,base_height=8,base_width=14)

legends_sq <- plot_grid(reslegend,extrlegend,
                        sq_avg+theme(legend.position='none'),
                        ncol=1,
                        rel_heights = c(1,1,3),
                        labels=c("","","B"))
legends_sq

heatmap_res <- plot_grid(heatmap_res_sputum+theme(legend.position='none'),
                         legends_sq, 
                         #labels=c("DNA","Sputum extracts"),
                         #label_y=1.05,
                         labels=c("A",""),
                         ncol=2)
heatmap_res
save_plot("Supp_heatmap_all.svg",heatmap_res,
          base_height=8,base_width=12)


#Now I'm going to reformat the metadata a bit so it works nicely in nextstrain
TB004_coverage <- read.delim('data/TB004_coverage.tsv')
TB005_coverage <-read.delim('data/TB005_coverage.tsv')
TB003_coverage <- read.delim('../3_read_coverage/data/TB.coverage.tsv')
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
         NGS_Run_ID,numreads,covbases,coverage,meandepth,susceptibility,lineage,Template_dilution)%>%
  filter(is.na(Template_dilution) | Template_dilution==.01) %>% #remove sample duplicates (don't include a bunch of dilutions)
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
context_metadata <- read.delim('../../pipelines/nextstrain/data/input_metadata.tsv') %>%
  mutate(across(all_of(newdrugs),as.factor))
merged_metadata_ns <- context_metadata %>%
  full_join(metadata_ns,by=c(strain='Seq_ID',date='collection_date','EMB','ETH','INH','KAN','OFL','PZA','RIF','STR',
                             'continent','region','country','study'))

#check these
glab_merged <- merged_metadata_ns %>%
  filter(country=="Peru")

write.csv(merged_metadata_ns,file='../../pipelines/nextstrain/data/metadata_merged.csv')


