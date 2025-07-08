rm(list = ls())

library(tidyverse)
library(readr)
library(cowplot)
library(rlang)
library(forcats)
library(readxl)
library(stringr)


all_mykrobe <- read.csv('data/combined_mykrobe.csv')
sample_metadata <-read.csv('data/M.tuberculosis.csv',
                           na.strings=c("","NaN","NA"))
sputum_extracts_summary <- read_excel('data/sputum_extraction_summary.xlsx',sheet=1) %>%
  mutate("Heat_inactivation"=as_factor(Heat_inactivation))%>%
  mutate("Chemical_inactivation"=as_factor(Chemical_inactivation))%>%
  mutate("liquefaction"=as_factor(liquefaction))%>%
  mutate("fastprep"=as_factor(fastprep))%>%
  mutate("pt_sample_num"=as_factor(pt_sample_num))

mykrobe_combined <- all_mykrobe %>%
  left_join(sample_metadata, by=c(sample='Seq_ID')) %>%
  filter(Sample_source %in% c("Cayetano","Moldova")) %>%
  mutate(drug=tolower(drug))


l1_2HRZE = c("isoniazid","rifampicin","pyrazinamide","ethambutol")
#l1_2HPZM = c("isoniazid","rifampicin","pyrazinamide","moxifloxacin")
l2_fluoro = c("moxifloxacin","levofloxacin")
l2_inject = c("amikacin","kanamycin","capreomycin")
others = c("ciprofloxacin","delamanid","ethionamide","linezolid","ofloxacin","streptomycin")
drugorder <- unique(c(l1_2HRZE,l2_fluoro,l2_inject,others))
mykrobe_combined$drug <- factor(mykrobe_combined$drug, levels=rev(drugorder),ordered=T)
mykrobe_combined <- mykrobe_combined %>% 
  mutate(drug_type = case_when(drug %in% l1_2HRZE ~ "1st",
                               drug %in% c(l2_fluoro,l2_inject) ~ "2nd",
                               drug %in% others ~ "other"))


# mykrobe_combined %>% 
#   filter(NGS_Run_ID %in% c("TB005","TB006")) %>% 
#   rename(ID=Original_ID) %>% 
#   select(ID,drug,susceptibility,drug_type,NGS_Prep_Method)
# 
# mykrobe_combined %>% 
#   filter(NGS_Run_ID %in% c("TB005","TB006")) %>% 
#   rename(ID=Original_ID) %>% 
#   select(ID,drug,susceptibility,NGS_Prep_Method)
# 
# 

drug_lookup <- read.csv("data/DataDictionaryMoldova_share_11Jan2021.csv",header=T) %>% 
  filter(grepl("dst",new_var)) %>%
  mutate(drug = gsub(" from.*","",gsub("dst from.*: ","",description))) %>%
  mutate(drug_short = gsub(".*dst_","",new_var)) %>%
  select(drug_short,drug) %>% 
  mutate(drug_type = case_when(drug %in% l1_2HRZE ~ "1st",
                               drug %in% c(l2_fluoro,l2_inject) ~ "2nd",
                               drug %in% others ~ "other"))


drug_lookup <- unique(drug_lookup)

moldova_dst <- read.table("data/moldova_dst_results.txt",header=T,sep="\t",na.strings = c("NA",""))

moldova_xpert <- moldova_dst %>% 
  rename(ID=tubeLabel) %>% 
  select(ID,uniqueID,age,MDR_phenotype,MDR_WGS,xpertRifR) %>%
  mutate(dst_method="xpert",
         drug_short="r") %>%
  rename(susceptibility = xpertRifR) %>%
  filter(!is.na(susceptibility))

moldova_dst <- moldova_dst %>% 
        pivot_longer(contains("dst_"),names_pattern = "(.*)dst_(.*)",names_to=c("dst_method","drug_short"),values_to="susceptibility") %>%
        rename(ID=tubeLabel) %>% 
        select(ID,uniqueID,age,MDR_phenotype,MDR_WGS,dst_method,drug_short,susceptibility)

moldova_dst <- rbind(moldova_dst,moldova_xpert)

moldova_dst <- merge(moldova_dst,unique(drug_lookup),by="drug_short",all.x=T) %>%
                        filter(drug %in% drugorder) %>%               
                        mutate(drug = factor(drug,levels=rev(drugorder)))              
moldova_dst$dst_method <- fct_recode(moldova_dst$dst_method,"liquid media"="mgit","solid media"="lj")





#######
#  Make all plots
######

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

best_ext_samples <- sputum_extracts_summary %>%
                          filter(liquefaction=="NALC 0.5% - NaOH"&fastprep=="40s") %>%
                          pull("ID")

mykrobe_best_peru <- mykrobe_combined %>% 
                        filter(Original_ID %in% best_ext_samples) %>% 
                        rename(ID=sample) %>% 
                        select(ID,Original_ID,drug,susceptibility,drug_type,NGS_Run_ID,NGS_Prep_Method)

heatmap_res_sputum_small <- mykrobe_best_peru %>%
  filter(!is.na(drug)) %>%
  filter(NGS_Run_ID == "TB005") %>%
  ggplot()+
  geom_tile(aes(x=Original_ID,y=drug,fill=susceptibility,height=0.95,width=0.95))+ 
  scale_fill_manual(
    values=c("r"="pink","R"="#B20026","S"="gray"),
    name="Susceptibility",
    na.value='white',
    labels=c("R"="Resistant","r"="Heteroresistant","S"="Susceptible")
  ) +
  labs(x='Sample',
       y='Drug') +
  theme_half_open()+
  theme(axis.title=element_blank(),
        axis.text.x=element_text(angle=45,vjust=1,hjust=1),
        panel.background = element_rect(fill='white')) 
plot(heatmap_res_sputum_small)

mykrobe_moldova <- mykrobe_combined %>% 
                        filter(!is.na(Template_dilution),!is.na(drug),NGS_Run_ID %in% c("TB001","TB003")) %>%
                        filter(!Template_dilution %in% c("2e-08")) %>%
                        rename(ID=Original_ID,
                                     genes=starts_with("genes"),
                               variants=starts_with("variants"),) %>% 
                        mutate(Template_dilution = factor(Template_dilution,levels=unique(as.numeric(Template_dilution)))) %>% 
                        select(ID,drug,susceptibility,drug_type,NGS_Prep_Method,Template_dilution,genes,variants)
mykrobe_moldova$NGS_Prep_Method <- fct_recode(mykrobe_moldova$NGS_Prep_Method,"TB-seq"="Mpox/Hybrid CovidSeq","direct-seq"="Hybrid CovidSeq without amplification")

full_heatmap_moldova <- mykrobe_moldova %>%
  #filter(!is.na(drug)) %>%
  ggplot()+
  geom_tile(aes(x=Template_dilution,y=drug,fill=susceptibility,height=0.95,width=0.95))+ 
  scale_fill_manual(
    values=c("pink","#B20026","gray"),
    name="Susceptibility",
    na.value='gray'
  ) +
  labs(x='Sample',
       y='Drug') +
  theme_half_open()+
  theme(axis.title=element_blank(),
        axis.text.x=element_text(angle=45,vjust=1,hjust=1),
        panel.background = element_rect(fill='white')) + 
  facet_grid(NGS_Prep_Method ~ ID)

plot(full_heatmap_moldova)

dst_cf <- rbind(moldova_dst %>% 
                  filter(!is.na(susceptibility)) %>%
                  select(ID,dst_method,drug,drug_type,susceptibility) ,
                mykrobe_moldova %>%
                  rename(dst_method=NGS_Prep_Method) %>% 
                  filter(Template_dilution=="0.01") %>%
                  select(ID,dst_method,drug,drug_type,susceptibility)
          ) %>% 
        mutate(dst_method = factor(dst_method,levels=c("xpert","solid media","liquid media","direct-seq","TB-seq"),ordered=T)) %>%
        mutate(ID = factor(ID,levels=c("111-10712-18","444-2261-18",
                                       "111-10300-18",
                                       "444-2281-18","444-3403-18", "444-2588-18", "111-2565-18","111-5264-18") ))

                           
dst_cf_plot <- ggplot(dst_cf,aes(x=dst_method,y=drug,fill=susceptibility))+
  geom_point(data=dst_cf,size=4,shape=21) + 
  geom_tile(data=subset(dst_cf,dst_method %in% c("TB-seq")),height=0.95,width=0.95)+ 
  scale_fill_manual(
    values=c("pink","#B20026","gray"),
    name="Susceptibility",
    na.value='white',
    labels=c("R"="Resistant","r"="Heteroresistant","S"="Susceptible")
  ) +
  theme_half_open()+
  theme(axis.title=element_blank(),
        axis.text.x=element_text(angle=45,vjust=1,hjust=1),
        panel.background = element_rect(fill='white')) + 
  facet_grid(. ~ ID)

dst_cf_plot
save_plot("dst_comparison_moldova.svg",dst_cf_plot,base_width=12.5,base_height=4)


dst_cf_plot_pw <- dst_cf_plot + 
                    facet_grid(drug_type ~ ID,scale="free_y",space="free_y",switch="y") +
                    scale_y_discrete(position = "right") + 
                    theme(axis.text.x=element_text(size=9),
                          strip.text.x=element_text(size=9),
                          legend.position = "bottom")
heatmap_res_sputum_pw <- heatmap_res_sputum_small + 
                            facet_grid(drug_type ~ .,scale="free_y",space="free_y") +
                            theme(legend.position="none",
                                  strip.background = element_blank(),
                                  strip.text = element_blank(),
                                axis.text.x=element_text(size=9),
                                axis.text.y=element_blank(),
                          )

dst_cf_plot_pw + heatmap_res_sputum_pw + guide_area() + 
  plot_layout(design="ab\ncc",widths=c(7,3),heights=c(9,0),guides = 'collect') + 
  plot_annotation(tag_levels = 'a',)

ggsave("Fig2_heatmap.png",width=330,height=120,units="mm")
ggsave("Fig2_heatmap.svg",width=330,height=120,units="mm")


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
save_plot("extr_summary.svg",extr_summ,base_height=8,base_width=14)

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



