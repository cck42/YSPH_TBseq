rm(list = ls())

library(tidyverse)


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
  filter(Sample_source %in% c("Moldova")) %>%
  mutate(drug=tolower(drug)) %>%
  filter(Template_dilution=="0.01") %>% 
  mutate(NGS_Prep_Method = fct_recode(NGS_Prep_Method,
                                                   "seq"="Hybrid CovidSeq without amplification",
                                                   "TB-seq"="Mpox/Hybrid CovidSeq"))


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
moldova_dst$dst_method <- fct_recode(moldova_dst$dst_method,"liquid"="mgit","solid"="lj")



dstcf <- moldova_dst %>% select(ID,drug,dst_method,drug_type,"phenotype"=susceptibility) %>% 
  pivot_wider(names_from = dst_method,values_from=c(phenotype)) %>% 
  merge(mykrobe_combined %>% 
          select(Original_ID,drug,NGS_Prep_Method,"predicted"=susceptibility) %>%
          pivot_wider(names_from = NGS_Prep_Method,values_from=c(predicted)),
        by.x=c("ID","drug"),by.y=c("Original_ID","drug")) %>%
  unique() %>% 
  arrange(ID,drug_type,drug)

write.table(dstcf,"moldova_dst_mykrobe_table.tsv",sep="\t",row.names=F,quote=F)
