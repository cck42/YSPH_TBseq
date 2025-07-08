rm(list = ls())

library(tidyverse)
library(forcats)
library(stringr)
library(jsonlite)

# read meta / dst ---------------------------------------------------------

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



# 
# drug_lookup <- read.csv("data/DataDictionaryMoldova_share_11Jan2021.csv",header=T) %>% 
#   filter(grepl("dst",new_var)) %>%
#   mutate(drug = gsub(" from.*","",gsub("dst from.*: ","",description))) %>%
#   mutate(drug_short = gsub(".*dst_","",new_var)) %>%
#   select(drug_short,drug) %>% 
#   mutate(drug_type = case_when(drug %in% l1_2HRZE ~ "1st",
#                                drug %in% c(l2_fluoro,l2_inject) ~ "2nd",
#                                drug %in% others ~ "other"))
# 
# 
# drug_lookup <- unique(drug_lookup)
# 
# moldova_dst <- read.table("data/moldova_dst_results.txt",header=T,sep="\t",na.strings = c("NA",""))
# 
# moldova_xpert <- moldova_dst %>% 
#   rename(ID=tubeLabel) %>% 
#   select(ID,uniqueID,age,MDR_phenotype,MDR_WGS,xpertRifR) %>%
#   mutate(dst_method="xpert",
#          drug_short="r") %>%
#   rename(susceptibility = xpertRifR) %>%
#   filter(!is.na(susceptibility))
# 
# moldova_dst <- moldova_dst %>% 
#   pivot_longer(contains("dst_"),names_pattern = "(.*)dst_(.*)",names_to=c("dst_method","drug_short"),values_to="susceptibility") %>%
#   rename(ID=tubeLabel) %>% 
#   select(ID,uniqueID,age,MDR_phenotype,MDR_WGS,dst_method,drug_short,susceptibility)
# 
# moldova_dst <- rbind(moldova_dst,moldova_xpert)
# rm(moldova_xpert)
# moldova_dst <- merge(moldova_dst,unique(drug_lookup),by="drug_short",all.x=T) %>%
#   filter(drug %in% drugorder) %>%               
#   mutate(drug = factor(drug,levels=rev(drugorder)))              
# moldova_dst$dst_method <- fct_recode(moldova_dst$dst_method,"liquid media"="mgit","solid media"="lj")
# 


# read depths -------------------------------------------------------------

name_lookup <- sample_metadata %>% 
  select(Seq_ID,Original_ID) %>% 
  rename(sample=Seq_ID,ID=Original_ID)

genedepths <- read.table("data/all_genedepth.txt",header=T) %>%
                    filter(name!="name") %>% 
                    merge(name_lookup,by='sample')


#d/l panel 20230928 from	https://figshare.com/ndownloader/files/42494211	20230928

gene2res = list()
jsonfiles = list.files("./data/mykrobe",pattern = "*.json",full.names = T)
for(jf in jsonfiles) {
  print(jf)
  resjson <- fromJSON(paste(readLines(jf), collapse=""),flatten=T)
  for(r in names(resjson) ){
    gene = strsplit(r,"_")[[1]][1]
    res <- resjson[r]
    #print(paste(res,gene))
    if(!gene %in% gene2res) {
      gene2res[[gene]] <- c()
    }
    gene2res[[gene]] <- unique(c(gene2res[[gene]],res))
  }
}

resgenes <- names(gene2res)
generestab <- data.frame('gene'=character(),
                         'resistance'=character())
for(gene in resgenes) {
  i <- nrow(generestab)+1
  generestab[i,"gene"] = gene
  generestab[i,"resistance"] = paste(gene2res[[gene]][[1]],collapse=",")
}
generestab

# get_relevant_samples ----------------------------------------------------


sputum_samples <- sputum_extracts_summary %>%
  filter(liquefaction=="NALC 0.5% - NaOH"&fastprep=="40s") %>%
  pull(`Yale-ID`)

md_samples <- mykrobe_combined %>%
  filter(!is.na(Template_dilution),
         drug=="isoniazid",
         NGS_Run_ID %in% c("TB001","TB003"),
         NGS_Prep_Method == 'Mpox/Hybrid CovidSeq',
         Sample_source == "Moldova") %>% 
  rename(dst_method=NGS_Prep_Method) %>% 
  filter(Template_dilution=="0.01") %>%
  pull(sample) %>% 
  unique()


# get depth table ---------------------------------------------------------


mdgenedepths <- genedepths %>% filter(sample %in% md_samples,
                      name %in% resgenes) %>%
  rename(gene=name) %>%
  mutate(across(all_of(c('depth','start','end')),as.numeric)) %>% 
  mutate(length=end-start) %>% 
  select(ID,gene,length,depth) %>%
  pivot_wider(names_from=ID,values_from=depth) %>% 
  merge(generestab,by="gene") %>%
  mutate(resorder = gsub("Quinolones","Moxifloxacin,Levofloxacin",resistance)) %>% 
  mutate(resorder = tolower(gsub("\\,.*","",resorder))) %>% 
  mutate(resorder = factor(resorder,levels=drugorder)) %>%
  arrange(resorder) %>% 
  select(-resorder,-length) %>% 
  relocate(resistance,.before=gene)

write.table(mdgenedepths,file="moldova_resistance_gene_depths.txt",sep="\t",col.names=T,quote=F,row.names = F)


pegenedepths <- genedepths %>% filter(sample %in% sputum_samples,
                                      name %in% resgenes) %>%
  rename(gene=name) %>%
  mutate(across(all_of(c('depth','start','end')),as.numeric)) %>% 
  mutate(length=end-start) %>% 
  select(ID,gene,length,depth) %>%
  pivot_wider(names_from=ID,values_from=depth) %>% 
  merge(generestab,by="gene") %>%
  mutate(resorder = gsub("Quinolones","Moxifloxacin,Levofloxacin",resistance)) %>% 
  mutate(resorder = tolower(gsub("\\,.*","",resorder))) %>% 
  mutate(resorder = factor(resorder,levels=drugorder)) %>%
  arrange(resorder) %>% 
  select(-resorder,-length) %>% 
  relocate(resistance,.before=gene)

write.table(pegenedepths,file="peru_resistance_gene_depths.txt",sep="\t",col.names=T,quote=F,row.names = F)


