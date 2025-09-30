rm(list = ls())

library(tidyverse)
library(forcats)
library(stringr)
library(jsonlite)
library(readxl)

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




# read depths (genes) -------------------------------------------------------

name_lookup <- sample_metadata %>% 
  select(Seq_ID,Original_ID) %>% 
  dplyr::rename(sample=Seq_ID,ID=Original_ID)

genedepths <- read.table("data/all_genedepth.txt",header=T) %>%
                    filter(name!="name") %>% 
                    merge(name_lookup,by='sample')

#get gene-to-locus mapping
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


# read depths (resist loci) -----------------------------------------------

resdepths <- read.table("data/all_resistdepth.txt",header=T) %>%
  filter(name!="name") %>% 
  merge(name_lookup,by='sample')




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
  dplyr::rename(dst_method=NGS_Prep_Method) %>% 
  filter(Template_dilution=="0.01") %>%
  pull(sample) %>% 
  unique()


# get depth table (genes) ---------------------------------------------------------


mdgenedepths <- genedepths %>% filter(sample %in% md_samples,
                      name %in% resgenes) %>%
  dplyr::rename(gene=name) %>%
  dplyr::mutate(across(all_of(c('depth','start','end')),as.numeric)) %>% 
  dplyr::mutate(length=end-start) %>% 
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
  dplyr::rename(gene=name) %>%
  dplyr::mutate(across(all_of(c('depth','start','end')),as.numeric)) %>% 
  dplyr::mutate(length=end-start) %>% 
  select(ID,gene,length,depth) %>%
  pivot_wider(names_from=ID,values_from=depth) %>% 
  merge(generestab,by="gene") %>%
  dplyr::mutate(resorder = gsub("Quinolones","Moxifloxacin,Levofloxacin",resistance)) %>% 
  dplyr::mutate(resorder = tolower(gsub("\\,.*","",resorder))) %>% 
  dplyr::mutate(resorder = factor(resorder,levels=drugorder)) %>%
  arrange(resorder) %>% 
  select(-resorder,-length) %>% 
  relocate(resistance,.before=gene)

write.table(pegenedepths,file="peru_resistance_gene_depths.txt",sep="\t",col.names=T,quote=F,row.names = F)


# get depth table (res loci) ---------------------------------------------------------


mdresdepths <- resdepths %>% filter(sample %in% md_samples) %>%
  dplyr::mutate(across(all_of(c('depth','start','end')),as.numeric)) %>% 
  dplyr::mutate(length=end-start) %>% 
  select(ID,name,length,depth) %>%
  pivot_wider(names_from=ID,values_from=depth)


#View(mdresdepths)
write.table(mdresdepths,file="moldova_resistance_locus_depths.txt",sep="\t",col.names=T,quote=F,row.names = F)


peresdepths <- resdepths %>% filter(sample %in% sputum_samples) %>%
  dplyr::mutate(across(all_of(c('depth','start','end')),as.numeric)) %>% 
  dplyr::mutate(length=end-start) %>% 
  select(ID,name,length,depth) %>%
  pivot_wider(names_from=ID,values_from=depth)

#View(peresdepths)

write.table(peresdepths,file="peru_resistance_res_depths.txt",sep="\t",col.names=T,quote=F,row.names = F)



