

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(argparse)
library(ggplot2)
library(cowplot)
library(ggthemes)

coverageinfo <- read_delim('/Users/chaneykalinich/Documents/PGCoE/github/combinedcoverage.tsv',delim='\t')
datasheet <- read_delim('/Users/chaneykalinich/Documents/PGCoE/github/pgcoe_pipeline/pgcoe_pipeline/data/PGCOE_Datasheet - M.tuberculosis_Seq.csv',delim=',')

data_combined <- coverageinfo %>%
  left_join(datasheet, by=c('sample'='Seq_ID'))

#statsheaders <- c('sample','subsample','rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')

coverageplot <- data_combined %>%
  mutate(Template_dilution=replace_na(Template_dilution,1)) %>%
  mutate(NGS_primers=replace_na(NGS_primers,'None'))%>%
  filter(Original_ID!='NTC') %>%
  ggplot() +
  geom_point(aes(x=Template_dilution,y=coverage,color=Original_ID, shape=NGS_primers)) +
  scale_x_log10() +
  #facet_wrap('subsample') +
  #theme_minimal_hgrid()+
  #theme_stata() +
  theme(axis.text.x = element_text(angle = 300, vjust = 0.5, hjust=0)) +
  labs(x="Fold dilution from template DNA",
       y="Coverage of H37Rv Reference (%)") +
  scale_color_discrete(name="Sample")
plot(coverageplot)

coverageplot_ct <-data_combined %>%
  #mutate(Template_dilution=replace_na(Template_dilution,1)) %>%
  filter(Original_ID!='NTC') %>%
  filter(subsample==1.0) %>%
  #filter(subsample %in% c(1,0.5,0.25,0.05)) %>%
  mutate(NGS_primers=replace_na(NGS_primers,'None'))%>%
  #mutate(ct2=replace_na(CT,15))
  #mutate(logdil=log10(Template_dilution)+15) %>%
  #mutate(ct2=log10(Template_dilution)+15) %>%
  ggplot() +
  geom_hline(aes(yintercept=80)) +
  geom_smooth(aes(x=CT,y=coverage),method='loess') +
  geom_point(aes(x=CT,y=coverage,color=Original_ID,shape=NGS_primers)) +
  #scale_x_log10() +
  #facet_grid(cols=vars(subsample), rows=vars(NGS_primers))  +
  facet_grid('NGS_primers') +
  labs(x='Ct',y="Coverage of H37Rv Reference (%)") 
  #geom_hline(aes(yintercept=80,linetype="dashed")) +
  #geom_smooth(aes(x=CT,y=coverage),method='loess')
plot(coverageplot_ct)


foldchange <- data_combined %>%
  mutate(NGS_primers=replace_na(NGS_primers,'None'))%>%
  pivot_wider(names_from="NGS_primers",values_from='coverage',id_cols=c(Original_ID,CT,subsample,Template_dilution)) %>%
  filter(Original_ID %in% c('111-10712-18','111-5264-18')) %>%
  ggplot() +
  geom_linerange(aes(ymin=None,ymax=TBv1.1,x=CT)) +
  geom_point(aes(x=CT,y=None,color=Original_ID)) +
  geom_point(aes(x=CT,y=TBv1.1,color=Original_ID)) +
  facet_wrap('subsample') +
  labs(x='Ct',y="Coverage of H37Rv Reference (%)") 
plot(foldchange)


foldchange_reads<- data_combined %>%
  mutate(NGS_primers=replace_na(NGS_primers,'None'))%>%
  pivot_wider(names_from="NGS_primers",values_from='numreads',id_cols=c(Original_ID,CT,subsample,Template_dilution)) %>%
  filter(Original_ID %in% c('111-10712-18','111-5264-18')) %>%
  mutate(changereads=TBv1.1/None)
  ggplot() +
  geom_linerange(aes(ymin=None,ymax=TBv1.1,x=CT)) +
  geom_point(aes(x=CT,y=None,color=Original_ID)) +
  geom_point(aes(x=CT,y=TBv1.1,color=Original_ID)) +
  facet_wrap('subsample') +
  labs(x='Ct',y="Number of reads aligned to H37Rv reference") +
  scale_y_log10()
plot(foldchange_reads)
