library(ggsci)
library(tidyverse)
library(magrittr)
library(egg)
#plotting the results from onek1k demultiplexing
setwd("~/Downloads/genotype_free/")
#freemuxlet 
ll <- system("ls pool*/genetic/gene_demulti/freemuxlet/freemuxlet_*/freemuxlet_out.clust1.samples.gz", intern = T)
p1 <- list()
for (fl in ll){
  pool<-sapply(strsplit(fl, "/", fixed=T),"[[",1)
  px<- read.table(fl, header = T)[,c("BARCODE" ,"DROPLET.TYPE","SNG.BEST.GUESS")]
  colnames(px) <- c("Barcode" ,"Freemuxlet","donor")
  p1[[pool]] <- px
}
p1<-bind_rows(p1,.id="pool")

#souporcell 
ll <- system("ls pool*/genetic/gene_demulti/souporcell/souporcell_pool*/souporcell_out/clusters.tsv", intern = T)
p2 <- list()
for (fl in ll){
  pool<-sapply(strsplit(fl, "/", fixed=T),"[[",1)
  px<- read.table(fl, header = T)[,c("barcode"  ,"status","assignment")]
  colnames(px) <- c("Barcode" ,"Souporcell","donor")
  p2[[pool]] <- px
}

p2<-bind_rows(p2,.id="pool")

#scsplit
ll <- system("ls pool*/genetic/gene_demulti/scSplit/scsplit_pool*/scsplit_out/scSplit_result.csv", intern = T)
p3 <- list()
for (fl in ll){
  pool<-sapply(strsplit(fl, "/", fixed=T),"[[",1)
  px<- read.table(fl, header = T)
  colnames(px) <- c("Barcode" ,"scSplit")
  px$donor <- px$scSplit
  p3[[pool]] <- px
}


p3<-bind_rows(p3,.id="pool")

#vireo

ll <- system("ls pool*/genetic/gene_demulti/vireo/vireo_pool*/vireo_out/donor_ids.tsv", intern = T)
p4 <- list()
for (fl in ll){
  pool<-sapply(strsplit(fl, "/", fixed=T),"[[",1)
  px<- read.table(fl, header = T)[,c("cell","donor_id")]
  colnames(px) <- c("Barcode" ,"vireo")
  px$donor <- px$vireo
  p4[[pool]] <- px
}


p4<-bind_rows(p4,.id="pool")


setwd("~/Downloads/genotype_aware/")
#demuxlet 
ll <- system("ls pool*/genetic/gene_demulti/demuxlet/demuxlet_*/demuxlet_res.best", intern = T)
p5 <- list()
for (fl in ll){
  pool<-sapply(strsplit(fl, "/", fixed=T),"[[",1)
  px<- read.table(fl, header = T)[,c("BARCODE" ,"DROPLET.TYPE","SNG.BEST.GUESS")]
  colnames(px) <- c("Barcode" ,"Demuxlet(GT)","donor")
  p5[[pool]] <- px
}
p5<-bind_rows(p5,.id="pool")

#vireo

ll <- system("ls pool*/genetic/gene_demulti/vireo/vireo_pool*/vireo_out/donor_ids.tsv", intern = T)
p6 <- list()
for (fl in ll){
  pool<-sapply(strsplit(fl, "/", fixed=T),"[[",1)
  px<- read.table(fl, header = T)[,c("cell","donor_id")]
  colnames(px) <- c("Barcode" ,"Vireo(GT)")
  px$donor <- px$vireo
  p6[[pool]] <- px
}

p6<-bind_rows(p6,.id="pool")


#count the cells in the pools

p1 %>% group_by(pool) %>% summarize(frenocell=n())
p2 %>% group_by(pool) %>% summarize(sounocell=n())
p3 %>% group_by(pool) %>% summarize(scsnocell=n())
p4 %>% group_by(pool) %>% summarize(virnocell=n())
#scsplit identifies less cells !

#order pools by increasing number of cells
p1 %>% group_by(pool) %>% summarize(nocell=n()) %>% arrange(nocell) %>% pull(pool) ->pool_order

p6 %>% 
  mutate(classification = ifelse(grepl("doublet",`Vireo(GT)`),"doublet","singlet")) %>% 
  group_by(pool,classification) %>% summarise(cellby=n()) %>% 
  group_by(pool) %>% mutate(nocell=sum(cellby), method = "Vireo(GT)") ->p6_summary

p5 %>% 
  mutate(classification = ifelse(grepl("SNG",`Demuxlet(GT)`),"singlet","doublet")) %>% 
  group_by(pool,classification) %>% summarise(cellby=n()) %>% 
  group_by(pool) %>% mutate(nocell=sum(cellby), method = "Demuxlet(GT)") ->p5_summary


p4 %>% 
  mutate(classification = ifelse(grepl("donor",vireo),"singlet","doublet")) %>% 
  group_by(pool,classification) %>% summarise(cellby=n()) %>% 
  group_by(pool) %>% mutate(nocell=sum(cellby), method = "Vireo") ->p4_summary


p3 %>%  
  mutate(classification = ifelse(grepl("SNG",scSplit),"singlet","doublet")) %>% 
  group_by(pool,classification) %>% summarise(cellby=n()) %>% 
  group_by(pool) %>% mutate(nocell=sum(cellby), method = "scSplit") ->p3_summary


p2 %>%  
  mutate(classification = ifelse(grepl("singlet",Souporcell),"singlet","doublet")) %>% 
  group_by(pool,classification) %>% summarise(cellby=n()) %>% 
  group_by(pool) %>% mutate(nocell=sum(cellby), method = "Souporcell") ->p2_summary

p1 %>% 
  mutate(classification = ifelse(grepl("SNG",Freemuxlet),"singlet","doublet")) %>% 
  group_by(pool,classification) %>% summarise(cellby=n()) %>% 
  group_by(pool) %>% mutate(nocell=sum(cellby), method = "Freemuxlet") ->p1_summary

sel <- c("pool",
  "cellby",
  "classification",
  "nocell",
  "method")

summ<-rbind(p1_summary[,sel] , p2_summary[,sel] ,p3_summary[,sel], p4_summary[,sel],p5_summary[,sel],p6_summary[,sel] )
summ %<>% mutate(percent = 100*cellby/nocell)

summ %<>% mutate(pool =factor(pool, levels = pool_order))
mo <- c("Freemuxlet","Vireo","Souporcell","scSplit","Demuxlet(GT)","Vireo(GT)")
summ %<>% mutate(method =factor(method, levels = mo),
                 GT = ifelse(grepl("GT", method),"genotype-aware","genotype-free"))


summ %>% dplyr::filter(classification=="singlet") %>% 
  ggplot(aes(x= pool,y=percent, color=method, group=method)) + 
  geom_point() + 
  geom_smooth(method="lm", se=FALSE, alpha=0.2) + theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13)) +
  scale_color_jco() + facet_grid(~GT) +
  xlab("Pools ordered by increasing no. of cells ") + ylab("Percent Singlets") ->pb

# check model freemuxlet
summ %>% dplyr::filter(classification=="singlet", method=="Freemuxlet") ->df
summary(lm(percent~nocell, data = df))$r.squared
summary(lm(percent~nocell, data = df))$coefficient[,4]

models <- summ %>% dplyr::filter(classification=="singlet") %>% 
  split(.$method) %>% 
  map(~lm(percent ~ nocell, data = .))

map_dfr(models, ~summary(.)$r.squared)
map_dfr(models, ~summary(.)$coefficient[,4]) %>% mutate(padj = p.adjust(nocell))
# Read onek1k single cell metadata for comparison
onek<- read.delim("~/Downloads/onek1k-cell_metadata.tsv", header=T, sep="\t")


# this number doesn't reflect the number of donors in the pools as deposited on https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196830
onek %>% group_by(pool_number,donor_id) %>% distinct(pool_number,donor_id) %>% group_by(pool_number) %>% summarize(ndon=n()) %>% arrange(ndon)
onek %<>% select(barcode, pool_number,donor_id)

onek %<>% mutate(pool = paste0("pool",pool_number),
                Barcode = paste0(sapply(strsplit(barcode,"-"),"[[",1),"-1"))
onek$ground_t <- onek$donor_id


dist_don<-function(method_data){
  method_data %>% distinct(pool,donor) %>% 
    group_by(pool) %>% summarize(ndon=n()) %>% arrange(ndon) -> kk
  return(kk)
}

onek %>% distinct(pool,donor_id) %>% 
  group_by(pool) %>% summarize(ndon=n()) %>% arrange(ndon) ->gt
p6 %>% dplyr::filter(grepl("donor",`Vireo(GT)`)) %>% distinct(pool, `Vireo(GT)`)

#parsed https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196830
dep <- read.table("~/Documents/projects/donor_deconvolution/poolgenotype_truth_onek1k", sep="\t", header=T)
head(dep)
dep$pool<-gsub(" ", "",dep$pool)
list(dep,gt, 
     dist_don(p1),
     dist_don(p2 %>% dplyr::filter(Souporcell=="singlet")),
     dist_don(p3 %>% dplyr::filter(grepl("SNG",scSplit))),
     dist_don(p4 %>% dplyr::filter(!grepl("doublet",vireo))),
     dist_don(p5 %>% dplyr::filter(grepl("SNG",`Demuxlet(GT)`))), 
     dist_don(p6 %>% mutate(donor = `Vireo(GT)`) %>% dplyr::filter(!grepl("doublet",donor)) )) %>% 
  reduce(left_join, by = "pool") ->pinfo

colnames(pinfo) <- c("pool","genotypes_true","anndata_baseline","freemuxlet","souporcell","scsplit","vireo","demuxlet","vireo_GT")


ll <- system("ls ~/Downloads/donor_gt_pools/*/*.csv", intern=T)

plist<-list()
for (fil in ll){
  fn <- gsub(".csv","",sapply(strsplit(basename(fil),"_"),"[[",4))
  plist[[fn]] <- read.table(fil,sep=',', header = T)[,c(1,2), drop=FALSE]
}
identity<-bind_rows(plist, .id="pool")

#there's 1015 unique donors,
length(unique(identity$donor_gt))
# but some donors appear in multiple pools
length(identity$donor_gt)
# the anndata obj doesn't have all the donors, so we use the baseline only
baseline_donors<-identity$donor_gt[identity$donor_gt %in% unique(onek$ground_t)]

#both genotype aware calls share the same donor ids.
any(!unique(p5$donor)%in%unique(p6$`Vireo(GT)`))

#generate a consensus call between demuxlet and vireo to deanonimize the genotype_free calls
merge(p5,p6,by=c("pool","Barcode")) %>% 
  mutate(callout= ifelse(`Vireo(GT)`==donor, "match","mismatch" )) ->adj_gt

adj_gt %>% group_by(callout) %>% 
  summarize(call = n()) %>% 
  mutate(tot=sum(call), percent = 100* call/tot)

adj_gt %>% dplyr::filter(callout =="match") ->common_gt

dim(common_gt)

common_gt %>% dplyr::filter(donor %in% baseline_donors) ->deanonimizer


deanonimizer %<>% select(pool,Barcode,donor)

colnames(deanonimizer)[3] <-"donor_id"

p1%>% dplyr::filter(pool=="pool1") %>% 
  left_join(deanonimizer, by=c("pool","Barcode")) %>% 
  group_by(donor, donor_id) %>% 
  summarize(pmax=n()) %>% 
  dplyr::filter(pmax == max(pmax)) %>% 
  arrange(donor,donor_id) %>% select(donor, donor_id)




dict_func <- function(method_data, pool_use){
   method_data %>% dplyr::filter( pool==!!pool_use) %>% 
    left_join(deanonimizer, by=c("pool","Barcode")) %>% 
    group_by(donor, donor_id) %>% 
    summarize(pmax=n()) %>% 
    dplyr::filter(pmax == max(pmax)) %>%
    arrange(donor,donor_id) %>% select(donor, donor_id) -> key
  key$pool <- rep(pool_use)
  return(key)
}  

pools<-unique(onek$pool)  

callout_func<-function(method_data){
  d<-list()
    for (pp in pools){
      d[[pp]]<-dict_func(method_data, pool_use=pp)
    }
  
  bind_rows(d) %>% as.data.frame() ->legend
  legend$donor <- as.character(legend$donor)
  method_data$donor <- as.character(method_data$donor)
  
  dm<-full_join(method_data, legend, by=c("pool","donor"))
  
  
  
  merge(dm, onek %>% select(-donor_id), by=c("pool","Barcode")) %>% 
    group_by(donor_id, ground_t,pool) %>% summarize(nocells=n()) %>%  
    rowwise() %>% 
    mutate(callout= ifelse(donor_id==ground_t, "match","mismatch" )) %>% 
    group_by(callout,pool) %>% 
    summarize(call = sum(nocells)) %>% 
    group_by(pool) %>% 
    mutate(tot=sum(call), percent = 100* call/tot) -> pcout
    return(pcout)

}


#freemuxlet
df1<-callout_func(p1)
df1$method = "Freemuxlet"

#souporcell
df2<-callout_func(p2)
df2$method = "Souporcell"

#scsplit  
df3<-callout_func(p3)
df3$method = "scSplit"

#vireo
df4<-callout_func(p4)
df4$method = "Vireo"

rbind(df1,df2,df3, df4) ->matched_sng

#for donor aware you just need to merge with onek1k 

merge(p6, onek %>% select(-donor_id), by=c("pool","Barcode")) %>% #first merge by pool and barcode so the cells are intersected based on their name
  mutate(donor = `Vireo(GT)`) %>%
  dplyr::filter(donor %in% baseline_donors) %>%   #now only use cells in the intersection to evaluate
  group_by(donor, ground_t,pool) %>% summarize(nocells=n()) %>%  
  rowwise() %>% 
  mutate(callout= ifelse(donor==ground_t, "match","mismatch" )) %>% 
  group_by(callout,pool) %>% 
  summarize(call = sum(nocells)) %>% 
  group_by(pool) %>% 
  mutate(tot=sum(call), percent = 100* call/tot) -> df6
df6$method = "Vireo(GT)"

unique(p5$`Demuxlet(GT)`)

merge(p5, onek %>% select(-donor_id), by=c("pool","Barcode")) %>% 
  dplyr::filter(donor %in% baseline_donors) %>% 
  group_by(donor, ground_t,pool) %>% summarize(nocells=n()) %>%  
  rowwise() %>% 
  mutate(callout= ifelse(donor==ground_t, "match","mismatch" )) %>% 
  group_by(callout,pool) %>% 
  summarize(call = sum(nocells)) %>% 
  group_by(pool) %>% 
  mutate(tot=sum(call), percent = 100* call/tot) -> df5
df5$method = "Demuxlet(GT)"

  


########
rbind(matched_sng,df5,df6) ->matched_sng

matched_sng %<>% mutate(pool =factor(pool, levels = pool_order))

matched_sng %<>% mutate(method =factor(method, levels = mo),
                 GT = ifelse(grepl("GT", method),"genotype-aware","genotype-free"))



matched_sng %>% dplyr::filter(callout=="match") %>% 
  ggplot(aes(pool, percent, color=method, group=method)) +
  geom_point()+ 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13)) +
  scale_color_jco() + facet_grid(~GT) +
  xlab("") + ylab("Percent Matched Singlets") ->pa



g1<-egg::ggarrange(plots = list(pa,pb), ncol=1, common.legend=TRUE)

ggsave(g1, filename = "~/Documents/papers/hadge/onek1k_genotype_free_consensus_call.pdf", height = 6, width = 14)



matched_sng %>% dplyr::filter(callout=="match") %>% 
  ggplot(aes(method, percent, color=method, group=method)) +
  geom_boxplot()+ 
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle=45,hjust=1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13)) +
  scale_color_jco() +
  xlab("") + ylab("Percent Matched Singlets") ->bp1 


summ %>% dplyr::filter(classification=="singlet") %>% 
  ggplot(aes(method, percent, color=method, group=method)) +
  geom_boxplot()+ 
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle=45,hjust=1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13)) +
  scale_color_jco() +
  xlab("") + ylab("Percent Singlets") ->bp2 

ggsave(bp1, filename ="~/Documents/papers/hadge/summary_matchedSNG.pdf", height = 4, width = 5)
ggsave(bp2, filename ="~/Documents/papers/hadge/summary_SNG.pdf", height = 4, width = 5)

models <- summ %>% dplyr::filter(classification=="singlet") %>% 
  split(.$method) %>% 
  map(~lm(percent ~ nocell, data = .))

map_dfr(models, ~summary(.)$r.squared)
map_dfr(models, ~summary(.)$coefficient[,4]) %>% mutate(padj = p.adjust(nocell))

dim(summ)

p1 %>% group_by(pool) %>% summarize(nocell=n()) %>% arrange(nocell) ->data_pools

models_SNG <- matched_sng %>% dplyr::filter(callout=="match") %>% 
  left_join(data_pools) %>% 
  split(.$method) %>% 
  map(~lm(percent ~ nocell, data = .))

map_dfr(models_SNG, ~summary(.)$r.squared)
map_dfr(models_SNG, ~summary(.)$coefficient[,4]) %>% mutate(padj = p.adjust(nocell))



identity %>% distinct(pool, donor_gt) %>% group_by(pool) %>% summarise(ndonors=n()) ->ndonors

ndonors %>% arrange(ndonors)
table(ndonors$ndonors)



matched_sng %>% dplyr::filter(callout=="match") %>% left_join(ndonors) %>% 
  ggplot(aes(reorder(pool,ndonors), percent, color=method, group=method)) +
  geom_point()+ 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13)) +
  scale_color_jco() + facet_grid(~GT) +
  xlab("") + ylab("Percent Matched Singlets") ->g3

summ %>% dplyr::filter(classification=="singlet") %>% left_join(ndonors) %>% 
  ggplot(aes(x= reorder(pool,ndonors),y=percent, color=method, group=method)) + 
  geom_point() + 
  geom_smooth(method="lm", se=FALSE, alpha=0.2) + theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13)) +
  scale_color_jco() + facet_grid(~GT) +
  xlab("Pools ordered by number of donors in pool") + ylab("Percent Singlets") ->g4

g2<-egg::ggarrange(plots = list(g3,g4), ncol=1, common.legend=TRUE)

ggsave(g2, filename = "~/Documents/papers/hadge/onek1k_Ndonors_common_calls.pdf", height = 6, width = 14)


matched_sng %>% dplyr::filter(callout=="match") %>% left_join(ndonors) %>% 
  mutate(ndonors = factor(ndonors, levels = c(9L, 11L, 12L, 13L, 14L, 15L))) %>% 
  ggplot(aes(method, percent, color=ndonors)) +
  geom_boxplot()+ 
  theme_bw() +
  theme(#axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13)) +
  scale_color_jco() + 
  xlab("") + ylab("Percent Matched Singlets") -> g7
ggsave(g7, filename = "~/Documents/papers/hadge/boxplot_number_of_donors_in_pool.pdf", height = 4, width = 7)

ndonors

matched_sng %>% dplyr::filter(callout=="match") %>% left_join(ndonors) %>% 
  split(.$method) %>% 
  map(~lm(percent ~ ndonors, data = .)) ->mod

map_dfr(mod, ~summary(.)$r.squared)
map_dfr(mod, ~summary(.)$coefficient[,4]) %>% mutate(padj = p.adjust(ndonors))
