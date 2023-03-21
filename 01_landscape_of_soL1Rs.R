######## Fig 1a (Experiment design : drawings)

######## Fig 1b (L1 count by tissue types)
df_count_total <- df_st1 %>% select(sampleID) %>% left_join(df_st2 %>% filter(MEI_TYPE %in% c("Solo-L1","Partnered_transduction","Orphan_transduction")) %>% count(sampleID)) %>% rename(L1_count = n) %>% 
  mutate(tissue = ifelse(str_detect(sampleID,"tumor"), "Colon carcinoma (19)",
                         ifelse(str_detect(sampleID,"FP"), "Colon adenoma (12)", 
                                ifelse(str_detect(sampleID,"BM[HP]"), "HSC (140)", 
                                       ifelse(str_detect(sampleID,"HC"), "Colon (406)", "Fibroblast (341)"))))) %>% 
  mutate(L1_count = ifelse(is.na(L1_count),0,L1_count)) %>% 
  mutate(group = ifelse(L1_count > 50, ">50", 
                        ifelse(L1_count > 10, "11-50",
                               ifelse(L1_count >= 6, "6-10", 
                                      ifelse(L1_count >= 3, "3-5", 
                                             ifelse(L1_count >=1, "1-2", "0"))))))

df_count_total %>% 
  select(sampleID,tissue,group) %>% 
  ggplot(aes(factor(tissue,levels=c("Colon carcinoma (19)","Colon adenoma (12)","Colon (406)","Fibroblast (341)","HSC (140)")),
             fill=factor(group,levels=c("0","1-2","3-5","6-10","11-50",">50"))))+
  geom_bar(position = "fill")+theme_classic()+coord_flip()+
  scale_fill_manual(values = c("#c5c5c5","#e5c53a","#e39139","#e93e27","#bb1d3e","#731124"))+
  xlab("")+ylab("Clones (%)")+
  theme(axis.line.y = element_blank(), axis.text = element_text(size=16), axis.title = element_text(size=16))+
  theme(legend.position = "top",legend.title=element_text(size=16),legend.text=element_text(size=16))+
  guides(fill=guide_legend(title="# of soL1Rs per clone",nrow=2))


######## Fig 1c (L1 count by individual)
df_count <- df_count_total %>% filter(tissue == "Colon (406)")

df_count_pt_info <- df_count %>% mutate(patientID = str_sub(sampleID,1,4)) %>% 
  group_by(patientID) %>% summarize(total_L1_count = sum(L1_count), n_organoid = n()) %>% ungroup()

pID_list <- df_count_pt_info %>%
  mutate(patientID = str_c(patientID," (",n_organoid,")")) %>%
  arrange(desc(total_L1_count)) %>% pull(patientID)

df_count %>% 
  mutate(patientID = str_sub(sampleID,1,4)) %>%
  left_join(df_count_pt_info) %>%
  mutate(patientID = str_c(patientID," (",n_organoid,")")) %>%
  mutate(group = ifelse(L1_count > 10, "11-20",
                        ifelse(L1_count >= 6, "6-10",
                               ifelse(L1_count >= 3, "3-5",
                                      ifelse(L1_count >=1, "1-2", "0"))))) %>%
  select(sampleID,patientID,group) %>%
  ggplot(aes(factor(patientID,levels=rev(pID_list)),fill=factor(group,levels=c("0","1-2","3-5","6-10","11-20"))))+
  geom_bar(position = "fill")+theme_classic()+coord_flip()+
  scale_fill_manual(values = c("#c5c5c5","#e5c53a","#e39139","#e93e27","#bb1d3e"))+
  xlab("")+ylab("Clones (%)")+
  theme(axis.line.y = element_blank(), axis.text = element_text(size=16), axis.title = element_text(size=16))+
  theme(legend.position = "top",legend.title = element_text(size=16), legend.text = element_text(size=16))+
  guides(fill=guide_legend(title="# of soL1Rs per clone",nrow=2))+
  scale_x_discrete(drop = FALSE,position="top")


####### Fig 1d (L1 count vs age)
df_count_pt <- df_count %>% 
  mutate(patientID = str_sub(sampleID,1,4)) %>% 
  group_by(patientID) %>% summarize(total_L1_count = sum(L1_count), min_L1_count = min(L1_count), max_L1_count = max(L1_count), N_organoid = n()) %>% ungroup() %>% 
  left_join(df_st1 %>% select(patientID,age) %>% distinct(), by=c("patientID")) %>% mutate(mean_L1_count = total_L1_count/N_organoid)

ggplot()+
  geom_point(data=df_count_pt, aes(age,mean_L1_count))+
  geom_text_repel(data=subset(df_count_pt, patientID %in% c("HC15","HC06")),aes(x=age,y=mean_L1_count,label = patientID),size=5)+
  geom_errorbar(data=df_count_pt,aes(x=age,y=mean_L1_count,ymin=min_L1_count, ymax=max_L1_count))+
  geom_smooth(data=df_count_pt %>% filter(!patientID %in% c("HC15","HC06")), aes(x=age,y=mean_L1_count),method = "lm", formula = y ~ x, level=0.95,fullrange=TRUE, color = "blue")+
  theme_classic()+
  theme(plot.title = element_text(size = 20))+
  theme(axis.text = element_text(size=16), axis.title = element_text(size=16))+
  xlab("Age (years)")+ylab("Ave. # of soL1Rs per clone")+
  xlim(0,100)+ylim(0,18)

######## Fig 1e, 1f (phylogenetic tree of HC14/HC19 : drawing)
######## Fig 1g (soL1R rate in various developmental stages and cell types : drawing) 
