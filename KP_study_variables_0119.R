###################################################
# 20/1/16
# Inyoung Jun
# Objective: Creating variables
###################################################
rm(list=ls())

##### Setting
setwd("S:/UF_IDR_AMR/KP")

library(tidyverse)
library(dplyr)
library(comorbidity)
library(table1)
library(reshape2)
library(tidyselect)

`%!in%` = Negate(`%in%`)

daydiff <- function(a,b,c,d){
  result<-(a-b)*365+c-d
  return(result)
} #a:Later Year, b:Former Year, c:Later Date, d:Former Date



##### Import
d<-readRDS("data/dnew.RDS") #diagnosis
a<-readRDS("data/a.RDS") #antibiogram
demo<-readRDS("data/demo.RDS") #demographics
#ad<-readRDS("data/ad.RDS") #admission
m<-readRDS("data/m.RDS") #medication



##### Start from study population dataset
study_pop<-readRDS("data_modified/study_pop.RDS")
ID <- unique(study_pop$PID)#1400
head(study_pop)


##### For this PHHP research day abstract, focuing on the first event
study_pop_first <- study_pop %>% filter(Event_num ==1) 
                  

table(study_pop_first$Outpatient)
168+1232#1400
str(study_pop_first)

table(study_pop_first$AD_Year,useNA = "ifany")#168 NA
table(study_pop_first$DC_Year,useNA = "ifany")#168 NA


##### Demographics
demo_pop <- demo %>% filter(PID %in% ID) %>% 
                 mutate(Race_group = case_when(Race == "WHITE" ~ "WHITE",
                                               Race == "BLACK" ~ "BLACK",
                                               Race == "ASIAN" ~ "ASIAN",
                                               TRUE ~ "OTHER")) %>% 
                 mutate(Ethnicity_group = case_when(Ethnicity == "HISPANIC" ~ "HISPANIC",
                                                                    TRUE ~  "NOT HISPANIC"))
saveRDS(demo_pop,"data_modified/demo_pop.RDS")

##### Comorbidities
d_comorbidity<- d %>% 
                rename(PID = Deidentified_Patient_ID, D_Year = Year, D_Date = Date_of_Service) %>%
                select(PID, D_Year, D_Date, ICD_Code, ICD_Code2, ICD_Description)%>% 
                filter(PID %in% ID) %>% 
                left_join(select(study_pop_first,c(PID,First_BSI_Year,First_BSI_Date)),by="PID") %>% 
                filter(daydiff(First_BSI_Year,D_Year,First_BSI_Date,D_Date)>=365)

length(unique(d_comorbidity$PID))#1399
#PID 98298 does not have any comorbidity records

CI_ID<-d_comorbidity[,c("PID", "ICD_Code2")]
CI<-CI_ID %>% 
  rename(
    id  = PID,
    code  = ICD_Code2
  )
str(CI)

CI_pop<-comorbidity(
  x=CI,
  id="id",
  code="code",
  score="charlson",
  assign0=TRUE,
  icd = "icd9",
  factorise = FALSE,
  labelled = TRUE,
  tidy.codes = FALSE
) %>% rename(PID = id)

head(CI_pop)
length(CI_pop$PID) #1399 obs 

saveRDS(CI_pop,"data_modified/CI_pop.RDS")
head(study_pop_first,20)

##### Nosocomial (infected >2 days after admission)
NOSO_pop <- study_pop_first %>% mutate(Noso = ifelse(daydiff(First_BSI_Year, AD_Year, First_BSI_Date, AD_Date)>2,1,0)) %>% 
                          select(PID, First_BSI_Year, First_BSI_Date, AD_Year, AD_Date, Noso)
head(NOSO_pop)
table(NOSO_pop$Noso)
313/(313+919)#25.4

##### ESBL phenotype 
########## considering
a_ESBL <- a %>%
  rename(PID = Deidentified_Patient_ID, A_Year = Year, A_Date = Date_of_Service) %>% 
  filter(PID %in% ID) %>% 
  filter(A_Year>=2011 & A_Year<=2018) %>%
  filter(antibiotic =="ESBL Screen") %>% 
  filter(str_detect(tolower(organism), "klebsiella")) %>% 
  select(PID, A_Year,A_Date,antibiotic,suscept) %>% 
  left_join(select(study_pop_first,c(PID,First_BSI_Year,First_BSI_Date,Last_BSI_Year,Last_BSI_Date,AD_Year,AD_Date,DC_Year,DC_Date,Outpatient))) %>% 
#  filter(is.na(BSI_Year)==FALSE) %>% 
  mutate(ESBL_upon = ifelse(
                           (Outpatient==0 & daydiff(A_Year, AD_Year, A_Date, AD_Date)>=0 & daydiff(DC_Year, A_Year, DC_Date, A_Date)>=0)|
                           (Outpatient==1 & daydiff(A_Year, First_BSI_Year, A_Date, AD_Date)>=0 & daydiff(Last_BSI_Year, A_Year, Last_BSI_Date, A_Date)>=0)
                           ,1,0)) %>%
  mutate(ESBL_before = ifelse(daydiff(A_Year,First_BSI_Year,A_Date,First_BSI_Date)<0, 1, 0)) %>% 
  mutate(ESBL_test = ifelse(is.na(ESBL_upon)==FALSE,ESBL_upon, ESBL_before)) %>% 
  filter(ESBL_test ==1) %>% 
  distinct(PID, A_Year, A_Date,.keep_all = TRUE)
  
head(a_ESBL)
table(a_ESBL$ESBL_test)#411
table(a_ESBL$ESBL_upon)#387
table(a_ESBL$ESBL_before)#145 or 266

table(a_ESBL$A_Year)#from 2016
table(a_ESBL$suscept,useNA="ifany") #moderately sensitive/Positive/Negative
length(unique(a_ESBL$PID))#363 ESBL test results that is upon diagnosis or before BSI

head(a_ESBL)
ESBL_pop<- as.data.frame(ID) %>% rename(PID = ID) %>% 
          left_join(select(a_ESBL,c("PID","A_Year","A_Date","suscept")),by ="PID") %>% 
          rename(ESBL_Year=A_Year,ESBL_Date=A_Date,ESBL_sus=suscept) %>% 
          mutate(ESBL_result=case_when(ESBL_sus=="Positive"~"Positive",
                                       ESBL_sus=="Negative"~"Negative",
                                       ESBL_sus=="Moderately sensitive"~"Moderately sensitive",
                                       TRUE ~ "Unknown")) %>% 
          distinct(PID,ESBL_result,.keep_all=TRUE) %>% 
          arrange(PID,-ESBL_Year,-ESBL_Date) %>% 
          distinct(PID,.keep_all = TRUE)

table(ESBL_pop$ESBL_result)#65 positive
288+65+1047#1400
length(unique(ESBL_pop$PID))#1400


##### Therapeutic medication (within the first five days after the first diagnosis)
m_pop<-m %>% 
       filter(PID %in% ID) %>% 
       select(PID, M_Year, M_Date, Simple_Generic_Med_Name) %>% 
       left_join(select(study_pop_first,c("PID","First_BSI_Year","First_BSI_Date","Last_BSI_Year","Last_BSI_Date",
                                          "AD_Year","AD_Date","DC_Year","DC_Date","Outpatient"))) %>% 
       mutate(trt_upon = ifelse((Outpatient==0 & 
                                  daydiff(M_Year, AD_Year, M_Date, AD_Date)>=0 & 
                                  daydiff(DC_Year, M_Year, DC_Date, M_Date)>=0)|
                                (Outpatient==1 & 
                                 daydiff(M_Year, First_BSI_Year, M_Date, First_BSI_Date)>=0 & 
                                 daydiff(Last_BSI_Year, M_Year, Last_BSI_Date, M_Date)>=0),1,0)) %>%
       mutate(trt_before = ifelse(daydiff(M_Year,First_BSI_Year,M_Date,First_BSI_Date)<0, 1, 0)) 

length(unique(m_pop$PID[m_pop$trt_upon==1]))#1342
length(unique(m_pop$PID[m_pop$trt_before==1]))#1252
length(unique(m_pop$PID))#1390


m_pop_upon <- m_pop %>% 
                   filter(trt_upon==1) %>% 
                   mutate(upon_acetaminophen=ifelse(str_detect(tolower(Simple_Generic_Med_Name), "acetaminophen"),1,0)) %>% 
                   mutate(upon_pantoprazole=ifelse(str_detect(tolower(Simple_Generic_Med_Name), "pantoprazole"),1,0)) %>% 
                   mutate(upon_aspirin=ifelse(str_detect(tolower(Simple_Generic_Med_Name), "aspirin"),1,0)) %>% 
                   mutate(upon_opium=ifelse(str_detect(tolower(Simple_Generic_Med_Name), c("oxycodone|fentanyl|morphine")),1,0)) %>% 
                   group_by(PID) %>% 
                   summarise_at(vars(upon_acetaminophen,upon_pantoprazole,upon_aspirin,upon_opium), funs(sum)) %>%
                   mutate(PID=as.character(PID)) %>% 
                   mutate_if(is.numeric, ~1 * (. > 0)) %>% 
                   mutate_if(is.numeric,as.factor)
                             
length(unique(m_pop_upon$PID))#1342               
summary(m_pop_upon)


m_pop_bf <- m_pop %>% 
                   filter(trt_before==1) %>% 
                   mutate(bf_acetaminophen=ifelse(str_detect(tolower(Simple_Generic_Med_Name), "acetaminophen"),1,0)) %>% 
                   mutate(bf_pantoprazole=ifelse(str_detect(tolower(Simple_Generic_Med_Name), "pantoprazole"),1,0)) %>% 
                   mutate(bf_aspirin=ifelse(str_detect(tolower(Simple_Generic_Med_Name), "aspirin"),1,0)) %>% 
                   mutate(bf_opium=ifelse(str_detect(tolower(Simple_Generic_Med_Name), c("oxycodone|fentanyl|morphine")),1,0)) %>% 
                   group_by(PID, bf_acetaminophen,bf_pantoprazole,bf_aspirin,bf_opium) %>% 
                   group_by(PID) %>% 
                   summarise_at(vars(bf_acetaminophen,bf_pantoprazole,bf_aspirin,bf_opium), funs(sum)) %>%
                   mutate(PID=as.character(PID)) %>% 
                   mutate_if(is.numeric, ~1 * (. > 0)) %>% 
                   mutate_if(is.numeric,as.factor)

summary(m_pop_bf)

trt_pop<- study_pop_first %>% 
          select(PID) %>%
          mutate(PID = as.character(PID)) %>% 
          left_join(m_pop_upon, by="PID") %>% 
          left_join(m_pop_bf, by="PID") %>% 
          mutate(PID = as.integer(PID))

saveRDS(trt_pop,"data_modified/trt_pop.RDS")

##### Antibiogram resistance
#KP by antibiogram records
a_pop <- a %>% 
  rename(PID = Deidentified_Patient_ID, A_Year = Year, A_Date = Date_of_Service) %>% 
  filter(PID %in% ID) %>% 
  filter(str_detect(tolower(organism), "klebsiella")) %>% 
  select(PID, A_Year, A_Date, specimen_type, specimen_source, antibiotic, suscept) %>% 
  filter(suscept=="Resistant") %>%
  left_join(select(study_pop_first,c("PID","First_BSI_Year","First_BSI_Date",
                                     "Last_BSI_Year","Last_BSI_Date",
                                     "AD_Year","AD_Date","DC_Year","DC_Date","Outpatient"))) %>%  
  mutate(ab_upon = ifelse((Outpatient==0 & 
                              daydiff(A_Year, AD_Year, A_Date, AD_Date)>=0 & 
                              daydiff(DC_Year, A_Year, DC_Date, A_Date)>=0)|
                             (Outpatient==1 & 
                                daydiff(A_Year, First_BSI_Year, A_Date, First_BSI_Date)>=0 & 
                                daydiff(Last_BSI_Year, A_Year, Last_BSI_Date, A_Date)>=0),1,0)) %>% 
  mutate(ab_before = ifelse(daydiff(A_Year,First_BSI_Year,A_Date,First_BSI_Date)<0, 1, 0)) 
  
length(unique(a_pop$PID)) #1356 showed resistance to at least one antibiotic

a_pop_upon <- a_pop %>% filter(ab_upon== 1) %>% 
          distinct(PID, antibiotic,.keep_all = TRUE) %>% 
          mutate(Var = 1) %>% dcast(PID~antibiotic) %>% 
          replace(is.na(.),0) %>% 
          mutate(mdrcheck=rowSums(.[-1])) %>% mutate(MDR = ifelse(mdrcheck>=2,1,0)) %>% 
          rename_at(vars(-PID),function(x) paste0("res_upon_",x))

head(a_pop_upon)

a_pop_bf <- a_pop %>% filter(ab_before == 1) %>% 
          distinct(PID, antibiotic,.keep_all = TRUE) %>% 
          mutate(Var = 1) %>% dcast(PID~antibiotic) %>% 
          replace(is.na(.),0) %>% 
          mutate(mdrcheck=rowSums(.[-1])) %>% mutate(MDR = ifelse(mdrcheck>=2,1,0)) %>% 
          rename_at(vars(-PID),function(x) paste0("res_bf_",x))
head(a_pop_bf)

mdr_pop<- study_pop_first %>% select(PID) %>% 
          left_join(a_pop_upon, by="PID") %>% 
          left_join(a_pop_bf, by="PID") %>% 
          replace(is.na(.),0)

head(mdr_pop)
glimpse(mdr_pop)
saveRDS(mdr_pop, "data_modified/mdr_pop.RDS")
head(mdr_pop)


#Attach class
top_res_rank_class<-read.csv("data_modified/top_res_rank_class.csv")

table(top_res_rank_class$Antibiotics_class)
aminoglycosides<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="aminoglycosides"]
betalactams<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="betalactams"]
fluoroquinolones<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="fluoroquinolones"]
glycopeptides<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="glycopeptides"]
sulfonamides<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="sulfonamides"]
tetracyclines<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="tetracyclines"]
others<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="others"]


mdr_pop_top_AB <- mdr_pop %>%
                  mutate(class_res_upon_aminoglycosides = ifelse(
                  sum(c_across(intersect(starts_with("res_upon"),
                                        contains(c(aminoglycosides)))))>0,1,0)) %>% 
                  mutate(class_res_upon_betalactams = ifelse(
                  sum(c_across(intersect(starts_with("res_upon"),
                                        contains(c(betalactams)))))>0,1,0)) %>% 
                  mutate(class_res_upon_fluoroquinolones = ifelse(
                    sum(c_across(intersect(starts_with("res_upon"),
                                           contains(c(fluoroquinolones)))))>0,1,0)) %>% 
                  mutate(class_res_upon_glycopeptides = ifelse(
                    sum(c_across(intersect(starts_with("res_upon"),
                                           contains(c(glycopeptides)))))>0,1,0)) %>%
                  mutate(class_res_upon_sulfonamides = ifelse(
                    sum(c_across(intersect(starts_with("res_upon"),
                                           contains(c(sulfonamides)))))>0,1,0)) %>%
                  mutate(class_res_upon_tetracyclines = ifelse(
                    sum(c_across(intersect(starts_with("res_upon"),
                                           contains(c(tetracyclines)))))>0,1,0)) %>%
                  mutate(class_res_upon_others = ifelse(
                    sum(c_across(intersect(starts_with("res_upon"),
                                           contains(c(others)))))>0,1,0)) %>% 
                  mutate(class_res_bf_aminoglycosides = ifelse(
                    sum(c_across(intersect(starts_with("res_bf"),
                                           contains(c(aminoglycosides)))))>0,1,0)) %>% 
                  mutate(class_res_bf_betalactams = ifelse(
                    sum(c_across(intersect(starts_with("res_bf"),
                                           contains(c(betalactams)))))>0,1,0)) %>% 
                  mutate(class_res_bf_fluoroquinolones = ifelse(
                    sum(c_across(intersect(starts_with("res_bf"),
                                           contains(c(fluoroquinolones)))))>0,1,0)) %>% 
                  mutate(class_res_bf_glycopeptides = ifelse(
                    sum(c_across(intersect(starts_with("res_bf"),
                                           contains(c(glycopeptides)))))>0,1,0)) %>%
                  mutate(class_res_bf_sulfonamides = ifelse(
                    sum(c_across(intersect(starts_with("res_bf"),
                                           contains(c(sulfonamides)))))>0,1,0)) %>%
                  mutate(class_res_bf_tetracyclines = ifelse(
                    sum(c_across(intersect(starts_with("res_bf"),
                                           contains(c(tetracyclines)))))>0,1,0)) %>%
                  mutate(class_res_bf_others = ifelse(
                    sum(c_across(intersect(starts_with("res_bf"),
                                           contains(c(others)))))>0,1,0)) %>% 
                  select(PID, class_res_bf_aminoglycosides,
                              class_res_bf_betalactams,
                              class_res_bf_fluoroquinolones,
                              class_res_bf_glycopeptides,
                              class_res_bf_sulfonamides,
                              class_res_bf_tetracyclines,
                              class_res_bf_others,
                              class_res_upon_aminoglycosides,
                              class_res_upon_betalactams,
                              class_res_upon_fluoroquinolones,
                              class_res_upon_glycopeptides,
                              class_res_upon_sulfonamides,
                              class_res_upon_tetracyclines,
                              class_res_upon_others,
                              res_bf_mdrcheck,
                              res_bf_MDR,
                              res_upon_mdrcheck,
                              res_upon_MDR)
                  
head(mdr_pop_top_AB)
saveRDS(mdr_pop_top_AB,"data_modified/mdr_pop_top_AB.RDS")
# #Check Top antibiotics that are resistant
# 
# # transpose all but the first column (name)
# top_res <- as.data.frame(t(mdr_pop[,-1])) %>% 
#            mutate(Res_sum = rowSums(.)) %>% 
#            mutate(Antibiotics = colnames(mdr_pop[,-1])) %>% 
#            select(Antibiotics,Res_sum) %>% 
#            arrange(-Res_sum)
#       
# top_res_upon <- top_res %>% 
#                 filter(str_detect(Antibiotics, "res_upon")) %>% 
#                 filter(!str_detect(Antibiotics,"mdr|MDR")) %>%  
#                 rename(Upon_Res_sum = Res_sum) %>% 
#                 mutate(Antibiotics = trimws(Antibiotics, which="both")) %>%  
#                 mutate(Antibiotics = substr(Antibiotics, 10, 90))
# 
# top_res_bf <- top_res %>% 
#               filter(str_detect(Antibiotics, "res_bf")) %>% 
#               filter(!str_detect(Antibiotics,"mdr|MDR")) %>% 
#               rename(Bf_Res_sum = Res_sum) %>% 
#               mutate(Antibiotics = trimws(Antibiotics, which="both")) %>%  
#               mutate(Antibiotics = substr(Antibiotics, 8, 90))
# 
# top_res_rank <- top_res_upon %>% 
#                 left_join(top_res_bf, by="Antibiotics") %>% 
#                 select(Antibiotics, Bf_Res_sum, Upon_Res_sum)
#   
# top_res_rank
# write.csv(top_res_rank,"data_modified/top_res_rank.csv",row.names=FALSE)
# saveRDS(top_res_rank,"data_modified/top_res_rank.RDS")

#####Site of infection
specimen_source_pop<-read.csv("data_modified/specimen_source_pop_1210.csv")

infectionsite<- a %>% 
                rename(PID = Deidentified_Patient_ID, A_Year = Year, A_Date = Date_of_Service) %>% 
                filter(PID %in% study_pop$PID) %>%
                filter(str_detect(tolower(organism), "klebsiella")) %>% 
                mutate(a_KP_Year = A_Year, a_KP_Date = A_Date) %>% 
                arrange(PID, a_KP_Year, a_KP_Date) %>% 
                select(PID, A_Year,A_Date, specimen_type, specimen_source) %>% 
                left_join(select(study_pop_first,c(PID|First_BSI_Year|First_BSI_Date|Last_BSI_Year|Last_BSI_Date|AD_Year|AD_Date|DC_Year|DC_Date)),by="PID") %>% 
                filter(abs(daydiff(A_Year,AD_Year,A_Date,AD_Date))<=7|abs(daydiff(A_Year,First_BSI_Year,A_Date,First_BSI_Date))<=7) %>% 
                filter(daydiff(DC_Year,A_Year,DC_Date,A_Date)>=0|daydiff(Last_BSI_Year,A_Year,Last_BSI_Date,A_Date)>=0) %>% 
                left_join(specimen_source_pop, by="specimen_source") %>% 
                mutate(inf_site1=ifelse(specimen_type!="", paste0(specimen_type),
                                        paste0(class_specimen_source))) %>% 
                distinct(PID,inf_site1,.keep_all = TRUE) %>% 
                mutate(infection_site = ifelse(str_detect(tolower(inf_site1),c("blood|venipuncture")),"Site_Blood",
                                 ifelse(str_detect(tolower(inf_site1),c("urine|catheter")),"Site_Urine",
                                        ifelse(str_detect(tolower(inf_site1),c("skin|wound|tissue")),"Site_Skin",
                                               ifelse(str_detect(tolower(inf_site1),c("sputum|throat|lung")),"Site_Throat_lung",
                                                      ifelse(str_detect(tolower(inf_site1),c("intra-abdominal")),"Site_Intra_abdominal",
                                                            ifelse(str_detect(tolower(inf_site1),c("other|fluid|bone|cardio")),"Site_Other","Site_Unknown"))))))) %>% 
  mutate(Var = 1) %>% dcast(PID~infection_site) %>%     
  mutate(Site_Throat_lung = replace(Site_Throat_lung, Site_Throat_lung >=1, 1)) %>% 
  mutate(Site_Other = replace(Site_Other, Site_Other >=1, 1)) %>% 
  mutate_if(is.numeric,as.factor) %>% 
  mutate(PID=as.integer(paste(PID)))


summary(infectionsite)                       
head(infectionsite)
str(infectionsite)
saveRDS(infectionsite,"data_modified/infectionsite.RDS")
                    

length(unique(infectionsite$PID))#1319                
####################################################################
#####Attach variable
data <- study_pop_first %>% 
        left_join(select(demo_pop,-c("Current_Age","Death_Year","Death_Date")), by="PID") %>% 
        left_join(CI_pop, by="PID") %>%
        left_join(select(NOSO_pop,c("PID","Noso")), by="PID") %>%
        left_join(trt_pop, by="PID") %>%
        left_join(mdr_pop_top_AB, by="PID") %>%
        left_join(select(ESBL_pop,c("PID","ESBL_result")),by="PID") %>%
        left_join(infectionsite,by="PID") %>% 
        mutate_all(~replace(., is.na(.), 0))

head(study_pop_first)
head(data)
glimpse(data)
write.csv(data,"data_modified/data.csv",row.names=FALSE)
saveRDS(data,"data_modified/data.RDS")
glimpse(data)

table(data$mortality30)
198+1202 #1400
write.csv(data,"data_modified/data.csv",row.names = FALSE)
