###################################################
# 20/1/15
# Inyoung Jun
# Objective: study population
###################################################


##### Setting
setwd("S:/UF_IDR_AMR/KP")

library(tidyverse)
library(dplyr)
library(comorbidity)
library(table1)

`%!in%` = Negate(`%in%`)

daydiff <- function(a,b,c,d){
  result<-(a-b)*365+c-d
  return(result)
} #a:Later Year, b:Former Year, c:Later Date, d:Former Date



##### Import
d<-readRDS("data/dnew.RDS") #diagnosis
a<-readRDS("data/a.RDS") #antibiogram
demo<-readRDS("data/demo.RDS") #demographics
ad<-readRDS("data/ad.RDS") #admission
#m<-readRDS("data/m.RDS") #medication

###################################################

##### blood stream infection
d_BSI <- d %>% 
  rename(PID = Deidentified_Patient_ID, D_Year = Year, D_Date = Date_of_Service) %>% 
  filter(D_Year>=2011 & D_Year<=2018) %>% # Year 2011 - 2018
  filter(ICD_Code %in% c("R78.81","R65.20","R65.21","A41","A41.59","A41.89","A41.9")|
           ICD_Code2 %in% c("7907","99591","99592")|
           ICD_Code3 %in% c("038")) %>% 
  mutate(BSI_Year = D_Year, BSI_Date = D_Date) %>% 
  arrange(PID, BSI_Year, BSI_Date) %>% 
  select(-c(D_Year, D_Date, ICD_Type, Relation_to_Last_Antibiogram,ICD_Code3))

length(unique(d_BSI$PID)) #14106

##### KP infection

#KP by diagnosis
d_KP <- d %>% 
  rename(PID = Deidentified_Patient_ID, D_Year = Year, D_Date = Date_of_Service) %>% 
  filter(D_Year>=2011 & D_Year<=2018) %>% # Year 2011 - 2018
  filter(ICD_Code %in% c("B96.1", "041.3")) %>% 
  mutate(d_KP_Year = D_Year, d_KP_Date = D_Date) %>% 
  arrange(PID, d_KP_Year, d_KP_Date) %>%  
  select(PID, d_KP_Year, d_KP_Date)

length(unique(d_KP$PID)) #2992
head(d_KP)

#KP by antibiogram records
a_KP <- a %>% 
  rename(PID = Deidentified_Patient_ID, A_Year = Year, A_Date = Date_of_Service) %>% 
  filter(A_Year>=2011 & A_Year<=2018) %>%
  filter(str_detect(tolower(organism), "klebsiella")) %>% 
  mutate(a_KP_Year = A_Year, a_KP_Date = A_Date) %>% 
  arrange(PID, a_KP_Year, a_KP_Date) %>% 
  select(PID, a_KP_Year, a_KP_Date)

length(unique(a_KP$PID)) #9704

#Construct data_KP which is integrated KP records
data_KP <- a_KP %>% filter(PID %in% d_BSI$PID) %>% 
           left_join(d_KP, by="PID") %>% 
           distinct(PID, a_KP_Year, a_KP_Date, d_KP_Year, d_KP_Date)
head(data_KP)
length(unique(data_KP$PID)) #3590

##### KP infection records after 0 ~ 7 days from BSI record
d_BSIKP<- full_join(d_BSI, data_KP, by = c("PID")) %>% 
  select(PID, BSI_Year, BSI_Date, a_KP_Year, a_KP_Date, d_KP_Year, d_KP_Date) %>%
  arrange(PID, BSI_Year, BSI_Date) %>% 
  filter(daydiff(BSI_Year,d_KP_Year,BSI_Date,d_KP_Date)<=7 & daydiff(BSI_Year,d_KP_Year,BSI_Date,d_KP_Date)>=0|
         daydiff(BSI_Year,a_KP_Year,BSI_Date,a_KP_Date)<=7 & daydiff(BSI_Year,a_KP_Year,BSI_Date,a_KP_Date)>=0) %>% 
  distinct(PID, BSI_Year, BSI_Date, d_KP_Year, d_KP_Date, a_KP_Year, a_KP_Date) %>% 
  mutate(KP_Year = BSI_Year) %>% 
  mutate(KP_Date = ifelse(is.na(d_KP_Date)==TRUE,a_KP_Date,ifelse(abs(BSI_Date-a_KP_Date)<=abs(BSI_Date-d_KP_Date), a_KP_Date, d_KP_Date))) %>% 
  arrange(PID, BSI_Year, BSI_Date, KP_Year, KP_Date) %>% 
  distinct(PID, BSI_Year, BSI_Date,.keep_all = TRUE) 

#saveRDS(d_BSIKP,"data_modified/d_BSIKP.RDS")
#d_BSIKP<-readRDS("data_modified/d_BSIKP.RDS")

head(d_BSIKP)
length(unique(d_BSIKP$PID))#1942

##### Exclude not adults
d_BSIKP_adult <- d_BSIKP %>% 
                 select(PID, BSI_Year,BSI_Date) %>% 
                 left_join(demo, by="PID") %>% 
                 mutate(BSI_Age = Current_Age - (2019 - BSI_Year)) %>% 
                 filter(BSI_Age>=18) %>% 
                 arrange(PID, BSI_Year, BSI_Date) %>% 
                 distinct(PID,BSI_Year,BSI_Date,.keep_all = TRUE)

length(unique(d_BSIKP_adult$PID))#1938


##### Exclude medical history <2

d_BSIKP_adult_mh2 <-d %>% 
                    rename(PID = Deidentified_Patient_ID, D_Year = Year, D_Date = Date_of_Service) %>% 
                    filter(PID %in% d_BSIKP_adult$PID) %>% 
                    arrange(PID, BSI_Year, BSI_Date) %>% 
                    mutate(DIFF_BSI_D = daydiff(BSI_Year,D_Year,BSI_Date,D_Date)) %>% 
                    mutate(Medical2yr=ifelse(DIFF_BSI_D >=730,1,0)) %>% 
                    filter(Medical2yr == 1)             

d_BSIKP_adult_mh2_BSI <- d_BSIKP_adult %>%
                         filter(PID %in% d_BSIKP_adult_mh2$PID) %>% 
                         group_by(PID) %>% 
                         mutate(P_BSI_Year = lag(BSI_Year, order_by = PID)) %>% 
                         mutate(P_BSI_Date = lag(BSI_Date, order_by = PID)) %>% 
                         mutate(Event_same = ifelse(daydiff(BSI_Year,P_BSI_Year,BSI_Date,P_BSI_Date)<=30,1,0)) %>% 
                         mutate(Event_same = ifelse(is.na(Event_same)==TRUE,1,Event_same)) %>% 
                         group_by(PID, Event_same) %>% 
                         arrange(PID, BSI_Year,BSI_Date) %>% 
                       #  distinct(PID,Event_same,.keep_all = TRUE) %>% 
                       #  arrange(PID, BSI_Year, BSI_Date) %>% 
                         group_by(PID) %>% 
                         mutate(Record_num = seq(1:length(PID)))
                         
length(unique(d_BSIKP_adult_mh2_BSI$PID))
summary(d_BSIKP_adult_mh2_BSI$BSI_Age)
head(d_BSIKP_adult_mh2_BSI,30)

#########Attach first date of each episode

first_bsi_date1 <- d_BSIKP_adult_mh2_BSI %>% filter(Event_same == 1) %>% 
                  arrange(PID, BSI_Year, BSI_Date, P_BSI_Year, P_BSI_Date) %>% 
                  group_by(PID) %>% 
                  mutate(Event_num = 1) %>% 
                  mutate(First_BSI_Year = BSI_Year, First_BSI_Date = BSI_Date) %>% 
                  distinct(PID,.keep_all = TRUE) %>% 
                  select(PID, Event_num, First_BSI_Year, First_BSI_Date)
table(first_bsi_date1$Event_num)

first_bsi_date2 <- d_BSIKP_adult_mh2_BSI %>% filter(Event_same != 1) %>% 
                   arrange(PID, BSI_Year, BSI_Date, P_BSI_Year, P_BSI_Date) %>% 
                   group_by(PID) %>% 
                   mutate(Event_num = seq(1:length(PID))+1) %>% 
                   mutate(First_BSI_Year = BSI_Year, First_BSI_Date = BSI_Date) %>% 
                   distinct(PID,Event_num, .keep_all = TRUE) %>% 
                   select(PID, Event_num, First_BSI_Year, First_BSI_Date)
table(first_bsi_date2$Event_num)

#United single episodes + multiple episodes
first_bsi_date <- rbind(first_bsi_date1, first_bsi_date2) %>% 
                  arrange(PID,First_BSI_Year, First_BSI_Date) %>% 
                  distinct(PID,Event_num, .keep_all = TRUE)
head(first_bsi_date)
table(first_bsi_date$Event_num)

first_bsi_date[first_bsi_date$PID==985,]

#########Attach last date of each episode

last_bsi_date <- d_BSIKP_adult_mh2_BSI %>% 
  mutate(Event_same2=ifelse(Event_same==1,0,1)) %>%
  group_by(PID) %>% 
  mutate(cumsum = cumsum(Event_same2)) %>% 
  arrange(PID,-BSI_Year,-BSI_Date,cumsum) %>% 
  distinct(PID,cumsum,.keep_all = TRUE) %>% 
  arrange(PID, BSI_Year,BSI_Date) %>% 
  mutate(Event_num =seq(1:length(PID))) %>% 
  mutate(Last_BSI_Year=BSI_Year, Last_BSI_Date=BSI_Date) %>% 
  select(PID, Event_num, Last_BSI_Year, Last_BSI_Date)

date <- first_bsi_date %>% left_join(last_bsi_date,by=c("PID","Event_num"))
head(date)

d_BSIKP_final <- d_BSIKP_adult_mh2_BSI %>% 
                 left_join(select(date,c("PID","First_BSI_Year","First_BSI_Date","Last_BSI_Year","Last_BSI_Date","Event_num")), by=c("PID")) %>% 
                 distinct(PID,Event_num,First_BSI_Year,First_BSI_Date,Last_BSI_Year,Last_BSI_Date,.keep_all = TRUE) 
head(d_BSIKP_final)

saveRDS(d_BSIKP_final,"data_modified/d_BSIKP_final.RDS")
d_BSIKP_final<-readRDS("data_modified/d_BSIKP_final.RDS")
head(d_BSIKP_final,30)
length(unique(d_BSIKP_final$PID))#1400
summary(d_BSIKP_final$BSI_Age)

##### Admission and Discharge
ad_BSIKP<- ad %>% 
           rename(PID = Deidentified_Patient_ID, AD_Year = Year_of_Admit, AD_Date = admit_date, 
           DC_Year = Year_of_Discharge, DC_Date = discharge_date, ICU = icu_y_n) %>% 
           select(PID, AD_Year, AD_Date, DC_Year, DC_Date, ICU) %>% 
           filter(PID %in% d_BSIKP_final$PID) 
head(ad_BSIKP)
length(unique(ad_BSIKP$PID))#1400

d_BSIKP_ad<-  d_BSIKP_final %>%  
              left_join(select(demo,c(PID,Current_Age)),by="PID") %>% 
              full_join(ad_BSIKP, by = c("PID")) %>% 
              arrange(PID, Event_num, AD_Year, AD_Date) %>% 
              filter(abs(daydiff(First_BSI_Year,AD_Year,First_BSI_Date,AD_Date))<=7 |
                     abs(daydiff(DC_Year,Last_BSI_Year,DC_Date,Last_BSI_Date))<=7 ) %>% 
              distinct(PID, Event_num,.keep_all = TRUE)
head(d_BSIKP_ad)

length(unique(d_BSIKP_ad$PID))#Admission records 1253


##### Readmission
Readmission <- d_BSIKP_final %>%  
               mutate(mortality30 = ifelse(daydiff(Death_Year,First_BSI_Year,Death_Date,First_BSI_Date)<=30,1,0))%>% 
               mutate_at(vars(mortality30), ~replace(., is.na(.), 0)) %>% 
               filter(mortality30==0) %>% 
               mutate(BSI_Age = Current_Age - (2019 - First_BSI_Year)) %>% 
               filter(BSI_Age>=18) %>% 
               full_join(ad_BSIKP, by = c("PID")) %>% 
               arrange(PID, Event_num, AD_Year, AD_Date) %>% 
               filter(abs(daydiff(First_BSI_Year,AD_Year,First_BSI_Date,AD_Date))<=7 |
                       abs(daydiff(DC_Year,Last_BSI_Year,DC_Date,Last_BSI_Date))<=7 ) %>% 
               filter(Event_num==1) %>% 
               mutate(P_DC_Year = lag(DC_Year, order_by = PID)) %>% 
               mutate(P_DC_Date = lag(DC_Date, order_by = PID)) %>% 
               mutate(readmission = ifelse(
               daydiff(AD_Year,P_DC_Year,AD_Date,P_DC_Date)<=30 & daydiff(AD_Year,P_DC_Year,AD_Date,P_DC_Date)>=0,1,0)) %>% 
               arrange(PID,-readmission) %>% 
               distinct(PID,.keep_all = TRUE) %>% 
               select(PID, AD_Year,AD_Date,DC_Year,DC_Date,P_DC_Year,P_DC_Date,Event_num,readmission,Death_Year,Death_Date,mortality30) %>% 
               filter(readmission==1)

length(unique(Readmission$PID))#42 labeled as one
head(Readmission)
table(Readmission$mortality30,useNA = "ifany")

#####Recurrence
Recurrence <- d_BSIKP_final %>% 
              mutate(mortality30 = ifelse(daydiff(Death_Year,First_BSI_Year,Death_Date,First_BSI_Date)<=30,1,0)) %>% 
              mutate_at(vars(mortality30), ~replace(., is.na(.), 0)) %>% 
              filter(mortality30==0) %>% 
              filter(Event_num>=2) %>% arrange(PID,-Event_num) %>%
              distinct(PID,.keep_all = TRUE) %>% 
              mutate(recurrence = 1) %>% 
              rename(Max_Event_num = Event_num)

length(unique(Recurrence$PID))#143 labeled as one
head(Recurrence)
table(Recurrence$mortality30,useNA = "ifany")


study_pop <- d_BSIKP_final %>% 
             select(PID,Event_num,First_BSI_Year,First_BSI_Date,Last_BSI_Year,Last_BSI_Date) %>%
             left_join(select(d_BSIKP_ad,c(PID,Event_num,AD_Year,AD_Date,DC_Year,DC_Date,ICU)),by=c("PID","Event_num"))%>%
             left_join(select(Readmission, c(PID,readmission)),by="PID") %>% 
             left_join(select(Recurrence,c(PID,recurrence,Max_Event_num)), by="PID") %>% 
             left_join(select(demo,c(PID,Current_Age,Death_Year,Death_Date)),by="PID") %>%
             mutate(BSI_Age = Current_Age - (2019 - First_BSI_Year)) %>% 
             mutate(mortality30 = ifelse(daydiff(Death_Year,First_BSI_Year,Death_Date,First_BSI_Date)<=30,1,0)) %>% 
             mutate_at(vars(mortality30,recurrence,readmission), ~replace(., is.na(.), 0)) %>% 
             mutate_at(vars(Max_Event_num), ~replace(., is.na(.), 1)) %>% 
             mutate(Outpatient=ifelse(is.na(AD_Year)==TRUE, 1,0)) #1400 

length(unique(study_pop$PID))#1400

table(study_pop$mortality30)
table(study_pop$recurrence,useNA="ifany")
table(study_pop$readmission,useNA="ifany")

head(study_pop)
saveRDS(study_pop,"data_modified/study_pop.RDS")
