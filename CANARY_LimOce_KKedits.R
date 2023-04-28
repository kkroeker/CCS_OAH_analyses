#Written by Emily Donham 11/1/2020, Edited by KK 2023
#Plots processed pH data

######################################################################################################
######################################################################################################
library(plyr); library(dplyr);library(broom);library(ggplot2); library(lubridate);library(LoLinR); library(stringr)
library(lme4); library(lmerTest); library(multcomp); library(phytotools); library(googledrive); library(Rmisc)
library(tibble); library(ggpubr); library(wesanderson); library(tidyverse);library(vegan);
library(lsmeans); library(RLRsim); library(pracma); library(chron); library(caTools); library(TTR);
library(zoo); library(datetimeutils); library(metR); library(cowplot); library(khroma); library(grid)
######################################################################################################
######################################################################################################

rm(list = ls())
######################################################################################################
pHDat = read.csv(paste('/Volumes/GoogleDrive/My Drive/UCSC/R/PA_all.csv'),
                 skip = 7)
pHDat$Lat = as.numeric(pHDat$Lat)
pHDat$Lon = as.numeric(pHDat$Lon)
pHDat$Depth = as.numeric(pHDat$Depth)
pHDat$pH = as.numeric(pHDat$pH)
pHDat$pH.Temp = as.numeric(pHDat$pH.Temp)
pHDat$QF_pH = as.numeric(pHDat$QF_pH )
pHDat$QF_TC = as.numeric(pHDat$QF_TC)
pHDat$DO..umol.kg. = as.numeric(pHDat$DO..umol.kg.)
pHDat$DO.mgL <- pHDat$DO..umol.kg./31.2512
pHDat$DO..Sat = as.numeric(pHDat$DO..Sat)
pHDat$DO.Temp = as.numeric(pHDat$DO.Temp)
pHDat$QF_DO = as.numeric(pHDat$QF_DO)
pHDat$year = as.numeric(str_sub(pHDat$pH_Date_UTC_ddmmyyyy,-4,-1))
pHDat$month = as.numeric(str_sub(pHDat$pH_Date_UTC_ddmmyyyy,-6,-5))
pHDat$day = as.numeric(str_sub(pHDat$pH_Date_UTC_ddmmyyyy,-8,-7))
pHDat$DT = as.POSIXct(paste(pHDat$day, pHDat$month, pHDat$year, 
                            pHDat$pH_Time_UTC_hhmmss), format = "%d%m%Y %H:%M:%S", tz = "UTC")

pHDat <- pHDat[complete.cases(pHDat$Lat),]
for (i in 1:nrow(pHDat)) {
  if(pHDat$DT[i] < '2018-05-03 6:00:00' & pHDat$DT[i] > '2017-11-02 22:00') {
    pHDat$DEP[i] = 1
  }else if(pHDat$DT[i] > '2018-05-03 22:00:00' & pHDat$DT[i] < '2018-7-10 6:00') {
    pHDat$DEP[i] = 2
  }else if(pHDat$DT[i] > '2018-07-10 22:00:00' & pHDat$DT[i] < '2018-9-26 6:00') {
    pHDat$DEP[i] = 3
  }else if(pHDat$DT[i] > '2018-09-27 22:00:00' & pHDat$DT[i] < '2018-12-20 6:00') {
    pHDat$DEP[i] = 4
  }else if(pHDat$DT[i] > '2019-04-01 22:00:00' & pHDat$DT[i] < '2019-11-07 6:00') {
    pHDat$DEP[i] = 5
  }else if(pHDat$DT[i] > '2019-11-07 22:00:00' & pHDat$DT[i] < '2020-5-5 6:00') {
    pHDat$DEP[i] = 6
  }else if(pHDat$DT[i] > '2020-07-16 22:00:00' & pHDat$DT[i] < '2020-10-07 6:00') {
    pHDat$DEP[i] = 7
  }else if(pHDat$DT[i] > '2020-10-07 7:30:00' & pHDat$DT[i] < '2021-07-27 8:00') {
    pHDat$DEP[i] = 8
  }else if(pHDat$DT[i] > '2021-07-27 22:00:00' & pHDat$DT[i] < '2021-11-18 6:00') {
    pHDat$DEP[i] = 9
  }else{pHDat$DEP[i] = NA
  }
}

#Average every hour
pH = pHDat %>%
  filter(QF_pH == "1") 

pH <- pH %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH))
pH <- pH[complete.cases(pH$avg),]
pH <- merge(x = pH, y = pHDat, by = "DT", all.x = TRUE)
pH <- pH[c(1,2,20)]

for (i in 1:nrow(pH)) {
  if(pH$DT[i] < '2018-05-03 6:00:00' & pH$DT[i] > '2017-11-02 22:00') {
    pH$DEP[i] = 1
  }else if(pH$DT[i] > '2018-05-03 22:00:00' & pH$DT[i] < '2018-7-10 6:00') {
    pH$DEP[i] = 2
  }else if(pH$DT[i] > '2018-07-10 22:00:00' & pH$DT[i] < '2018-9-26 6:00') {
    pH$DEP[i] = 3
  }else if(pH$DT[i] > '2018-09-27 22:00:00' & pH$DT[i] < '2018-12-20 6:00') {
    pH$DEP[i] = 4
  }else if(pH$DT[i] > '2019-04-01 22:00:00' & pH$DT[i] < '2019-11-07 6:00') {
    pH$DEP[i] = 5
  }else if(pH$DT[i] > '2019-11-07 22:00:00' & pH$DT[i] < '2020-5-5 6:00') {
    pH$DEP[i] = 6
  }else if(pH$DT[i] > '2020-07-16 22:00:00' & pH$DT[i] < '2020-10-07 6:00') {
    pH$DEP[i] = 7
  }else if(pH$DT[i] > '2020-10-07 7:30:00' & pH$DT[i] < '2021-07-27 8:00') {
    pH$DEP[i] = 8
  }else if(pH$DT[i] > '2021-07-27 22:00:00' & pH$DT[i] < '2021-11-18 6:00') {
    pH$DEP[i] = 9
  }else{pH$DEP[i] = NA
  }
}
pH <- pH[complete.cases(pH$DEP),]


DO <- pHDat[complete.cases(pHDat$DO.mgL),]
DO = DO %>%
  filter(QF_DO == "1")
DO <- DO %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(DO.mgL))
DO <- DO[complete.cases(DO$avg),]
DO <- merge(x = DO, y = pHDat, by = "DT", all.x = TRUE)
DO <- DO[c(1,2,20)]

DO <- DO[complete.cases(DO$DEP),]

for (i in 1:nrow(pHDat)) {
  if(pHDat$QF_TC[i]==1) {
    pHDat$TC_Both[i] = pHDat$pH.Temp[i]
  }else{
    pHDat$TC_Both[i] = pHDat$DO.Temp[i]
  }
}

TC <- pHDat[complete.cases(pHDat$TC_Both),]

TC <- TC %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(TC_Both))
TC <- merge(x = TC, y = pHDat, by = "DT", all.x = TRUE)
TC <- TC[c(1,2,20)]

for (i in 1:nrow(TC)) {
  if(TC$DT[i] < '2018-05-03 6:00:00') {
    TC$DEP[i] = 1
  }else{
    TC$DEP[i] = TC$DEP[i]
  }
}
TC <- TC[complete.cases(TC$DEP),]

AllPA <- merge(x = pH, y = DO, by = "DT", all = TRUE)
AllPA <- merge(x = AllPA, y = TC, by = "DT", all = TRUE)
names(AllPA) <- c('DT', 'pH','Dep1','DO','Dep2','TC','Dep3')

AllPA_Good <- AllPA
AllPA <- AllPA[complete.cases(AllPA$pH),]
AllPA <- AllPA[complete.cases(AllPA$DO),]
AllPA <- AllPA[complete.cases(AllPA$TC),]

#Turn all dates to same year in order to separate by month
AllPA <- AllPA %>%
  mutate(date=ymd_hm(format(AllPA$DT, "2017-%m-%d-%H:%M")))

AllPA_Good <- AllPA_Good %>%
  mutate(date=ymd_hm(format(AllPA_Good$DT, "2017-%m-%d-%H:%M")))

pal = wes_palette("Zissou1",6,type = "continuous")


pHPA <- ggplot(AllPA_Good, aes(x = DT, y = pH, group = Dep1)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)", y = "pH", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
  pHPA 

TCPA <- ggplot(AllPA_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)",element_text(size = 20)) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
TCPA 

DOPA <- ggplot(AllPA_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_line() +
  theme_classic() +
  labs(y="Dissolved Oxygen (mg/L)", x="Date (MM/YY)", element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0))
DOPA 

PA <- ggarrange(pHPA, TCPA, DOPA, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3, 
                legend = "bottom", common.legend = FALSE, align = "hv")
PA

#Regressions of time series 
pHDO <- lm(AllPA$DO ~ AllPA$pH)
summary(pHDO)
pHTC <- lm(AllPA$TC ~ AllPA$pH)
summary(pHTC)
DOTC <- lm(AllPA$TC ~ AllPA$DO)
summary(DOTC)
TCpH <- lm(AllPA$pH ~AllPA$TC)
summary(TCpH)

pHQ <- quantile(AllPA_Good$DO, na.rm = TRUE)
pHave <- mean(AllPA_Good$DO, na.rm = TRUE)
pHstd <- sd(AllPA_Good$DO, na.rm = TRUE)

PAscat <- ggplot(AllPA, aes(x = pH, y = DO, color = TC)) + 
  geom_point(size = 4, shape = 1, position = "jitter") +
  scale_color_gradientn(colours = c("salmon2", "red4", "gold2", "forestgreen", "turquoise2", "royalblue4", "blueviolet"),
                        values = c(1.0,0.9,0.7,0.5,0.3,0.1,0), name = "Temp (\u00B0C)", limits = c(7, 24)) +
  labs(y=element_blank(), x=element_blank()) +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  xlim(7.3, 8.4) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'))
PAscat <- PAscat + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5, barwidth = 1, barheight = 20))+ scale_y_continuous(limits = c(0.5,11), breaks = c(2,4,6,8,10))
PAscat

PAscatpHtemp <- ggplot(AllPA, aes(x = TC, y = pH, color = DO)) + 
  geom_point(size = 4, shape = 1, position = "jitter") +
  scale_color_gradientn(colours = c("salmon2", "red4", "gold2", "forestgreen", "turquoise2", "royalblue4", "blueviolet"),
                        values = c(1.0,0.9,0.7,0.5,0.3,0.1,0), name = "DO mg/L", limits = c(1, 10)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(y=element_blank(), x=element_blank()) +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  xlim(4, 20) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'))
PAscatpHtemp <- PAscatpHtemp + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5, barwidth = 1, barheight = 20))+ scale_y_continuous(limits = c(7.5,8.2), breaks = c(7.5,7.6,7.7, 7.8,7.9,8.0,8.1,8.2))
PAscatpHtemp
PA1<-PAscatpHtemp + annotate(geom="text", x=16, y=7.6, label="y = 0.069 + 7.05, Adj R2 = 0.58",
         color="black")
PA1

pHPAmodel = read.csv("phdot_pa.csv", header = FALSE)
colnames(pHPAmodel) <- c("pH", "DO", "T")
pHPAmodel$DOconverted <-pHPAmodel$DO* 1.42903 

PAmodelscatpHtemp<- ggplot(pHPAmodel, aes(x = T, y = pH, color = DOconverted)) + 
  geom_point(size = 4, shape = 1, position = "jitter") +
  scale_color_gradientn(colours = c("salmon2", "red4", "gold2", "forestgreen", "turquoise2", "royalblue4", "blueviolet"),
                        values = c(1.0,0.9,0.7,0.5,0.3,0.1,0), name = "DO mg/L", limits = c(1, 10)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(y=element_blank(), x=element_blank()) +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  xlim(4, 20) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'))
PAmodelscatpHtemp <- PAmodelscatpHtemp + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5, barwidth = 1, barheight = 20))+ scale_y_continuous(limits = c(7.5,8.2), breaks = c(7.5,7.6,7.7, 7.8,7.9,8.0,8.1,8.2))
PAmodelscatpHtemp

PA2<-PAmodelscatpHtemp + annotate(geom="text", x=16, y=7.6, label="y = 0.054 + 7.27, Adj R2 = 0.54",
                             color="black")
PA2

#Regressions of time series 
pHDOmodel <- lm(pHPAmodel$DOconverted ~ pHPAmodel$pH)
summary(pHDOmodel)
pHTCmodel <- lm(pHPAmodel$T ~ pHPAmodel$pH)
summary(pHTCmodel)
DOTCmodel <- lm(pHPAmodel$T ~ pHPAmodel$DOconverted)
summary(DOTCmodel)
TCpHmodel <- lm(pHPAmodel$pH ~ pHPAmodel$T)
summary(TCpHmodel)
######################################################################################################
pHDat2 = read.csv(paste('/Volumes/GoogleDrive/My Drive/UCSC/R/VD_all.csv'),
                 skip = 7)
pHDat2$Lat = as.numeric(pHDat2$Lat)
pHDat2$Lon = as.numeric(pHDat2$Lon)
pHDat2$Depth = as.numeric(pHDat2$Depth)
pHDat2$pH = as.numeric(pHDat2$pH)
pHDat2$pH.Temp = as.numeric(pHDat2$pH.Temp)
pHDat2$QF_pH = as.numeric(pHDat2$QF_pH )
pHDat2$QF_TC = as.numeric(pHDat2$QF_TC)
pHDat2$DO..umol.kg. = as.numeric(pHDat2$DO..umol.kg.)
pHDat2$DO.mgL <- pHDat2$DO..umol.kg./31.2512
pHDat2$DO..Sat = as.numeric(pHDat2$DO..Sat)
pHDat2$DO.Temp = as.numeric(pHDat2$DO.Temp)
pHDat2$QF_DO = as.numeric(pHDat2$QF_DO)
pHDat2$year = as.numeric(str_sub(pHDat2$pH_Date_UTC_ddmmyyyy,-4,-1))
pHDat2$month = as.numeric(str_sub(pHDat2$pH_Date_UTC_ddmmyyyy,-6,-5))
pHDat2$day = as.numeric(str_sub(pHDat2$pH_Date_UTC_ddmmyyyy,-8,-7))
pHDat2$DT = as.POSIXct(paste(pHDat2$day, pHDat2$month, pHDat2$year, 
                            pHDat2$pH_Time_UTC_hhmmss), format = "%d%m%Y %H:%M:%S", tz = "UTC")

pHDat2 <- pHDat2[complete.cases(pHDat2$Lat),]

for (i in 1:nrow(pHDat2)) {
  if(pHDat2$DT[i] > '2017-11-03 22:00:00' & pHDat2$DT[i] < '2018-04-26 8:00') {
    pHDat2$DEP[i] = 1
  }else if(pHDat2$DT[i] > '2018-04-26 22:00:00' & pHDat2$DT[i] < '2018-07-08 8:00') {
    pHDat2$DEP[i] = 2
  }else if(pHDat2$DT[i] > '2018-07-08 22:00:00' & pHDat2$DT[i] < '2019-01-25 8:00') {
    pHDat2$DEP[i] = 3
  }else if(pHDat2$DT[i] > '2019-01-30 22:00:00' & pHDat2$DT[i] < '2019-04-04 8:00') {
    pHDat2$DEP[i] = 4
  }else if(pHDat2$DT[i] > '2019-04-04 22:00:00' & pHDat2$DT[i] < '2019-11-07 8:00') {
    pHDat2$DEP[i] = 5
  }else if(pHDat2$DT[i] > '2019-11-07 22:00:00' & pHDat2$DT[i] < '2020-05-05 8:00') {
    pHDat2$DEP[i] = 6
  }else if(pHDat2$DT[i] > '2020-07-16 22:00:00' & pHDat2$DT[i] < '2020-10-19 8:00') {
    pHDat2$DEP[i] = 7
  }else if(pHDat2$DT[i] > '2020-10-19 22:00:00' & pHDat2$DT[i] < '2021-07-27 8:00') {
    pHDat2$DEP[i] = 8
  }else if(pHDat2$DT[i] > '2021-07-27 12:30:00' & pHDat2$DT[i] < '2022-04-07 11:00') {
    pHDat2$DEP[i] = 9
  }else{pHDat2$DEP[i] = NA
  }
}

#Average every day
pH <- subset(pHDat2, pHDat2$pH< 8.5)
pH = pH %>%
  filter(QF_pH == "1") 

pH <- pH %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH))
pH <- pH[complete.cases(pH$avg),]
pH <- merge(x = pH, y = pHDat2, by = "DT", all.x = TRUE)
pH <- pH[c(1,2,20)]

for (i in 1:nrow(pH)) {
  if(pH$DT[i] > '2017-11-03 22:00:00' & pH$DT[i] < '2018-04-26 6:00') {
    pH$DEP[i] = 1
  }else if(pH$DT[i] > '2018-04-26 22:00:00' & pH$DT[i] < '2018-07-08 6:00') {
    pH$DEP[i] = 2
  }else if(pH$DT[i] > '2018-07-08 22:00:00' & pH$DT[i] < '2019-01-25 6:00') {
    pH$DEP[i] = 3
  }else if(pH$DT[i] > '2019-01-30 22:00:00' & pH$DT[i] < '2019-04-04 6:00') {
    pH$DEP[i] = 4
  }else if(pH$DT[i] > '2019-04-04 22:00:00' & pH$DT[i] < '2019-11-07 6:00') {
    pH$DEP[i] = 5
  }else if(pH$DT[i] > '2019-11-07 22:00:00' & pH$DT[i] < '2020-05-05 6:00') {
    pH$DEP[i] = 6
  }else if(pH$DT[i] > '2020-07-16 22:00:00' & pH$DT[i] < '2020-10-19 6:00') {
    pH$DEP[i] = 7
  }else if(pH$DT[i] > '2020-10-19 22:00:00' & pH$DT[i] < '2021-07-27 6:00') {
    pH$DEP[i] = 8
  }else if(pH$DT[i] > '2021-07-27 12:30:00' & pH$DT[i] < '2022-04-07 11:00') {
    pH$DEP[i] = 9
  }else{pH$DEP[i] = NA
  }
}
pH <- pH[complete.cases(pH$DEP),]

DO <- pHDat2[complete.cases(pHDat2$DO.mgL),]
DO = DO %>%
  filter(QF_DO == "1")
DO <- DO %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(DO.mgL))
DO <- DO[complete.cases(DO$avg),]
DO <- merge(x = DO, y = pHDat2, by = "DT", all.x = TRUE)
DO <- DO[c(1,2,20)]

for (i in 1:nrow(DO)) {
  if(DO$DT[i] > '2017-11-03 22:00:00' & DO$DT[i] < '2018-04-26 6:00') {
    DO$DEP[i] = 1
  }else if(DO$DT[i] > '2018-04-26 22:00:00' & DO$DT[i] < '2018-07-08 6:00') {
    DO$DEP[i] = 2
  }else if(DO$DT[i] > '2018-07-08 22:00:00' & DO$DT[i] < '2019-01-25 6:00') {
    DO$DEP[i] = 3
  }else if(DO$DT[i] > '2019-01-30 22:00:00' & DO$DT[i] < '2019-04-04 6:00') {
    DO$DEP[i] = 4
  }else if(DO$DT[i] > '2019-04-04 22:00:00' & DO$DT[i] < '2019-11-07 6:00') {
    DO$DEP[i] = 5
  }else if(DO$DT[i] > '2019-11-07 22:00:00' & DO$DT[i] < '2020-05-05 6:00') {
    DO$DEP[i] = 6
  }else if(DO$DT[i] > '2020-07-16 22:00:00' & DO$DT[i] < '2020-10-19 6:00') {
    DO$DEP[i] = 7
  }else if(DO$DT[i] > '2020-10-19 22:00:00' & DO$DT[i] < '2021-07-27 6:00') {
    DO$DEP[i] = 8
  }else if(DO$DT[i] > '2021-07-27 12:30:00' & DO$DT[i] < '2022-04-07 11:00') {
    DO$DEP[i] = 9
  }else{DO$DEP[i] = NA
  }
}
DO <- DO[complete.cases(DO$DEP),]

TC <- pHDat2[complete.cases(pHDat2$pH.Temp),]

TC = TC %>%
  filter(QF_TC == "1") 

TC <- TC %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH.Temp))
TC <- TC[complete.cases(TC$avg),]
TC <- merge(x = TC, y = pHDat2, by = "DT", all.x = TRUE)
TC <- TC[c(1,2,20)]

for (i in 1:nrow(TC)) {
  if(TC$DT[i] > '2017-11-03 22:00:00' & TC$DT[i] < '2018-04-26 6:00') {
    TC$DEP[i] = 1
  }else if(TC$DT[i] > '2018-04-26 22:00:00' & TC$DT[i] < '2018-07-08 6:00') {
    TC$DEP[i] = 2
  }else if(TC$DT[i] > '2018-07-08 22:00:00' & TC$DT[i] < '2019-01-25 6:00') {
    TC$DEP[i] = 3
  }else if(TC$DT[i] > '2019-01-30 22:00:00' & TC$DT[i] < '2019-04-04 6:00') {
    TC$DEP[i] = 4
  }else if(TC$DT[i] > '2019-04-04 22:00:00' & TC$DT[i] < '2019-11-07 6:00') {
    TC$DEP[i] = 5
  }else if(TC$DT[i] > '2019-11-07 22:00:00' & TC$DT[i] < '2020-05-05 6:00') {
    TC$DEP[i] = 6
  }else if(TC$DT[i] > '2020-07-16 22:00:00' & TC$DT[i] < '2020-10-19 6:00') {
    TC$DEP[i] = 7
  }else if(TC$DT[i] > '2020-10-19 22:00:00' & TC$DT[i] < '2021-07-27 6:00') {
    TC$DEP[i] = 8
  }else if(TC$DT[i] > '2021-07-27 12:30:00' & TC$DT[i] < '2022-04-07 11:00') {
    TC$DEP[i] = 9  
  }else{TC$DEP[i] = NA
  }
}
TC <- TC[complete.cases(TC$DEP),]

AllVD <- merge(x = pH, y = DO, by = "DT", all = TRUE)
AllVD <- merge(x = AllVD, y = TC, by = "DT", all = TRUE)

names(AllVD) <- c('DT', 'pH','Dep1','DO','Dep2','TC','Dep3')
AllVD_Good <- AllVD

AllVD <- AllVD[complete.cases(AllVD$pH),]
AllVD <- AllVD[complete.cases(AllVD$DO),]
AllVD <- AllVD[complete.cases(AllVD$TC),]

#Turn all dates to same year in order to separate by month
AllVD <- AllVD %>%
  mutate(date=ymd_hm(format(AllVD$DT, "2017-%m-%d-%H:%M")))

AllVD_Good <- AllVD_Good %>%
  mutate(date=ymd_hm(format(AllVD_Good$DT, "2017-%m-%d-%H:%M")))

pHVD <- ggplot(AllVD_Good, aes(x = DT, y = pH, group = Dep1)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)", y = "pH", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
pHVD

TCVD <- ggplot(AllVD_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)",element_text(size = 20)) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
TCVD 

DOVD <- ggplot(AllVD_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_line() +
  theme_classic() +
  labs(y="Dissolved Oxygen (mg/L)", x="Date (MM/YY)", element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0))
DOVD 

VD <- ggarrange(pHVD, TCVD, DOVD, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3, 
                legend = "bottom", common.legend = FALSE, align = "hv")
VD

#Regressions of time series 
pHDO <- lm(AllVD$DO ~ AllVD$pH)
summary(pHDO)
pHTC <- lm(AllVD$TC ~ AllVD$pH)
summary(pHTC)
DOTC <- lm(AllVD$TC ~ AllVD$DO)
summary(DOTC)

pHQ <- quantile(AllVD_Good$DO, na.rm = TRUE)
pHave <- mean(AllVD_Good$DO, na.rm = TRUE)
pHstd <- sd(AllVD_Good$DO, na.rm = TRUE)

######################################################################################################
pHDat3 = read.csv(paste('/Volumes/GoogleDrive/My Drive/UCSC/R/LB_all.csv'),
                 skip = 7)
pHDat3$Lat = as.numeric(pHDat3$Lat)
pHDat3$Lon = as.numeric(pHDat3$Lon)
pHDat3$Depth = as.numeric(pHDat3$Depth)
pHDat3$pH = as.numeric(pHDat3$pH)
pHDat3$pH.Temp = as.numeric(pHDat3$pH.Temp)
pHDat3$QF_pH = as.numeric(pHDat3$QF_pH )
pHDat3$QF_TC = as.numeric(pHDat3$QF_TC)
pHDat3$DO..umol.kg. = as.numeric(pHDat3$DO..umol.kg.)
pHDat3$DO.mgL <- pHDat3$DO..umol.kg./31.2512
pHDat3$DO..Sat = as.numeric(pHDat3$DO..Sat)
pHDat3$DO.Temp = as.numeric(pHDat3$DO.Temp)
pHDat3$QF_DO = as.numeric(pHDat3$QF_DO)
pHDat3$year = as.numeric(str_sub(pHDat3$pH_Date_UTC_ddmmyyyy,-4,-1))
pHDat3$month = as.numeric(str_sub(pHDat3$pH_Date_UTC_ddmmyyyy,-6,-5))
pHDat3$day = as.numeric(str_sub(pHDat3$pH_Date_UTC_ddmmyyyy,-8,-7))
pHDat3$DT = as.POSIXct(paste(pHDat3$day, pHDat3$month, pHDat3$year, 
                             pHDat3$pH_Time_UTC_hhmmss), format = "%d%m%Y %H:%M:%S", tz = "UTC")

pHDat3 <- pHDat3[complete.cases(pHDat3$Lat),]

for (i in 1:nrow(pHDat3)) {
  if(pHDat3$DT[i] > '2017-12-12 24:00:00' & pHDat3$DT[i] < '2018-01-09 6:00') {
    pHDat3$DEP[i] = 1
  }else if(pHDat3$DT[i] > '2018-06-12 24:00:00' & pHDat3$DT[i] < '2018-09-05 6:00') {
    pHDat3$DEP[i] = 2
  }else if(pHDat3$DT[i] > '2018-09-05 24:00:00' & pHDat3$DT[i] < '2018-12-10 6:00') {
    pHDat3$DEP[i] = 3
  }else if(pHDat3$DT[i] > '2018-12-10 24:00:00' & pHDat3$DT[i] < '2019-03-14 6:00') {
    pHDat3$DEP[i] = 4
  }else if(pHDat3$DT[i] > '2019-03-14 24:00:00' & pHDat3$DT[i] < '2019-06-10 6:00') {
    pHDat3$DEP[i] = 5
  }else if(pHDat3$DT[i] > '2019-06-10 24:00:00' & pHDat3$DT[i] < '2019-12-18 6:00') {
    pHDat3$DEP[i] = 6
  }else if(pHDat3$DT[i] > '2019-12-18 24:00:00' & pHDat3$DT[i] < '2020-03-04 6:00') {
    pHDat3$DEP[i] = 7
  }else if(pHDat3$DT[i] > '2020-03-04 24:00:00' & pHDat3$DT[i] < '2020-9-27 6:00') {
    pHDat3$DEP[i] = 8
  }else if(pHDat3$DT[i] > '2020-11-19 24:00:00' & pHDat3$DT[i] < '2021-04-02 6:00') {
    pHDat3$DEP[i] = 9
  }else if(pHDat3$DT[i] > '2021-04-02 24:00:00' & pHDat3$DT[i] < '2021-10-29 6:00') {
    pHDat3$DEP[i] = 10
  }else if(pHDat3$DT[i] > '2021-10-29 10:30:00' & pHDat3$DT[i] < '2022-06-26 6:00') {
    pHDat3$DEP[i] = 11
  }else{pHDat3$DEP[i] = NA
  }
}


#Average every hour
pH <- subset(pHDat3, pHDat3$pH< 8.5)
pH = pH %>%
  filter(QF_pH == "1") 

pH <- pH %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH))
pH <- pH[complete.cases(pH$avg),]
pH <- merge(x = pH, y = pHDat3, by = "DT", all.x = TRUE)
pH <- pH[c(1,2,20)]

for (i in 1:nrow(pH)) {
  if(pH$DT[i] > '2017-12-12 24:00:00' & pH$DT[i] < '2018-01-09 6:00') {
    pH$DEP[i] = 1
  }else if(pH$DT[i] > '2018-06-12 24:00:00' & pH$DT[i] < '2018-09-05 6:00') {
    pH$DEP[i] = 2
  }else if(pH$DT[i] > '2018-09-05 24:00:00' & pH$DT[i] < '2018-12-10 6:00') {
    pH$DEP[i] = 3
  }else if(pH$DT[i] > '2018-12-10 24:00:00' & pH$DT[i] < '2019-03-14 6:00') {
    pH$DEP[i] = 4
  }else if(pH$DT[i] > '2019-03-14 24:00:00' & pH$DT[i] < '2019-06-10 6:00') {
    pH$DEP[i] = 5
  }else if(pH$DT[i] > '2019-06-10 24:00:00' & pH$DT[i] < '2019-12-18 6:00') {
    pH$DEP[i] = 6
  }else if(pH$DT[i] > '2019-12-18 24:00:00' & pH$DT[i] < '2020-03-04 6:00') {
    pH$DEP[i] = 7
  }else if(pH$DT[i] > '2020-03-04 24:00:00' & pH$DT[i] < '2020-4-27 6:00') {
    pH$DEP[i] = 8
  }else if(pH$DT[i] > '2020-11-19 24:00:00' & pH$DT[i] < '2021-03-02 6:00') {
    pH$DEP[i] = 9
  }else if(pH$DT[i] > '2021-04-02 24:00:00' & pH$DT[i] < '2021-10-29 6:00') {
    pH$DEP[i] = 10
  }else if(pH$DT[i] > '2021-10-29 10:30:00' & pH$DT[i] < '2022-06-26 6:00') {
    pH$DEP[i] = 11
  }else{pH$DEP[i] = NA
  }
}
pH <- pH[complete.cases(pH$DEP),]

DO <- pHDat3[complete.cases(pHDat3$DO.mgL),]
DO = DO %>%
  filter(QF_DO == "1")
DO <- DO %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(DO.mgL))
DO <- DO[complete.cases(DO$avg),]
DO <- merge(x = DO, y = pHDat3, by = "DT", all.x = TRUE)
DO <- DO[c(1,2,20)]

for (i in 1:nrow(DO)) {
  if(DO$DT[i] > '2017-12-12 24:00:00' & DO$DT[i] < '2018-01-09 6:00') {
    DO$DEP[i] = 1
  }else if(DO$DT[i] > '2018-06-12 24:00:00' & DO$DT[i] < '2018-09-05 6:00') {
    DO$DEP[i] = 2
  }else if(DO$DT[i] > '2018-09-05 24:00:00' & DO$DT[i] < '2018-12-10 6:00') {
    DO$DEP[i] = 3
  }else if(DO$DT[i] > '2018-12-10 24:00:00' & DO$DT[i] < '2019-03-01 6:00') {
    DO$DEP[i] = 4
  }else if(DO$DT[i] > '2019-03-14 24:00:00' & DO$DT[i] < '2019-06-10 6:00') {
    DO$DEP[i] = 5
  }else if(DO$DT[i] > '2019-06-10 24:00:00' & DO$DT[i] < '2019-12-01 6:00') {
    DO$DEP[i] = 6
  }else if(DO$DT[i] > '2019-12-18 24:00:00' & DO$DT[i] < '2020-03-04 6:00') {
    DO$DEP[i] = 7
  }else if(DO$DT[i] > '2020-03-04 24:00:00' & DO$DT[i] < '2020-9-27 6:00') {
    DO$DEP[i] = 8
  }else if(DO$DT[i] > '2020-11-19 24:00:00' & DO$DT[i] < '2021-03-02 6:00') {
    DO$DEP[i] = 9
  }else if(DO$DT[i] > '2021-04-02 24:00:00' & DO$DT[i] < '2021-10-29 6:00') {
    DO$DEP[i] = 10
  }else if(DO$DT[i] > '2021-10-29 10:30:00' & DO$DT[i] < '2022-06-26 6:00') {
    DO$DEP[i] = 11
  }else{DO$DEP[i] = NA
  }
}
DO <- DO[complete.cases(DO$DEP),]

TC <- pHDat3[complete.cases(pHDat3$pH.Temp),]

TC = TC %>%
  filter(QF_TC == "1") 

TC <- TC %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH.Temp))
TC <- TC[complete.cases(TC$avg),]
TC <- merge(x = TC, y = pHDat3, by = "DT", all.x = TRUE)
TC <- TC[c(1,2,20)]

for (i in 1:nrow(TC)) {
  if(TC$DT[i] > '2017-12-12 24:00:00' & TC$DT[i] < '2018-01-09 6:00') {
    TC$DEP[i] = 1
  }else if(TC$DT[i] > '2018-06-12 24:00:00' & TC$DT[i] < '2018-09-05 6:00') {
    TC$DEP[i] = 2
  }else if(TC$DT[i] > '2018-09-05 24:00:00' & TC$DT[i] < '2018-12-10 6:00') {
    TC$DEP[i] = 3
  }else if(TC$DT[i] > '2018-12-10 24:00:00' & TC$DT[i] < '2019-03-01 6:00') {
    TC$DEP[i] = 4
  }else if(TC$DT[i] > '2019-03-14 24:00:00' & TC$DT[i] < '2019-06-10 6:00') {
    TC$DEP[i] = 5
  }else if(TC$DT[i] > '2019-06-10 24:00:00' & TC$DT[i] < '2019-12-01 6:00') {
    TC$DEP[i] = 6
  }else if(TC$DT[i] > '2019-12-18 24:00:00' & TC$DT[i] < '2020-03-04 6:00') {
    TC$DEP[i] = 7
  }else if(TC$DT[i] > '2020-03-04 24:00:00' & TC$DT[i] < '2020-9-27 6:00') {
    TC$DEP[i] = 8
  }else if(TC$DT[i] > '2020-11-19 24:00:00' & TC$DT[i] < '2021-03-02 6:00') {
    TC$DEP[i] = 9
  }else if(TC$DT[i] > '2021-04-02 24:00:00' & TC$DT[i] < '2021-10-29 6:00') {
    TC$DEP[i] = 10
  }else if(TC$DT[i] > '2021-10-29 10:30:00' & TC$DT[i] < '2022-06-26 6:00') {
    TC$DEP[i] = 11  
  }else{TC$DEP[i] = NA
  }
}
TC <- TC[complete.cases(TC$DEP),]

AllLB <- merge(x = pH, y = DO, by = "DT", all = TRUE)
AllLB <- merge(x = AllLB, y = TC, by = "DT", all = TRUE)

names(AllLB) <- c('DT', 'pH','Dep1','DO','Dep2','TC','Dep3')
AllLB_Good <- AllLB

AllLB <- AllLB[complete.cases(AllLB$pH),]
AllLB <- AllLB[complete.cases(AllLB$DO),]
AllLB <- AllLB[complete.cases(AllLB$TC),]

#Turn all dates to same year in order to separate by month
AllLB <- AllLB %>%
  mutate(date=ymd_hm(format(AllLB$DT, "2017-%m-%d-%H:%M")))

AllLB_Good <- AllLB_Good %>%
  mutate(date=ymd_hm(format(AllLB_Good$DT, "2017-%m-%d-%H:%M")))

pHLB <- ggplot(AllLB_Good, aes(x = DT, y = pH, group = Dep1)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)", y = "pH", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
pHLB

TCLB <- ggplot(AllLB_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)",element_text(size = 20)) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
TCLB 

DOLB <- ggplot(AllLB_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_line() +
  theme_classic() +
  labs(y="Dissolved Oxygen (mg/L)", x="Date (MM/YY)", element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0))
DOLB 

LB <- ggarrange(pHLB, TCLB, DOLB, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3, 
                legend = "bottom", common.legend = FALSE, align = "hv")
LB

#Regressions of time series 
pHDO <- lm(AllLB$DO ~ AllLB$pH)
summary(pHDO)
pHTC <- lm(AllLB$TC ~ AllLB$pH)
summary(pHTC)
DOTC <- lm(AllLB$TC ~ AllLB$DO)
summary(DOTC)

pHQ <- quantile(AllLB_Good$TC, na.rm = TRUE)
pHave <- mean(AllLB_Good$TC, na.rm = TRUE)
pHstd <- sd(AllLB_Good$TC, na.rm = TRUE)

######################################################################################################
pHDat4 = read.csv(paste('/Volumes/GoogleDrive/My Drive/UCSC/R/CI_all.csv'),
                  skip = 7)
pHDat4 = read.csv(paste("CI_all.csv"),
                  skip = 7)
pHDat4$Lat = as.numeric(pHDat4$Lat)
pHDat4$Lon = as.numeric(pHDat4$Lon)
pHDat4$Depth = as.numeric(pHDat4$Depth)
pHDat4$pH = as.numeric(pHDat4$pH)
pHDat4$pH.Temp = as.numeric(pHDat4$pH.Temp)
pHDat4$QF_pH = as.numeric(pHDat4$QF_pH )
pHDat4$QF_TC = as.numeric(pHDat4$QF_TC)
pHDat4$DO..umol.kg. = as.numeric(pHDat4$DO..umol.kg.)
pHDat4$DO.mgL <- pHDat4$DO..umol.kg./31.2512
pHDat4$DO..Sat = as.numeric(pHDat4$DO..Sat)
pHDat4$DO.Temp = as.numeric(pHDat4$DO.Temp)
pHDat4$QF_DO = as.numeric(pHDat4$QF_DO)
pHDat4$year = as.numeric(str_sub(pHDat4$pH_Date_UTC_ddmmyyyy,-4,-1))
pHDat4$month = as.numeric(str_sub(pHDat4$pH_Date_UTC_ddmmyyyy,-6,-5))
pHDat4$day = as.numeric(str_sub(pHDat4$pH_Date_UTC_ddmmyyyy,-8,-7))
pHDat4$DT = as.POSIXct(paste(pHDat4$day, pHDat4$month, pHDat4$year, 
                             pHDat4$pH_Time_UTC_hhmmss), format = "%d%m%Y %H:%M:%S", tz = "UTC")

pHDat4 <- pHDat4[complete.cases(pHDat4$Lat),]


for (i in 1:nrow(pHDat4)) {
  if(pHDat4$DT[i] > '2017-12-11 22:00:00' & pHDat4$DT[i] < '2018-03-01 6:00') {
    pHDat4$DEP[i] = 1
  }else if(pHDat4$DT[i] > '2018-03-09 22:00:00' & pHDat4$DT[i] < '2018-06-11 6:00') {
    pHDat4$DEP[i] = 2
  }else if(pHDat4$DT[i] > '2018-06-11 22:00:00' & pHDat4$DT[i] < '2018-09-05 6:00') {
    pHDat4$DEP[i] = 3
  }else if(pHDat4$DT[i] > '2018-12-10 22:00:00' & pHDat4$DT[i] < '2019-03-14 6:00') {
    pHDat4$DEP[i] = 4
  }else if(pHDat4$DT[i] > '2019-03-14 22:00:00' & pHDat4$DT[i] < '2019-06-10 6:00') {
    pHDat4$DEP[i] = 5
  }else if(pHDat4$DT[i] > '2019-06-10 22:00:00' & pHDat4$DT[i] < '2019-12-18 6:00') {
    pHDat4$DEP[i] = 6
  }else if(pHDat4$DT[i] > '2019-12-18 22:00:00' & pHDat4$DT[i] < '2020-03-03 6:00') {
    pHDat4$DEP[i] = 7
  }else if(pHDat4$DT[i] > '2020-03-03 22:00:00' & pHDat4$DT[i] < '2020-10-14 6:00') {
    pHDat4$DEP[i] = 8
  }else if(pHDat4$DT[i] > '2020-10-14 22:00:00' & pHDat4$DT[i] < '2021-06-26 6:00') {
    pHDat4$DEP[i] = 9
  }else{pHDat4$DEP[i] = NA
  }
}

#Average every hour
pH <- subset(pHDat4, pHDat4$pH< 8.5)
pH = pH %>%
  filter(QF_pH == "1") 

pH <- pH %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH))
pH <- pH[complete.cases(pH$avg),]
pH <- merge(x = pH, y = pHDat4, by = "DT", all.x = TRUE)
pH <- pH[c(1,2,20)]

for (i in 1:nrow(pH)) {
  if(pH$DT[i] > '2017-12-11 22:00:00' & pH$DT[i] < '2018-03-01 6:00') {
    pH$DEP[i] = 1
  }else if(pH$DT[i] > '2018-03-09 22:00:00' & pH$DT[i] < '2018-06-11 6:00') {
    pH$DEP[i] = 2
  }else if(pH$DT[i] > '2018-06-11 22:00:00' & pH$DT[i] < '2018-09-05 6:00') {
    pH$DEP[i] = 3
  }else if(pH$DT[i] > '2018-12-10 22:00:00' & pH$DT[i] < '2019-03-14 6:00') {
    pH$DEP[i] = 4
  }else if(pH$DT[i] > '2019-03-14 22:00:00' & pH$DT[i] < '2019-06-10 6:00') {
    pH$DEP[i] = 5
  }else if(pH$DT[i] > '2019-06-10 22:00:00' & pH$DT[i] < '2019-12-18 6:00') {
    pH$DEP[i] = 6
  }else if(pH$DT[i] > '2019-12-18 22:00:00' & pH$DT[i] < '2020-03-03 6:00') {
    pH$DEP[i] = 7
  }else if(pH$DT[i] > '2020-03-03 22:00:00' & pH$DT[i] < '2020-10-14 6:00') {
    pH$DEP[i] = 8
  }else if(pH$DT[i] > '2020-10-14 22:00:00' & pH$DT[i] < '2021-06-26 6:00') {
    pH$DEP[i] = 9
  }else{pH$DEP[i] = NA
  }
}
pH <- pH[complete.cases(pH$DEP),]

DO <- pHDat4[complete.cases(pHDat4$DO.mgL),]
DO = DO %>%
  filter(QF_DO == "1")
DO <- DO %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(DO.mgL))
DO <- DO[complete.cases(DO$avg),]
DO <- merge(x = DO, y = pHDat4, by = "DT", all.x = TRUE)
DO <- DO[c(1,2,20)]

for (i in 1:nrow(DO)) {
  if(DO$DT[i] > '2017-12-11 22:00:00' & DO$DT[i] < '2018-03-01 6:00') {
    DO$DEP[i] = 1
  }else if(DO$DT[i] > '2018-03-09 22:00:00' & DO$DT[i] < '2018-06-11 6:00') {
    DO$DEP[i] = 2
  }else if(DO$DT[i] > '2018-06-11 22:00:00' & DO$DT[i] < '2018-09-05 6:00') {
    DO$DEP[i] = 3
  }else if(DO$DT[i] > '2018-12-10 22:00:00' & DO$DT[i] < '2019-03-14 6:00') {
    DO$DEP[i] = 4
  }else if(DO$DT[i] > '2019-03-14 22:00:00' & DO$DT[i] < '2019-06-10 6:00') {
    DO$DEP[i] = 5
  }else if(DO$DT[i] > '2019-06-10 22:00:00' & DO$DT[i] < '2019-12-18 6:00') {
    DO$DEP[i] = 6
  }else if(DO$DT[i] > '2019-12-18 22:00:00' & DO$DT[i] < '2020-03-03 6:00') {
    DO$DEP[i] = 7
  }else if(DO$DT[i] > '2020-03-03 22:00:00' & DO$DT[i] < '2020-10-14 6:00') {
    DO$DEP[i] = 8
  }else if(DO$DT[i] > '2020-10-14 22:00:00' & DO$DT[i] < '2021-06-26 6:00') {
    DO$DEP[i] = 9
  }else{DO$DEP[i] = NA
  }
}
DO <- DO[complete.cases(DO$DEP),]


for (i in 1:nrow(pHDat4)) {
  if(pHDat4$QF_TC[i]==1) {
    pHDat4$TC_Both[i] = pHDat4$pH.Temp[i]
  }else{
    pHDat4$TC_Both[i] = pHDat4$DO.Temp[i]
  }
}

TC <- pHDat4[complete.cases(pHDat4$TC_Both),]

TC <- TC %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(TC_Both))
TC <- merge(x = TC, y = pHDat4, by = "DT", all.x = TRUE)
TC <- TC[c(1,2,20)]

for (i in 1:nrow(TC)) {
  if(TC$DT[i] > '2017-12-11 22:00:00' & TC$DT[i] < '2018-03-01 6:00') {
    TC$DEP[i] = 1
  }else if(TC$DT[i] > '2018-03-09 22:00:00' & TC$DT[i] < '2018-06-11 6:00') {
    TC$DEP[i] = 2
  }else if(TC$DT[i] > '2018-06-11 22:00:00' & TC$DT[i] < '2018-09-05 6:00') {
    TC$DEP[i] = 3
  }else if(TC$DT[i] > '2018-12-10 22:00:00' & TC$DT[i] < '2019-03-14 6:00') {
    TC$DEP[i] = 4
  }else if(TC$DT[i] > '2019-03-14 22:00:00' & TC$DT[i] < '2019-06-10 6:00') {
    TC$DEP[i] = 5
  }else if(TC$DT[i] > '2019-06-10 22:00:00' & TC$DT[i] < '2019-12-18 6:00') {
    TC$DEP[i] = 6
  }else if(TC$DT[i] > '2019-12-18 22:00:00' & TC$DT[i] < '2020-03-03 6:00') {
    TC$DEP[i] = 7
  }else if(TC$DT[i] > '2020-03-03 22:00:00' & TC$DT[i] < '2020-10-14 6:00') {
    TC$DEP[i] = 8
  }else if(TC$DT[i] > '2020-10-14 22:00:00' & TC$DT[i] < '2021-06-26 6:00') {
    TC$DEP[i] = 9
  }else{TC$DEP[i] = NA
  }
}
TC <- TC[complete.cases(TC$DEP),]

AllCI <- merge(x = pH, y = DO, by = "DT", all = TRUE)
AllCI <- merge(x = AllCI, y = TC, by = "DT", all = TRUE)

names(AllCI) <- c('DT', 'pH','Dep1','DO','Dep2','TC','Dep3')
AllCI_Good <- AllCI

AllCI <- AllCI[complete.cases(AllCI$pH),]
AllCI <- AllCI[complete.cases(AllCI$DO),]
AllCI <- AllCI[complete.cases(AllCI$TC),]

#Turn all dates to same year in order to separate by month
AllCI <- AllCI %>%
  mutate(date=ymd_hm(format(AllCI$DT, "2017-%m-%d-%H:%M")))

AllCI_Good <- AllCI_Good %>%
  mutate(date=ymd_hm(format(AllCI_Good$DT, "2017-%m-%d-%H:%M")))

pHCI <- ggplot(AllCI_Good, aes(x = DT, y = pH, group = Dep1)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)", y = "pH", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
pHCI

TCCI <- ggplot(AllCI_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)",element_text(size = 20)) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
TCCI 

DOCI <- ggplot(AllCI_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_line() +
  theme_classic() +
  labs(y="Dissolved Oxygen (mg/L)", x="Date (MM/YY)", element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0))
DOCI 

CI <- ggarrange(pHCI, TCCI, DOCI, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3, 
                legend = "bottom", common.legend = FALSE, align = "hv")
CI

#Regressions of time series 
pHDO <- lm(AllCI$DO ~ AllCI$pH)
summary(pHDO)
pHTC <- lm(AllCI$TC ~ AllCI$pH)
summary(pHTC)
DOTC <- lm(AllCI$TC ~ AllCI$DO)
summary(DOTC)
TCpH <- lm(AllCI$pH ~ AllCI$TC)
summary(TCpH)

pHQ <- quantile(AllCI_Good$pH, na.rm = TRUE)
pHave <- mean(AllCI_Good$pH, na.rm = TRUE)
pHstd <- sd(AllCI_Good$pH, na.rm = TRUE)

CIscat <- ggplot(AllCI, aes(x = pH, y = DO, color = TC)) + 
  geom_point(size = 4, shape = 1, position = "jitter") +
  scale_color_gradientn(colours = c("salmon2", "red4", "gold2", "forestgreen", "turquoise2", "royalblue4", "blueviolet"),
                        values = c(1.0,0.9,0.7,0.5,0.3,0.1,0), name = "Temp (\u00B0C)", limits = c(7, 24)) +
  labs(y=element_blank(), x="pH") +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  xlim(7.3, 8.4) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'), axis.title=element_text(size=16))
CIscat <- CIscat + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5, barwidth = 1, barheight = 20))+ scale_y_continuous(limits = c(0.5,11), breaks = c(2,4,6,8,10))
CIscat

CIscatpHtemp <- ggplot(AllCI, aes(x = TC, y = pH, color = DO)) + 
  geom_point(size = 4, shape = 1, position = "jitter") +
  scale_color_gradientn(colours = c("salmon2", "red4", "gold2", "forestgreen", "turquoise2", "royalblue4", "blueviolet"),
                        values = c(1.0,0.9,0.7,0.5,0.3,0.1,0), name = "DO mg/L", limits = c(1, 10)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(y=element_blank(), x=element_blank()) +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  xlim(4, 20) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'))
CIscatpHtemp <- CIscatpHtemp + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5, barwidth = 1, barheight = 20))+ scale_y_continuous(limits = c(7.5,8.2), breaks = c(7.5,7.6,7.7, 7.8,7.9,8.0,8.1,8.2))
CIscatpHtemp
CI1<-CIscatpHtemp + annotate(geom="text", x=16, y=7.6, label="y = -0.004 + 8.09, Adj R2 = 0.04",
                             color="black")
CI1

pHCImodel = read.csv("phdot_ci.csv", header = FALSE)
colnames(pHCImodel) <- c("pH", "DO", "T")
pHCImodel$DOconverted <-pHCImodel$DO* 1.42903 

CImodelscatpHtemp<- ggplot(pHCImodel, aes(x = T, y = pH, color = DOconverted)) + 
  geom_point(size = 4, shape = 1, position = "jitter") +
  scale_color_gradientn(colours = c("salmon2", "red4", "gold2", "forestgreen", "turquoise2", "royalblue4", "blueviolet"),
                        values = c(1.0,0.9,0.7,0.5,0.3,0.1,0), name = "DO mg/L", limits = c(1, 10)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(y=element_blank(), x=element_blank()) +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  xlim(4, 20) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'))
CImodelscatpHtemp <- CImodelscatpHtemp + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5, barwidth = 1, barheight = 20))+ scale_y_continuous(limits = c(7.5,8.2), breaks = c(7.5,7.6,7.7, 7.8,7.9,8.0,8.1,8.2))
CImodelscatpHtemp

TCpHmodel <- lm(pHCImodel$pH ~ pHCImodel$T)
summary(TCpHmodel)

CI2<-CImodelscatpHtemp + annotate(geom="text", x=16, y=7.6, label="y = -0.003 + 8.10, Adj R2 = 0.07",
                             color="black")
CI2


######################################################################################################
pHDat5 = read.csv(paste('/Volumes/GoogleDrive/My Drive/UCSC/R/PB_all.csv'),
                  skip = 7)
pHDat5$Lat = as.numeric(pHDat5$Lat)
pHDat5$Lon = as.numeric(pHDat5$Lon)
pHDat5$Depth = as.numeric(pHDat5$Depth)
pHDat5$pH = as.numeric(pHDat5$pH)
pHDat5$pH.Temp = as.numeric(pHDat5$pH.Temp)
pHDat5$QF_pH = as.numeric(pHDat5$QF_pH )
pHDat5$QF_TC = as.numeric(pHDat5$QF_TC)
pHDat5$DO..umol.kg. = as.numeric(pHDat5$DO..umol.kg.)
pHDat5$DO.mgL <- pHDat5$DO..umol.kg./31.2512
pHDat5$DO..Sat = as.numeric(pHDat5$DO..Sat)
pHDat5$DO.Temp = as.numeric(pHDat5$DO.Temp)
pHDat5$QF_DO = as.numeric(pHDat5$QF_DO)
pHDat5$year = as.numeric(str_sub(pHDat5$pH_Date_UTC_ddmmyyyy,-4,-1))
pHDat5$month = as.numeric(str_sub(pHDat5$pH_Date_UTC_ddmmyyyy,-6,-5))
pHDat5$day = as.numeric(str_sub(pHDat5$pH_Date_UTC_ddmmyyyy,-8,-7))
pHDat5$DT = as.POSIXct(paste(pHDat5$day, pHDat5$month, pHDat5$year, 
                             pHDat5$pH_Time_UTC_hhmmss), format = "%d%m%Y %H:%M:%S", tz = "UTC")

pHDat5 <- pHDat5[complete.cases(pHDat5$Lat),]

for (i in 1:nrow(pHDat5)) {
  if(pHDat5$DT[i] > '2017-12-07 12:00:00' & pHDat5$DT[i] < '2018-02-15 6:00') {
    pHDat5$DEP[i] = 1
  }else if(pHDat5$DT[i] > '2018-02-15 11:00:00' & pHDat5$DT[i] < '2018-07-21 9:00') {
    pHDat5$DEP[i] = 2
  }else if(pHDat5$DT[i] > '2018-07-21 11:00:00' & pHDat5$DT[i] < '2018-09-22 11:00') {
    pHDat5$DEP[i] = 3
  }else if(pHDat5$DT[i] > '2018-09-22 11:00:00' & pHDat5$DT[i] < '2019-03-04 10:00') {
    pHDat5$DEP[i] = 4
  }else if(pHDat5$DT[i] > '2019-03-04 11:00:00' & pHDat5$DT[i] < '2019-07-20 9:00') {
    pHDat5$DEP[i] = 5
  }else if(pHDat5$DT[i] > '2019-07-20 11:00:00' & pHDat5$DT[i] < '2019-11-22 9:00') {
    pHDat5$DEP[i] = 6
  }else if(pHDat5$DT[i] > '2019-11-22 11:00:00' & pHDat5$DT[i] < '2020-02-20 9:00') {
    pHDat5$DEP[i] = 7
  }else if(pHDat5$DT[i] > '2020-02-20 10:00:00' & pHDat5$DT[i] < '2020-08-12 10:00') {
    pHDat5$DEP[i] = 8
  }else if(pHDat5$DT[i] > '2020-08-12 10:00:00' & pHDat5$DT[i] < '2020-11-12 9:00') {
    pHDat5$DEP[i] = 9
  }else if(pHDat5$DT[i] > '2020-11-12 10:00:00' & pHDat5$DT[i] < '2021-06-01 7:00') {
    pHDat5$DEP[i] = 10
  }else if(pHDat5$DT[i] > '2021-06-02 13:00:00' & pHDat5$DT[i] < '2021-09-12 7:00') {
    pHDat5$DEP[i] = 11
  }else{pHDat5$DEP[i] = NA
  }
}

#Average every hour
pH <- subset(pHDat5, pHDat5$pH< 8.5)
pH = pH %>%
  filter(QF_pH == "1") 

pH <- pH %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH))
pH <- pH[complete.cases(pH$avg),]
pH <- merge(x = pH, y = pHDat5, by = "DT", all.x = TRUE)
pH <- pH[c(1,2,20)]

for (i in 1:nrow(pH)) {
  if(pH$DT[i] > '2017-12-07 12:00:00' & pH$DT[i] < '2018-02-15 6:00') {
    pH$DEP[i] = 1
  }else if(pH$DT[i] > '2018-02-15 11:00:00' & pH$DT[i] < '2018-05-19 9:00') {
    pH$DEP[i] = 2
  }else if(pH$DT[i] > '2018-09-22 11:00:00' & pH$DT[i] < '2019-03-04 10:00') {
    pH$DEP[i] = 4
  }else if(pH$DT[i] > '2019-03-04 11:00:00' & pH$DT[i] < '2019-07-20 9:00') {
    pH$DEP[i] = 5
  }else if(pH$DT[i] > '2019-07-20 11:00:00' & pH$DT[i] < '2019-11-22 9:00') {
    pH$DEP[i] = 6
  }else if(pH$DT[i] > '2019-11-22 11:00:00' & pH$DT[i] < '2020-02-20 9:00') {
    pH$DEP[i] = 7
  }else if(pH$DT[i] > '2020-02-20 10:00:00' & pH$DT[i] < '2020-07-07 10:00') {
    pH$DEP[i] = 8
  }else if(pH$DT[i] > '2020-08-12 10:00:00' & pH$DT[i] < '2020-08-14 9:00') {
    pH$DEP[i] = 9
  }else if(pH$DT[i] > '2020-11-12 10:00:00' & pH$DT[i] < '2021-06-02 7:00') {
    pH$DEP[i] = 10
  }else if(pH$DT[i] > '2021-06-02 13:00:00' & pH$DT[i] < '2021-09-12 7:00') {
    pH$DEP[i] = 11
  }else{pH$DEP[i] = NA
  }
}
pH <- pH[complete.cases(pH$DEP),]

DO <- pHDat5[complete.cases(pHDat5$DO.mgL),]
DO = DO %>%
  filter(QF_DO == "1")
DO <- DO %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(DO.mgL))
DO <- DO[complete.cases(DO$avg),]
DO <- merge(x = DO, y = pHDat5, by = "DT", all.x = TRUE)
DO <- DO[c(1,2,20)]

for (i in 1:nrow(DO)) {
  if(DO$DT[i] > '2017-12-07 12:00:00' & DO$DT[i] < '2018-02-15 6:00') {
    DO$DEP[i] = 1
  }else if(DO$DT[i] > '2018-02-15 11:00:00' & DO$DT[i] < '2018-07-21 9:00') {
    DO$DEP[i] = 2
  }else if(DO$DT[i] > '2018-07-21 11:00:00' & DO$DT[i] < '2018-09-22 11:00') {
    DO$DEP[i] = 3
  }else if(DO$DT[i] > '2018-09-22 11:00:00' & DO$DT[i] < '2019-03-04 10:00') {
    DO$DEP[i] = 4
  }else if(DO$DT[i] > '2019-03-04 11:00:00' & DO$DT[i] < '2019-07-20 9:00') {
    DO$DEP[i] = 5
  }else if(DO$DT[i] > '2019-07-20 11:00:00' & DO$DT[i] < '2019-11-22 9:00') {
    DO$DEP[i] = 6
  }else if(DO$DT[i] > '2019-11-22 11:00:00' & DO$DT[i] < '2020-02-20 9:00') {
    DO$DEP[i] = 7
  }else if(DO$DT[i] > '2020-02-20 10:00:00' & DO$DT[i] < '2020-08-12 10:00') {
    DO$DEP[i] = 8
  }else if(DO$DT[i] > '2020-08-12 10:00:00' & DO$DT[i] < '2020-11-12 9:00') {
    DO$DEP[i] = 9
  }else if(DO$DT[i] > '2020-11-12 10:00:00' & DO$DT[i] < '2021-06-02 7:00') {
    DO$DEP[i] = 10
  }else if(DO$DT[i] > '2021-06-02 13:00:00' & DO$DT[i] < '2021-09-12 7:00') {
    DO$DEP[i] = 11
  }else{DO$DEP[i] = NA
  }
}
DO <- DO[complete.cases(DO$DEP),]

TC <- pHDat5[complete.cases(pHDat5$pH.Temp),]

TC = TC %>%
  filter(QF_TC == "1") 

TC <- TC %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH.Temp))
TC <- TC[complete.cases(TC$avg),]
TC <- merge(x = TC, y = pHDat5, by = "DT", all.x = TRUE)
TC <- TC[c(1,2,20)]

for (i in 1:nrow(TC)) {
  if(TC$DT[i] > '2017-12-07 12:00:00' & TC$DT[i] < '2018-02-15 6:00') {
    TC$DEP[i] = 1
  }else if(TC$DT[i] > '2018-02-15 11:00:00' & TC$DT[i] < '2018-07-21 9:00') {
    TC$DEP[i] = 2
  }else if(TC$DT[i] > '2018-07-21 11:00:00' & TC$DT[i] < '2018-09-22 11:00') {
    TC$DEP[i] = 3
  }else if(TC$DT[i] > '2018-09-22 11:00:00' & TC$DT[i] < '2019-03-04 10:00') {
    TC$DEP[i] = 4
  }else if(TC$DT[i] > '2019-03-04 11:00:00' & TC$DT[i] < '2019-07-20 9:00') {
    TC$DEP[i] = 5
  }else if(TC$DT[i] > '2019-07-20 11:00:00' & TC$DT[i] < '2019-11-22 9:00') {
    TC$DEP[i] = 6
  }else if(TC$DT[i] > '2019-11-22 11:00:00' & TC$DT[i] < '2020-02-20 9:00') {
    TC$DEP[i] = 7
  }else if(TC$DT[i] > '2020-02-20 10:00:00' & TC$DT[i] < '2020-07-07 10:00') {
    TC$DEP[i] = 8
  }else if(TC$DT[i] > '2020-08-12 10:00:00' & TC$DT[i] < '2020-11-12 9:00') {
    TC$DEP[i] = 9
  }else if(TC$DT[i] > '2020-11-12 10:00:00' & TC$DT[i] < '2021-05-01 7:00') {
    TC$DEP[i] = 10
  }else if(TC$DT[i] > '2021-06-02 13:00:00' & TC$DT[i] < '2021-09-12 7:00') {
    TC$DEP[i] = 11
  }else{TC$DEP[i] = NA
  }
}
TC <- TC[complete.cases(TC$DEP),]

AllPB <- merge(x = pH, y = DO, by = "DT", all = TRUE)
AllPB <- merge(x = AllPB, y = TC, by = "DT", all = TRUE)

names(AllPB) <- c('DT', 'pH','Dep1','DO','Dep2','TC','Dep3')
AllPB_Good <- AllPB


AllPB <- AllPB[complete.cases(AllPB$pH),]
AllPB <- AllPB[complete.cases(AllPB$DO),]
AllPB <- AllPB[complete.cases(AllPB$TC),]

#Turn all dates to same year in order to separate by month
AllPB <- AllPB %>%
  mutate(date=ymd_hm(format(AllPB$DT, "2017-%m-%d-%H:%M")))

AllPB_Good <- AllPB_Good %>%
  mutate(date=ymd_hm(format(AllPB_Good$DT, "2017-%m-%d-%H:%M")))

pHPB <- ggplot(AllPB_Good, aes(x = DT, y = pH, group = Dep1)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)", y = "pH", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
pHPB

TCPB <- ggplot(AllPB_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)",element_text(size = 20)) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
TCPB 

DOPB <- ggplot(AllPB_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_line() +
  theme_classic() +
  labs(y="Dissolved Oxygen (mg/L)", x="Date (MM/YY)", element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0))
DOPB 

PB <- ggarrange(pHPB, TCPB, DOPB, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3, 
                legend = "bottom", common.legend = FALSE, align = "hv")
PB

PBscat <- ggplot(AllPB, aes(x = pH, y = DO, color = TC)) + 
  geom_point(size = 4, shape = 1, position = "jitter") +
  scale_color_gradientn(colours = c("salmon2", "red4", "gold2", "forestgreen", "turquoise2", "royalblue4", "blueviolet"),
                        values = c(1.0,0.9,0.7,0.5,0.3,0.1,0), name = "Temp (\u00B0C)", limits = c(7, 24)) +
  labs(y=element_blank(), x=element_blank()) +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  xlim(7.3, 8.4) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'))
PBscat <- PBscat + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5, barwidth = 1, barheight = 20))+ scale_y_continuous(limits = c(0.5,11), breaks = c(2,4,6,8,10))
PBscat

PBscatpHtemp <- ggplot(AllPB, aes(x = TC, y = pH, color = DO)) + 
  geom_point(size = 4, shape = 1, position = "jitter") +
  scale_color_gradientn(colours = c("salmon2", "red4", "gold2", "forestgreen", "turquoise2", "royalblue4", "blueviolet"),
                        values = c(1.0,0.9,0.7,0.5,0.3,0.1,0), name = "DO mg/L", limits = c(1, 10)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(y=element_blank(), x=element_blank()) +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  xlim(4, 20) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'))
PBscatpHtemp <- PBscatpHtemp + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5, barwidth = 1, barheight = 20))+ scale_y_continuous(limits = c(7.5,8.2), breaks = c(7.5,7.6,7.7,7.8,7.9,8.0,8.1,8.2))
PBscatpHtemp

PB1<-
  PBscatpHtemp + annotate(geom="text", x=15, y=7.6, label="y = 0.04 + 7.39, Adj R2 = 0.31",
                               color="black")
PB1

#Regressions of time series 
pHDO <- lm(AllPB$DO ~ AllPB$pH)
summary(pHDO)
pHTC <- lm(AllPB$TC ~ AllPB$pH)
summary(pHTC)
DOTC <- lm(AllPB$TC ~ AllPB$DO)
summary(DOTC)
TCpH <- lm(AllPB$pH ~ AllPB$TC)
summary(TCpH)

pHQ <- quantile(AllPB_Good$TC, na.rm = TRUE)
pHave <- mean(AllPB_Good$TC, na.rm = TRUE)
pHstd <- sd(AllPB_Good$TC, na.rm = TRUE)

pHPBmodel = read.csv("phdot_pb.csv", header = FALSE)
colnames(pHPBmodel) <- c("pH", "DO", "T")
pHPBmodel$DOconverted <-pHPBmodel$DO* 1.42903 

PBmodelscatpHtemp<- ggplot(pHPBmodel, aes(x = T, y = pH, color = DOconverted)) + 
  geom_point(size = 4, shape = 1, position = "jitter") +
  scale_color_gradientn(colours = c("salmon2", "red4", "gold2", "forestgreen", "turquoise2", "royalblue4", "blueviolet"),
                        values = c(1.0,0.9,0.7,0.5,0.3,0.1,0), name = "DO mg/L", limits = c(1, 10)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(y=element_blank(), x=element_blank()) +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  xlim(4, 20) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'))
PBmodelscatpHtemp <- PBmodelscatpHtemp + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5, barwidth = 1, barheight = 20))+ scale_y_continuous(limits = c(7.5,8.2), breaks = c(7.5,7.6,7.7, 7.8,7.9,8.0,8.1,8.2))
PBmodelscatpHtemp

PB2<-
  PBmodelscatpHtemp + annotate(geom="text", x=15, y=7.6, label="y = 0.05 + 7.32, Adj R2 = 0.61",
                                   color="black")
PB2

#Regressions of time series 
pHDOmodel <- lm(pHPBmodel$DOconverted ~ pHPBmodel$pH)
summary(pHDOmodel)
pHTCmodel <- lm(pHPBmodel$T ~ pHPBmodel$pH)
summary(pHTCmodel)
DOTCmodel <- lm(pHPBmodel$T ~ pHPBmodel$DOconverted)
summary(DOTCmodel)
TCpHmodel <- lm(pHPBmodel$pH ~ pHPBmodel$T)
summary(TCpHmodel)

################################################################################################################
pHDat6 = read.csv(paste('/Volumes/GoogleDrive/My Drive/UCSC/R/BC_all.csv'),
                  skip = 7)
pHDat6$Lat = as.numeric(pHDat6$Lat)
pHDat6$Lon = as.numeric(pHDat6$Lon)
pHDat6$Depth = as.numeric(pHDat6$Depth)
pHDat6$pH = as.numeric(pHDat6$pH)
pHDat6$pH.Temp = as.numeric(pHDat6$pH.Temp)
pHDat6$QF_pH = as.numeric(pHDat6$QF_pH )
pHDat6$QF_TC = as.numeric(pHDat6$QF_TC)
pHDat6$DO..umol.kg. = as.numeric(pHDat6$DO..umol.kg.)
pHDat6$DO.mgL <- pHDat6$DO..umol.kg./31.2512
pHDat6$DO..Sat = as.numeric(pHDat6$DO..Sat)
pHDat6$DO.Temp = as.numeric(pHDat6$DO.Temp)
pHDat6$QF_DO = as.numeric(pHDat6$QF_DO)
pHDat6$year = as.numeric(str_sub(pHDat6$pH_Date_UTC_ddmmyyyy,-4,-1))
pHDat6$month = as.numeric(str_sub(pHDat6$pH_Date_UTC_ddmmyyyy,-6,-5))
pHDat6$day = as.numeric(str_sub(pHDat6$pH_Date_UTC_ddmmyyyy,-8,-7))
pHDat6$DT = as.POSIXct(paste(pHDat6$day, pHDat6$month, pHDat6$year, 
                             pHDat6$pH_Time_UTC_hhmmss), format = "%d%m%Y %H:%M:%S", tz = "UTC")

pHDat6 <- pHDat6[complete.cases(pHDat6$Lat),]

for (i in 1:nrow(pHDat6)) {
  if(pHDat6$DT[i] > '2017-10-26 14:00:00' & pHDat6$DT[i] < '2018-01-30 6:00') {
    pHDat6$DEP[i] = 1
  }else if(pHDat6$DT[i] > '2018-01-30 10:00:00' & pHDat6$DT[i] < '2018-05-18 11:00') {
    pHDat6$DEP[i] = 2
  }else if(pHDat6$DT[i] > '2018-05-18 10:00:00' & pHDat6$DT[i] < '2018-08-18 10:00') {
    pHDat6$DEP[i] = 3
  }else if(pHDat6$DT[i] > '2018-08-18 10:00:00' & pHDat6$DT[i] < '2018-11-11 10:00') {
    pHDat6$DEP[i] = 4
  }else if(pHDat6$DT[i] > '2018-11-11 12:00:00' & pHDat6$DT[i] < '2019-02-01 11:00') {
    pHDat6$DEP[i] = 5
  }else if(pHDat6$DT[i] > '2019-02-24 15:00:00' & pHDat6$DT[i] < '2019-07-26 12:00') {
    pHDat6$DEP[i] = 6
  }else if(pHDat6$DT[i] > '2019-07-26 12:00:00' & pHDat6$DT[i] < '2019-10-08 9:00') {
    pHDat6$DEP[i] = 7
  }else if(pHDat6$DT[i] > '2019-10-08 10:00:00' & pHDat6$DT[i] < '2020-02-14 11:00') {
    pHDat6$DEP[i] = 8
  }else if(pHDat6$DT[i] > '2020-02-14 11:00:00' & pHDat6$DT[i] < '2020-07-23 12:00') {
    pHDat6$DEP[i] = 9
  }else if(pHDat6$DT[i] > '2020-07-23 12:00:00' & pHDat6$DT[i] < '2021-03-28 8:00') {
    pHDat6$DEP[i] = 10
  }else if(pHDat6$DT[i] > '2021-03-28 12:00:00' & pHDat6$DT[i] < '2021-08-13 10:00') {
    pHDat6$DEP[i] = 11
  }else if(pHDat6$DT[i] > '2021-08-13 12:00:00' & pHDat6$DT[i] < '2022-04-18 8:00') {
    pHDat6$DEP[i] = 12
  }else{pHDat6$DEP[i] = NA
  }
}

#Average every hour
pH <- subset(pHDat6, pHDat6$pH< 8.5)
pH = pH %>%
  filter(QF_pH == "1") 

pH <- pH %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH))
pH <- pH[complete.cases(pH$avg),]
pH <- merge(x = pH, y = pHDat6, by = "DT", all.x = TRUE)
pH <- pH[c(1,2,20)]

for (i in 1:nrow(pH)) {
  if(pH$DT[i] > '2017-10-26 14:00:00' & pH$DT[i] < '2018-01-30 6:00') {
    pH$DEP[i] = 1
  }else if(pH$DT[i] > '2018-01-30 10:00:00' & pH$DT[i] < '2018-05-18 11:00') {
    pH$DEP[i] = 2
  }else if(pH$DT[i] > '2018-05-18 10:00:00' & pH$DT[i] < '2018-08-18 10:00') {
    pH$DEP[i] = 3
  }else if(pH$DT[i] > '2018-08-18 10:00:00' & pH$DT[i] < '2018-11-11 10:00') {
    pH$DEP[i] = 4
  }else if(pH$DT[i] > '2018-11-11 12:00:00' & pH$DT[i] < '2019-02-01 11:00') {
    pH$DEP[i] = 5
  }else if(pH$DT[i] > '2019-02-24 15:00:00' & pH$DT[i] < '2019-07-26 12:00') {
    pH$DEP[i] = 6
  }else if(pH$DT[i] > '2019-07-26 12:00:00' & pH$DT[i] < '2019-10-08 9:00') {
    pH$DEP[i] = 7
  }else if(pH$DT[i] > '2019-10-08 10:00:00' & pH$DT[i] < '2020-02-14 11:00') {
    pH$DEP[i] = 8
  }else if(pH$DT[i] > '2020-02-14 11:00:00' & pH$DT[i] < '2020-06-27 12:00') {
    pH$DEP[i] = 9
  }else if(pH$DT[i] > '2020-07-23 12:00:00' & pH$DT[i] < '2020-08-08 8:00') {
    pH$DEP[i] = 10
  }else if(pH$DT[i] > '2021-03-28 12:00:00' & pH$DT[i] < '2021-08-13 10:00') {
    pH$DEP[i] = 11
  }else if(pH$DT[i] > '2021-08-13 12:00:00' & pH$DT[i] < '2022-04-18 8:00') {
    pH$DEP[i] = 12
  }else{pH$DEP[i] = NA
  }
}
pH <- pH[complete.cases(pH$DEP),]

DO <- pHDat6[complete.cases(pHDat6$DO.mgL),]
DO = DO %>%
  filter(QF_DO == "1")
DO <- DO %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(DO.mgL))
DO <- DO[complete.cases(DO$avg),]
DO <- merge(x = DO, y = pHDat6, by = "DT", all.x = TRUE)
DO <- DO[c(1,2,20)]

for (i in 1:nrow(DO)) {
  if(DO$DT[i] > '2017-10-26 14:00:00' & DO$DT[i] < '2018-01-30 6:00') {
    DO$DEP[i] = 1
  }else if(DO$DT[i] > '2018-01-30 10:00:00' & DO$DT[i] < '2018-05-18 11:00') {
    DO$DEP[i] = 2
  }else if(DO$DT[i] > '2018-05-18 10:00:00' & DO$DT[i] < '2018-08-18 10:00') {
    DO$DEP[i] = 3
  }else if(DO$DT[i] > '2018-08-18 10:00:00' & DO$DT[i] < '2018-11-11 10:00') {
    DO$DEP[i] = 4
  }else if(DO$DT[i] > '2018-11-11 12:00:00' & DO$DT[i] < '2019-02-01 11:00') {
    DO$DEP[i] = 5
  }else if(DO$DT[i] > '2019-02-24 15:00:00' & DO$DT[i] < '2019-07-26 12:00') {
    DO$DEP[i] = 6
  }else if(DO$DT[i] > '2019-07-26 12:00:00' & DO$DT[i] < '2019-10-08 9:00') {
    DO$DEP[i] = 7
  }else if(DO$DT[i] > '2019-10-08 10:00:00' & DO$DT[i] < '2020-02-14 11:00') {
    DO$DEP[i] = 8
  }else if(DO$DT[i] > '2020-02-14 11:00:00' & DO$DT[i] < '2020-07-23 12:00') {
    DO$DEP[i] = 9
  }else if(DO$DT[i] > '2020-07-23 12:00:00' & DO$DT[i] < '2021-03-28 8:00') {
    DO$DEP[i] = 10
  }else if(DO$DT[i] > '2021-03-28 12:00:00' & DO$DT[i] < '2021-08-13 10:00') {
    DO$DEP[i] = 11
  }else if(DO$DT[i] > '2021-08-13 12:00:00' & DO$DT[i] < '2022-04-18 8:00') {
    DO$DEP[i] = 12
  }else{DO$DEP[i] = NA
  }
}
DO <- DO[complete.cases(DO$DEP),]

for (i in 1:nrow(pHDat6)) {
  if(pHDat6$QF_TC[i]==1) {
    pHDat6$TC_Both[i] = pHDat6$pH.Temp[i]
  }else{
    pHDat6$TC_Both[i] = pHDat6$DO.Temp[i]
  }
}

TC <- pHDat6[complete.cases(pHDat6$TC_Both),]

TC <- TC %>%
  mutate(DT = floor_date(DT,"days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(TC_Both))

for (i in 1:nrow(TC)) {
  if(TC$DT[i] > '2017-10-26 14:00:00' & TC$DT[i] < '2018-01-30 6:00') {
    TC$DEP[i] = 1
  }else if(TC$DT[i] > '2018-01-30 10:00:00' & TC$DT[i] < '2018-05-18 11:00') {
    TC$DEP[i] = 2
  }else if(TC$DT[i] > '2018-05-18 10:00:00' & TC$DT[i] < '2018-08-18 10:00') {
    TC$DEP[i] = 3
  }else if(TC$DT[i] > '2018-08-18 10:00:00' & TC$DT[i] < '2018-11-11 10:00') {
    TC$DEP[i] = 4
  }else if(TC$DT[i] > '2018-11-11 12:00:00' & TC$DT[i] < '2019-02-01 11:00') {
    TC$DEP[i] = 5
  }else if(TC$DT[i] > '2019-02-24 15:00:00' & TC$DT[i] < '2019-07-26 12:00') {
    TC$DEP[i] = 6
  }else if(TC$DT[i] > '2019-07-26 12:00:00' & TC$DT[i] < '2019-10-08 9:00') {
    TC$DEP[i] = 7
  }else if(TC$DT[i] > '2019-10-08 10:00:00' & TC$DT[i] < '2020-02-14 11:00') {
    TC$DEP[i] = 8
  }else if(TC$DT[i] > '2020-02-14 11:00:00' & TC$DT[i] < '2020-06-27 12:00') {
    TC$DEP[i] = 9
  }else if(TC$DT[i] > '2020-07-23 12:00:00' & TC$DT[i] < '2021-03-28 8:00') {
    TC$DEP[i] = 10
  }else if(TC$DT[i] > '2021-03-28 12:00:00' & TC$DT[i] < '2021-08-13 10:00') {
    TC$DEP[i] = 11
  }else if(TC$DT[i] > '2021-08-13 12:00:00' & TC$DT[i] < '2022-04-18 8:00') {
    TC$DEP[i] = 12  
  }else{TC$DEP[i] = NA
  }
}
TC <- TC[complete.cases(TC$DEP),]

AllBC <- merge(x = pH, y = DO, by = "DT", all = TRUE)
AllBC <- merge(x = AllBC, y = TC, by = "DT", all = TRUE)

names(AllBC) <- c('DT', 'pH','Dep1','DO','Dep2','TC','Dep3')
AllBC_Good <- AllBC


AllBC <- AllBC[complete.cases(AllBC$pH),]
AllBC <- AllBC[complete.cases(AllBC$DO),]
AllBC <- AllBC[complete.cases(AllBC$TC),]

#Turn all dates to same year in order to separate by month
AllBC <- AllBC %>%
  mutate(date=ymd_hm(format(AllBC$DT, "2017-%m-%d-%H:%M")))
AllBC_Good <- AllBC_Good %>%
  mutate(date=ymd_hm(format(AllBC_Good$DT, "2017-%m-%d-%H:%M")))

pHBC <- ggplot(AllBC_Good, aes(x = DT, y = pH, group = Dep1)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)", y = "pH", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
pHBC

TCBC <- ggplot(AllBC_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)",element_text(size = 20)) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
TCBC 

DOBC <- ggplot(AllBC_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_line() +
  theme_classic() +
  labs(y="Dissolved Oxygen (mg/L)", x="Date (MM/YY)", element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0))
DOBC

BC <- ggarrange(pHBC, TCBC, DOBC, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3, 
                legend = "bottom", common.legend = FALSE, align = "hv")
BC

#Regressions of time series 
pHDO <- lm(AllBC$DO ~ AllBC$pH)
summary(pHDO)
pHTC <- lm(AllBC$TC ~ AllBC$pH)
summary(pHTC)
DOTC <- lm(AllBC$TC ~ AllBC$DO)
summary(DOTC)

pHQ <- quantile(AllBC_Good$DO, na.rm = TRUE)
pHave <- mean(AllBC_Good$DO, na.rm = TRUE)
pHstd <- sd(AllBC_Good$DO, na.rm = TRUE)

#################################################################################################################
##FigS2 Sub region######

AllPB_Good$Site <- "Point Buchon"
AllBC_Good$Site <- "Big Creek"
AllCen <- rbind(AllPB_Good, AllBC_Good)
AllCen$Group1 <- paste(AllCen$Dep1, AllCen$Site)
AllCen$Site <- factor(AllCen$Site, levels=c('Point Buchon','Big Creek'))

pHCen <- ggplot(AllCen, aes(x = DT, y = pH, group = Group1, color = Site)) + 
  geom_line() +
  theme_classic() +
  geom_hline(yintercept = 7.7, linetype = "dashed", color = "black", size = 0.75) +
  scale_color_manual(values=c("#33a02c","#b2df8a" )) +
  labs(x=element_blank(), y = element_blank(), element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank(),
        legend.position = c(0.85,0.85), legend.title = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0)) +
  theme(axis.text.y=element_blank()) +
  scale_y_continuous(limits = c(7.4, 8.4), breaks = c(7.4,7.6,7.8,8.0,8.2))
pHCen

TCCen <- ggplot(AllPB_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line(color = "#33a02c") +
  geom_line(data = AllBC_Good, aes(x = DT, y = TC, group = Dep3), color = "#b2df8a") +
  theme_classic() +
  labs(x=element_blank(), y = element_blank(), element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  theme(axis.text.y=element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0)) +
  ylim(7, 23)
TCCen 

DOCen <- ggplot(AllPB_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_line(color = "#33a02c") +
  geom_line(data = AllBC_Good, aes(x = DT, y = DO, group = Dep2), color = "#b2df8a") +
  theme_classic() +
  geom_hline(yintercept = 4.6, linetype = "dashed", color = "black", size = 0.75) +
  labs(x=element_blank(), y = element_blank(), element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 10), breaks = c(3,5,7,9)) +
  theme(axis.text.y=element_blank())
DOCen

Cen <- ggarrange(pHCen, TCCen, DOCen, labels = c("(b)", "(e)", "(h)"), ncol = 1, hjust = 0.5, nrow = 3, 
                 common.legend = FALSE, align = "hv")
Cen


AllLB_Good$Site <- "Laguna Beach"
AllCI_Good$Site <- "Catalina Island"
AllSo <- rbind(AllLB_Good, AllCI_Good)
AllSo$Group1 <- paste(AllSo$Dep1, AllSo$Site)
AllSo$Site <- factor(AllSo$Site, levels=c('Catalina Island','Laguna Beach'))

pHSo <- ggplot(AllSo, aes(x = DT, y = pH, group = Group1, color = Site)) + 
  geom_line() +
  theme_classic() +
  geom_hline(yintercept = 7.7, linetype = "dashed", color = "black", size = 0.75) +
  scale_color_manual(values=c("#e31a1c","#fb9a99")) +
  labs(x=element_blank(), y = element_blank(), element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank(),
        legend.position = c(0.85,0.85), legend.title = element_blank()) +
  theme(axis.text.y=element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0)) +
  scale_y_continuous(limits = c(7.4, 8.4), breaks = c(7.4,7.6,7.8,8.0,8.2))
pHSo

TCSo <- ggplot(AllCI_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line(color = "#e31a1c") +
  geom_line(data = AllLB_Good, aes(x = DT, y = TC, group = Dep3), color = "#fb9a99") +
  theme_classic() +
  labs(x=element_blank(), y = element_blank(), element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank(),
        axis.text.y=element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0)) +
  ylim(7, 23)
TCSo 

DOSo <- ggplot(AllCI_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_line(color = "#e31a1c") +
  geom_line(data = AllLB_Good, aes(x = DT, y = DO, group = Dep2), color = "#fb9a99") +
  theme_classic()+
  theme(axis.text.y=element_blank()) +
  geom_hline(yintercept = 4.6, linetype = "dashed", color = "black", size = 0.75) +
  labs(x=element_blank(), y = element_blank(), element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 10), breaks = c(3,5,7,9))
DOSo

So <- ggarrange(pHSo, TCSo, DOSo, labels = c("(c)", "(f)", "(i)"), ncol = 1, hjust = 0.5, nrow = 3, 
                  common.legend = FALSE, align = "hv")
So


AllVD_Good$Site <- "Van Damme"
AllPA_Good$Site <- "Point Arena"
AllNo <- rbind(AllVD_Good, AllPA_Good)
AllNo$Group1 <- paste(AllNo$Dep1, AllNo$Site)
AllNo$Site <- factor(AllNo$Site, levels=c('Point Arena','Van Damme'))

pHNo <- ggplot(AllNo, aes(x = DT, y = pH, group = Group1, color = Site)) + 
  geom_line() +
  theme_classic() +
  geom_hline(yintercept = 7.7, linetype = "dashed", color = "black", size = 0.75) +
  scale_color_manual(values=c("#1f78b4","#a6cee3")) +
  labs(x= element_blank(), y = "pH", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank(),
        legend.position = c(0.85,0.85), legend.title = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0)) +
  scale_y_continuous(limits = c(7.4, 8.4), breaks = c(7.4,7.6,7.8,8.0,8.2))
pHNo

TCNo <- ggplot(AllVD_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line(color = "#a6cee3") +
  geom_line(data = AllPA_Good, aes(x = DT, y = TC, group = Dep3), color = "#1f78b4") +
  theme_classic() +
  labs(x=element_blank(),element_text(size = 20)) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0)) +
  ylim(7, 23)
TCNo 

DONo <- ggplot(AllVD_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_line(color = "#a6cee3") +
  geom_line(data = AllPA_Good, aes(x = DT, y = DO, group = Dep2), color = "#1f78b4") +
  geom_hline(yintercept = 4.6, linetype = "dashed", color = "black", size = 0.75) +
  theme_classic() +
  labs(y="Dissolved Oxygen (mg/L)", x=element_blank(), element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 10), breaks = c(3,5,7,9))
DONo

No <- ggarrange(pHNo, TCNo, DONo, labels = c("(a)", "(d)", "(g)"), ncol = 1, hjust = -0.8, nrow = 3, 
                 common.legend = FALSE, align = "hv")
No

All_FT <- ggarrange(No, Cen, So, ncol = 3, nrow = 1, 
                 common.legend = FALSE, align = "hv")
All_FT

All_FT <- annotate_figure(All_FT,bottom = text_grob("Date (MM/YY)", size = 10)) 
All_FT

ggsave(plot = All_FT, file = "SuppFig.png", 
       type = "cairo-png",  bg = "white",
       width = 40, height = 25, units = "cm", dpi = 300)

#################################################################################################################
##Fig2 Subset plotted together###############################################################################

SubRegion <- rbind(AllPB_Good, AllPA_Good, AllCI_Good)
SubRegion$Group1 <- paste(SubRegion$Dep1, SubRegion$Site)

SubRegion$Site <- factor(SubRegion$Site, levels=c('Point Arena','Point Buchon', 'Catalina Island'))
pHSub <- ggplot(SubRegion, aes(x = DT, y = pH, group = Group1, color = Site)) + 
  geom_hline(yintercept = 7.7, linetype = "dashed", color = "black", size = 0.75) +
  geom_line() +
  theme_classic() +
  scale_color_manual(values=c("#1f78b4","#33a02c","#e31a1c")) +
  labs(x= element_blank(), y = "pH", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank(),
        legend.position = c(0.9,0.9), legend.title = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0)) +
  scale_y_continuous(limits = c(7.4, 8.4), breaks = c(7.4,7.6,7.8,8.0,8.2))
pHSub

TCSub <- ggplot(AllPA_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line(color = "#1f78b4") +
  geom_line(data = AllPB_Good, aes(x = DT, y = TC, group = Dep3), color = "#33a02c") +
  geom_line(data = AllCI_Good, aes(x = DT, y = TC, group = Dep3), color = "#e31a1c") +
  theme_classic() +
  labs(x=element_blank(),element_text(size = 20)) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0)) +
  ylim(7, 23)
TCSub 

DOSub <- ggplot(AllPA_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_hline(yintercept = 4.6, linetype = "dashed", color = "black", size = 0.75) +
  geom_line(color = "#1f78b4") +
  geom_line(data = AllPB_Good, aes(x = DT, y = DO, group = Dep2), color = "#33a02c") +
  geom_line(data = AllCI_Good, aes(x = DT, y = DO, group = Dep2), color = "#e31a1c") +
  theme_classic() +
  labs(y="Dissolved Oxygen (mg/L)", x=element_blank(), element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 10), breaks = c(3,5,7,9))
DOSub

SubReg <- ggarrange(pHSub, TCSub, DOSub, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3, 
                common.legend = FALSE, align = "hv")
SubReg
SubReg <- annotate_figure(SubReg,bottom = text_grob("Date (MM/YY)", size = 10)) 
SubReg

ggsave(plot = SubReg, file = "Fig2.png", 
       type = "cairo-png",  bg = "white",
       width = 30, height = 20, units = "cm", dpi = 300)

##########################################################################################
#####Fig3 Correlations#############################################
Scat <- ggarrange(PA1, PB1, CI1, PA2, PB2, CI2, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), vjust = 1.2, hjust = 0.2, ncol = 3, nrow = 2, 
                    common.legend = TRUE, legend = "right", align = "hv")
Scat

Scat <- annotate_figure(Scat, left = textGrob("Dissolved oxygen (mg/L)", rot = 90, vjust = 0.2, gp = gpar(cex = 1.5)))
Scat
ggsave(plot = Scat, file = "Scatter.png", 
       bg = "white",
       width = 34, height = 27, units = "cm", dpi = 300)


##########################################################################################
#####FigS1 Climatology##############
AllCI_Good$Year <- year(AllCI_Good$DT)
AllCI_Good$Group1 <- paste(AllCI_Good$Dep1,AllCI_Good$Year)
AllCI_Good$Group3 <- paste(AllCI_Good$Dep3,AllCI_Good$Year)
AllCI_Good$Group2 <- paste(AllCI_Good$Dep2,AllCI_Good$Year)

AllPA_Good$Year <- year(AllPA_Good$DT)
AllPA_Good$Group3 <- paste(AllPA_Good$Dep3,AllPA_Good$Year)
AllPA_Good$Group2 <- paste(AllPA_Good$Dep2,AllPA_Good$Year)

AllPB_Good$Year <- year(AllPB_Good$DT)
AllPB_Good$Group1 <- paste(AllPB_Good$Dep1,AllPB_Good$Year)
AllPB_Good$Group3 <- paste(AllPB_Good$Dep3,AllPB_Good$Year)
AllPB_Good$Group2 <- paste(AllPB_Good$Dep2,AllPB_Good$Year)

pHYr <- ggplot(AllPA_Good, aes(x = date, y = pH, group = factor(Year), color = factor(Year))) + 
  geom_hline(yintercept = 7.7, linetype = "dashed", color = "black", size = 0.75) +
  geom_line() +
  scale_color_manual(values=c("#1f78b4","#33a02c","#e31a1c","#984ea3", "#ff7f00")) +
  theme_classic() +
  guides(color=guide_legend(title="Year")) +
  labs(x= element_blank(), y = "pH", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "2 months", expand = c(0, 0)) +
  scale_y_continuous(limits = c(7.4, 8.4), breaks = c(7.4,7.6,7.8,8.0,8.2))
pHYr

TCYr <- ggplot(AllPA_Good, aes(x = date, y = TC, group = Group3, color = factor(Year))) + 
  geom_line() +
  scale_color_manual(values=c("#1f78b4","#33a02c","#e31a1c","#984ea3", "#ff7f00")) +
  theme_classic() +
  labs(x=element_blank(),element_text(size = 20)) +
  guides(color=guide_legend(title="Year")) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "2 months", expand = c(0, 0)) +
  ylim(7, 23)
TCYr 

DOYr <- ggplot(AllPA_Good, aes(x = date, y = DO, group = Group2, color = factor(Year))) + 
  geom_hline(yintercept = 4.6, linetype = "dashed", color = "black", size = 0.75) +
  geom_line() +
  scale_color_manual(values=c("#1f78b4","#33a02c","#e31a1c","#984ea3", "#ff7f00")) +
  theme_classic() +
  guides(color=guide_legend(title="Year")) +
  labs(y="Dissolved Oxygen (mg/L)", x=element_blank(), element_text(size = 20)) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b", expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 10), breaks = c(3,5,7,9))
DOYr

PAClim <- ggarrange(pHYr, TCYr, DOYr, labels = c("(a)", "(d)", "(g)"), ncol = 1, nrow = 3, 
                    common.legend = TRUE, align = "hv", hjust = -0.8)
PAClim


pHYr <- ggplot(AllPB_Good, aes(x = date, y = pH, group = Group1, color = factor(Year))) + 
  geom_hline(yintercept = 7.7, linetype = "dashed", color = "black", size = 0.75) +
  geom_line() +
  scale_color_manual(values=c("#1f78b4","#33a02c","#e31a1c","#984ea3", "#ff7f00")) +
  theme_classic() +
  guides(color=guide_legend(title="Year")) +
  theme(axis.text.y=element_blank()) +
  labs(x= element_blank(), y = element_blank(), element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "2 months", expand = c(0, 0)) +
  scale_y_continuous(limits = c(7.4, 8.4), breaks = c(7.4,7.6,7.8,8.0,8.2))
pHYr

TCYr <- ggplot(AllPB_Good, aes(x = date, y = TC, group = Group3, color = factor(Year))) + 
  geom_line() +
  scale_color_manual(values=c("#1f78b4","#33a02c","#e31a1c","#984ea3", "#ff7f00")) +
  theme_classic() +
  theme(axis.text.y=element_blank()) +
  labs(x= element_blank(), y = element_blank(), element_text(size = 20)) +
  guides(color=guide_legend(title="Year")) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "2 months", expand = c(0, 0)) +
  ylim(7, 23)
TCYr 

DOYr <- ggplot(AllPB_Good, aes(x = date, y = DO, group = Group2, color = factor(Year))) + 
  geom_hline(yintercept = 4.6, linetype = "dashed", color = "black", size = 0.75) +
  geom_line() +
  scale_color_manual(values=c("#1f78b4","#33a02c","#e31a1c","#984ea3", "#ff7f00")) +
  theme_classic() +
  guides(color=guide_legend(title="Year")) +
  theme(axis.text.y=element_blank()) +
  labs(x= element_blank(), y = element_blank(), element_text(size = 20)) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b", expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 10), breaks = c(3,5,7,9))
DOYr

PBClim <- ggarrange(pHYr, TCYr, DOYr, labels = c("(b)", "(e)", "(h)"), ncol = 1, nrow = 3, 
                    common.legend = TRUE, align = "hv", hjust = 0.5)
PBClim

pHYr <- ggplot(AllCI_Good, aes(x = date, y = pH, group = Group1, color = factor(Year))) + 
  geom_hline(yintercept = 7.7, linetype = "dashed", color = "black", size = 0.75) +
  geom_line() +
  scale_color_manual(values=c("#1f78b4","#33a02c","#e31a1c","#984ea3", "#ff7f00")) +
  theme_classic() +
  guides(color=guide_legend(title="Year")) +
  theme(axis.text.y=element_blank()) +
  labs(x= element_blank(), y = element_blank(), element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "2 months", expand = c(0, 0)) +
  scale_y_continuous(limits = c(7.4, 8.4), breaks = c(7.4,7.6,7.8,8.0,8.2))
pHYr

TCYr <- ggplot(AllCI_Good, aes(x = date, y = TC, group = Group3, color = factor(Year))) + 
  geom_line() +
  scale_color_manual(values=c("#1f78b4","#33a02c","#e31a1c","#984ea3", "#ff7f00")) +
  theme_classic() +
  theme(axis.text.y=element_blank()) +
  labs(x= element_blank(), y = element_blank(), element_text(size = 20)) +
  guides(color=guide_legend(title="Year")) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "2 months", expand = c(0, 0)) +
  ylim(7, 23)
TCYr 

DOYr <- ggplot(AllCI_Good, aes(x = date, y = DO, group = Group2, color = factor(Year))) + 
  geom_hline(yintercept = 4.6, linetype = "dashed", color = "black", size = 0.75) +
  geom_line() +
  scale_color_manual(values=c("#1f78b4","#33a02c","#e31a1c","#984ea3", "#ff7f00")) +
  theme_classic() +
  guides(color=guide_legend(title="Year")) +
  theme(axis.text.y=element_blank()) +
  labs(x= element_blank(), y = element_blank(), element_text(size = 20)) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b", expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 10), breaks = c(3,5,7,9))
DOYr

CIClim <- ggarrange(pHYr, TCYr, DOYr, labels = c("(c)", "(f)", "(i)"), ncol = 1, nrow = 3, 
                    common.legend = TRUE, align = "hv", hjust = 0.5)
CIClim

All_Clim <- ggarrange(PAClim, PBClim, CIClim, ncol = 3, nrow = 1, 
                    common.legend = FALSE, align = "hv")
All_Clim

All_Clim <- annotate_figure(All_Clim,bottom = text_grob("Date (Month)", size = 10)) 
All_Clim

ggsave(plot = All_Clim, file = "SuppFig1.png", 
       type = "cairo-png",  bg = "white",
       width = 40, height = 25, units = "cm", dpi = 300)

##########################################################################################

