
path <- ""
setwd(path)

library(tidyverse)
library(dplyr)
library(vcd)
library(car)
library(MASS)
library(faraway)
library(jtools)
library(ggcorrplot)
library(AER)
library(devtools)
library(Rtools)
library(merTools)
library(gtsummary)


library(lme4)
library(glmmTMB)
#############################
#########READ ALL DATA#######
#############################

#main
ongoing.waits <- read.csv("sot_performance_ongoing_waits_mar22.csv")
#step 1
data <- ongoing.waits[ongoing.waits["PatientType"]=="Inpatient/Day case",]
hb.demography <- read.csv("HBdemography.csv")


#the rest
simd <- read.csv("simd2020v2_22062020.csv", quote="")
hb.region <- read.csv("HB_region.csv")
speciality.aggregates <- read.csv("speciality_aggregates.csv")
hb.14.19 <- read.csv("hb14_hb19.csv")
specialty.ref <- read.csv("specialty-reference.csv")

######################
###Data Preparation###
######################

#checking for aggregate specialty 
#step 14 
sum(data[data$Specialty=="Z9","NumberWaiting"], na.rm = TRUE) #14818720
sum(data[data$Specialty!="Z9","NumberWaiting"], na.rm = TRUE) #14686927
##the difference between is 131793
nrow(data[data$Specialty=="Z9",]) #1826
nrow(data[data$Specialty!="Z9",]) #36598
#trying to find what Z9 might have used for
unique(data[data$Specialty=="Z9","HBT"]) #for all HBT types including null


#get rid of unnecessary columns
cols.of.int <- c("MonthEnding", "HBT", "Specialty", "NumberWaiting", "NumberWaitingQF")
data.reduced <- data[cols.of.int]

#split date into month and year 
data.reduced["Year"] <- substring(data.reduced$MonthEnding, first=1, last=4)
data.reduced["Month"] <- substring(data.reduced$MonthEnding, first=5, last=6)
data.reduced["Year"] <- strtoi(data.reduced$Year)
data.month.year <- data.reduced[,-1]

#get rid of scotland and jubilee HBT 
#scotland S92000003
#jubilee HBT SB0801
data.HBT.updated <- data.month.year[!(data.month.year$HBT=="S92000003"),]
data.HBT.updated <- data.HBT.updated[!(data.HBT.updated$HBT=="SB0801"),]

#get rid of aggregate specialty 
data.z9.omit <- data.HBT.updated[data.HBT.updated$Specialty!="Z9",]

#get rid of "null" entry in HBT (there are total of 4)
data.HBT.14 <-data.z9.omit[!data.z9.omit$HBT=="null",]

#check number of NAs 
data.HBT.14 %>%
  summarise_all(funs(sum(is.na(.))))
#Only NumberWaiting has 3530 NAs 

#get rid of all NA rows
#step 5
data.complete <- na.omit(data.HBT.14)



##group months into four groups
data.seasonal <- data.complete
data.seasonal["Season"] <- NA
data.seasonal[data.seasonal$Month %in% c("01","02","03"), "Season"] = 1
data.seasonal[data.seasonal$Month %in% c("04","05","06"), "Season"] = 2
data.seasonal[data.seasonal$Month %in% c("07","08","09"), "Season"] = 3
data.seasonal[data.seasonal$Month %in% c("10","11","12"), "Season"] = 4

#add demography covariates
data.pop <- data.seasonal
hb.demography.reduced <- hb.demography[hb.demography$Sex=="All", c(2,3,7)]

#add average age 
age.cols <- 8:98
age.levels <- 0:90
matrix.demography <- data.matrix(hb.demography[hb.demography$Sex=="All", age.cols])
m.total.age <-  t(t(matrix.demography) * age.levels)
age.sum <- rowSums(m.total.age)
hb.demography.reduced["AgeSum"] <-  age.sum
hb.demography.reduced["AvgAge"] <- hb.demography.reduced$AgeSum/hb.demography.reduced$AllAges
hb.demography.reduced$AgeSum <- NULL

#add female to male ratio
females <- hb.demography[hb.demography$Sex=="Female", "AllAges"]
males <- hb.demography[hb.demography$Sex=="Male", "AllAges"]
hb.demography.reduced["FMRatio"] <- females / males

#join demographic covariates to main dataframe
data.pop <- left_join(data.pop, hb.demography.reduced, 
                      by=c('Year'='Year', 'HBT'='HB'))


#insert average simd for each hbt on dataset
simd.by.hb = simd %>% 
  group_by(HB)  %>%
  summarise(AvgSimd20 = mean(SIMD2020V2Rank),
            .groups = 'drop')
data.simd <- left_join(data.pop, simd.by.hb, by=c("HBT"="HB"))

#Add regions to dataframe
##aggregate hb into regions 
hb.region.df <- left_join(hb.14.19, hb.region,
                          by=c("HBName"="Board.Name"))
hb.region.df <- hb.region.df[!is.na(hb.region.df$Region),]
hb.region.reduced <- hb.region.df[c("HB", "Region","HBName")]
##add it to data frame
data.region <- left_join(data.simd, hb.region.reduced, by=c("HBT"="HB"))
##drop HB codes


#Add aggregate specialties to dataframe
spec.list <- unique(data.region$Specialty)
spec.df <- data.frame(Specialty = spec.list)
spec.df <- left_join(spec.df, specialty.ref, by=("Specialty"="Specialty"))
spec.df <- left_join(spec.df, speciality.aggregates, 
                     by=c("SpecialtyName"="Speciality"))
data.spec.agg <- left_join(data.region, spec.df, by=Specialty)

#exclude last two years to use demographic data
#last two years excluded to use population data
data.pop.not.na <- data.spec.agg[!is.na(data.spec.agg$AllAges),]

#factorize categorical variables 
#step 9
data.cat <- data.pop.not.na
data.cat$SpecialtyName <- factor(data.cat$SpecialtyName)
data.cat$HBName <- factor(data.cat$HBName)
data.cat$Aggregated.speciality <- factor(data.cat$Aggregated.speciality)
data.cat$Month <- factor(data.cat$Month)
data.cat$Season <- factor(data.cat$Season)
data.cat$Region <- factor(data.cat$Region)
data.cat$Specialty <- as.factor(data.cat$Specialty)
data.cat$HBT <- as.factor(data.cat$HBT)
data.cat["CatYear"] <- factor(data.cat$Year)

#normalize continuous vars
#step 10 
data.scaled <- data.cat
data.scaled["ScaledYear"] <- data.scaled$Year - min(data.scaled$Year)
data.scaled["ScaledPop"] <- scale(data.scaled$AllAges)
data.scaled["ScaledSimd"] <- scale(data.scaled$AvgSimd20)
data.scaled["ScaledAvgAge"] <- scale(data.scaled$AvgAge)
data.scaled["ScaledFMRatio"] <- scale(data.scaled$FMRatio)


#check if there is any unreliable data left
n.unrel <- sum(data.scaled[data.scaled$NumberWaitingQF==":u"])
data.scaled$NumberWaitingQF <- NULL

#just change the name of population column so that it's more clear
colnames(data.scaled)[7] <- "Population"

#add a standardised NumberWaiting column by population 
data.scaled.NW <- data.scaled
data.scaled.NW["ScaledNW"] <- data.scaled.NW$NumberWaiting / data.scaled.NW$Population



##########################################################
#last dataset to work on without transformation###########
#if any transformation is needed do it before this section
##########################################################,
saveRDS(data.scaled.NW, file = "data_model.RDS")
# rm(list=ls())
# gc()

#######################
### Final Dataframe ###
#######################

#main dataframe with all possible covariates
data.model.pop <- readRDS("data_model.RDS")



#################################################
### Exploratory Data Analysis on cleaned data####
#################################################


#max value and minimum interval for 90% of the data
max(data.model.pop$NumberWaiting) #8749
y <- data.model.pop$NumberWaiting
y <- sort(y)
n <- length(y) #23530
y.90 <- y[round(n*0.90)] #768
y.mean <- mean(y) #260.4
y.var <- var(y) #293520

#plot hist of NumberWaiting
# Histogram with density plot
ggplot(data.model.pop, aes(x=NumberWaiting)) + 
  geom_histogram(binwidth=50, color="black", fill="gray") +
  geom_vline(xintercept = y.mean, linetype="dashed", color = "red", size=0.5)


#plot bar plot for HB totals
by.HB <- data.model.pop %>% 
  group_by(HBName)  %>%
  summarise(Waiting = sum(NumberWaiting),
            Population = sum(Population),
            .groups = 'drop')
by.HB <- by.HB[order(by.HB$Waiting),]
by.HB <- left_join(by.HB, hb.region, by=c("HBName"="Board.Name"))
ggplot(by.HB,  aes(x = reorder(HBName, +Waiting), y=Waiting, fill=Region)) + 
  geom_bar(stat="identity", position="dodge") +
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) +
  labs(y= "patient count", x = "health board") +
  scale_fill_brewer(palette="Greens") +
  theme(text = element_text(size=9)) 

#plot population ratio 
by.HB["Ratio"] <- by.HB$Waiting/by.HB$Population
ggplot(by.HB, aes(x=reorder(HBName, +Population), y=Ratio, group=1)) +
  geom_line(stat="identity", color="black") +
  ylim(0, 0.00145)+
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) +
  labs(y= "patient count/population", x = "health board") +
  theme(text = element_text(size=9)) 

#region totals
by.region <- by.HB %>%
  group_by(Region) %>%
  summarise(Waiting = sum(Waiting),
            .groups = 'drop') 

#box plot of HB counts
ggplot(data.model.pop,  aes(x = reorder(HBName, +NumberWaiting), y=NumberWaiting, fill=Region)) + 
  geom_boxplot(position="dodge") +
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) +
  labs(y= "patient count", x = "health board") +
  scale_fill_brewer(palette="Greens") +
  theme(text = element_text(size=9)) +
  scale_x_discrete(limits=by.HB$HBName)

#bar plot for spec totals
by.spec <- data.model.pop %>% 
  group_by(SpecialtyName)  %>%
  summarise(Waiting = sum(NumberWaiting),
            .groups = 'drop')
by.spec <- by.spec[order(by.spec$Waiting),]
by.spec <- left_join(by.spec, speciality.aggregates, by=c("SpecialtyName"="Speciality"))
by.agg.spec <- by.spec  %>% 
  group_by(Aggregated.speciality)  %>%
  summarise(Waiting = sum(Waiting),
            .groups = 'drop')
by.agg.spec <- by.agg.spec[order(by.agg.spec$Waiting),]
ggplot(by.spec,  aes(x = reorder(SpecialtyName, +Waiting), y=Waiting, fill=Aggregated.speciality)) + 
  geom_bar(stat="identity", position="dodge") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1)) +
  labs(y= "patient count", x = "speciality") +
  theme(text = element_text(size=9)) +
  scale_fill_brewer(palette="Greens") 


#box plot of spec counts
ggplot(data.model.pop,  aes(x = reorder(Aggregated.speciality, +NumberWaiting), y=NumberWaiting, fill=Aggregated.speciality)) + 
  geom_boxplot(position="dodge") +
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) +
  labs(y= "patient count", x = "speciality") +
  theme(text = element_text(size=9)) +
  scale_fill_brewer(palette="Greens") 



#check for correlation between continuous variables 
corr.data <- data.model.pop[c(3,7,8,9,10)]
corr <- round(cor(corr.data), 1)
ggcorrplot(corr, hc.order = TRUE, type = "lower", lab = TRUE, lab_size = 4, 
           method="circle", colors = c("#006633","white"), 
           outline.color = "gray", show.legend = TRUE, show.diag = FALSE, 
           title="Correlogram of continuous variables")+
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust=1),
        axis.text.y=element_text(size=9),
        panel.grid.major=element_blank()) 





#group by months to use quarter or seasonal aggregates
#2012 is excluded because it only has november and december data
by.month = data.model.pop[!data.model.pop$Year==2012,] %>% 
  group_by(Month)  %>%
  summarise(TotalWaiting = sum(NumberWaiting),
            .groups = 'drop')
plot(by.month$TotalWaiting, ylim=c(470000, 540000), xlab="months", ylab="NumberWaiting")
abline(v=c(3.5, 6.5, 9.5), col="darkgreen", lty=2)
ratio.month <- max(by.month$TotalWaiting)/ min(by.month$TotalWaiting)


#check the relationship between NumberWaiting and
#population with respect to time
by.year.month <- data.model.pop %>%
  group_by(Year, Month) %>%
  summarise(sum.waiting = sum(NumberWaiting),
            sum.pop = sum(Population))
plot(by.year.month$sum.waiting)
plot(by.year.month$sum.pop, ylim=c(100000000, 135000000))



# check homogeneity of variances in between groups
leveneTest(data.model.pop$NumberWaiting, data.model.pop$HBT, center = mean)

#checking for collinearity 
m.lm <- lm(NumberWaiting ~  ScaledPop +
             ScaledSimd + ScaledAvgAge + ScaledFMRatio, data=data.model.pop)
vif(m.lm) #around 5 indicated multicollinearrity but everything is below

#######################
######## MODEL ########
#######################


pseudoR <- function(resdev, nulldev){
  return(1-round(resdev/nulldev, 2))
}

###model assumption checks###

#overdispertion test
m.poisson <- glm(NumberWaiting ~ Aggregated.speciality + Region + Season + ScaledYear +
                   ScaledPop + ScaledAvgAge + ScaledFMRatio,  
                 data=data.model.pop, family = poisson)
dispersiontest(m.poisson, trafo=1)

#fitted values(mu) vs variance
plot(log(fitted(m.poisson)),log((data.model.pop$NumberWaiting-fitted(m.poisson))^2),
     xlab=expression(hat(mu)),ylab="")
abline(0,1, col="red", lwd=2)
title(ylab=expression((y-hat(mu))^2), line=2)




#####################
### Quasi-poisson####
#####################


m.qp <- glm(NumberWaiting ~  Aggregated.speciality +  Region + ScaledAvgAge + Season +
              ScaledPop + ScaledFMRatio + ScaledSimd + ScaledYear,
            data=data.model.pop, family = quasipoisson)

#check the significance of variables with F test
drop1(m.qp,test="F")
#drop season 
m.qp.2 <- glm(NumberWaiting ~  Aggregated.speciality +  Region + ScaledAvgAge +
                ScaledPop + ScaledFMRatio + ScaledSimd + ScaledYear,
              data=data.model.pop, family = quasipoisson)
#check again
drop1(m.qp.2,test="F")

#summary of the final model
summary(m.qp.2)
pseudo.m.qp <- pseudoR(6485334,13624445)
# (Dispersion parameter for quasipoisson family taken to be 350.3138)
# Null deviance: 13624445  on 23529  degrees of freedom
# Residual deviance:  6485334  on 23515  degrees of freedom

#mean vs variance 
plot(log(fitted(m.qp.2)),log((data.model.pop$NumberWaiting-fitted(m.qp.2))^2),
     xlab=expression(hat(mu)),ylab=expression((y-hat(mu))^2))
abline(0,1, col="red")

#plot summary 

tbl_regression(m.gp.2, exponentiate = TRUE)

##############################################
####Poisson Regression with Random Effects####
##############################################

#did not converge
m.pois.re <- glmer(NumberWaiting ~ (1 | HBName) + (1|SpecialtyName)  + Season +
                     ScaledSimd + ScaledFMRatio + ScaledAvgAge + ScaledPop + ScaledYear, 
                   data = data.model.pop, family = poisson(link = "log"))


#did not converge
m.pois.re.2 <- glmer(NumberWaiting ~ (1 | HBName) + (1|SpecialtyName)  +
                       ScaledSimd + ScaledFMRatio + ScaledAvgAge + ScaledPop + ScaledYear, 
                     data = data.model.pop, family = poisson(link = "log"))

#ScaledPop is insignificant
m.pois.re.3 <- glmer(NumberWaiting ~ (1 | HBName) + (1|SpecialtyName)  +
                       ScaledSimd + ScaledFMRatio + ScaledAvgAge + ScaledPop, 
                     data = data.model.pop, family = poisson(link = "log"))

#final model
m.pois.re.4 <- glmer(NumberWaiting ~  (1 | HBName) + (1|SpecialtyName) +
                       ScaledSimd + ScaledFMRatio + ScaledAvgAge , 
                     data = data.model.pop, family = poisson(link = "log"))

###model diagnostics 

#summary 
summary(m.pois.re.4)

## QQ plot
plot(ranef(m.pois.re.4))

## Caterpillar plot
lattice::dotplot(ranef(m.pois.re.4, condVar = TRUE))






