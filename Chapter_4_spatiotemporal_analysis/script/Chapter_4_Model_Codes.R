#====================================================================================================
#  LOAD IMPORTANT LIBRARIES
#====================================================================================================

library(scales)
library(mgcv)
library(tidyr)
library(data.table)
library(grid)
library(animation)
library(MASS)
library(nlme)
library('plyr')
library(sjmisc)
library(rlang)
library(texreg)
library(lubridate)
library(cartography)
library(sf)
library(tidyverse)
library(anytime)
library(zoo)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(AICcmodavg)
library(bbmle)
library(MuMIn)
library(patchwork)
library(purrr)
library(plotly)
library(hrbrthemes)
library(ggpubr)
library(gridExtra)
library(grid)
library(lattice)
library(cowplot)
library(xlsx)
library(spdep)
library(mapproj)
library(plyr)
library(RColorBrewer)
library(INLA)
library(viridis)
library(ggthemes)
library(spdep)
library(foreign)
library(sp)
library(mapview)
library(car)
library(ggrepel)
require(tmap)
require(tidyverse)
library(directlabels)
library(pacman)
p_load(tidyverse, sf, tmap, INLA, car, ggplot2, readxl, spdep)
devtools::install_github("gfalbery/ggregplot")
library(spdep)
library(sf)
library(INLA)




################################################################################

###    M  O  D  E  L  S       

################################################################################

data <- read.csv("------------------/all_data_with_envi.csv")

data$PropFemales<- 1-data$PropMales
  
library(boot)    
data <- data %>% mutate(PovertyProp = logit(PovertyProp),
                        PropMales = logit(PropMales),
                        PropFemales = logit(PropFemales),
                        PropVacc = logit(PropVacc),
                        datetime = as.Date(Date, format = "%d/%m/%Y"), 
                        date_month = lubridate::month(datetime),
                        date_month_poverty = date_month,
                        date_month_pop = date_month) 


#==============================================================================
# GLM
#==============================================================================

fit <- glm(Cases ~ PovertyProp + PopDensity+Agegroup+ PropMales + PropVacc +
             offset(log(ExpectedCases)) + RH + temp + prep,
           data = data, family=poisson(link="log"))

saveRDS(fit, "fit_result.rds")
# fit <- readRDS("fit_result.rds")
summary(fit)



#==============================================================================
# Model I: Age, poverty, population density and precipitation  
#==============================================================================
formula <- Cases ~ Agegroup +  prep +
                   f(time3ID, PovertyProp, model = "rw1" ) +
                   f(time4ID, PopDensity, model = "rw1" )
res.S1 <- inla(formula,family = "poisson", data = data, E = ExpectedCases,
             control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))

saveRDS(res.S1, "res.S1_result.rds")
# res.S1 <- readRDS("res.S1_result.rds")
summary(res.S1) 
 



#Create adjacency matrix
shp <- st_read("C:/Users/mwand/OneDrive - The University of Manchester/aaPhD Manchetster/Project2/PHIM_COVID-19 Data_SpatialAnalysis/Risk Maps/Manuscript2Materials/gadm41_MWI_shp/gadm41_MWI_1.shp")
mal2 <- shp[-10, ]  # remove Likoma district
mal2
library(spdep)

nb <- poly2nb(mal2)
nb2INLA("map.adj", nb)
library(INLA)

Adj.Mat <- inla.read.graph(filename = "map.adj")
 


#==============================================================================
# Model II: Age, poverty, population density, precipitation and convolution 
#==============================================================================
formula <- Cases ~ Agegroup +
                  prep +
                f(time3ID, PovertyProp, model = "rw1" ) +
                f(time4ID, PopDensity, model = "rw1" ) +
                f(areaID, model = "besag", graph=Adj.Mat) +f(area1ID, model = "iid")
res.ST0 <- inla(formula,family = "poisson", data = data, E = ExpectedCases,
            control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE)) #control.predictor = list(compute = TRUE),
saveRDS(res.ST0, "resST0_result.rds")
res.ST0 <- readRDS("resST0_result.rds")
summary(res.ST0)  
 


#==============================================================================
# Model III: Age, poverty, population density, precipitation, convolution and general time trend                                          
#==============================================================================
formula <- Cases ~ Agegroup+timeID +
                   prep +
                f(time3ID, PovertyProp, model = "rw1" ) +
                f(time4ID, PopDensity, model = "rw1" ) +
                f(areaID, model = "besag", graph=Adj.Mat) + f(area1ID, model = "iid")
res.ST1 <- inla(formula,family = "poisson", data = data, E = ExpectedCases,
            control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))

saveRDS(res.ST1, "resST1_result.rds")
#res.ST1 <- readRDS("resST1_result.rds")
summary(res.ST1)
 


#==============================================================================
# Model IV: Age, poverty, population density, precipitation, convolution and random walk 1  
#==============================================================================
set.seed(120)
formula <- Cases ~ Agegroup+ #PropMales + PropVacc +
                    prep +
                  f(time3ID, PovertyProp, model = "rw1" ) +
                  f(time4ID, PopDensity, model = "rw1" ) +
                  f(areaID, model = "besag", graph=Adj.Mat) + f(area1ID, model = "iid") +
                  f(timeID, model = "rw1") +f(time1ID, model = "iid")
res.ST2 <- inla(formula,family = "poisson", data = data, E = ExpectedCases,
            control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))
saveRDS(res.ST2, "resST2_result.rds")
res.ST2 <- readRDS("resST2_result.rds")
summary(res.ST2)
 


#==============================================================================
# Model V: Age, poverty, population density, precipitation, convolution and random walk 2  
#==============================================================================
set.seed(120)
formula <- Cases ~ Agegroup+ #PropMales + PropVacc +
                    prep +
                f(time3ID, PovertyProp, model = "rw2" ) +
                f(time4ID, PopDensity, model = "rw2" ) +
                f(areaID, model = "besag", graph=Adj.Mat ) + f(area1ID, model = "iid") +
                f(timeID, model = "rw2" ) + f(time1ID, model = "iid")
res.ST3 <- inla(formula,family = "poisson", data = data, E = ExpectedCases,
            control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))

 saveRDS(res.ST3, "resST3_result.rds")
#res.ST3 <- readRDS("resST3_result.rds")
summary(res.ST3) 
 



#============================================================================== 
# Model VI: Age, poverty, population density, precipitation, convolution, random walk 1 and ST-interactionType 1  
#==============================================================================
set.seed(120)
formula <- Cases ~ Agegroup + #PropMales + PropVacc +
                      prep +
                    f(time3ID, PovertyProp, model = "rw1" ) +
                    f(time4ID, PopDensity, model = "rw1" ) +
                    f(areaID, model = "besag", graph=Adj.Mat ) + f(area1ID, model = "iid") +
                    f(timeID, model = "rw1" ) + f(time1ID, model = "iid") +
                    f(STinteractionIndex, model = "iid" )
res.ST4 <- inla(formula,family = "poisson", data = data, E = ExpectedCases,
             control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))

saveRDS(res.ST4, "resST4_result.rds")
#res.ST4 <- readRDS("resST4_result.rds")
summary(res.ST4)
 




#==============================================================================
# Model VII: Age, poverty, population density, precipitation, convolution, random walk 1 and ST-interactionType 2
#==============================================================================
set.seed(120)
formula <- Cases ~ Agegroup + #PropMales + PropVacc +
                    prep +
                f(time3ID, PovertyProp, model = "rw1" ) +
                f(time4ID, PopDensity, model = "rw1" ) +
              f(areaID, model = "besag", graph=Adj.Mat ) + f(area1ID, model = "iid" ) +
              f(timeID, model = "rw1" ) + f(time1ID, model = "iid" ) +
              f(area2ID, model = "iid" , group=time2ID, control.group= list(model = "rw1"  )    )
res.ST5 <- inla(formula,family = "poisson", data = data, E = ExpectedCases,
             control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))

saveRDS(res.ST5, "resST5_result.rds")
#res.ST5 <- readRDS("resST5_result.rds")
summary(res.ST5)





#==============================================================================
# Model VIII: Age, poverty, population density, precipitation, convolution, random walk 1 and ST-interactionType 3
#==============================================================================
set.seed(120)
formula <- Cases ~ Agegroup +
                    prep +
                f(time3ID, PovertyProp, model = "rw1" ) +
                f(time4ID, PopDensity, model = "rw1" ) +
            f(areaID, model = "besag", graph=Adj.Mat) + f(area1ID, model = "iid" ) +
            f(timeID, model = "rw2" ) + f(time1ID, model = "iid" ) +
            f(time2ID, model = "iid", group=area2ID, control.group= list(model = "besag", graph=Adj.Mat))
res.ST6 <- inla(formula,family = "poisson", data = data, E = ExpectedCases,
             control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))

saveRDS(res.ST6, "resST6_result.rds")
#res.ST6 <- readRDS("resST6_result.rds")
summary(res.ST6) 




#============================================================================== 
# Model IX: Age, poverty, population density, precipitation, convolution, random walk 1 and ST-interactionType 4 
#==============================================================================
set.seed(120)
 formula <- Cases ~ Agegroup +
                       prep +
                   f(time3ID, PovertyProp, model = "rw1" ) +
                   f(time4ID, PopDensity, model = "rw1" ) +
                 f(areaID, model = "bym2", graph=Adj.Mat ) +
                 f(timeID, model = "rw1" )  +  f(time1ID, model = "iid" ) +
                 f(area2ID, model="besag", graph=Adj.Mat, group=time2ID, control.group= list(model = "rw1"),  constr = TRUE)#+

res.ST7 <- inla(formula,family = "poisson", data = data, E = ExpectedCases,
            control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE, return.marginals.predictor=TRUE),
            control.predictor = list(compute = TRUE))

saveRDS(res.ST7, "resST7_result.rds")
#res.ST7 <- readRDS("resST7_result.rds")
summary(res.ST7)
 




#==============================================================================
#Model X: Age, poverty, population density, precipitation, convolution, random walk 1, ST-interactionType 4 and STA-interaction
#==============================================================================
set.seed(120)
formula <- Cases ~ Agegroup + #PropMales + PropVacc +
                   prep +
                f(time3ID, PovertyProp, model = "rw1" ) +
                f(time4ID, PopDensity, model = "rw1" ) +
              f(areaID, model = "bym2", graph=Adj.Mat ) +
              f(timeID, model = "rw1" )  +  f(time1ID, model = "iid" ) +
              f(area2ID, model="besag", graph=Adj.Mat, group=time2ID, control.group= list(model = "rw1"),  constr = TRUE)+
              f(SpaceTimeAgeIndex, model="iid" )
res.STA0 <- inla(formula,family = "poisson", data = data, E = ExpectedCases,
              control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE, return.marginals.predictor=TRUE),
              control.predictor = list(compute = TRUE))
saveRDS(res.STA0, "resSTA0_result.rds")
#res.STA0 <- readRDS("resSTA0_result.rds")
summary(res.STA0)






















###############################################################################

# MAIN FIGURES

###############################################################################




RR <- res.STA0$summary.fitted.values[, "mean"] 
min(RR)
max(RR)
mean(RR)



#===============================================================================
# Temporal variation of Poverty   
#===============================================================================
pv<-res.STA0$summary.random$time3ID
pov<-data.frame(pv)
pov$date<-data$Date[1:105]
pov<-pov%>% mutate(date =dmy(date))

ggplot(pov, aes(x = date, y = mean)) + 
  geom_line() +
  geom_ribbon(aes(ymin = X0.025quant, ymax = X0.975quant), alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Red line at y = 0
  ggtitle("Variation of poverty over time") +
  xlab("Date") +
  ylab("Effect of poverty")



#================================================================================
# Temporal variation of Population density   
#================================================================================

pp<-res.STA0$summary.random$time4ID
pop<-data.frame(pp)
pop$date<-data$Date[1:105]
pop<-pop%>% mutate(date =dmy(date))
ggplot(pop, aes(x = date, y = mean)) + 
  geom_line() +
  geom_ribbon(aes(ymin = X0.025quant, ymax = X0.975quant), alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Red line at y = 0
  ggtitle("Variation of population density over time") +
  xlab("Date") +
  ylab("Effect of population density")



#================================================================================
# Maps of theta 
#================================================================================

RR <- res.STA0$summary.fitted.values[, "mean"] 
data.tj <- data.frame(RR = RR, Week = data$timeID, 
                      area = data$areaID, agegroup =data$Agegroup)
dat<-data.tj

 
lake_malawi_shapefile<- st_read("C:/Users/mwand/OneDrive - The University of Manchester/aaPhD Manchetster/Project2/PHIM_COVID-19 Data_SpatialAnalysis/Risk Maps/Manuscript2Materials/gadm41_MWI_shp/Lake Malawi.shp")

 
if (st_crs(mal2) != st_crs(lake_malawi_shapefile)) {
  lake_malawi_shapefile <- st_transform(lake_malawi_shapefile, st_crs(malawi_shapefile))
}

mal2$area <- 1:27
shp2 <- mal2 %>% left_join(dat, by="area")


shp3 <- shp2 %>% filter(agegroup =="0-19", Week %in% c(5,15,25,36,42,47,61,69,77,86,91,98))   
tm_shape(shp3) + 
  tm_polygons("RR", title= "Relative \nRisk") +
  tm_facets("Week", free.coords=FALSE, ncol=4) +
  tm_scale_bar(position = c("LEFT", "BOTTOM")) +
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1, panel.label.bg.color="white",
            panel.labels=paste0(c("26 April - 2 May 2020", "5 - 11 Jul. 2020", "13 - 19 Sept. 2020", "29 Nov. - 5 Dec. 2020", "10 - 16 Jan. 2021", "14 - 20 Feb. 2021",
                                  "23 - 29 May 2021", "18 - 24 Jul. 2021", "12 - 18 Sept. 2021", "14 - 20 Nov. 2021", "19 - 25 Dec. 2021", "6 - 12 Feb. 2022")))+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7)
 

  
shp3 <- shp2 %>% filter(agegroup =="20-29", Week %in% c(5,15,25,36,42,47,61,69,77,86,91,98))   
tm_shape(shp3) + 
  tm_polygons("RR", title= "Relative \nRisk") +
  tm_facets("Week", free.coords=FALSE, ncol=4) +
  tm_scale_bar(position = c("LEFT", "BOTTOM")) +
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1, panel.label.bg.color="white",
            panel.labels=paste0(c("26 April - 2 May 2020", "5 - 11 Jul. 2020", "13 - 19 Sept. 2020", "29 Nov. - 5 Dec. 2020", "10 - 16 Jan. 2021", "14 - 20 Feb. 2021",
                                  "23 - 29 May 2021", "18 - 24 Jul. 2021", "12 - 18 Sept. 2021", "14 - 20 Nov. 2021", "19 - 25 Dec. 2021", "6 - 12 Feb. 2022")))+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7) 



shp3 <- shp2 %>% filter(agegroup =="30-39", Week %in% c(5,15,25,36,42,47,61,69,77,86,91,98))     
tm_shape(shp3) + 
  tm_polygons("RR", title= "Relative \nRisk") +
  tm_facets("Week", free.coords=FALSE, ncol=4) +
  tm_scale_bar(position = c("LEFT", "BOTTOM")) +
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1, panel.label.bg.color="white",
            panel.labels=paste0(c("26 April - 2 May 2020", "5 - 11 Jul. 2020", "13 - 19 Sept. 2020", "29 Nov. - 5 Dec. 2020", "10 - 16 Jan. 2021", "14 - 20 Feb. 2021",
                                  "23 - 29 May 2021", "18 - 24 Jul. 2021", "12 - 18 Sept. 2021", "14 - 20 Nov. 2021", "19 - 25 Dec. 2021", "6 - 12 Feb. 2022")))+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7) 



shp3 <- shp2 %>% filter(agegroup =="40-49", Week %in% c(5,15,25,36,42,47,61,69,77,86,91,98))     
tm_shape(shp3) + 
  tm_polygons("RR", title= "Relative \nRisk") +
  tm_facets("Week", free.coords=FALSE, ncol=4) +
  tm_scale_bar(position = c("LEFT", "BOTTOM")) +
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1, panel.label.bg.color="white",
            panel.labels=paste0(c("26 April - 2 May 2020", "5 - 11 Jul. 2020", "13 - 19 Sept. 2020", "29 Nov. - 5 Dec. 2020", "10 - 16 Jan. 2021", "14 - 20 Feb. 2021",
                                  "23 - 29 May 2021", "18 - 24 Jul. 2021", "12 - 18 Sept. 2021", "14 - 20 Nov. 2021", "19 - 25 Dec. 2021", "6 - 12 Feb. 2022")))+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7) 



shp3 <- shp2 %>% filter(agegroup =="50+", Week %in% c(5,15,25,36,42,47,61,69,77,86,91,98))     
tm_shape(shp3) + 
  tm_polygons("RR", title= "Relative \nRisk") +
  tm_facets("Week", free.coords=FALSE, ncol=4) +
  tm_scale_bar(position = c("LEFT", "BOTTOM")) +
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1, panel.label.bg.color="white",
            panel.labels=paste0(c("26 April - 2 May 2020", "5 - 11 Jul. 2020", "13 - 19 Sept. 2020", "29 Nov. - 5 Dec. 2020", "10 - 16 Jan. 2021", "14 - 20 Feb. 2021",
                                  "23 - 29 May 2021", "18 - 24 Jul. 2021", "12 - 18 Sept. 2021", "14 - 20 Nov. 2021", "19 - 25 Dec. 2021", "6 - 12 Feb. 2022")))+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7) 



#================================================================================
# Risk in week between 19 and 25 Dec for all age groups     
#================================================================================

shpDec <- shp2 %>% filter(Week ==91, agegroup %in% c("0-19", "20-29", "30-39", "40-49", "50+"))     

tm_shape(shpDec) + 
  tm_polygons("RR", title= "Relative \nRisk") +
  tm_facets("agegroup", free.coords=FALSE, ncol=3) +
  tm_scale_bar(position = c("LEFT", "BOTTOM")) +
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1, panel.label.bg.color="white",
            panel.labels=paste0(c("Age 0-19", "Age 20-29", "Age 30-39", "Age 40-49", "Age 50 and above")))+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7) 

head(shp2)



#================================================================================
# Exceedance probability     
#================================================================================

# first get the thresholds in selected weeks
RR <- res.STA0$summary.fitted.values[, "mean"]
plot.week <- c(5,15,25,36,42,47,61,69,77,86,91,98)
data.tj <- data.frame(RR = RR, Week = data$timeID,
                      area = data$areaID, agegroup =data$Agegroup, ID = 1:14175)
# threshold
data.tj %>% filter(Week %in% plot.week, agegroup == "40-49") %>%
  group_by(Week) %>% summarise(thres = mean(RR))



#=========================
# 1. Threshold in week 5
#=========================
IDD <- data.tj  %>% filter(Week == 5, agegroup == "40-49") %>% pull(ID)
exc <- sapply(lapply(IDD[c(1:8, 10:27)], function(x)res.STA0$marginals.fitted.values[[x]]),
              FUN = function(marg){
                1- inla.pmarginal(q=  0.00111, marginal = marg)
              })
exc1 <- append(exc, 0, after = 8)#c(exc[1:8], 0, exc[10:26])
area=c(1:27)
shp <- mal2
lake_malawi_shapefile<- st_read("C:/Users/mwand/OneDrive - The University of Manchester/aaPhD Manchetster/Project2/PHIM_COVID-19 Data_SpatialAnalysis/Risk Maps/Manuscript2Materials/gadm41_MWI_shp/Lake Malawi.shp")
Exceedance1<-data.frame(exc1, area)
shp$area <- 1:27
shp2 <- shp %>% left_join(Exceedance1, by="area")
shp2$label_name <- ifelse(shp2$area %in% c(2, 10, 16), shp2$NAME_1, NA)


create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }

  return(labs)
}

brks <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0)

labs <- create_labels(brks, greater = F)

tm_shape(shp2) +
  tm_polygons("exc1", title= "Exceedance \nProbability \n(Threshold=0.0011)", labels = labs, breaks = brks, style = "fixed") +
  tm_text("label_name", size = 0.7, fontface = "bold", col = "black",
          bg.color = "white", bg.alpha = 0.5)+

  tm_scale_bar(position = c("LEFT", "BOTTOM"),
               text.size = 1.0, width = 0.3) +
  tm_legend(position = c(c(0.7, 0.15)))+
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1.6, panel.label.bg.color="white",title.size = 2.5,
            legend.title.size = 1.5, legend.title.fontface = "bold",
            legend.text.size = 0.9, panel.label.size = 1.3,
            panel.labels=paste0("26 April - 2 May 2020" )  )+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7)




#=========================
# 2. Threshold in week 15
#=========================
IDD <- data.tj  %>% filter(Week == 15, agegroup == "40-49") %>% pull(ID)
exc <- sapply(lapply(IDD[c(1:8, 10:27)], function(x)res.STA0$marginals.fitted.values[[x]]),
              FUN = function(marg){
                1- inla.pmarginal(q=  0.0256, marginal = marg)
              })
exc1 <- append(exc, 0, after = 8)
area=c(1:27)
shp <- mal2
Exceedance1<-data.frame(exc1, area)
shp$area <- 1:27
shp2 <- shp %>% left_join(Exceedance1, by="area")

create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }

  return(labs)
}

brks <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0)

labs <- create_labels(brks, greater = F)

tm_shape(shp2) +
  tm_polygons("exc1", title= "Exceedance \nProbability \n(Threshold=0.0256)",  labels = labs, breaks = brks, style = "fixed") +
  tm_scale_bar(position = c("LEFT", "BOTTOM"),
               text.size = 1.0, width = 0.3) +
  tm_legend(position = c(c(0.7, 0.15)))+
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1.6, panel.label.bg.color="white",title.size = 2.5,
            legend.title.size = 1.5, legend.title.fontface = "bold",
            legend.text.size = 0.9, panel.label.size = 1.3,
            panel.labels=paste0("5 - 11 Jul. 2020" )  )+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7)



#=========================
# 3. Threshold in week 25
#=========================
IDD <- data.tj  %>% filter(Week == 25, agegroup == "40-49") %>% pull(ID)
exc <- sapply(lapply(IDD[c(1:8, 10:27)], function(x)res.STA0$marginals.fitted.values[[x]]),
              FUN = function(marg){
                1- inla.pmarginal(q=  0.000853, marginal = marg)
              })
exc1 <- append(exc, 0, after = 8)
area=c(1:27)
shp <- mal2
Exceedance1<-data.frame(exc1, area)
shp$area <- 1:27
shp2 <- shp %>% left_join(Exceedance1, by="area")

create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }

  return(labs)
}

brks <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0)

labs <- create_labels(brks, greater = F)

tm_shape(shp2) +
  tm_polygons("exc1", title= "Exceedance \nProbability \n(Threshold=0.0009)",   labels = labs, breaks = brks, style = "fixed") +
  tm_scale_bar(position = c("LEFT", "BOTTOM"),
               text.size = 1.0, width = 0.3) +
  tm_legend(position = c(c(0.7, 0.15)))+
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1.6, panel.label.bg.color="white",title.size = 2.5,
            legend.title.size = 1.5, legend.title.fontface = "bold",
            legend.text.size = 0.9, panel.label.size = 1.3,
            panel.labels=paste0("13 - 19 Sept. 2020" )  )+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7)



#=========================
# 4. Threshold in week 36
#=========================
IDD <- data.tj  %>% filter(Week == 36, agegroup == "40-49") %>% pull(ID)
exc <- sapply(lapply(IDD[c(1:8, 10:27)], function(x)res.STA0$marginals.fitted.values[[x]]),
              FUN = function(marg){
                1- inla.pmarginal(q=  0.000937, marginal = marg)
              })
exc1 <- append(exc, 0, after = 8)
area=c(1:27)
shp <- mal2
Exceedance1<-data.frame(exc1, area)
shp$area <- 1:27
shp2 <- shp %>% left_join(Exceedance1, by="area")

create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }

  return(labs)
}

brks <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0)

labs <- create_labels(brks, greater = F)

tm_shape(shp2) +
  tm_polygons("exc1", title= "Exceedance \nProbability \n(Threshold=0.0009)",   labels = labs, breaks = brks, style = "fixed") +
  tm_scale_bar(position = c("LEFT", "BOTTOM"),
               text.size = 1.0, width = 0.3) +
  tm_legend(position = c(c(0.7, 0.15)))+
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1.6, panel.label.bg.color="white",title.size = 2.5,
            legend.title.size = 1.5, legend.title.fontface = "bold",
            legend.text.size = 0.9, panel.label.size = 1.3,
            panel.labels=paste0("29 Nov. - 5 Dec. 2020" )  )+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7)



#=========================
# 5. Threshold in week 42
#=========================
IDD <- data.tj  %>% filter(Week == 42, agegroup == "40-49") %>% pull(ID)

exc <- sapply(lapply(IDD[c(1:8, 10:27)], function(x)res.STA0$marginals.fitted.values[[x]]),
              FUN = function(marg){
                1- inla.pmarginal(q=  0.0941 , marginal = marg)
              })
exc1 <- append(exc, 0, after = 8)
area=c(1:27)
shp <- mal2
Exceedance1<-data.frame(exc1, area)
shp$area <- 1:27
shp2 <- shp %>% left_join(Exceedance1, by="area")

create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }

  return(labs)
}

brks <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0)

labs <- create_labels(brks, greater = F)

tm_shape(shp2) +
  tm_polygons("exc1", title= "Exceedance \nProbability \n(Threshold=0.0941)",  labels = labs, breaks = brks, style = "fixed") +
  tm_scale_bar(position = c("LEFT", "BOTTOM"),
               text.size = 1.0, width = 0.3) +
  tm_legend(position = c(c(0.7, 0.15)))+
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1.6, panel.label.bg.color="white",title.size = 2.5,
            legend.title.size = 1.5, legend.title.fontface = "bold",
            legend.text.size = 0.9, panel.label.size = 1.3,
            panel.labels=paste0("10 - 16 Jan. 2021" )  )+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7)



#=========================
# 6. Threshold in week 47
#=========================
IDD <- data.tj  %>% filter(Week == 47, agegroup == "40-49") %>% pull(ID)
exc <- sapply(lapply(IDD[c(1:8, 10:27)], function(x)res.STA0$marginals.fitted.values[[x]]),
              FUN = function(marg){
                1- inla.pmarginal(q=  0.0309, marginal = marg)
              })
exc1 <- append(exc, 0, after = 8)
area=c(1:27)
shp <- mal2
Exceedance1<-data.frame(exc1, area)
shp$area <- 1:27
shp2 <- shp %>% left_join(Exceedance1, by="area")

create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }

  return(labs)
}

brks <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0)

labs <- create_labels(brks, greater = F)

tm_shape(shp2) +
  tm_polygons("exc1", title= "Exceedance \nProbability \n(Threshold=0.0309)",   labels = labs, breaks = brks, style = "fixed") +
  tm_scale_bar(position = c("LEFT", "BOTTOM"),
               text.size = 1.0, width = 0.3) +
  tm_legend(position = c(c(0.7, 0.15)))+
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1.6, panel.label.bg.color="white",title.size = 2.5,
            legend.title.size = 1.5, legend.title.fontface = "bold",
            legend.text.size = 0.9, panel.label.size = 1.3,
            panel.labels=paste0("14 - 20 Feb. 2021" )  )+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7)



#=========================
# 7. Threshold in week 61
#=========================
IDD <- data.tj  %>% filter(Week == 61, agegroup == "40-49") %>% pull(ID)
exc <- sapply(lapply(IDD[c(1:8, 10:27)], function(x)res.STA0$marginals.fitted.values[[x]]),
              FUN = function(marg){
                1- inla.pmarginal(q=  0.00194, marginal = marg)
              })
exc1 <- append(exc, 0, after = 8)
area=c(1:27)
shp <- mal2
Exceedance1<-data.frame(exc1, area)
shp$area <- 1:27
shp2 <- shp %>% left_join(Exceedance1, by="area")

create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }

  return(labs)
}

brks <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0)

labs <- create_labels(brks, greater = F)

tm_shape(shp2) +
  tm_polygons("exc1", title= "Exceedance \nProbability \n(Threshold=0.0019)",  labels = labs, breaks = brks, style = "fixed") +
  tm_scale_bar(position = c("LEFT", "BOTTOM"),
               text.size = 1.0, width = 0.3) +
  tm_legend(position = c(c(0.7, 0.15)))+
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1.6, panel.label.bg.color="white",title.size = 2.5,
            legend.title.size = 1.5, legend.title.fontface = "bold",
            legend.text.size = 0.9, panel.label.size = 1.3,
            panel.labels=paste0("23 - 29 May 2021" )  )+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7)



#=========================
# 8. Threshold in week 69
#=========================
IDD <- data.tj  %>% filter(Week == 69, agegroup == "40-49") %>% pull(ID)
exc <- sapply(lapply(IDD[c(1:8, 10:27)], function(x)res.STA0$marginals.fitted.values[[x]]),
              FUN = function(marg){
                1- inla.pmarginal(q=  0.137, marginal = marg)
              })
exc1 <- append(exc, 0, after = 8)
area=c(1:27)
shp <- mal2
Exceedance1<-data.frame(exc1, area)
shp$area <- 1:27
shp2 <- shp %>% left_join(Exceedance1, by="area")

create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }

  return(labs)
}

brks <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0)

labs <- create_labels(brks, greater = F)

tm_shape(shp2) +
  tm_polygons("exc1", title= "Exceedance \nProbability \n(Threshold=0.1370)",  labels = labs, breaks = brks, style = "fixed") +
  tm_scale_bar(position = c("LEFT", "BOTTOM"),
               text.size = 1.0, width = 0.3) +
  tm_legend(position = c(c(0.7, 0.15)))+
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1.6, panel.label.bg.color="white",title.size = 2.5,
            legend.title.size = 1.5, legend.title.fontface = "bold",
            legend.text.size = 0.9, panel.label.size = 1.3,
            panel.labels=paste0("18 - 24 Jul. 2021" )  )+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7)


 

#=========================
# 9. Threshold in week 77
#=========================
IDD <- data.tj  %>% filter(Week == 77, agegroup == "40-49") %>% pull(ID)
exc <- sapply(lapply(IDD[c(1:8, 10:27)], function(x)res.STA0$marginals.fitted.values[[x]]),
              FUN = function(marg){
                1- inla.pmarginal(q=  0.00927, marginal = marg)
              })
exc1 <- append(exc, 0, after = 8)
area=c(1:27)
shp <- mal2
Exceedance1<-data.frame(exc1, area)
shp$area <- 1:27
shp2 <- shp %>% left_join(Exceedance1, by="area")

create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }

  return(labs)
}

brks <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0)

labs <- create_labels(brks, greater = F)

tm_shape(shp2) +
  tm_polygons("exc1", title= "Exceedance \nProbability \n(Threshold=0.0093)", labels = labs, breaks = brks, style = "fixed") +
  tm_scale_bar(position = c("LEFT", "BOTTOM"),
               text.size = 1.0, width = 0.3) +
  tm_legend(position = c(c(0.7, 0.15)))+
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1.6, panel.label.bg.color="white",title.size = 2.5,
            legend.title.size = 1.5, legend.title.fontface = "bold",
            legend.text.size = 0.9, panel.label.size = 1.3,
            panel.labels=paste0("12 - 18 Sept. 2021" )  )+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7)



#=========================
# 10. Threshold in week 86
#=========================
IDD <- data.tj  %>% filter(Week == 86, agegroup == "40-49") %>% pull(ID)
exc <- sapply(lapply(IDD[c(1:8, 10:27)], function(x)res.STA0$marginals.fitted.values[[x]]),
              FUN = function(marg){
                1- inla.pmarginal(q=  0.000900, marginal = marg)
              })
exc1 <- append(exc, 0, after = 8)
area=c(1:27)
shp <- mal2
Exceedance1<-data.frame(exc1, area)
shp$area <- 1:27
shp2 <- shp %>% left_join(Exceedance1, by="area")

create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }

  return(labs)
}

brks <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0)

labs <- create_labels(brks, greater = F)

tm_shape(shp2) +
  tm_polygons("exc1", title= "Exceedance \nProbability \n(Threshold=0.0.0009)", labels = labs, breaks = brks, style = "fixed") +
  tm_scale_bar(position = c("LEFT", "BOTTOM"),
               text.size = 1.0, width = 0.3) +
  tm_legend(position = c(c(0.7, 0.15)))+
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1.6, panel.label.bg.color="white",title.size = 2.5,
            legend.title.size = 1.5, legend.title.fontface = "bold",
            legend.text.size = 0.9, panel.label.size = 1.3,
            panel.labels=paste0("14 - 20 Nov. 2021" )  )+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7)



#=========================
# 11. Threshold in week 91
#=========================
IDD <- data.tj  %>% filter(Week == 91, agegroup == "40-49") %>% pull(ID)
exc <- sapply(lapply(IDD[c(1:8, 10:27)], function(x)res.STA0$marginals.fitted.values[[x]]),
              FUN = function(marg){
                1- inla.pmarginal(q=  0.118, marginal = marg)
              })
exc1 <- append(exc, 0, after = 8)
area=c(1:27)
shp <- mal2
Exceedance1<-data.frame(exc1, area)
shp$area <- 1:27
shp2 <- shp %>% left_join(Exceedance1, by="area")

create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }

  return(labs)
}

brks <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0)

labs <- create_labels(brks, greater = F)

tm_shape(shp2) +
  tm_polygons("exc1", title= "Exceedance \nProbability \n(Threshold=0.1180)",labels = labs, breaks = brks, style = "fixed") +
  tm_scale_bar(position = c("LEFT", "BOTTOM"),
               text.size = 1.0, width = 0.3) +
  tm_legend(position = c(c(0.7, 0.15)))+
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1.6, panel.label.bg.color="white",title.size = 2.5,
            legend.title.size = 1.5, legend.title.fontface = "bold",
            legend.text.size = 0.9, panel.label.size = 1.3,
            panel.labels=paste0("19 - 25 Dec. 2021" )  )+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7)



#=========================
# 12. Threshold in week 98
#=========================
IDD <- data.tj  %>% filter(Week == 98, agegroup == "40-49") %>% pull(ID)
exc <- sapply(lapply(IDD[c(1:8, 10:27)], function(x)res.STA0$marginals.fitted.values[[x]]),
              FUN = function(marg){
                1- inla.pmarginal(q=  0.00459, marginal = marg)
              })
exc1 <- append(exc, 0, after = 8)
area=c(1:27)
shp <- mal2
Exceedance1<-data.frame(exc1, area)
shp$area <- 1:27
shp2 <- shp %>% left_join(Exceedance1, by="area")

create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }

  return(labs)
}

brks <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0)

labs <- create_labels(brks, greater = F)

tm_shape(shp2) +
  tm_polygons("exc1", title= "Exceedance \nProbability \n(Threshold=0.0046)", labels = labs, breaks = brks, style = "fixed") +
  tm_scale_bar(position = c("LEFT", "BOTTOM"),
               text.size = 1.0, width = 0.3) +
  tm_legend(position = c(c(0.7, 0.15)))+
  tm_compass(position = c("RIGHT", "TOP"), size = 0.5) +
  tm_layout(scale = 1, asp = 1.6, panel.label.bg.color="white",title.size = 2.5,
            legend.title.size = 1.5, legend.title.fontface = "bold",
            legend.text.size = 0.9, panel.label.size = 1.3,
            panel.labels=paste0("6 - 12 Feb. 2022" )  )+
  tm_shape(lake_malawi_shapefile) +
  tm_polygons(col = "lightblue",alpha = 0.7)










#






