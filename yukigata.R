###
#
#   Legacy effect of snow cover on alpine landscapes
#   Insights from the Lautaret area
#   First Version 20191021 
#   Last update 20241114
#   Ph. Choler
#
###

# SHIFT+ALT+O to unfold all
# ALT+O to fold all workbench (use capital letter)

# RATIONALE ----
  # a. Objective of this study
# flat areas where vegetation cover is sparse and soil has low SOM correspond to permanently snow-covered sites during the little ice age
# legacy effect of long snow cover duration on current biodiversity distribution
# a follow-up to Choler & al. Waning snowfields have transformed into greening hotspots (in revision). Nature Climate Change

  # b. Challenges in data analysis
# 1. assemble 30m resolution maps of SMOD using Sentinel-2, SPOT4-Take5, Landsat for the period 1984-2023
# 2. implement a phenomenological model of fractional Snow Cover Area driven by degree days
# 3. locate areas that were permanently snow-covered during the little ice age using historical met records

  # c. Other ideas
    # (i) article sur le dispositif suivi habitat et CBNA data - capacité de l'imagerie a capturer les niches de nivation et SDMs

    # (ii) legacy effect resulting from long snow covered sites avec diagnostic des anomalies de NDVImax à très haute résolution 

  # (iii) composition floristique des emerging snowbeds (nvlle niche thermique) over the last 4 decades using -> secteur Thabor only ?

# LOAD LIBRARIES ----
.libPaths()

library(mcp)          # breakpoint detection
library(rjags)        # breakpoint detection - need to install jags before !
library(patchwork)    # plot breakpoint detection
library(signal)

library(terra)
library(sf)           # for st_graticule
library(sp)
 # library(maptools)  " no longer on Cran 20240718
library(maptree)
library(PBSmapping)
library(ade4)
library(maps)
library(fields)
library(smatr)

library(RNetCDF)
library(gplots)
library(ggplot2)

source("C:/Users/cholerp/Documents/Rissue.r")

library(coda)
# library(rgeos)        # will be retired during 2023 -> sf functions
library(mblm)
library(RgoogleMaps)
library(gtools)
library(ggpubr)
library(cowplot)
library(palmerpenguins)
library(ggdensity)
library(gdata)
library(httr)
library(scales)       # alpha
# library(hypervolume)  # hypervolume - ISSUE (20250118)
library(multcompView)   # Tukey post-hoc
library(spatialEco)

library(reshape)        # cast & melt
library(vegan)          # vegdist

source("C:/Users/cholerp/Documents/Rissue.r")

library(lubridate)    # yday and year
library(zyp)          # zip.sen
library(plotrix)      # draw.circle
require(snowfall)     # for parallelization
library(Kendall)      # for MannKendall test
library(trend)        # for MannKendall test with continuity

library(grid)         # to combine maps using the function viewport()
library(tabularaster) # extract data from multiple polygons
library(ggformula)    # pour ??
library(corrplot)

library(Hmisc)        # smedian.hollow
library(tmap)

library(tidyverse) 
library(lqr)          # Logistic Linear Quantile Regression
library(data.table)   # rollmean
library(corrplot)
library(quantreg)     # quantile regression
library(rworldmap)
library(rworldxtra)
library(spData)       # example datasets
library(extrafont)
library(SingleCaseES) # Log Response Ratio (LRRd)
library(metafor)      # meta-analysis of LRR
library(RColorBrewer)
graphics.off();windows(15,15);display.brewer.all()


font_import()
yloadfonts(device = "win")
windowsFonts()

library(ggh4x)

library(MASS)
library(zoo)
library(RStoolbox)

show_line_types()

library(ggplot2)
library(nlraa)          # bootstrap & CI after nls
library(nlme)
library(mgcv)
library(nlstools)       # evaluate nls models

library(smoother)       # smth.gaussian for time series

library(plyr)
library(dplyr)
library(tidyr)
library(tidyterra)      # hypso.colors2
library(quantreg)
library(Metrics)

library(qmap)           # quantile mapping

# library(tidyquant)    # SOURCE OF THE ISSUE geom_ma moving average

source("C:/Users/cholerp/Documents/Rissue.r") # NO MORE ISSUE


# START (20240718) ----
  rm(list=ls())
  WD <- "~/MANUSCRITSmc/ONGOING_FIRSTAUTHOR/MFP3_Choler&al_Yukigata/DATA_ANALYSIS"
  setwd(WD); load("YUKIGATA.RData")    
  WD <- "~/MANUSCRITSmc/ONGOING_FIRSTAUTHOR/MFP3_Choler&al_Yukigata/DATA_ANALYSIS"
  DATA_WD <- "D:/DATA_YUKIGATA"
  
  source('~/SCRIPTSmc/palette.r')
  
  # load functions
  sigmoid = function(params, x) {
    params[1] / (1 + exp(-params[2] * (x - params[3])))
  }
  
  invsig  = function(params,x){
    -(log(params[1]/x -1)-params[3]*params[2])/params[2]
  }
  
  # log logistic with three parameters
  LL3 <- function(params,x){
    params[1]/(1+exp(params[2]*(log(x)-log(params[3]))))
  }
  
  # model evaluation
  m.eval  = function(meas,pred){
    # A. Function m.eval MODEL evaluation
    # INPUT: 
    # meas : observed/measured
    # pred :(dataframe or matrix) qui peut contient n series of predicted values
    # attention predicetd values in  y !!! 
    if (is.vector(pred)==TRUE) {p<-1} else {p<-dim(pred)[2]} # number of temptative predictings
    n<-length(meas) # number of observations
    RES<-matrix(0,p,17)
    colnames(RES)<-c("n","r2","elev.SMA","pvalue.elev=0","slope.SMA","pvalue.slop=1","MSE","MSEs","MSEu","CVMSE","CVMSEs","CVMSEu","RMSE","NRMSE","CVRMSE","MAE","CVMAE")
    rownames(RES)<-colnames(pred)
    for (i in 1:p) {
      if (p>1) {tmp.pred<-pred[,i]} else {tmp.pred<-pred} # number of temptative predictings
      
      SMA<-line.cis(tmp.pred,meas,method="SMA")
      SMA<-line.cis(tmp.pred,meas,method="OLS")
      
      RES[i,1]<-length(meas)
      RES[i,2]<-cor(meas,tmp.pred,use='complete.obs') # PB ici
      RES[i,3]<-smatr::elev.test(tmp.pred,meas,test.value=0,method="SMA")$a
      RES[i,4]<-smatr::elev.test(tmp.pred,meas,test.value=0,method="SMA")$p
      RES[i,5]<-smatr::slope.test(tmp.pred,meas,test.value=1,method="SMA")$b
      RES[i,6]<-smatr::slope.test(tmp.pred,meas,test.value=1,method="SMA")$p
      
      # Mean predicted from the smatr using OLS 
      tmp.predSM= SMA[2,1] * meas + SMA[1,1]
      
      # three nxt to divide by /n in the original file
      RES[i,7]<-sum((tmp.pred - meas)^(2))/n # MSE mean square error
      RES[i,8]<-sum((tmp.predSM - meas)^(2))/n # MSEs systematic MSE
      RES[i,9]<-sum((tmp.pred - tmp.predSM)^(2))/n # MSEu unsystematic MSE
      
      RES[i,10]<-RES[i,7]/mean(meas) # CVMSE
      RES[i,11]<-RES[i,8]/mean(meas) # CVMSEs
      RES[i,12]<-RES[i,9]/mean(meas) # CVMSEu
      
      RES[i,13]<-sqrt(RES[i,7]) # RMSE
      RES[i,14]<-RES[i,10]/(max(meas)-min(meas)) # NRMSE
      RES[i,15]<-RES[i,13]/mean(meas) # (CV)RMSE TO SHOW !!!
      RES[i,16]<-sum(abs(tmp.pred-meas))/n # MAE
      RES[i,17]<-RES[i,16]/mean(meas) # CVMAE
    }
    RES<-round(RES,5)
    
    return(RES)
    
  }
  f       = function(x) {
    r <- quantile(x, probs = c(0.1, 0.33, 0.5, 0.66, 0.9))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  # rouding to the nearest multiple of a base number
  mround <- function(x,base){base*round(x/base)}
  # will be applied centrally, a boxcar filter by default, wxith percentage of weight created
  filter_with_NA <- function(x, window_length=12,myfilter=rep(1/window_length,window_length), max_percentage_NA=25){
    # make the signal longer at both sides
    signal <- c(rep(NA,window_length),x,rep(NA,window_length))
    # see where data are present and not NA
    present <- is.finite(signal)
    
    # replace the NA values by zero
    signal[!is.finite(signal)] <- 0
    # apply the filter
    filtered_signal <- as.numeric(stats::filter(signal,myfilter, sides=2))
    
    # find out which percentage of the filtered signal was created by non-NA values
    # this is easy because the filter is linear
    original_weight <- as.numeric(stats::filter(present,myfilter, sides=2))
    # where this is lower than one, the signal is now artificially smaller 
    # because we added zeros - compensate that
    filtered_signal <- filtered_signal / original_weight
    # but where there are too few values present, discard the signal
    filtered_signal[100*(1-original_weight) > max_percentage_NA] <- NA
    
    # cut away the padding to left and right which we previously inserted
    filtered_signal <- filtered_signal[((window_length+1):(window_length+length(x)))]
    return(filtered_signal)
  }
  
  # load rasters 11340 x 9920 at 0.5 m res
  setwd(DATA_WD) 
 
  myDEM     <- terra::rast("myDEM.tif") 
  myDAH     <- terra::rast("myDAH.tif") 
  myHILL    <- terra::rast("myHILL.tif") 
  mySLOPE   <- terra::rast("mySLOPE.tif")
  myALPHA   <- terra::rast("myALPHA.tif")
  myALPHA.HR<- terra::rast("myALPHA.HR.tif")
  
  my3HILL    <- terra::rast("my3HILL30.tif")
  my3ALPHA   <- terra::rast("my3ALPHA.tif")
  HILL       <- terra::rast("myHILL30.tif")   
  
  # load vectors
  OBS.PTS   <- terra::vect("OBS.PTS.shp")
  SOIL.PTS  <- terra::vect("SOIL.PTS.shp")
  my1CONT   <- terra::vect("myCONT.shp")  # only RNPI 
  my2CONT   <- terra::vect("my2CONT.shp") # RNPI + LAUZ 
  my3CONT   <- terra::vect("my3CONT.shp") # RNPI + LAUZ + MAND
  values(my3CONT)$Name <- c("RNPI","LAUZ","MAND")
  
  SAG3lines <- terra::vect("SAG3lines.shp")
  
  RNPI.CONT <- terra::vect("RNPI.CONT.shp")
  LAUZ.CONT <- terra::vect("LAUZ.CONT.shp")
  MAND.CONT <- terra::vect("MAND.CONT.shp")
  
  myROUT    <- terra::vect("~/SIG/BD_CARTO/BDCARTO_5-0_TOUSTHEMES_SHP_LAMB93_D005_2024-03-15/BDCARTO/1_DONNEES_LIVRAISON_2024-07-00018/BDC_5-0_SHP_LAMB93_D005-ED2024-03-15/TRANSPORT/ROUTE_NUMEROTEE_OU_NOMMEE.shp")
  myHYDRO   <- terra::vect("~/SIG/BD_CARTO/BDCARTO_5-0_TOUSTHEMES_SHP_LAMB93_D005_2024-03-15/BDCARTO/1_DONNEES_LIVRAISON_2024-07-00018/BDC_5-0_SHP_LAMB93_D005-ED2024-03-15/HYDROGRAPHIE/COURS_D_EAU.shp")
  
  # reference grids for RN/PISSOU and WS3
  myREF01   <- terra::rast("myREF01.tif")
  myREF10   <- terra::rast("myREF10.tif")
  myREF30   <- terra::rast("myREF30.tif")
  my3REF01  <- terra::rast("my3REF01.tif")
  my3REF10  <- terra::rast("my3REF10.tif")
  my3REF30  <- terra::rast("my3REF30.tif")
  
  my3EXT    <- terra::ext(964500, 969450, 6443020, 6448690)
  
  myfont = "MS Reference Sans Serif"
  
  myWS <- c("RNPI","LAUZ","MAND","TOT")
  my3WS <- c("2.Roche Noire", "3.Lauzette", "1.Mandette")
  
  # graphical checking - issue with legend !
  myY = 23 # pick a year
  SCAtmp.L3 <- terra::rast(paste0("SCA",myY,".L3.tif"))
  SCAtmp.L3[is.na(SCAtmp.L3)]<-(-1)
  for (i in 1:dim(SCAtmp.L3)[3]) terra::coltab(SCAtmp.L3[[i]]) <- coltb
  
  graphics.off();windows(14,16)
  terra::plot(SCAtmp.L3,mar=c(0.1, 0.1, 0.1, 0.1),loc.main="bottomleft",col.main="black",cex.main=1.5,axes=F,alpha=my3ALPHA,maxnl=16,fun=function() {
    terra::plot(my3HILL,add=T,alpha=0.5,col=grey(0:100 / 100),legend=F);
    terra::plot(my3CONT,add=T,border="black"); 
    terra::plot(SAG3lines,add=T,col="gray");
    terra::plot(OBS.PTS[2],add=T,pch=22,bg="blue")
  })
  
  # Default colors in R for habitats (x7)
  # OR
  HABcol <- rev(c('#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#0c2c84'))
  HABcol <- rev(c('#d73027','#fc8d59','#fee090','#ffffbf','#e0f3f8','#91bfdb','#4575b4'))
  HABcol <- c('#762a83','#af8dc3','#e7d4e8','#f7f7f7','#d9f0d3','#7fbf7b','#1b7837')
  # https://colorbrewer2.org/#type=diverging&scheme=PRGn&n=7
  # rgb(118,42,131)','rgb(175,141,195)','rgb(231,212,232)','rgb(247,247,247)','rgb(217,240,211)','rgb(127,191,123)','rgb(27,120,55)
  
  library(scales) 
  DEF.COL <- hue_pal()(7)
  tmp <- hex2RGB(DEF.COL)@coords
  round(255*tmp)
  
  # Color of snowcolor pixels (from the DIVcbs_BGpal)
  SNOWcol <- '#4575b4'
  
  # Color of snowfree pixels (from the DIVcbs_BGpal)
  NOSNOWcol <- '#fee090'
  
# (NOT TO RUN) UTILITIES -> see SNOW_ROCHENOIRE.RData ----
  # does not include the spatial data ! 
  SMOD_pal <- colorRampPalette(
    # c(colors()[81],colors()[50],gray(0.9),'#67a9cf','#2166ac'), 
    c("palegreen4","palegreen1",gray(0.9),"deepskyblue","deepskyblue4"), 
    #c(colors()[81],colors()[50],colors()[86],gray(0.9),'#67a9cf','#2166ac'), 
    space = "Lab",interpolate="linear"
  )
  graphics.off();windows(5,1);par(mar=c(0,0,0,0));colorstrip(SMOD_pal(100), "SMOD_pal",ShowAxis=FALSE)
  
  SRTM_pal2<- colorRampPalette(
    c( rgb(128,184,134,max=255), # added
       rgb(148,191,139,max=255),
       rgb(168,198,143,max=255),
       rgb(189,204,150,max=255),
       rgb(209,215,171,max=255),
       rgb(239,235,192,max=255),
       rgb(222,214,163,max=255),
       rgb(202,185,130,max=255),
       rgb(192,154, 83,max=255),
       rgb(182,124, 53,max=255)  # added
    ),
    space = "Lab"
  )
  colorstrip(SRTM_pal2(100),"SRTM_pal2",ShowAxis=FALSE)
  
  GRAY_pal1 <- colorRampPalette(
    c(gray(0),gray(1)), 
    space = "Lab",
    interpolate='linear'
  )
  colorstrip(GRAY_pal1(100),"GRAY_pal2",ShowAxis=FALSE)
  
  # color for binary map
  coltb <- data.frame(value=c(0,1,-1), col=c("palegreen","deepskyblue","red"))
  
  # strip for facet_wrap
  COL.WS  <- SMOD_pal(10)[c(1,7,9)]
  STRIP   <- ggh4x::strip_themed(background_x = elem_list_rect(fill = COL.WS))
  STRIP4  <- ggh4x::strip_themed(background_x = elem_list_rect(fill = c(COL.WS,"gray")))
  
  # Number of cells for RN/PISSOU and WS3
  myYEAR   <- 13:23
  WS       <- c("RNPI","LAUZ","MAND")
  
  NCELL01  <- sum(terra::extract(my3REF01,my3CONT)[,2])
  NCELL10  <- sum(terra::extract(my3REF10,my3CONT)[,2])
  NCELL30  <- sum(terra::extract(my3REF30,my3CONT)[,2])
  
  NCELL01.WS3 <- NCELL10.WS3 <- NCELL30.WS3 <- rep(NA,3)
  for (i in 1:3) {
    NCELL01.WS3[i] <- sum(terra::extract(my3REF01,my3CONT[i,])[,2])
    NCELL10.WS3[i] <- sum(terra::extract(my3REF10,my3CONT[i,])[,2])
    NCELL30.WS3[i] <- sum(terra::extract(my3REF30,my3CONT[i,])[,2])
  }
  
  NCELL30tot.WS3 <- sum(NCELL30.WS3)
  
  # SMOD estimates
  SMOD.MODIS <- read.csv("SMOD.MODIS.csv")
  SMOD.CROC  <- read.csv("SMOD.CROC.csv")
  SMOD.FLUX  <- read.csv("SMOD.FLUX.csv")
  SMOD.NIVO  <- read.csv("SMOD.NIVO.csv")
  SMOD.SOIL  <- read.csv("SMOD.SOIL.csv")
  
  # save start
  setwd(WD); save.image("SNOW_ROCHENOIRE.RData")
# I. (NOT TO RUN) SITES, AREA OF INTEREST & CONTOURS (NOT TO RUN) ----
  # I.A. AREA and SITES of INTEREST (20240301) ----
    # only RN/Pissou-> my1EXT, my1CONT, my1REF30 ----
  # includes 128 row x 150 col
  # Selected projection is Lambert 93 - epsg:2154
  # Contours of ROCHE NOIRE & PISSOU watershed  -> my1CONT & my1EXT
  terra::plot(RNcont <- terra::vect("LAST_ROCHE_NOIRE.kml"))
  RNcontL93     <- terra::project(RNcont,"epsg:2154") # FINE
  setwd(WD); terra::writeVector(RNcontL93,"RNPI.CONT.shp",overwrite=T)

  # myCONT <- rbind(PIcontL93,RNcontL93,BCcontL93)
  my1CONT <- RNcontL93
  graphics.off() ; terra::plot(my1CONT)
  setwd(WD); terra::writeVector(my1CONT,"my1CONT.shp",overwrite=T)
  myEXT <- terra::ext(my1CONT)
  # adjust to multiples of 30 (& 10)
  myEXT     <- terra::ext(964500, 969000, 6443020, 6446860)
  myREF01   <- terra::rast(myEXT,res=1,crs="epsg:2154",vals=1)
  myREF10   <- terra::rast(myEXT,res=10,crs="epsg:2154",vals=1)
  myREF30   <- terra::rast(myEXT,res=30,crs="epsg:2154",vals=1)  
  
  tmp <- terra::project(myREF30,"epsg:32631")
  myEXT_UTM31 <- terra::ext(tmp)  
  
  setwd(WD); terra::writeRaster(myREF10,"myREF10.tif",overwrite=T)
  setwd(WD); terra::writeRaster(myREF30,"myREF30.tif",overwrite=T)
  setwd(WD); terra::writeRaster(myREF01,"myREF01.tif",overwrite=T)  
  
  graphics.off();windows(10,10)
  terra::plot(my1CONT)
  terra::plot(my1EXT,border="blue",add=T)
  
  # create a raster of alpha values
  myALPHA <- terra::mask(myREF30,my1CONT)
  myALPHA[is.na(myALPHA)]<-0.3
  setwd(WD); terra::writeRaster(myALPHA,"myALPHA.tif",overwrite=T)
  
    # RN/Pissou + Lauzette -> my2CONT and my2REF30  ----
    # includes with 189 row x 150 col
  setwd(WD);   terra::plot(LZcont <- terra::vect("LAST_LAUZETTE.kml"))
  LZcontL93     <- terra::project(LZcont,"epsg:2154") # FINE
  
  setwd(WD); terra::writeVector(LZcontL93,"LAUZ.CONT.shp",overwrite=T)
  
  my2CONT <- rbind(RNcontL93,LZcontL93)
  
  terra::plot(my2CONT) # only need to change the ymax border

  setwd(WD); terra::writeVector(my2CONT,"my2CONT.shp",overwrite=T)
  
  my2EXT <- terra::ext(my2CONT)
  
  # adjust to multiples of 30 (& 10) -> 189 row x 150 col
  my2EXT     <- terra::ext(964500, 969000, 6443020, 6448690)
  my2REF01   <- terra::rast(my2EXT,res=1,crs="epsg:2154",vals=1)
  my2REF10   <- terra::rast(my2EXT,res=10,crs="epsg:2154",vals=1)
  my2REF30   <- terra::rast(my2EXT,res=30,crs="epsg:2154",vals=1)

  tmp <- terra::project(my2REF30,"epsg:32631")
  my2EXT_UTM31 <- terra::ext(tmp)

  setwd(WD); terra::writeRaster(my2REF10,"my2REF10.tif",overwrite=T)
  setwd(WD); terra::writeRaster(my2REF30,"my2REF30.tif",overwrite=T)
  setwd(WD); terra::writeRaster(my2REF01,"my2REF01.tif",overwrite=T)
  
  graphics.off();windows(10,10)
  terra::plot(my2CONT)
  terra::plot(my2EXT,border="blue",add=T)
  
  # create a raster of alpha values
  my2ALPHA <- terra::mask(my2REF30,my2CONT)
  my2ALPHA[is.na(my2ALPHA)]<-0.3
  setwd(WD); terra::writeRaster(my2ALPHA,"my2ALPHA.tif",overwrite=T)
  
  # Fluxalp and Nivose
  tmp  <- data.frame(X_wgs84=c(6.38994,6.41061),Y_wgs84=c(45.05495,45.04128))
  tmp1 <- terra::vect(tmp,geom=c("X_wgs84","Y_wgs84"),crs="epsg:4326")
  OBS.PTS <- terra::project(tmp1,"epsg:2154")
  setwd(WD); terra::writeVector(OBS.PTS,"OBS.PTS.shp",overwrite=T)
  
  
    # RN/Pissou + Lauzette + Mandette -> my3CONT and my3REF30 ----
  # includes 189 row x 165 col
  setwd(WD);   terra::plot(MAcont <- terra::vect("LAST_MANDETTE.kml"))
  MAcontL93  <- terra::project(MAcont,"epsg:2154") # FINE
  
  setwd(WD); terra::writeVector(MAcontL93,"MAND.CONT.shp",overwrite=T)
  
  my3CONT    <- rbind(RNcontL93,LZcontL93,MAcontL93)
  
  terra::plot(my3CONT) # only need to change the ymax border
  
  setwd(WD); terra::writeVector(my3CONT,"my3CONT.shp",overwrite=T)
  
  my3EXT <- terra::ext(my3CONT)
  # only change the xmax
  # adjust to multiples of 30 (& 10) -> 189 row x 165 col
  my3EXT     <- terra::ext(964500, 969460, 6443020, 6448690)
  my3REF01   <- terra::rast(my3EXT,res=1,crs="epsg:2154",vals=1)
  my3REF10   <- terra::rast(my3EXT,res=10,crs="epsg:2154",vals=1)
  my3REF30   <- terra::rast(my3EXT,res=30,crs="epsg:2154",vals=1)
  
  tmp <- terra::project(my3REF30,"epsg:32631")
  my3EXT_UTM31 <- terra::ext(tmp)
  
  setwd(WD); terra::writeRaster(my3REF10,"my3REF10.tif",overwrite=T)
  setwd(WD); terra::writeRaster(my3REF30,"my3REF30.tif",overwrite=T)
  setwd(WD); terra::writeRaster(my3REF01,"my3REF01.tif",overwrite=T)
  
  graphics.off();windows(10,10)
  terra::plot(my3CONT)
  terra::plot(my3EXT,border="blue",add=T)
  
  # create a raster of alpha values
  my3ALPHA <- terra::mask(my3REF30,my3CONT)
  my3ALPHA[is.na(my3ALPHA)]<-0.3
  setwd(WD); terra::writeRaster(my3ALPHA,"my3ALPHA.tif",overwrite=T)
  
    # Fluxalp and Nivose ----
  tmp  <- data.frame(X_wgs84=c(6.38994,6.41061),Y_wgs84=c(45.05495,45.04128))
  tmp1 <- terra::vect(tmp,geom=c("X_wgs84","Y_wgs84"),crs="epsg:4326")
  OBS.PTS <- terra::project(tmp1,"epsg:2154")
  setwd(WD); terra::writeVector(OBS.PTS,"OBS.PTS.shp",overwrite=T)
  
    # Emprise plateforme et Lautaret ----
  setwd("~/COORDINATIONmc/eLTER/LTSER_LAUTARET_OISANS/202403_EMPRISES")
  graphics.off(); LTSER <- terra::vect("Projet_Plateforme_LTSER_LAUTARET_OISANS_2024.shp")
  LTSER.l93 <- terra::project(LTSER,"epsg:2154")
  LTSERsat <- terra::vect("satellites.shp")
  LTSER.l93 <- terra::project(LTSER,"epsg:2154")
  LTSER.l93 <- terra::project(LTSER,"epsg:2154")
  writeVector(LTSER,"LTSER.kml")
  writeVector(LTSERsat,"LTSERsat.kml")
  # I.B. PLOT DATA -> ALL.PTS (20250121) ------
  # using all plots of inventaires capteurs microclimat
  # this includes the GAL and Aravo sectors
  ALL.PTS <- read.table("~/CLIMATOmc/20250121_PLOTS.txt",header=T,sep="\t")
  # select area of interest
  ALL.PTS <- ALL.PTS[which(ALL.PTS$LOC=="LAU" | ALL.PTS$LOC=="DAC" | ALL.PTS$LOC=="RON" | ALL.PTS$LOC=="ARA" | ALL.PTS$LOC=="GAL"),]
  # add Lambert93 coordinates
  ALL.PTS.vec <- terra::vect(ALL.PTS,geom=c("LONG_DD","LAT_DD"),crs="epsg:4326")
  ALL.PTS.vec <- terra::project(ALL.PTS.vec,"epsg:2154")
  tmp <- terra::crds(ALL.PTS.vec)
  colnames(tmp) <- c("X_L93","Y_L93")
  ALL.PTS <- cbind(ALL.PTS,tmp)
  
  # PLOTS with soil data [1-50]
  setwd(DATA_WD)
  SOIL.PTS <- read.csv("all_plot_soil.csv",sep=";")
  colnames(SOIL.PTS)[2]<-"ID"
  table(SOIL.PTS$id_community)
  which(is.na(match(SOIL.PTS$ID,ALL.PTS$ID))) # should return integer(0)
  SOIL.PTS$CNratio <- SOIL.PTS$X.C/SOIL.PTS$X.N
  
  ALL.PTS <- merge(x=ALL.PTS,y=SOIL.PTS[,-c(5:6)], by="ID", all.x=TRUE)
  
  # DO NOusing RNED plots (not necessary as included in ALL.PTS)
  # PLOTS  <- read.csv("~/PROJETSmc/ROCHENOIRE/DATA_ANALYSIS/RNED_PLOTS.csv",header=T,sep=";")[1:134,] 
  # RNED_PLOTS.csv is the csv export of the sheet PLOTS of RNED_WorkingDoc.xlsx
  # PLOTS$ID
  # PLOTS$IDs <- gsub("_","",gsub("RON_","",PLOTS$ID))
  # SOIL.PTS <- terra::vect(PLOTS,geom=c("X_L93","Y_L93"),crs="epsg:2154")

  # PLOTS with soil temperature SOIL.PTS
  # Files have been processed in climato.r

  SMOD.PTS   <- read.csv("~/CLIMATOmc/DATA_ANALYSIS/20250121_BIOCLIM_SOIL_TEMP_GH_SMOD.csv")  
  colnames(SMOD.PTS)<-c("ID",paste0("SMOD",1999:2023))
  ALL.PTS   <- merge(x=ALL.PTS,y=SMOD.PTS, by="ID", all.x=TRUE)
  
  SOD.PTS   <- read.csv("~/CLIMATOmc/DATA_ANALYSIS/20250121_BIOCLIM_SOIL_TEMP_GH_SOD.csv")  
  colnames(SOD.PTS)<-c("ID",paste0("SOD",1999:2023))
  ALL.PTS   <- merge(x=ALL.PTS,y=SOD.PTS, by="ID", all.x=TRUE)

  ALL.PTS$SMODmean <- apply(ALL.PTS[,grep("SMOD",colnames(ALL.PTS))],1,mean,na.rm=T) 
  ALL.PTS$SODmean <- apply(ALL.PTS[,grep("SOD",colnames(ALL.PTS))],1,mean,na.rm=T) 
  
  # selection of habitats
  # adjust and select COM for consistency
  ALL.PTS$COMcor <- ALL.PTS$COM
  # Combine AP and CF
  ALL.PTS$COMcor <- ifelse(ALL.PTS$COMcor=="AP","CF",ALL.PTS$COMcor)
  # Combine AR and RK and PA
  ALL.PTS$COMcor <- ifelse(ALL.PTS$COMcor=="AR","PA",ALL.PTS$COMcor)
  ALL.PTS$COMcor <- ifelse(ALL.PTS$COMcor=="RK","PA",ALL.PTS$COMcor)
  # Combine KD and KS
  ALL.PTS$COMcor <- ifelse(ALL.PTS$COMcor=="KD","KS",ALL.PTS$COMcor)
  ALL.PTS$COMcor <- ifelse(ALL.PTS$COMcor=="EY","EN",ALL.PTS$COMcor)
  
  # select the North exposed flatty sequence EN, CF, CT, VM, FP, KS
  ALL.PTS <- ALL.PTS[ALL.PTS$COMcor %in% c("EN","CF","CT","VM","FP","KS","PA"),]
  table(ALL.PTS$COMcor)
  graphics.off();windows(10,10);boxplot(SMODmean~COMcor,data=ALL.PTS)
  graphics.off();windows(10,10);boxplot(SODmean~COMcor,data=ALL.PTS)
  graphics.off();windows(10,10);boxplot(CNratio~COMcor,data=ALL.PTS)
  
  setwd(DATA_WD) ; write.csv(ALL.PTS,"ALL.PTS.csv",row.names=F)
  
# OLD - NOT TO RUN PLOTS with simulated soil temperatures (Cesar Berger - CROCUS/SAFRAN) 
# Files have been processed in climato.r
  # tmp   <- read.csv("~/CLIMATOmc/DATA_ANALYSIS/OLD/BIOCLIM_SOIL_TEMP_S_FSFD.csv",row.names=1)
  # sel      <- PLOTS$ID%in%rownames(tmp)
  # mySITsim <- PLOTS[sel,]  
  
  # I.C. HYDROGRAPHY ----
  # http://sandre.eaufrance.fr/geonetwork/srv/fr/main.home
  tmp   <- terra::vect("D:/HYDRO/ZRPEHY/ZRPEHY.shp")
  tmp   <- terra::vect("D:/HYDRO/LWBODY/LWBODY.shp")
  
  tmp      <- terra::vect("D:/HYDRO/RWBODY/RWBODY.shp")
  myRWBODY <- terra::crop(terra::project(tmp,"epsg:2154"),myEXT)
  # tmp <- terra::crop(terra::project(SH,"epsg:2154"),myEXT)
  # terra::plot(tmp)
  
  # I.D. DEM ----
    # I.D.1. DEM LIDAR IGN aggregated at 1m -> myDEM  (20250315) ----
    # crop the Lidar data prepared by Arthur on March 2023 -> myDEM
    setwd(DATA_WD)
    tmp  <- terra::rast("mnt_50cm_lautaret.tif")
    tmp  <- terra::crop(tmp, my3EXT)
    terra::writeRaster(tmp,"myDEM.tif",overwrite=T)
    
    DEM <- terra::rast("myDEM.tif")
    
    # compute topographical indices at 0.5 m resolution
    setwd(DATA_WD)
    terra::terrain(DEM, v='slope',unit="radians",filename="mySLOPE.tif",overwrite=T)
    terra::terrain(DEM, v='aspect',unit="radians",filename="myASPECT.tif",overwrite=T)
    SLO <- terra::rast("mySLOPE.tif")
    ASP <- terra::rast("myASPECT.tif")
    DAH    <- cos(202.5*pi/180-ASP)*atan(SLO)
    terra::writeRaster(DAH,"myDAH.tif",overwrite=T)
    
    # A better hill shade may be achieved by combining different angles and directions
    hh   <- shade(SLO, ASP, angle = c(45, 45, 45, 80), direction = c(225, 270, 315, 135))
    h1   <- Reduce(mean, hh)
    HILL <- mean(hh)
    # HILL   <- terra::shade(SLO, ASP, 30, 202,normalize=T)
    terra::writeRaster(HILL,filename="myHILL.tif",overwrite=T)
    
    myALPHA.HR <- terra::mask(myDEM,my3EXT)
    myALPHA.HR[!is.na(myALPHA.HR)] <- 0.9
    myALPHA.HR[is.na(myALPHA.HR)]  <- 0.3
    terra::writeRaster(myALPHA.HR,"myALPHA.HR.tif",overwrite=T)    
    
    
    # aggregate to 1M and 2M resolution    
    setwd(DATA_WD)
    tmp2m <- terra::aggregate(DEM,fact=4,fun="mean")
    terra::writeRaster(tmp2m,"DEM_2M.tif",overwrite=T)
    tmp1m <- terra::aggregate(DEM,fact=2,fun="mean")
    terra::writeRaster(tmp1m,"DEM_1M.tif",overwrite=T)
    tmp2m <- terra::aggregate(DEM,fact=20,fun="mean")
    terra::writeRaster(tmp2m,"DEM_10M.tif",overwrite=T)
    
    tmp2m <- terra::aggregate(SLO,fact=4,fun="mean")
    terra::writeRaster(tmp2m,"SLO_2M.tif",overwrite=T)
    tmp1m <- terra::aggregate(SLO,fact=2,fun="mean")
    setwd(DATA_WD) ;terra::writeRaster(tmp1m,"SLO_1M.tif",overwrite=T)
    
    tmp2m <- terra::aggregate(HILL,fact=4,fun="mean")
    terra::writeRaster(tmp2m,"HILL_2M.tif",overwrite=T)
    tmp1m <- terra::aggregate(HILL,fact=2,fun="mean")
    setwd(DATA_WD) ;terra::writeRaster(tmp1m,"HILL_1M.tif",overwrite=T)
    
    tmp1m <- terra::aggregate(DAH,fact=2,fun="mean")
    terra::writeRaster(tmp1m,"DAH_1M.tif",overwrite=T)
    tmp2m <- terra::aggregate(DAH,fact=4,fun="mean")
    terra::writeRaster(tmp2m,"DAH_2M.tif",overwrite=T)
    tmp1m <- terra::aggregate(DAH,fact=20,fun="mean")
    terra::writeRaster(tmp1m,"DAH_10M.tif",overwrite=T)
    

    # I.D.2. DEM SINTEGRA at 50CM resolution UNUSED ! (20240205) ----
  # # see rochenoire.r section I.D.1 for the production of the DEM
  # setwd("~/PROJETSmc/ROCHENOIRE/DATA_ANALYSIS")
  # tmp <- terra::rast("RN_dem_50CM.tif")
  # terra::crs(tmp) <- "epsg:2154"
  # myDEM <- terra::crop(tmp,myEXT)
  # terra::crs(myDEM) <- "epsg:2154"
  
    # I.D.3. NOT TO RUN SAGA-derived variables on myDEM ---- 
  setwd(WD)
  RSAGA::rsaga.import.gdal("myDEM.tif","myDEM.sgrd")

  RSAGA::rsaga.contour("myDEM.sgrd", "SAGlinesHR", 50, 1500, 3500)
  graphics.off();windows(15,15);par(mar=c(2,2,0,0),oma=c(0,0,0,0))
  
  # Wetness index - quite similat to TPI
  rsaga.wetness.index("DEM10.tif","SWI.sgrd")
  SWI <- raster("SWI.sdat") ; SWI@crs <- LAMB93
  SWI25M <- projectRaster(SWI,REF25)
  graphics.off();windows(10,10);plot(SWI25M)
  
  plot((SWI25M*POI)[],(DAH*POI)[],pch=".")

  # topdown processing and ridge/hollow delineation - quite long
  RSAGA::rsaga.topdown.processing("myDEM.tif", out.carea = "CAREA")
  
  graphics.off();windows(20,20)
  terra::plot(TD<-terra::rast("CAREA.sdat"),range=c(0,100))

  RIDGE <- TD; RIDGE[]<-0;RIDGE[which(TD[]>=5)]<-1
  hist(TD[])
  
  graphics.off();windows(20,20);plot(RIDGE,col=c("black","gray"))
  
  RIDpt <- raster::extract(RIDGE,PLOTS[,15:16])
  names(RIDpt)<-PLOTS[,"ID"]
  table(RIDpt,PLOTS[,"COM"])
  
  # route from sink - useless
  RSAGA::rsaga.sink.route("myDEM.tif", out.sinkroute="myDRAIN")
  
  graphics.off();windows(20,20)
  plot(DR <- terra::rast("myDRAIN.sdat"),range=c(0,8))
  
  # Watershed delineation after filling sinks
  
  RSAGA::rsaga.fill.sinks("myDEM.tif", "myWS.sgrd",out.wshed="myWSb", method="wang.liu.2006")
  
  graphics.off();windows(20,20);
  plot(WS <- terra::rast("myWSb.sdat"),col=rainbow(15))
  terra::plot(my3CONT,add=T)
  
  # WIND SHELTER index  -> SHELTt 
  # used as a proxy for snow accumulation on the lee side
  # first trial with a degraded dem model at 10m resolution

  # direction = direction from which the wind originates
  # for example -45?=-pi/4 north-west means 
  # the more negative the most exposed to prevailing wind conditions

  writeRaster(raster("DEM10.tif"),"DEM10.asc",format="ascii",overwrite=T)
  
  RADIUS   <- 10  # in number of cells -> 60 m
  CELLSIZE <- 10 # DEM resolution
  TOL      <- pi/10 # pour simul b (pi/12 pour sinon)
  ctrl   <- wind.shelter.prep(radius=RADIUS, direction= -pi/2, tolerance = pi/10,cellsize = CELLSIZE)
  focal.function("DEM10.asc", fun=wind.shelter, control=ctrl,radius=RADIUS, search.mode="circle")
  writeRaster(raster("windshelter.asc"),"SHELT10r100b.tif",format="GTiff")
  
  graphics.off()
  SHELTR60 <- raster("SHELT10r60.tif") ; SHELTR60@crs <- LAMB93
  SHELT25R60 <- projectRaster(SHELTR60,REF25)
  windows(10,8) ; plot(crop(SHELT25R60,RNse),col=rev(GRAY_pal1(100)))
  
  SHELTR100 <- raster("SHELT10r100.tif") ; SHELTR100@crs <- LAMB93
  SHELT25R100 <- projectRaster(SHELTR100,REF25)
  windows(10,8) ; plot(crop(SHELT25R100,RNse),col=rev(GRAY_pal1(100)))
  
  SHELTR100b <- raster("SHELT10r100b.tif") ; SHELTR100b@crs <- LAMB93
  SHELT25R100b <- projectRaster(SHELTR100b,REF25)
  windows(10,8) ; plot(crop(SHELT25R100b,RNse),col=rev(GRAY_pal1(100)))
  
  SHELTt <- SHELT25R100*POI 
  
  # solar radiation
  TEST <- rsaga.pisr2("SAGDEM10m.sgrd",out.total.grid="TOTSRAD10m",latitude=45,time.step=1,start.date=list(day=15,month=11,year=2013),end.date=list(day=16,month=11,year=2013))
  
  # solar radiation
  TEST <- rsaga.pisr2("SAGDEM10m.sgrd",out.total.grid="TOTSRAD10m",latitude=45,time.step=1,start.date=list(day=15,month=2,year=2015),end.date=list(day=17,month=2,year=2015))
  # use month 2 for January !
  plot(raster("TOTSRAD10m.sdat"))
  
  
    # I.D.4. DEM IGN reprojected at 30m -> myDEM30 (20240217) ----
        # DEM at 25m from BD_ALTI
  DEM1      <- terra::rast("D:/DEM/DEPTS25m/Dpt_05_asc.asc")
  DEM2      <- terra::rast("D:/DEM/DEPTS25m/Dpt_38_asc.asc")
  DEM3      <- terra::rast("D:/DEM/DEPTS25m/Dpt_73_asc.asc")
  TMP       <- terra::merge(DEM1,DEM2,DEM3)
  terra::crs(TMP) <- "epsg:2154"
 
      # for myREF ---- 
  DEM30     <- terra::project(TMP,myREF30)
  setwd(WD)
  terra::writeRaster(DEM30,"myDEM30.tif",overwrite=T)
  terra::terrain(DEM30, v='slope',unit="radians",filename="mySLOPE30.tif",overwrite=T)
  terra::terrain(DEM30, v='aspect',unit="radians",filename="myASPECT30.tif",overwrite=T)
  slo <- terra::rast("mySLOPE30.tif"); asp <- terra::rast("myASPECT30.tif")
  DAH    <- cos(202.5*pi/180-asp)*atan(slo)
  HILL   <- terra::shade(slo, asp, 40, 202)
  terra::writeRaster(HILL,filename="myHILL30.tif",overwrite=T)
  terra::writeRaster(DAH,filename="myDAH30.tif",overwrite=T)

  tm_shape(HILL) + tm_raster(alpha = 0.4,legend.show=F,palette="-Greys",n=100)
  
  terra::plot(HILL,alpha=0.4,col=grey(0:100 / 100),legend=F,axes=F)
  plot(my3CONT,add=T)
  
      # same for my3REF30 ----
  DEM30     <- terra::project(TMP,my3REF30)
  setwd(WD)
  terra::writeRaster(DEM30,"my3DEM30.tif",overwrite=T)
  terra::terrain(DEM30, v='slope',unit="radians",filename="my3SLOPE30.tif",overwrite=T)
  terra::terrain(DEM30, v='aspect',unit="radians",filename="my3ASPECT30.tif",overwrite=T)
  slo <- terra::rast("my3SLOPE30.tif"); asp <- terra::rast("my3ASPECT30.tif")
  DAH    <- cos(202.5*pi/180-asp)*atan(slo)
  HILL   <- terra::shade(slo, asp, 40, 202)
  terra::writeRaster(HILL,filename="my3HILL30.tif",overwrite=T)
  terra::writeRaster(DAH,filename="my3DAH30.tif",overwrite=T)
  
  tm_shape(HILL) + tm_raster(alpha = 0.4,legend.show=F,palette="-Greys",n=100)
  
  terra::plot(HILL,alpha=0.4,col=grey(0:100 / 100),legend=F,axes=F)
  plot(my3CONT,add=T)
  
    # I.D.4. NOT TO RUN SAGA-derived variables on myDEM30 ---- 
  setwd(WD)
  RSAGA::rsaga.import.gdal("myDEM30.tif","myDEM.sgrd")
  RSAGA::rsaga.import.gdal("my3DEM30.tif","my3DEM.sgrd")
  
  RSAGA::rsaga.contour("myDEM.sgrd", "SAGlines", 100, 1500, 3500)
  RSAGA::rsaga.contour("my3DEM.sgrd", "SAG3lines", 100, 1400, 3400)
  
  graphics.off();windows(15,15);par(mar=c(2,2,0,0),oma=c(0,0,0,0))
  terra::plot(terra::vect("SAGlines.shp"),col="gray")
  
  graphics.off();windows(15,15);par(mar=c(2,2,0,0),oma=c(0,0,0,0))
  terra::plot(terra::vect("SAG3lines.shp"),col="gray")
  
  # Wetness index - quite similat to TPI
  rsaga.wetness.index("DEM10.tif","SWI.sgrd")
  SWI <- raster("SWI.sdat") ; SWI@crs <- LAMB93
  SWI25M <- projectRaster(SWI,REF25)
  graphics.off();windows(10,10);plot(SWI25M)
  
  plot((SWI25M*POI)[],(DAH*POI)[],pch=".")
  
  # topdown processing and ridge/hollow delineation - quite long
  RSAGA::rsaga.topdown.processing("myDEM.tif", out.carea = "CAREA")
  
  graphics.off();windows(20,20)
  terra::plot(TD<-terra::rast("CAREA.sdat"),range=c(0,100))
  
  RIDGE <- TD; RIDGE[]<-0;RIDGE[which(TD[]>=5)]<-1
  hist(TD[])
  
  graphics.off();windows(20,20);plot(RIDGE,col=c("black","gray"))
  
  RIDpt <- raster::extract(RIDGE,PLOTS[,15:16])
  names(RIDpt)<-PLOTS[,"ID"]
  table(RIDpt,PLOTS[,"COM"])
  
  # route from sink - useless
  RSAGA::rsaga.sink.route("myDEM.tif", out.sinkroute="myDRAIN")
  
  graphics.off();windows(20,20)
  plot(DR <- terra::rast("myDRAIN.sdat"),range=c(0,8))
  
  # Watershed delineation after filling sinks
  
  RSAGA::rsaga.fill.sinks("myDEM.tif", "myWS.sgrd",out.wshed="myWSb", method="wang.liu.2006")
  
  graphics.off();windows(20,20);
  plot(WS <- terra::rast("myWSb.sdat"),col=rainbow(15))
  terra::plot(my3CONT,add=T)
  
  # WIND SHELTER index  -> SHELTt 
  # used as a proxy for snow accumulation on the lee side
  # first trial with a degraded dem model at 10m resolution
  
  # direction = direction from which the wind originates
  # for example -45?=-pi/4 north-west means 
  # the more negative the most exposed to prevailing wind conditions
  
  writeRaster(raster("DEM10.tif"),"DEM10.asc",format="ascii",overwrite=T)
  
  RADIUS   <- 10  # in number of cells -> 60 m
  CELLSIZE <- 10 # DEM resolution
  TOL      <- pi/10 # pour simul b (pi/12 pour sinon)
  ctrl   <- wind.shelter.prep(radius=RADIUS, direction= -pi/2, tolerance = pi/10,cellsize = CELLSIZE)
  focal.function("DEM10.asc", fun=wind.shelter, control=ctrl,radius=RADIUS, search.mode="circle")
  writeRaster(raster("windshelter.asc"),"SHELT10r100b.tif",format="GTiff")
  
  graphics.off()
  SHELTR60 <- raster("SHELT10r60.tif") ; SHELTR60@crs <- LAMB93
  SHELT25R60 <- projectRaster(SHELTR60,REF25)
  windows(10,8) ; plot(crop(SHELT25R60,RNse),col=rev(GRAY_pal1(100)))
  
  SHELTR100 <- raster("SHELT10r100.tif") ; SHELTR100@crs <- LAMB93
  SHELT25R100 <- projectRaster(SHELTR100,REF25)
  windows(10,8) ; plot(crop(SHELT25R100,RNse),col=rev(GRAY_pal1(100)))
  
  SHELTR100b <- raster("SHELT10r100b.tif") ; SHELTR100b@crs <- LAMB93
  SHELT25R100b <- projectRaster(SHELTR100b,REF25)
  windows(10,8) ; plot(crop(SHELT25R100b,RNse),col=rev(GRAY_pal1(100)))
  
  SHELTt <- SHELT25R100*POI 
  
  # solar radiation
  TEST <- rsaga.pisr2("SAGDEM10m.sgrd",out.total.grid="TOTSRAD10m",latitude=45,time.step=1,start.date=list(day=15,month=11,year=2013),end.date=list(day=16,month=11,year=2013))
  
  # solar radiation
  TEST <- rsaga.pisr2("SAGDEM10m.sgrd",out.total.grid="TOTSRAD10m",latitude=45,time.step=1,start.date=list(day=15,month=2,year=2015),end.date=list(day=17,month=2,year=2015))
  # use month 2 for January !
  plot(raster("TOTSRAD10m.sdat"))
  
  
  # I.E. S2M reanalyses -> SMOD.CROC.csv (20240209) ----
  SAFSRC <-"D:/DATA_S2M/CROCUS_flat_alpes_DAILY_CSV/"
  LF  <- list.files(SAFSRC,full=T)
  
  # tmp0 <- read.csv(LF[180])  # Oisans 1650 1950 is file 182
  # tmp2 <- read.csv(LF[182])  # Oisans 2250 2550 is file 182  
  
  tmp0 <- read.csv(LF[225])  # Thabor 1950 2250 is file 182  
  tmp1 <- read.csv(LF[181])  # Oisans 1950 2250 is file 182
  
  dfC <- data.frame(
    DATE      = as.Date(tmp0$DATE,format="%Y-%m-%d"),
    YEAR      = lubridate::year(tmp0$DATE),
    MONTH     = lubridate::month(tmp0$DATE),
    SNOWdepth = 0.5*(tmp0$SNOWdepth + tmp1$SNOWdepth)
  )
  
  SAFSRC <-"D:/DATA_S2M/SAFRAN_flat_alpes_DAILY_CSV/"
  LF  <- list.files(SAFSRC,full=T)
  tmpS0 <- read.csv(LF[225]) # Thabor 1950 2250 is file 182  
  tmpS1 <- read.csv(LF[181])  # Oisans 1950 2250 is file 182
 dim(tmpS0)
  
  # full time series
  dfS <- data.frame(
    DATE      = as.Date(tmpS0$DATE,format="%Y-%m-%d"),
    YEAR      = lubridate::year(tmpS0$DATE),
    MONTH     = lubridate::month(tmpS0$DATE),
    Tmean     = 0.5*(tmpS0$Tmean + tmpS1$Tmean),
    SWdown    = 0.5*(tmpS0$DIR_SWdown + tmpS1$DIR_SWdown + tmpS0$DIF_SWdown + tmpS1$DIF_SWdown)
  )

  graphics.off()
  windows(20,5)
  par(mar=c(3,3,0,0))
  plot(dfC$DATE,dfC$SNOWdepth,type="h",ylim=c(0,2.5),col="darkblue")
  
  
  # extract the mean snowdepth for January-March -> S2M.SnowDepth.csv
  tmp  <- dfC[which(dfC$MONTH<=3),c("DATE","YEAR","MONTH","SNOWdepth")]
  SNOW.DEPTH <- aggregate(tmp$SNOWdepth,list(tmp$YEAR),mean)
  colnames(SNOW.DEPTH) <- c("Year","SnowDepth")
  
  # extract the temp and SW radiation for April-June  
  tmp  <- dfS[which(dfS$MONTH>=4 & dfS$MONTH<=6),c("DATE","YEAR","MONTH","Tmean","SWdown")]
  SNOW.METEO1 <- aggregate(tmp$Tmean-273.15,list(tmp$YEAR),mean)
  SNOW.METEO2 <- aggregate(tmp$SWdown,list(tmp$YEAR),sum)
  SNOW.METEO  <- cbind(SNOW.METEO1,SNOW.METEO2[,2])
  colnames(SNOW.METEO) <- c("Year","Tmean_AMJ","SWdown_AMJ")
  
  setwd(WD)
  DATA <- read.csv("FLUXALP.daily.csv")
  DATA$DATE <- as.Date(DATA$DATE,format="%Y-%m-%d")

  graphics.off()
  windows(20,5)
  par(mar=c(3,3,0,0))
  plot(DATA$DATE,DATA$SNOW,type="h",ylim=c(-0.5,2.5),col="darkblue")
  lines(DATA[,"DATE"],DATA[,"NDVI"],col="green")
  abline(h=0.08,col="white")
  abline(h=-0.075,col="darkgreen")
  lines(tmp1$DATE,tmp1$SNOWdepth,type="l",ylim=c(0,2.5),col="red",lwd=2)  
  
  graphics.off()
  windows(15,5)
  plot(SNOW.DEPTH$Year,SNOW.DEPTH$SnowDepth,type="h",pch=21,bg="lightgray",ylab="JFM Mean Snow Depth (m)",lwd=3)
  abline(h=Hmisc::smedian.hilow(SNOW.DEPTH$SnowDepth,0.75),lty=c(1,2,2))
  
  tmp <- tmp1
  tmp1 <- tmp[which(lubridate::year(tmp$DATE)>=2013),c("DATE","SNOWdepth")]  
  
  graphics.off()
  windows(20,5)
  par(mar=c(3,3,0,0))
  plot(tmp1$YEAR,tmp1$SNOWdepth,type="h",ylim=c(0,2.5),col="darkblue")
  lines(DATA$DATE,DATA$SNOW,type="l",ylim=c(0,2.5),col="red",lwd=2)  
  
  # import snowcover duration from Crocus model
  SAFDES <-"D:/DATA_S2M/CROCUS_flat_alpes_SNOWCOVER_CSV/"
  LF  <- list.files(SAFDES,full=T)
  # tmp0  <- read.csv(LF[180])[,2:3]  # Oisans 1650 1950 is file 182
  # tmp1  <- read.csv(LF[181])[,2:3]   # Oisans 1950 2250 is file 182
  # tmp2  <- read.csv(LF[182])[,2:3]   # Oisans 2250 2550 is file 182
  
  tmp0  <- read.csv(LF[225])[,2:3]    # Thabor 1950 2250 is file 182
  tmp1  <- read.csv(LF[181])[,2:3]   # Oisans 1950 2250 is file 182
  
  # full time series
  df.smod <- data.frame(
    YEAR      = tmp0$YEAR,
    SMOD      =  0.5*(tmp0$SMOD + tmp1$SMOD)
  )
  
  SMOD.CROC <- cbind(SNOW.DEPTH,df.smod[,2],SNOW.METEO[,2:3])
  colnames(SMOD.CROC)[1:3] <- c("Year","SnowDepth_JFM","SMOD")
  setwd(WD); write.csv(SMOD.CROC,"SMOD.CROC.csv",row.names=F)
  
  # Compute a DELTA SMOD based on S2M 
  SAFDES <-"D:/DATA_S2M/CROCUS_flat_alpes_SNOWCOVER_CSV/"
  LF  <- list.files(SAFDES,full=T,pattern="Thabor")[c(4,6)]
  # compare the bands 1950-2250 & 2550-2850 (i.e. 2100 - 2700)
 
  SMOD.S2M <- data.frame(YEAR     = 1959:2023,
                         SMOD2100 = read.csv(LF[1])[,3],
                         SMOD2700 = read.csv(LF[2])[,3]
                         )
  SMOD.S2M$DELTA <- SMOD.S2M$SMOD2700-SMOD.S2M$SMOD2100
  
  plot(SMOD.S2M$YEAR,SMOD.S2M$DELTA,type="b",pch=21,bg="gray")
  plot(SMOD.S2M$YEAR[26:65],SMOD.S2M$DELTA[26:65],type="b",pch=21,bg="red",cex=2)
  
  # I.F. FLUXALP data   -> SMOD.FLUX.csv (20240209) ----
  # clean data of hydrological years are available on gricad
  # https://gricad-gitlab.univ-grenoble-alpes.fr/lautaret/fluxalp/-/tree/master/data_clean_hydro_year
  # see climato.r part III.2 for the update of data production
  
  setwd("~/PROJETSmc/FLUXALP/DATA"); DATA <- read.csv("FLUXALP.daily.csv")
  DATA$DATE <- as.Date(DATA$DATE,format="%Y-%m-%d")
  
  SEL <- which(lubridate::year(as.Date(DATA$DATE))>=2013)
  graphics.off()
  windows(20,14)
  par(mar=c(3,3,0,0))
  plot(DATA$DATE[SEL],DATA$SNOW[SEL],type="h",ylim=c(-0.5,2.5),col="darkblue")
  lines(DATA[SEL,"DATE"],DATA[SEL,"NDVI"],col="green")
  lines(DATA$DATE[SEL],tmp <- data.table::frollapply(DATA$NDVI[SEL], n=15,FUN="mean",align="center"),type="l",lwd=2)
  abline(h=0.08,col="black")
  abline(h=0,col="red")
  abline(h=0,col="darkgreen")
  
    # OPT FLUXALP - retrieve SMOD and optimize with imagery ----
  # SMODf <- function(x) 50+which.min(x[51:180])
  THRndvi<- 0.004
  SMODf <- function(x) 50+which(x[51:180]<0.09)[1]
  YEAR  <- lubridate::year(DATA$DATE)
  FLUX.SMOD1 <- aggregate(DATA$SNOW,list(year=YEAR),SMODf)[-1,]
  
  NDVIf <- function(x) 50+which(x[51:180]>(THRndvi))[1]
  FLUX.SMOD2 <- aggregate(DATA$NDVI,list(year=YEAR),NDVIf)[-1,]
  
  SMOD.FLUX <- cbind(FLUX.SMOD1,FLUX.SMOD2[,2])
  colnames(SMOD.FLUX) <- c("YEAR","SMOD_snow","SMOD_ndvi")
  
  YEAR <- 2013:2023; SMOD.FLUX$SMOD_esti <- NA
  for (i in 1:length(YEAR)){
    print(i)
    SEL <- which(lubridate::year(as.Date(DATA$DATE))==YEAR[i])
    tmp <- data.table::frollapply(DATA$NDVI[SEL], n=15,FUN="mean",align="center")
    SMOD.FLUX$SMOD_esti[i] <- which(tmp>0.1)[1]
  }
  
  # compare with imagery
  setwd(WD); LF <- list.files(pattern=".SMODs.tif")
  SMODs <- terra::rast(LF);names(SMODs) <- 2013:2023
  tmp <- t(terra::extract(terra::focal(SMODs,w=3,"mean"),OBS.PTS))[-1,]
  
  graphics.off();windows(10,10)
  plot(SMOD.FLUX$SMOD_ndvi,tmp[,2],xlim=c(120,160),ylim=c(120,160),pch=21,bg="green",cex=2,xlab="SMOD_ndvi",ylab="SMOD satellite")
  text(SMOD.FLUX$SMOD_ndvi,tmp[,2],2013:2023)
  abline(0,1); abline(-5,1,lty=2);abline(5,1,lty=2)
  cor(SMOD.FLUX$SMOD_ndvi,tmp[,2]) # optim at 0.897
  
  setwd(WD); write.csv(SMOD.FLUX,"SMOD.FLUX.csv",row.names = F)
  
  # underestimation of SMOD_snow in 2021 due to gaps in data
  # Conclusion: target the years 2013 (snowiest) and 2022 (driest) and 2015 (intermediate)
  
  # I.G. NIVOSE GALIBIER data    -> SMOD.NIVO.csv (20240209) ----
  
  setwd("~/CLIMATOmc/DATA/CLIMATHEQUE")
  NIVOd <- read.csv("NIVOSE_GALIBIER.txt",sep="\t")
  NIVOd$DATE <- as.Date(NIVOd$DATE,format="%d/%m/%Y")
  
  NIVOd <- NIVOd[,c("DATE","TM","NEIGETOT06")]
  colnames(NIVOd) <- c("DATE","TM","SNOW")
  
  setwd(WD);write.csv(NIVOd,"NIVOSE.daily.csv",row.names=F)
  
  graphics.off()
  windows(20,5)
  par(mar=c(3,3,0,0))
  plot(NIVOd$DATE,NIVOd$SNOW,type="h",ylim=c(-0.5,200),col="darkblue")

    # OPT NIVO - function to retrieve SMOD and optimize with imagery ----
  setwd(WD); NIVOd <- read.csv("NIVOSE.daily.csv")
  THRsno <- 5
  SMODf <- function(x) 180-length(which(x[101:180]<=THRsno))
  YEAR <- lubridate::year(NIVOd$DATE)
  SMOD.NIVO <- aggregate(NIVOd$SNOW,list(YEAR),SMODf)[-9,]
  colnames(SMOD.NIVO)<- c("YEAR","SMOD_nivose")
  
  # compare with imagery
  tmp <- t(terra::extract(terra::focal(SMODs,w=3,"mean"),OBS.PTS))[-1,]
  plot(SMOD.NIVO$SMOD_nivose[2:8],tmp[4:10,1],xlim=c(120,160),ylim=c(120,160),pch=21,bg="gray",cex=2)
  abline(0,1); abline(-5,1,lty=2); abline(5,1,lty=2)
  cor(SMOD.NIVO$SMOD_nivose[2:8],tmp[4:10,1]) # maximisé à 0.793
  
  setwd(WD); write.csv(SMOD.NIVO,"SMOD.NIVO.csv",row.names = F)
  
  # I.H. SoilTemp data  -> SMOD.SOIL.csv (20240209) ----
  setwd("~/CLIMATOmc/DATA_ANALYSIS")
  tmp1        <- read.csv("20231122_BIOCLIM_SOIL_TEMP_GH_FSFD.csv",row.names=1)
  colnames(tmp1) <- 1999:2023
  # select the recent period
  tmp1 <- tmp1[,which(as.numeric(colnames(tmp1))>=2013)]
  # select sites
  tmp1 <- tmp1[which(substr(rownames(tmp1),1,3)=="RON" & substr(rownames(tmp1),5,6)!="BE"),]

  head(tmp1)
  
  na_count <-sapply(tmp1, function(x) sum(length(which(!is.na(x)))))
  
  setwd(WD)
  write.csv(tmp1,"SMOD.SOIL.csv",row.names = T)
  
  # I.I. MOD09A1 data   -> SMOD.MODIS.csv (20240211)  ----
  # I used the Global subset tools at https://modis.ornl.gov/history/
  # Citation 
  # [1] ORNL DAAC. 2018. Terrestrial Ecology Subsetting & Visualization Services (TESViS) Global Subsets Tool. ORNL DAAC, Oak Ridge, Tennessee, USA. Accessed February 11, 2024. Subset obtained for MOD09A1 product at 45.055N,6.385E, time period: 1999-12-31 to 2024-02-03, and subset size: 20.5 x 20.5 km. https://doi.org/10.3334/ORNLDAAC/1379
  
  # [2] E. Vermote. 2021. MOD09A1 MODIS/Terra Surface Reflectance 8-Day L3 Global 500m SIN Grid V061. NASA EOSDIS Land Processes DAAC. https://doi.org/10.5067/MODIS/MOD09A1.061
  
  # tif files  are stored in 
  setwd("D:/DATA_ROCHENOIRE/MOD09A1")
  RED   <- terra::rast(list.files(pattern="b01"))
  NIR   <- terra::rast(list.files(pattern="b02"))
  GREEN <- terra::rast(list.files(pattern="b04"))
  MIR6  <- terra::rast(list.files(pattern="b06"))
  
  NDVI.tmp <- (NIR-RED) / (NIR+RED) 
  NDSI.tmp <- (GREEN-MIR6)/(GREEN+MIR6)
  
  tmp.date        <- substr(names(RED),10,16)
  names(NDVI.tmp) <- names(NDSI.tmp) <- tmp.date
  MODIS.DATE      <- as.Date(tmp.date, format="%Y%j") 
  
  my3CONTwgs <- terra::project(my3CONT,"epsg:4326")

  COO <- as.data.frame(terra::crds(NDVI.tmp))
  MOD09.PTSwgs <- terra::vect(COO,geom=c("x","y"),crs="epsg:4326")
  MOD09.PTS <- terra::project(MOD09.PTSwgs, "epsg:2154")
  
  # select points that fall within my3CONT
  tmp <- terra::is.related(MOD09.PTSwgs,my3CONTwgs,"within")
  ID  <- which(tmp==TRUE)
  COO <- terra::crds(MOD09.PTSwgs[ID])
  SEL <- terra::cellFromXY(NDVI.tmp,COO)
  
  # graphical checking
  graphics.off()
  windows(10,10)
  plot(myCONTwgs)
  plot(MOD09.PTSwgs[ID],add=T,cex=.2)
  text(COO,labels=SEL)
  
  # compare with
  PT       <- c(6.4106,45.04128)
  tmp      <- matrix(c(6.38994,6.41061,45.05495,45.04128),2,2)
  (terra::cellFromXY(NDVI.tmp,tmp))

  ID.FXA = 1550; ID.NIV = 1293
  
  df.ndsi <- data.frame(
    DATE = MODIS.DATE,
    NDSI = t(terra::extract(NDSI.tmp,SEL))
  )
  df.ndvi <- data.frame(
    DATE = MODIS.DATE,
    NDVI = t(terra::extract(NDVI.tmp,SEL))
  )
  colnames(df.ndsi)<-colnames(df.ndsi)<-c("DATE",SEL)
  
  # Smoothing time series - inspired from NDVImet.r
  # set parameters & gap fill first missing values
  YEAR         <- 2000:2023
  nyear        <- length(YEAR)
  TMP          <- expand.grid(YEAR,substr(as.character(seq(1001,1365,8)),2,4))
  TIME.REF.jd  <- sort(paste(TMP[,1],TMP[,2],sep=''))
  TIME.REF.dat <- as.Date(TIME.REF.jd,format="%Y%j") 
  TIME.REF.jday<- substr(TIME.REF.jd,5,7) 
  TIME.REF.yr  <- as.numeric(substr(TIME.REF.jd,1,4))
  ndate        <- length(TIME.REF.dat)
  TIME.REF.jdX <- paste("X",TIME.REF.jd,sep="")
  SEQ.DAY      <- seq(as.Date(paste0("2000","-01-01")),as.Date(paste0("2023","-12-31")),1)
  npoint       <- length(SEL) 
  
  NDVI.tsf     <- NDSI.tsf <- data.frame(matrix(NA,ndate,npoint,dimnames=list(TIME.REF.jdX,SEL)))
  MATCH        <- na.omit(match(paste0("X",rownames(df.ndvi)),TIME.REF.jdX))
  NDVI.tsf[MATCH,]<- df.ndvi[1:1096,-1] # what is the issue here ?
  NDSI.tsf[MATCH,]<- df.ndsi[1:1096,-1]
  
  NDVI.tsf[1:10,] 
  
  for (i in 1:length(SEL)){
    NDVI.tsf[1:7,i]    <- mean(NDVI.tsf[8:10,i],na.rm=T) 
    NDSI.tsf[1:7,i]    <- mean(NDSI.tsf[8:10,i],na.rm=T)
  }
  plot(NDVI.tsf$"X1550",type="l")
  
  # CORRECTION
  DATA <- as.matrix(NDVI.tsf)
  SMOOTHf <- function(DATA){
    M         <- DATA
    NTIMES    <- nrow(DATA)
    thr=0.2; forw=3
    # B. Bise1
    M.DIFF      <- M[-1,]-M[-NTIMES,]
    M.DIFF2     <- rbind(NA,M.DIFF)
    M.DEC.ID    <- which(M.DIFF2<0)
    M.DEC.TH    <- rbind(NA,M[-1,] - thr*M.DIFF)
    
    M.VAL.SP    <- caTools::runmax(M,k=forw,align="left",endrule="NA") 
    
    M.LAW       <- M.VAL.SP[M.DEC.ID] - M.DEC.TH[M.DEC.ID]
    M.REJ.DEC   <- M.DEC.ID[which(M.LAW>0)]
    
    M[M.REJ.DEC]<- NA
    # M           <- zoo::na.spline(M,method="natural") 
    M           <- zoo::na.approx(M) 
    
    # B. Bise2
    M.DIFF      <- M[-1,]-M[-NTIMES,]
    M.DIFF2     <- rbind(NA,M.DIFF)
    M.DEC.ID    <- which(M.DIFF2<0)
    M.DEC.TH    <- rbind(NA,M[-1,] - thr*M.DIFF)
    
    M.VAL.SP    <- caTools::runmax(M,k=forw,align="left",endrule="NA") 
    
    M.LAW       <- M.VAL.SP[M.DEC.ID]- M.DEC.TH[M.DEC.ID]
    M.REJ.DEC   <- M.DEC.ID[which(M.LAW>0)]
    
    M[M.REJ.DEC]<- NA
    # M           <- zoo::na.spline(M,method="natural") 
    M           <- zoo::na.approx(M) 
    NDVI_MOD09f <- t(M)
    

  tmp1        <-  split(NDVI_MOD09f, seq(nrow(NDVI_MOD09f)))
  names(tmp1) <-  "FLUXALP"
  tmp2        <- lapply(tmp1,function(x) {x[ndate]<-x[ndate-1]; return(x)})
  tmp3        <- lapply(tmp2,sgolayfilt, p=3, n=7, m=0)
  
  DAYINT <- function(x,seq.day,match) {
    NDAY        <- length(seq.day)
    NDATE       <- length(x)
    tmp1        <- rep(NA,NDAY)
    tmp1[match] <- x
    Z           <- stats::spline(tmp1,n=NDAY)$y
    return(Z)
  }
  
  MATCH     <- match(TIME.REF.dat,SEQ.DAY)
  
  tmp4        <- lapply(tmp3, FUN=DAYINT, seq.day=SEQ.DAY,match=MATCH)
   return(tmp4)
  }
  
  NDVIcor <- SMOOTHf(as.matrix(NDVI.tsf))
  NDSIcor <- SMOOTHf(as.matrix(NDSI.tsf))
  
  SMODf  <- function(x) 50+which.min(abs(x[51:210]-0.2))
  # SMODf  <- function(x) 50+which.min(abs(x[51:210]-0.4))
  SMODest <- matrix(NA,50,24)
  for (i in 1:length(SEL)){
    print(i)
    SMODest[i,] <- aggregate(NDVIcor[[i]]-NDSIcor[[i]],list(lubridate::year(SEQ.DAY)),SMODf)$x
  }
  
  colnames(SMODest) <- 2000:2023
  rownames(SMODest) <- SEL
  
  setwd(WD); write.csv(SMODest,"SMOD.MODIS.csv",row.names=T)
  
  SMOD.MODIS <- read.csv("SMOD.MODIS.csv",row.names=1)
  
  setwd(WD); SMOD.FLUX <- read.csv("SMOD.FLUX.csv")
  
  cor(as.numeric(SMOD.MODIS["1550",14:24]),SMOD.FLUX$SMOD_ndvi) # 0.864

  graphics.off(); windows(10,10)
  plot(as.numeric(SMOD.MODIS["1550",14:24]),SMOD.FLUX$SMOD_ndvi,type="p",pch=21,xlim=c(110,160),ylim=c(110,160),bg="red",cex=.2) 
  abline(0,1)
  text(as.numeric(SMOD.MODIS["1550",14:24]),SMOD.FLUX$SMOD_ndvi,labels=2013:2022)

  graphics.off(); windows(15,7)
  plot(2000:2023,as.numeric(SMOD.MODIS["1293",])-as.numeric(SMOD.MODIS["1550",]),type="b",pch=21,bg="gray",cex=2,ylim=c(15,50))
  abline(h=Hmisc::smedian.hilow(as.numeric(SMOD.MODIS["1293",])-as.numeric(SMOD.MODIS["1550",]),0.75),lty=c(1,2,2))
  
  graphics.off(); windows(15,7)
  plot(2000:2023,as.numeric(SMOD.MODIS["1550",]),type="b",pch=21,bg="gray",cex=2,ylab="MODIS SMOD at Fluxalp") 
  abline(h=Hmisc::smedian.hilow(as.numeric(SMOD.MODIS["1550",]),0.75),lty=c(1,2,2))

  
  # I.J. IGN HISTORIQUE - REMONTER LE TEMPS -> IGN.fSCA.csv (20240227) ----
    # I.J.1. Process orthorectified images (20250315) ----
  LF      <- list.files("E:/DATA_ROCHENOIRE/BD_HISTORIQUE/OR_DATA/BINARY",pattern=".tif",full=T)
  LFs     <- list.files("E:/DATA_ROCHENOIRE/BD_HISTORIQUE/OR_DATA/BINARY",pattern=".tif")
  OR.DATE <- (as.Date(substr(LFs,15,24)))
  
  HIST.fSCA <- data.frame(
    DATE = OR.DATE, 
    YEAR = lubridate::year(OR.DATE),
    DOY  = lubridate::yday(OR.DATE)
  )
  
  # Compute the thermal time in DD at Fluxalp
  HIST.fSCA$DD <- NA
  setwd(WD); DD  <- read.csv("FLUXALP.DD.csv")
  
  for (i in 1:nrow(HIST.fSCA)){
    tmpYEAR <- which(colnames(DD)==paste0("X",HIST.fSCA$YEAR[i]))
    HIST.fSCA$DD[i] <- DD[HIST.fSCA$DOY[i],tmpYEAR]
  }
  
  # compute the fSCA, including the TOT area
  HIST.fSCA$fSCA.MAND <- HIST.fSCA$fSCA.LAUZ <- HIST.fSCA$fSCA.RNPI <- HIST.fSCA$fSCA.TOT<- NA
  
  for (i in 1:length(LF)){
    print(i)
    tmp      <- terra::rast(LF[i])
    crs(tmp) <-"EPSG:2154"
    
    # keep the full area to check
    # reproject over each polygon of my3CONT at 1m
    # assume that no snowy pixels in the truncated part
    tmp1     <- terra::mask(terra::project(tmp,my3REF01,method="near"),my3CONT)
    m        <- c(-Inf, 254, 0, 254, +Inf, 1)
    rcl      <- matrix(m, ncol=3, byrow=TRUE)
    tmp2     <- terra::classify(tmp1,rcl)
    setwd(WD); writeRaster(tmp2,paste0(HIST.fSCA$DATE[i],"_snowbin_1M.tif"),overwrite=T) 
    
    HIST.fSCA$fSCA.TOT[i] <- sum(tmp2[],na.rm=T)/sum(NCELL01.WS3)  
    
    for (j in 1:3){
      tmp1     <- terra::mask(terra::project(tmp,my3REF01,method="near"),my3CONT[j,])
      m        <- c(-Inf, 254, 0, 254, +Inf, 1)
      rcl      <- matrix(m, ncol=3, byrow=TRUE)
      tmp2     <- terra::classify(tmp1,rcl)
      HIST.fSCA[i,5+j] <- sum(tmp2[],na.rm=T)/NCELL01.WS3[j]   
    }
    
  }
  
  HIST.fSCA
  
  # save as longitudinal data to match other datasets
  HIST.fSCA$SRC<-"IGN"
  
  IGN <- reshape2::melt(HIST.fSCA[,-1],id.vars=c("YEAR","DOY","DD","SRC"),variable.name="WS",value.name="fSCA")
  
  # reorganize columns 
  IGN$WS <- substr(IGN$WS,6,9)
  IGN$fNA <- 0 # to be checked
  IGN <- IGN[,c("YEAR","DOY","DD","WS","fSCA","fNA","SRC")]
  
  setwd(WD);write.csv(IGN,"IGN.fSCA.csv",row.names=F)

  
    # I.J.2. Other images to orthorectify ----
  # photos à orthorectifier en complement du travail d'Arthur Provost
  LF <- list.files("D:/DATA_ROCHENOIRE/BD_HISTORIQUE/SRC_DATA/",full=F,pattern="jp2")
  HIST <- unique(as.Date(substr(LF,15,24)))
  HIST.DB <- data.frame(
                  DATE = HIST, 
                  YEAR = lubridate::year(HIST),
                  DOY  = lubridate::yday(HIST)
                          )
  setwd(WD);write.csv(HIST.DB,"HIST.DB.csv",row.names=F)
  
  # thermal time of aerial photographs
  setwd(WD); TMP <- read.csv("HIST.DB.csv")
  setwd(WD); DD  <- read.csv("FLUXALP.DD.csv")
  
  TMP$DD <- NA
  for (i in 1:nrow(TMP)){
    tmpYEAR <- which(colnames(DD)==paste0("X",TMP$YEAR[i]))
    TMP$DD[i] <- DD[TMP$DOY[i],tmpYEAR]
  }
  
  # test with 23-06-1945
  tmp <- terra::rast("D:/DATA_ROCHENOIRE/BD_HISTORIQUE/OR_DATA/IGNF_PVA_1-0__1945-06-23_orthorect.tif")
  tmp[tmp==0]<-NA
  graphics.off();windows(10,5);terra::plot(tmp)
  terra::plot(my3CONT,add=T)
  
  # transform into binary
  hist(tmp,breaks=seq(0,255,1),ylim=c(0,15000),xlim=c(100,180))
  
  THR = 122
  m <- c(-Inf, THR, 0,THR, +Inf, 1)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  tmpb <- classify(tmp, rclmat, include.lowest=TRUE)
  graphics.off();windows(15,15);terra::plot(tmpb)
  terra::plot(my3CONT,add=T)
  
  BIN45 <- terra::project(tmpb,my3REF01,method="near")
  BIN45_2m <- terra::aggregate(BIN45,fact=2,fun="modal")
  plot(BIN45_2m)
  
  # count the fSCA
  for (j in 1:3){
    tmp1     <- terra::mask(terra::project(tmpb,my3REF01,method="near"),my3CONT[j,])
    table(tmp1[],useNA="always")
    sum(tmp1[],na.rm=T)/NCELL01.WS3[j] # 30%    
  }
  
    # I.J.3. Estimate DD  ----
  setwd(WD); DD  <- read.csv("FLUXALP.DD.csv")
  DD[lubridate::yday(as.Date("1945-06-23")),"X1945"]
  DD[lubridate::yday(as.Date("1981-06-15")),"X1981"]
  DD[lubridate::yday(as.Date("1990-07-21")),"X1990"]
  DD[lubridate::yday(as.Date("2022-05-20")),"X2022"]
  DD[lubridate::yday(as.Date("2022-06-06")),"X2023"]
  DD[lubridate::yday(as.Date("2010-07-10")),"X2010"]
  DD[lubridate::yday(as.Date("2011-06-25")),"X2011"]
  DD[lubridate::yday(as.Date("2011-05-29")),"X2011"]
  
  # I.K. HR RapidEye & SPOT6,7 on data terra / dinamis catalog (20250120) ----
    # I.K.1. SEARCH IMAGES on dinamis / data terra ----
     # trial 1. Available 2011 RapidEye/pleiades images on Summer and ids-dinamis ----
     # cloudless images on 07-04 ; 29-05 ; 25-06 ; 
     # 5m resolution RapiEye product in Lambert93
     # can also be downloaded at DINAMIS https://ids-dinamis.data-terra.org/web/guest/catalog1
     # trial 2. Find one additional RE images on Summer ----
  # only one cloudless images matching crop -> on 07-07-2010
 LF <- list.files("S:/LECA/PLATEAU-PASTIS/GIS_DATA/Alpes/IMAGERY/RAPIDEYE/_05",pattern=".tif",full=T)
 TMP <- terra::crop(terra::rast(LF[13]),my3CONT)
 TMP1 <- terra::subset(TMP,3:1)
 MAX  <- max(terra::global(TMP1,max))
 TMP2 <- 255*terra::subset(TMP1,3:1)/MAX # need to rescale
 terra::plotRGB(TMP2)
 RGB <- terra::RGB(TMP2, value=1:3,type="rgb")
 RGBsum <- terra::app(RGB,sum)
 terra::plot(RGBsum,legend=T)
 
 THR = 300
 m <- c(-Inf, THR, 0,THR, +Inf, 1)
 rclmat <- matrix(m, ncol=3, byrow=TRUE)
 tmpb2 <- terra::classify(RGBsum, rclmat);
 plot(tmpb2) # working

 setwd(DATA_WD);terra::writeRaster(tmpb2,"REbin_20100707.tif",overwrite=T)
  
     # trial 3. Additionnal images on dinamis / data terra ----
     # download at https://catalogue-dinamis.data-terra.org/
     # profile mis à jour sur https://sso.theia-land.fr/theia/common/success.xhtml
    # RapidEye bands correspond to 
      # Blue	  440 - 510
      # Green	  520 - 590
      # Red	    630 - 685
      # RedEdge	690 - 730
      # NIR	    760 - 850
  
     # I.K.2. Prepare binary snow cover maps (20250120) ----
  LF   <- list.files("E:/DATA_RAPIDEYE",pattern=".TIF",full=T)
  TMP  <- terra::crop(terra::rast(LF[1]),my3CONT)
  TMP1 <- terra::subset(TMP,3:1)
  MAX  <- max(terra::global(TMP1,max))
  TMP2 <- 255*terra::subset(TMP1,3:1)/MAX # need to rescale
  terra::plotRGB(TMP2)
  RGB <- terra::RGB(TMP2, value=1:3,type="rgb")
  # tmp <- terra::colorize(RGB, to="col")
  graphics.off();windows(10,10);terra::plotRGB(RGB)
  graphics.off();windows(10,10);terra::plot(subset(RGB,1))
  graphics.off();windows(10,10);terra::plot(subset(RGB,2))
  
  RGBsum <- terra::app(RGB,sum)
  windows(10,10);plotRGB(RGB)
  hist(RGBsum,breaks=seq(0,255*3,5))

  THR = 320
  m <- c(-Inf, THR, 0,THR, +Inf, 1)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  tmpb <- classify(RGBsum, rclmat, include.lowest=TRUE)
  windows(15,15);terra::plot(tmpb,col=c("white","blue"))
  terra::plot(my3CONT,add=T)
  
  setwd(DATA_WD);terra::writeRaster(tmpb,"REbin_20110407.tif")
  
  # Calcul du NDVI de RapidEye avec les bandes 3 (red) et 5 (NIR)
  NDVI_re <- (TMP[[5]]-TMP[[3]])/(TMP[[5]]+TMP[[3]])
  graphics.off();windows(10,10);plot(NDVI_re);plot(my3CONT,add=T)
  
    # I.K.3. Estimate DD (FluxAlp) & fSCA (20250120) ----
  DATE.HR <- as.Date(c("2010-07-07","2011-04-07","2011-05-29","2011-06-25","2018-06-28","2018-06-23","2016-07-07", "2016-07-19","2019-07-18","2019-06-26"))
  # NOTE : add 2024 DD on Fluxalp to process "2024-07-17" !
  
  df <- data.frame(DATE   = DATE.HR,
                   YEAR   = lubridate::year(DATE.HR),
                   DOY    = lubridate::yday(DATE.HR),
                   SOURCE = c("RE","RE","RE","RE","SPOT7","SPOT7","SPOT7","SPOT7","SPOT7","SPOT6"),
                   DD     = NA,   # degree days estimated at Fluxalp
                   SCA    = NA
  )
  
  setwd(WD); DD  <- read.csv("FLUXALP.DD.csv")
  
  for (i in 1:nrow(df)){
    print(i)
    tmpYEAR   <- which(colnames(DD)==paste0("X",df$YEAR[i]))
    df$DD[i]  <- round(DD[df$DOY[i],tmpYEAR])
    
    FILEbin      <- paste0(df$SOURCE[i],"bin_",gsub("-","",df$DATE[i]),".tif")
    setwd(DATA_WD) 
      tmp        <- try(terra::rast(FILEbin)) 
      if(!inherits(tmp, "try-error")){
        crs(tmp)   <-"epsg:2154"
        tmp1       <- terra::mask(tmp,my3CONT)
        land       <- length(which(terra::extract(tmp,my3CONT)$sum==0))
        snow       <- length(which(terra::extract(tmp1,my3CONT)$sum==1))
        df$SCA[i]  <- round(100*snow/(land+snow),1)
      }
      
  }

  setwd(DATA_WD);write.csv(df,"METADATA.RE.SPOT.csv")
  setwd(DATA_WD);tmp <- read.csv("METADATA.RE.SPOT.csv")
 
  # I.L. SPOT IRC for illustration ----
    tmp <- terra::rast("D:/DATA_ROCHENOIRE/SPOT/50512591109131049392J0.TIF")
    tmp1 <- terra::crop(tmp,my2CONT)
    graphics.off();windows(8,10)
    terra::plotRGB(tmp1)
    plot(my2CONT,add=T,border="yellow")
 
  # I.M. BD ORTHO at 50cm resampled at 2m ----
    # see rochenoire.r script part I.E
    # need to solve the issue of histogram matching between the two dpts
    setwd('D:/DATA_ROCHENOIRE/MY_BD_ORTHO/IRC_2m/')
    TMP05 <- terra::crop(terra::rast('RNED_ORTHO_D05_LA93-IRC_2m.tif'),my3EXT)
    TMP73 <- terra::crop(terra::rast('RNED_ORTHO_D73_LA93-IRC_2m.tif'),my3EXT)
    
    graphics.off();windows(10,15);terra::plotRGB(TMP73,axes=T)
    plot(my3CONT,add=T,border="yellow",lwd=2)
    
    graphics.off();windows(10,10);terra::plotRGB(TMP05,axes=T)
    plot(my3CONT,add=T,border="yellow",lwd=2)

    
    # define the refmask and the xmask by defining areas of NA
    MASK73 <- TMP73
    coo.tmp <- terra::crds(MASK73)
    CELL <- cellFromXY(MASK73,coo.tmp) 
    COND <- which(coo.tmp[,2]>6446000 & coo.tmp[,2]<6447000 & coo.tmp[,1]>966000 & coo.tmp[,1]<967000)  # remove the cloudy area in TMP73
    MASK73[CELL[COND]]<-NA
    MASK73[CELL[which(coo.tmp[,2]<6445000)]]<-NA
    MASK73[CELL[which(coo.tmp[,2]<6446000 & coo.tmp[,1]<966000)]]<-NA
    terra::plotRGB(MASK73,axes=T)
    
    UNI73 <- TMP73
    UNI73[CELL[which(coo.tmp[,2]<6447000 | coo.tmp[,1]<966000)]]<-NA
    graphics.off();terra::plotRGB(UNI73,axes=T)
    
    MASK05 <- TMP05
    MASK05u <- TMP05
    coo.tmp <- terra::crds(MASK05)
    CELL <- cellFromXY(MASK05,coo.tmp) 
    MASK05[CELL[which(coo.tmp[,2]>6447000 & coo.tmp[,1]>966000)]]<-NA
    MASK05u[CELL[which(coo.tmp[,2]<6447000 | coo.tmp[,1]<966000)]]<-NA
    
    terra::plotRGB(MASK05,axes=T)
    terra::plotRGB(MASK05u,axes=T)
    
    # good one ! hist match on the two areas of TMP05
    TMP <- RStoolbox::histMatch(x=MASK05,ref=MASK73,intersectOnly=T,nSamples=10^6)
    terra::plotRGB(TMP)
    
    TMPu <- RStoolbox::histMatch(x=MASK05u,ref=MASK73,intersectOnly=T,nSamples=10^5)
    terra::plotRGB(TMPu)
    
    # or all domain
    TMPt <- RStoolbox::histMatch(x=TMP05,ref=MASK73,intersectOnly=T,nSamples=10^6)
    terra::plotRGB(TMPt)
    
    RES <- terra::mosaic(TMP,TMPu,TMPt,MASK73,fun="mean")
    terra::plotRGB(RES)
    setwd(WD); terra::writeRaster(RES,"BD_ORTHO_IRC.tif")
    
    # NDVI x FLAT (20230331) -----
    setwd(WD); RES <- terra::rast("BD_ORTHO_IRC.tif")
    crs(RES)<-"epsg:2154"
    egi <- (2*RES[[1]]-RES[[2]]-RES[[3]])/(RES[[1]]+RES[[2]]+RES[[3]])
    # IRC-R / IRC+R
    # Bande 1 : infrarouge. Bande 2 : rouge. Bande 3 : vert
    NDVI <- (RES[[1]]-RES[[2]])/(RES[[1]]+RES[[2]])
    terra::plot(NDVI)
    plot(NDVI,breaks=seq(-0.6,0.6,0.01),legend=F) 
    plot(my3CONT,add=T,border="black",lwd=2)
    
    setwd(WD); SLO <- terra::rast("SLO_2M.tif")

    THR = 0.4
    m <- c(-Inf, THR, 1, THR, +Inf, NA)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    FLATb <- classify(SLO, rclmat, include.lowest=TRUE)
    terra::plot(FLATb)
    
    THR = 0.02
    m <- c(-Inf, THR, 1,THR, +Inf, NA)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    NDVIb <- classify(NDVI, rclmat, include.lowest=TRUE)
    terra::plot(NDVIb)
    
    plot(NDVIb*FLATb*BIN45_2m)
    plot(my3CONT,add=T,border="black",lwd=2)


    # other trials 
    # # still a difference in upper right
    # RES1 <- terra::mosaic(TMP,TMP73,fun="first")
    # plotRGB(RES1)
    # # still a difference in upper right
    # RES1 <- terra::mosaic(TMP,TMP73,fun="first")
    # plotRGB(RES1)
    # TMQ <- RStoolbox::histMatch(x=TMP05,ref=MASK73,intersectOnly=T)
    # terra::plotRGB(TMQ)
    # RES2 <- terra::mosaic(TMQ,TMP73,fun="first")
    # plotRGB(RES2)
    
  # SAVE ----
  setwd(WD); save(list=c("mySIT","mySITsim","myEXT","myRWBODY","colorstrip"),file="SNOW_ROCHENOIRE.RData") 
  

# II. FORCAGE METEO & WDD ESTIMATES  (20250320) ----
    # II.A. Prepare data sources XXXraw.csv files (20250316) ----
      # YEAR + MONTH + DOY + TEMP + WDD 
      # including na.spline and daily desaggregation if relevant 
      # HISraw : Histalp (20250125) ----
      # Ref is Auer up to 2010 only - delete the last four years
      # Absolute temperatures is given at 5'
    setwd("~/CLIMATOmc/DATA/DATA_HISTALP/")
    tmp   <- terra::rast("HISTALP_temperature_1780-2014.nc")
    tmp0  <- data.frame(X_wgs84=c(6.38994,6.41061),Y_wgs84=c(45.05495,45.04128))
    tmp1  <- terra::vect(tmp0,geom=c("X_wgs84","Y_wgs84"),crs="epsg:4326")
    tmp2  <- t(terra::extract(tmp,tmp1))[-(1:3),1]
    
    length(1780:2014)*12 #  extract monthly values
    
    st <- ymd("1780-01-15")
    en <- ymd("2014-12-15")
    YM <- st %m+% months(seq(0, round(interval(st, en) / months(1)), 1))
    
    M.HIST <- data.frame(YEAR  = lubridate::year(YM),
                         MONTH = lubridate::month(YM),
                         DOY   = lubridate::yday(YM),
                         TEMP  = tmp2
    )
    
    # delete years after 2010
    M.HIST <- M.HIST[-which(M.HIST$YEAR>2010),]
    
    # Disaggregate to daily data 
    # using Harmonic or Discrete Fourier Transform
    
    which(M.HIST$YEAR==1780)[1]
    nrow(M.HIST)
    
    # prepare the data.frame
    df <- M.HIST[which(M.HIST$YEAR==1780)[1]:nrow(M.HIST),]
    
    # add the corresponding date
    df$DATE <- seq(as.Date('1780-01-15'),as.Date('2010-12-15'),'month')
    # order of day
    seq.REF <- seq(as.Date('1780-01-01'),as.Date('2010-12-31'),'day')
    length(seq.REF)
    
    df$ORD <- match(df$DATE,seq.REF)
    
    graphics.off();windows(20,10);plot(df$DATE,df$TEMP,type="l")
    
    x=100+df$TEMP; n=nrow(df)/2; up=30;plot=TRUE;add=FALSE;main=NULL
    # using Fast Fourier Transform
    
    # The direct transformation
    # The first frequency is DC, the rest are duplicated
    
    dff = stats::fft(x)
    # The time
    # t = seq(from = 1, to = length(x))
    t = df$ORD
    
    # Upsampled time
    # nt = seq(from = 1, to = length(x)+1-1/up, by = 1/up)
    nt = df$ORD[1]:df$ORD[nrow(df)]
    
    #New spectrum
    ndff    = array(data = 0, dim = c(length(nt), 1L))
    ndff[1] = dff[1] #Always, it's the DC component
    
    if(n != 0){
      #The positive frequencies always come first
      ndff[2:(n+1)] = dff[2:(n+1)] 
      
      #The negative ones are trickier
      ndff[length(ndff):(length(ndff) - n + 1)] = dff[length(x):(length(x) - n + 1)]
    }
    
    #The inverses
    indff = fft(ndff/nrow(df), inverse = TRUE)
    idff  = fft(dff/nrow(df), inverse = TRUE)
    
    if(plot){
      if(!add){
        plot(x = t, y = x, pch = 16L, xlab = "Time", ylab = "Measurement",
             main = ifelse(is.null(main), paste(n, "harmonics"), main))
        lines(y = Mod(idff), x = t, col = "green")
      }
      lines(y = Mod(indff), x = nt,col="red")
    }
    
    SEQ <- as.Date(seq.REF[df$ORD[1]:df$ORD[nrow(df)]])
    HISraw = data.frame(DATE  = SEQ,
                        YEAR  = lubridate::year(SEQ),
                        MONTH = lubridate::month(SEQ),
                        DOY   = lubridate::yday(SEQ),
                        TEMP  = Mod(indff)-100 
                        )
    HISraw$WDD <- ifelse(HISraw$TEMP>0,HISraw$TEMP,0)
    head(HISraw)
    # was named ret before
    
    # graphical checking
    plot(SEQ[1:2000],HISraw$TEMP[1:2000],type='l',col="red")
    points(df$DATE,df$TEMP,pch=21,bg="green")
    abline(h=0)
    
    setwd(WD); write.csv(HISraw,"Histalp.Daily.csv")
    

      # MONraw : weather station La Monêtier les Bains ----
      # see climato.r section III.1 for the updated production of data
  setwd("~/CLIMATOmc/DATA/CLIMATHEQUE/DAILY_MEANS/")

  # BESSE en OISANS -> BESraw
  BES      <- read.csv("BESSE_METEO.csv")    # TM2 depuis 1959; RR since 2024
  BES$DATE <- as.Date(BES$DATE); 
  BES$TM2  <- 0.5*(BES$TX+BES$TN)  
  BESraw   <- data.frame(DATE  = BES$DATE,
                         YEAR  = lubridate::year(BES$DATE),
                         MONTH = lubridate::month(BES$DATE),
                         DOY   = lubridate::yday(BES$DATE),
                         TEMP  = BES$TM2
  )
  BESraw$WDD <- ifelse(BESraw$TEMP>0,BESraw$TEMP,0)
  head(BESraw)
  
  # check the number of missing data up to DOY212 and gap-fill if less than 10 NA
  ID <- which(BESraw$DOY<=212)
  tmp <- aggregate(BESraw$TEMP[ID],list(YEAR=BESraw$YEAR[ID]),function(x) length(which(is.na(x)==TRUE)))
  GAPf <- tmp[which(tmp$x<=10 & tmp$x>0),"YEAR"]
  # use NA spline for years 1939,1981,1982,1984
  for (i in GAPf){
    print(i)
    VEC <- zoo::na.spline(BESraw[which(BESraw$YEAR==i),"TEMP"][1:212])
    BESraw[which(BESraw$YEAR==i),"TEMP"][1:212] <- VEC
  }   
  
  CHK <- aggregate(BESraw$TEMP[ID],list(YEAR=BESraw$YEAR[ID]),function(x) length(which(is.na(x)==TRUE)))
  barplot(CHK$x,names=CHK$YEAR,las=2)

  # CHAMONIX -> CHAraw
  CHA       <- read.csv("CHAMONIX_METEO.csv") 
  CHA$DATE  <- as.Date(CHA$DATE)
  CHA$TM2   <- 0.5*(CHA$TX+CHA$TN)
  plot(CHA$DATE,CHA$TM2,type="l")
  CHAraw    <- data.frame(DATE  = CHA$DATE,
                       YEAR  = lubridate::year(CHA$DATE),
                       MONTH = lubridate::month(CHA$DATE),
                       DOY   = lubridate::yday(CHA$DATE),
                       TEMP  = CHA$TM2
  )
  CHAraw$WDD <- ifelse(CHAraw$TEMP>0,CHAraw$TEMP,0)
  head(CHAraw)
  
  # check the number of missing data up to DOY212 and gap-fill if less than 10 NA
  ID <- which(CHAraw$DOY<=212)
  tmp <- aggregate(CHAraw$TEMP[ID],list(YEAR=CHAraw$YEAR[ID]),function(x) length(which(is.na(x)==TRUE)))
  GAPf <- tmp[which(tmp$x<=10 & tmp$x>0),"YEAR"]
  # use NA spline for years 1939,1981,1982,1984
  for (i in GAPf){
    print(i)
    VEC <- zoo::na.spline(CHAraw[which(CHAraw$YEAR==i),"TEMP"][1:212])
    CHAraw[which(CHAraw$YEAR==i),"TEMP"][1:212] <- VEC
  }   
  
  CHK <- aggregate(CHAraw$TEMP[ID],list(YEAR=CHAraw$YEAR[ID]),function(x) length(which(is.na(x)==TRUE)))
  barplot(CHK$x,names=CHK$YEAR,las=2)
  
  # LE MONETIER les BAINS -> MONraw
  MON <- read.csv("MONETIER_METEO.csv") # TM2 depuis 1936; RR since 2024
  MON$DATE <- as.Date(MON$DATE); 
  MON$TM2 <- 0.5*(MON$TX+MON$TN)
  
    MONraw <- data.frame(DATE  = MON$DATE,
                         YEAR  = lubridate::year(MON$DATE),
                         MONTH = lubridate::month(MON$DATE),
                         DOY   = lubridate::yday(MON$DATE),
                         TEMP  = MON$TM2
  )
  MONraw$WDD <- ifelse(MONraw$TEMP>0,MONraw$TEMP,0)
  head(MONraw)
  
  # check the number of missing data up to DOY212 and gap-fill if less than 10 NA
  ID <- which(MONraw$DOY<=212)
  tmp <- aggregate(MONraw$TEMP[ID],list(YEAR=MONraw$YEAR[ID]),function(x) length(which(is.na(x)==TRUE)))
  GAPf <- tmp[which(tmp$x<=10 & tmp$x>0),"YEAR"]
  # use NA spline for years 1939,1981,1982,1984
  for (i in GAPf){
    print(i)
    VEC <- zoo::na.spline(MONraw[which(MONraw$YEAR==i),"TEMP"][1:212])
    MONraw[which(MONraw$YEAR==i),"TEMP"][1:212] <- VEC
  }    
  
  CHK <- aggregate(MONraw$TEMP[ID],list(YEAR=MONraw$YEAR[ID]),function(x) length(which(is.na(x)==TRUE)))
  barplot(CHK$x,names=CHK$YEAR,las=2)
  
  # use another weather station / of S2M to gap-fill years 1983 and 1986 ??
  
  # Discarded weather stations
  VALLOIRE <- read.csv("VALLOIRE_METEO.csv") # TM2 depuis 1992 only -> do not use
  VALLOIRE$DATE <- as.Date(VALLOIRE$DATE); VALLOIRE$TM2 <- 0.5*(VALLOIRE$TX+VALLOIRE$TN)
  
  # OTHER weather stations / a bit farther but probably useful to gap fill monetier ?
  ANNECY   <- read.csv("ANNECY_METEO.csv")  
  ANNECY$TM2 <- 0.5*(ANNECY$TX+ANNECY$TN) ; ANNECY$DATE <- as.Date(ANNECY$DATE)
  plot(ANNECY$DATE,ANNECY$TM2,type="l") # no recent data
  
      # FXAraw : FluxAlp climate data ----
      # need to complete with 2024 !
  setwd("~/PROJETSmc/FLUXALP/DATA")
  FXA       <- read.csv("FLUXALP.daily.csv")
  FXA$DATE  <- as.Date(FXA$DATE)
  FXA$TM2   <- 0.5*(FXA$TX+FXA$TN)
  
  FXAraw   <- data.frame(DATE  = FXA$DATE,
                         YEAR  = lubridate::year(FXA$DATE),
                         MONTH = lubridate::month(FXA$DATE),
                         DOY   = lubridate::yday(FXA$DATE),
                         TEMP  = FXA$TM2
  )
  FXAraw$WDD <- ifelse(FXAraw$TEMP>0,FXAraw$TEMP,0)
  head(FXAraw)
  
  # check the number of missing data up to DOY212 and gap-fill if less than 10 NA
  ID  <- which(FXAraw$DOY<=212)
  tmp <- aggregate(FXAraw$TEMP[ID],list(YEAR=FXAraw$YEAR[ID]),function(x) length(which(is.na(x)==TRUE)))
  GAPf <- tmp[which(tmp$x<=10 & tmp$x>0),"YEAR"] # no gap

  CHK <- aggregate(FXAraw$TEMP[ID],list(YEAR=FXAraw$YEAR[ID]),function(x) length(which(is.na(x)==TRUE)))
  barplot(CHK$x,names=CHK$YEAR,las=2)
  
  
      # GALraw : Nivose Galibier climate data ----
  setwd("~/CLIMATOmc/DATA/CLIMATHEQUE/DAILY_MEANS/")
  GAL <- read.csv("NIVOGALI_METEO.csv")  
  GAL$DATE  <- as.Date(GAL$DATE)
  GAL$TM2   <- 0.5*(GAL$TX+GAL$TN)
  
  GALraw   <- data.frame(DATE  = GAL$DATE,
                         YEAR  = lubridate::year(GAL$DATE),
                         MONTH = lubridate::month(GAL$DATE),
                         DOY   = lubridate::yday(GAL$DATE),
                         TEMP  = GAL$TM2
  )
  GALraw$WDD <- ifelse(GALraw$TEMP>0,GALraw$TEMP,0)
  head(GALraw)
  
  # check the number of missing data up to DOY212 and gap-fill if less than 10 NA
  ID  <- which(GALraw$DOY<=212)
  tmp <- aggregate(GALraw$TEMP[ID],list(YEAR=GALraw$YEAR[ID]),function(x) length(which(is.na(x)==TRUE)))
  GAPf <- tmp[which(tmp$x<=10 & tmp$x>0),"YEAR"] # no gap
  
  for (i in GAPf){
    print(i)
    VEC <- zoo::na.spline(GALraw[which(GALraw$YEAR==i),"TEMP"][1:212])
    GALraw[which(GALraw$YEAR==i),"TEMP"][1:212] <- VEC
  }    
  
  
  CHK <- aggregate(GALraw$TEMP[ID],list(YEAR=GALraw$YEAR[ID]),function(x) length(which(is.na(x)==TRUE)))
  barplot(CHK$x,names=CHK$YEAR,las=2)
  
  # Other NIVOSE for SNOW
  NIVOROCH <- read.csv("NIVOROCH_METEO.csv")
  NIVOECRI <- read.csv("NIVOECRI_METEO.csv")

    
      # SAFraw : S2M for the best correlated combination with MONraw
    
    # II.B. Quantile Regression / Mapping Analysis (20250316) ----
     # The initial script was developed by Alix Reverdy, PhD, quantile regression on historical temperatures
  # then we realized it should be done using quantile mapping
     # The quantile function of a scalar random variable Y is the inverse of its distribution function. Like the distribution function, the quantile function provides a complete description of the statistical properties of the random variable.
  # finally quantile mappinggive worst results than regression quantile
      # Adjusting MONraw to FXAraw -> MONirq ----
    MONirq <- MONraw%>%left_join(FXAraw,by=c("DATE","YEAR","MONTH","DOY"),suffix=c(".MON",".FXA"))
    MONirq$TEMP.FXAoverMON.rq <- NA
    MONirq$TEMP.FXAoverMON.qm <- NA
    MONirq$TEMP.FXAoverMON.lm <- NA
    
    qrfit <- list(length=length(unique(MONirq$MONTH)))
    qmfit <- list(length=length(unique(MONirq$MONTH)))
    lmfit <- list(length=length(unique(MONirq$MONTH)))
    
    for(i in 1:length(unique(MONirq$MONTH))){
      print(i)
      tmp <- MONirq[MONirq$MONTH==i,]
      
      # Quantile regression
      qrfit[[i]] <- rq(TEMP.FXA~TEMP.MON,data = tmp)
      MONirq[MONirq$MONTH==i,"TEMP.FXAoverMON.rq"] <- predict(qrfit[[i]],data.frame(TEMP.MON=MONirq[MONirq$MONTH==i,"TEMP.MON"]))
      
      # Quantile mapping
      TBF  <- which(!is.na(tmp$TEMP.MON)==TRUE)
      tmp1 <- tmp[TBF,]
      qmfit[[i]] <- fitQmap(obs=tmp1$TEMP.MON,mod=tmp1$TEMP.FXA,qstep=0.1,method="QUANT",nboot=1,wet.day=FALSE)
      tmp1$qm.a <- doQmap(tmp1$TEMP.MON,qmfit[[i]],type="tricub",method="QUANT")
      
      MONirq[MONirq$MONTH==i,"TEMP.FXAoverMON.qm"][TBF] <- tmp1$qm.a
      
      # REG <- 2100:2635
      # plot(REG,tmp$qm.a[REG],type="l")
      # lines(REG,tmp$TEMP.FXA[REG],col="red")
      
      # Linear regression (delete)
      lmfit[[i]]=lm(TEMP.FXA~TEMP.MON,data = tmp)
      MONirq[MONirq$MONTH==i,"TEMP.FXAoverMON.lm"]<- predict(lmfit[[i]],data.frame(TEMP.MON=MONirq[MONirq$MONTH==i,"TEMP.MON"]))
    }
    
    # checking and perforamce of the different model
    plot(MONirq$TEMP.FXAoverMON.rq,MONirq$TEMP.FXA)
    plot(MONirq$TEMP.FXAoverMON.qm,MONirq$TEMP.FXA)
    plot(MONirq$TEMP.FXAoverMON.lm,MONirq$TEMP.FXA)
    
    TMP.rq <- MONirq[,c("TEMP.FXAoverMON.rq","TEMP.FXA")]%>%drop_na()
    m.eval(TMP.rq$TEMP.FXAoverMON.rq,TMP.rq$TEMP.FXA)
    TMP.qm <- MONirq[,c("TEMP.FXAoverMON.qm","TEMP.FXA")]%>%drop_na()
    m.eval(TMP.qm$TEMP.FXAoverMON.qm,TMP.qm$TEMP.FXA)
    TMP.lm <- MONirq[,c("TEMP.FXAoverMON.lm","TEMP.FXA")]%>%drop_na()
    m.eval(TMP.lm$TEMP.FXAoverMON.lm,TMP.lm$TEMP.FXA)
    
      # Adjusting HISraw to MONirq -> HISirq ----
    HISirq <- HISraw%>%left_join(MONirq,by=c("DATE","YEAR","MONTH","DOY"),suffix=c(".MON",".FXA"))
    
    HISirq$TEMP.MONoverHIS.rq <- NA
    HISirq$TEMP.MONoverHIS.qm <- NA
    HISirq$TEMP.MONoverHIS.lm <- NA
    
    qrfit <- list(length=length(unique(MONirq$MONTH)))
    qmfit <- list(length=length(unique(MONirq$MONTH)))
    lmfit <- list(length=length(unique(MONirq$MONTH)))
    
    for(i in 1:length(unique(HISirq$MONTH))){
      print(i)
      tmp <- HISirq[HISirq$MONTH==i,]
      
      # Quantile regression
      qrfit[[i]] <- rq(TEMP.FXAoverMON.rq~TEMP,data = tmp)
      HISirq[HISirq$MONTH==i,"TEMP.MONoverHIS.rq"] <- predict(qrfit[[i]],data.frame(TEMP=HISirq[HISirq$MONTH==i,"TEMP"]))
      
      # Quantile mapping
      TBF  <- which(!is.na(tmp$TEMP)==TRUE)
      tmp1 <- tmp[TBF,]
      qmfit[[i]] <- fitQmap(obs=tmp1$TEMP,mod=tmp1$TEMP.FXAoverMON.rq,qstep=0.1,method="QUANT",nboot=1,wet.day=FALSE)
      tmp1$qm.a <- doQmap(tmp1$TEMP,qmfit[[i]],type="tricub",method="QUANT")
      
      HISirq[HISirq$MONTH==i,"TEMP.MONoverHIS.qm"][TBF] <- tmp1$qm.a
      
      # Linear regression
      lmfit[[i]]=lm(TEMP.FXAoverMON.rq~TEMP,data = tmp)
      HISirq[HISirq$MONTH==i,"TEMP.MONoverHIS.lm"]<- predict(lmfit[[i]],data.frame(TEMP=HISirq[HISirq$MONTH==i,"TEMP"]))
    }
    
    # checking and perforamce of the model
    plot(HISirq$TEMP.MONoverHIS.rq,HISirq$TEMP.FXAoverMON.rq)
    abline(0,1)
    
    TMP.rq <- HISirq[,c("TEMP.MONoverHIS.rq","TEMP.FXAoverMON.rq")]%>%drop_na()
    m.eval(TMP.rq$TEMP.MONoverHIS.rq,TMP.rq$TEMP.FXAoverMON.rq)
    
    TMP.qm <- HISirq[,c("TEMP.MONoverHIS.qm","TEMP.FXAoverMON.qm")]%>%drop_na()
    m.eval(TMP.qm$TEMP.MONoverHIS.qm,TMP.qm$TEMP.FXAoverMON.qm)
    
    TMP.lm <- HISirq[,c("TEMP.MONoverHIS.lm","TEMP.FXAoverMON.lm")]%>%drop_na()
    m.eval(TMP.lm$TEMP.MONoverHIS.lm,TMP.lm$TEMP.FXAoverMON.lm)
    
    # II.C. Estimate Warming Degree Days (20250316) ----
    # FXAwdd
    myret  <- FXAraw[,c("YEAR","DOY","WDD")]
    tmp    <- t(reshape::cast(myret,YEAR~DOY,value="WDD",add.missing=T))
    FXAwdd <- collapse::fcumsum(tmp)
    plot(FXAwdd[,"2016"])
    setwd(WD); write.csv(FXAwdd,"FXAwdd.csv",row.names=F)
    
    # MONwdd
    myret     <- MONirq[,c("YEAR","DOY","TEMP.FXAoverMON.rq")]
    myret$WDD <- ifelse(myret$TEMP.FXAoverMON.rq<=0,0,myret$TEMP.FXAoverMON.rq) 
    myret     <- myret[,c("YEAR","DOY","WDD")]
    tmp       <- t(reshape::cast(myret,YEAR~DOY,value="WDD",add.missing=T))
    MONwdd    <- collapse::fcumsum(tmp)
    plot(MONwdd[,"2016"])
    setwd(WD); write.csv(MONwdd,"MONwdd.csv",row.names=F)
    
    # HISwdd
    myret     <- HISirq[,c("YEAR","DOY","TEMP.MONoverHIS.rq")]
    myret$WDD <- ifelse(myret$TEMP.MONoverHIS.rq<=0,0,myret$TEMP.MONoverHIS.rq) 
    myret     <- myret[,c("YEAR","DOY","WDD")]
    # add the 14 missing days in 1780
    myret <-  rbind(data.frame(YEAR=1780,DOY=1:14,WDD=0),myret)
    tmp    <- t(reshape::cast(myret,YEAR~DOY,value="WDD",add.missing=T))
    HISwdd <- collapse::fcumsum(tmp)
    
    # Complete series for recent years using MONdd 2010->2023
    HISwdd <- HISwdd[,-which(colnames(HISwdd)=="2010")]
    BEG    <- which(colnames(MONwdd)=="2010")
    END    <- which(colnames(MONwdd)=="2024")
    HISwdd <- cbind(HISwdd,MONwdd[,BEG:END])
    setwd(WD); write.csv(HISwdd,"HISwdd.csv",row.names=F)
    
    test1  <- read.csv("FXAwdd.csv")
    colnames(test1) <- gsub("X", "",  colnames(test1))
    
    test2  <- read.csv("MONwdd.csv")
    colnames(test2) <- gsub("X", "",  colnames(test2))
    
    test3  <- read.csv("HISwdd.csv")
    colnames(test3) <- gsub("X", "",  colnames(test3))
    
    plot(test2[,"2006"],type="l")
    lines(test3[,"2006"],type="l",col="red")
    
    # DEPRECATED Predict FLUXALP TM using weather stations (20240226) ----
      # NAIVE LINEAR MODEL -> FORCAGE_METEO.csv ----
      # adjust to the longest time series which is MONETIER
  
  library(dplyr); library(tidyr)
  TMP1 <- MONETIER %>%left_join(BESSE,by="DATE")
  TMP2 <- MONETIER %>%left_join(VALLOIRE,by="DATE")
  TMP3 <- MONETIER %>%left_join(FLUXALP,by="DATE")
  TMP4 <- MONETIER %>%left_join(NIVOGALI,by="DATE")
  
  TMP <- cbind(TMP1[,c("DATE","TM2.x","TM2.y")],TMP2[,"TM2.y"],TMP3[,"TM2.y"],TMP4[,"TM2.y"])
  colnames(TMP)<- c("DATE","TM_MONETIER","TM_BESSE","TM_VALLOIRE","TM_FLUXALP","TM_NIVOGALI")
  TMP$MONTH <- lubridate::month(TMP$DATE)
  TMP$YEAR  <- lubridate::year(TMP$DATE)
  
  cor(TMP[,2:6],use="complete.obs")
  graphics.off()
  windows(10,10)
  plot(TMP[,2],TMP[,3],cex=2,pch=".")  
  
  
  # monthly correlation
  COR.l <- list()
  for (i in 1:12) {
    COR.l[[i]] <- cor(TMP[which(TMP$MONTH==i),2:6],use="complete.obs")
  }
  
  # FLUXALP ~ MONETIER
  FLUX2MONE <- unlist(lapply(COR.l,function(x) x[4,1]))
  FLUX2VALL <- unlist(lapply(COR.l,function(x) x[3,1]))
  FLUX2BESS <- unlist(lapply(COR.l,function(x) x[2,1]))
  
  graphics.off();windows(10,10)
  plot(0,0,xlim=c(1,12),ylim=c(0.6,1),type="n",xlab="month")
    lines(FLUX2MONE,type="b",pch=21,bg="gray",cex=2)
    lines(FLUX2VALL,type="b",pch=21,bg="black")
    lines(FLUX2BESS,type="b",pch=21,bg="red",cex=2)  

    # develop a predictive model per month
    LM0.l <- LM1.l <- LM2.l <- LM3.l <- list()
    for (i in 1:12) {
    LM0.l[[i]] <- lm(TM_FLUXALP~TM_MONETIER,data=TMP[which(TMP$MONTH==i),])
    LM1.l[[i]] <- lm(TM_FLUXALP~TM_BESSE,data=TMP[which(TMP$MONTH==i),])
    LM2.l[[i]] <- lm(TM_FLUXALP~TM_MONETIER+TM_BESSE,data=TMP[which(TMP$MONTH==i),])
    LM3.l[[i]] <- lm(TM_FLUXALP~TM_MONETIER+TM_BESSE+TM_VALLOIRE,data=TMP[which(TMP$MONTH==i),])
    }
    
    LM0.R2 <- unlist(lapply(LM0.l,function(x) summary(x)$r.squared))
    LM1.R2 <- unlist(lapply(LM1.l,function(x) summary(x)$r.squared))
    LM2.R2 <- unlist(lapply(LM2.l,function(x) summary(x)$r.squared))
    LM3.R2 <- unlist(lapply(LM3.l,function(x) summary(x)$r.squared))
    
    graphics.off();windows(10,10)
    plot(0,0,xlim=c(1,12),ylim=c(0.4,1),type="n",xlab="month")
    lines(LM0.R2,type="b",pch=21,bg="blue",cex=2)
    lines(LM1.R2,type="b",pch=21,bg="lightblue",cex=2)
    lines(LM2.R2,type="b",pch=21,bg="green",cex=2)
    lines(LM3.R2,type="b",pch=21,bg="darkgreen",cex=3)
    # Conclusion : the model is better constrained with LM3 or LM2 with BESSE
    # For a single station Besse is outperforming MONETIER !
  
  # Predict a full time series at FLUXALP
    TMP1.pred <- data.frame(matrix(NA,nrow(TMP),12))
    colnames(TMP1.pred) <- 1:12
    TMP0.pred <- TMP2.pred <- TMP3.pred <- TMP1.pred
    for (i in 1:12){
      TMP0.pred[which(TMP$MONTH==i),i] <- predict(LM0.l[[i]], newdata=TMP[which(TMP$MONTH==i),])
      TMP1.pred[which(TMP$MONTH==i),i] <- predict(LM1.l[[i]], newdata=TMP[which(TMP$MONTH==i),])
      TMP2.pred[which(TMP$MONTH==i),i] <- predict(LM2.l[[i]], newdata=TMP[which(TMP$MONTH==i),])
      TMP3.pred[which(TMP$MONTH==i),i] <- predict(LM3.l[[i]], newdata=TMP[which(TMP$MONTH==i),])
    }
    LM3.pred <- apply(TMP3.pred,1,mean,na.rm=T)
    plot(LM3.pred,type="l")
    LM2.pred <- apply(TMP2.pred,1,mean,na.rm=T)
    lines(LM2.pred,type="l",col="blue")
    LM1.pred <- apply(TMP1.pred,1,mean,na.rm=T)
    plot(LM1.pred,type="l",col="green")
    LM0.pred <- apply(TMP0.pred,1,mean,na.rm=T)
    lines(LM1.pred,type="l",col="green")
    
    TMP$FLUXALP.pred<-NA
    # RULE : observation > LM3 mod > LM2 mod > LM1 mod > LM0 mod
    TMP$FLUXALP.pred <- ifelse(!is.na(TMP$TM_FLUXALP),TMP$TM_FLUXALP,NA)
    plot(TMP$DATE,TMP$FLUXALP.pred,type="l")
    TMP$FLUXALP.pred <- ifelse(is.na(TMP$TM_FLUXALP)& !is.na(LM3.pred),LM3.pred,TMP$FLUXALP.pred)
    plot(TMP$DATE,TMP$FLUXALP.pred,type="l",col="blue")
    TMP$FLUXALP.pred <- ifelse(is.na(TMP$TM_FLUXALP) & is.na(LM3.pred) & !is.na(LM2.pred),LM2.pred,TMP$FLUXALP.pred)
    plot(TMP$DATE,TMP$FLUXALP.pred,type="l",col="darkgreen")
    TMP$FLUXALP.pred <- ifelse(is.na(TMP$TM_FLUXALP) & is.na(LM3.pred) & is.na(LM2.pred)& !is.na(LM1.pred),LM1.pred,TMP$FLUXALP.pred)
    plot(TMP$DATE,TMP$FLUXALP.pred,type="l",col="green")    
    TMP$FLUXALP.pred <- ifelse(is.na(TMP$TM_FLUXALP) & is.na(LM3.pred) & is.na(LM2.pred) & is.na(LM1.pred) & !is.na(LM0.pred),LM0.pred,TMP$FLUXALP.pred)
    plot(TMP$DATE,TMP$FLUXALP.pred,type="l",col="red")
    
    setwd(WD);write.csv(TMP,"FORCAGE_METEO.csv",row.names=F) 

   # Estimate of ELR (not to be used)  
    ELR1 <- coef(LM1)[1]/(ALTI[2]-ALTI[1])
    ELR2 <- coef(LM2)[2]/(ALTI[3]-ALTI[2])
    ELR3 <- coef(LM1)[3]/DIFF[2]
    
    mean(TMP[,4]-TMP[,2],na.rm=T)/(ALTI[2]-ALTI[1])
    mean(TMP[,3]-TMP[,2],na.rm=T)/(ALTI[3]-ALTI[2])
    mean(TMP[,5]-TMP[,2],na.rm=T)/(ALTI[4]-ALTI[2])
    
    # Monthly ELR along the gradient
    tmp <- aggregate(TMP[,4]-TMP[,2],list(MONTH=TMP$MO),mean,na.rm=T)
    tmp$ELR1 <- tmp$x/(ALTI[2]-ALTI[1]) 
    tmp1 <- aggregate(TMP[,3]-TMP[,2],list(MONTH=TMP$MO),mean,na.rm=T)
    tmp$ELR2 <- tmp1$x/(ALTI[2]-ALTI[3]) 
    tmp2 <- aggregate(TMP[,3]-TMP[,4],list(MONTH=TMP$MO),mean,na.rm=T)
    tmp$ELR3 <- tmp2$x/(ALTI[4]-ALTI[3]) 
    
    plot(100*tmp$ELR1,type="b",pch=21,ylim=c(0,1),bg="white",xlab="month",cex=2)
    lines(100*tmp$ELR2,type="b",pch=21,bg="gray",cex=2)
    lines(100*tmp$ELR3,type="b",pch=21,bg="black",cex=2)
  

      # OLD QUANTILE REGRESSION Alix -> FORCAGE_METEO_ALIX.csv ----
    
    # Requires a YEAR, DOY, val dataframe
    find_first_day_above_thresh=function(df,thresh=1200){
      df=pivot_wider(df,names_from = "YEAR",values_from = "val")
      df=df %>% arrange(DOY)
      df[is.na(df)]=0
      df[apply(df[,c(-1)],MARGIN=2,function(x) which(x>=thresh)[1]),]$DOY
    }
    
    setwd(WD) 
    temp <- read.csv("FORCAGE_METEO.csv")
    temp <- temp[!is.na(temp$TM_MONETIER),]
  
    ## Initializing quantile regression storage
    qrfit_monet_F <- list(length=length(unique(temp$MONTH)))
    qrfit_monet_N <- list(length=length(unique(temp$MONTH)))
    qrfit_bes_F   <- list(length=length(unique(temp$MONTH)))
    qrfit_bes_N   <- list(length=length(unique(temp$MONTH)))
    qrfit_val_F   <- list(length=length(unique(temp$MONTH)))
    qrfit_val_N   <- list(length=length(unique(temp$MONTH)))
    
    temp$TM_FLUXALP_predMONETIER    <- NA
    temp$TM_NIVOGALI_predMONETIER   <- NA
    temp$TM_FLUXALP_predBESSE       <- NA
    temp$TM_NIVOGALI_predBESSE      <- NA
    temp$TM_FLUXALP_predVALLOIRE    <- NA
    temp$TM_NIVOGALI_predVALLOIRE   <- NA
    
    ## Initializing linear regression storage for FluxAlp
    lmfit_monet_F <- list(length=length(unique(temp$MONTH)))
    lmfit_bes_F   <- list(length=length(unique(temp$MONTH)))
    lmfit_val_F   <- list(length=length(unique(temp$MONTH)))
    
    temp$TM_FLUXALP_predMONETIER_lm <- NA
    temp$TM_FLUXALP_predBESSE_lm    <- NA
    temp$TM_FLUXALP_predVALLOIRE_lm <- NA
    
    ## Quantile and linear regression by month for the local stations, each from the 3 distant stations
    for(i in 1:length(unique(temp$MONTH))){
      print(i)
      tmp=temp[temp$MONTH==i,]
      
      # Quantile regression
      qrfit_monet_F[[i]] <- rq(TM_FLUXALP~TM_MONETIER,data = tmp)
      temp[temp$MONTH==i,"TM_FLUXALP_predMONETIER"] <- predict(qrfit_monet_F[[i]],data.frame(TM_MONETIER=temp[temp$MONTH==i,"TM_MONETIER"]))
      
      qrfit_monet_N[[i]]=rq(TM_NIVOGALI~TM_MONETIER,data = tmp)
      temp[temp$MONTH==i,]$TM_NIVOGALI_predMONETIER=predict(qrfit_monet_N[[i]],data.frame(TM_MONETIER=temp[temp$MONTH==i,]$TM_MONETIER))
      
      qrfit_bes_F[[i]]=rq(TM_FLUXALP~TM_BESSE,data = tmp)
      temp[temp$MONTH==i,]$TM_FLUXALP_predBESSE=predict(qrfit_bes_F[[i]],data.frame(TM_BESSE=temp[temp$MONTH==i,]$TM_BESSE))
      
      qrfit_bes_N[[i]]=rq(TM_NIVOGALI~TM_BESSE,data = tmp)
      temp[temp$MONTH==i,]$TM_NIVOGALI_predBESSE=predict(qrfit_bes_N[[i]],data.frame(TM_BESSE=temp[temp$MONTH==i,]$TM_BESSE))
      
      qrfit_val_F[[i]]=rq(TM_FLUXALP~TM_VALLOIRE,data = tmp)
      temp[temp$MONTH==i,]$TM_FLUXALP_predVALLOIRE=predict(qrfit_val_F[[i]],data.frame(TM_VALLOIRE=temp[temp$MONTH==i,]$TM_VALLOIRE))
      
      qrfit_val_N[[i]]=rq(TM_NIVOGALI~TM_VALLOIRE,data = tmp)
      temp[temp$MONTH==i,]$TM_NIVOGALI_predVALLOIRE=predict(qrfit_val_N[[i]],data.frame(TM_VALLOIRE=temp[temp$MONTH==i,]$TM_VALLOIRE))
      
      # Linear regression
      lmfit_monet_F[[i]]=lm(TM_FLUXALP~TM_MONETIER,data = tmp)
      temp[temp$MONTH==i,]$TM_FLUXALP_predMONETIER_lm=predict(lmfit_monet_F[[i]],data.frame(TM_MONETIER=temp[temp$MONTH==i,]$TM_MONETIER))
      
      lmfit_bes_F[[i]]=lm(TM_FLUXALP~TM_BESSE,data = tmp)
      temp[temp$MONTH==i,]$TM_FLUXALP_predBESSE_lm=predict(lmfit_bes_F[[i]],data.frame(TM_BESSE=temp[temp$MONTH==i,]$TM_BESSE))
      
      lmfit_val_F[[i]]=lm(TM_FLUXALP~TM_VALLOIRE,data = tmp)
      temp[temp$MONTH==i,]$TM_FLUXALP_predVALLOIRE_lm=predict(lmfit_val_F[[i]],data.frame(TM_VALLOIRE=temp[temp$MONTH==i,]$TM_VALLOIRE))
    }
    
    # Quantile regression analysis
    temp$TM_FLUXALP_predMEAN=rowMeans(temp[,c("TM_FLUXALP_predMONETIER","TM_FLUXALP_predBESSE","TM_FLUXALP_predVALLOIRE")],na.rm=T)
    temp$TM_FLUXALP_predMIN=apply(temp[,c("TM_FLUXALP_predMONETIER","TM_FLUXALP_predBESSE","TM_FLUXALP_predVALLOIRE")],MARGIN = 1,FUN=min,na.rm=T)
    temp$TM_FLUXALP_predMAX=apply(temp[,c("TM_FLUXALP_predMONETIER","TM_FLUXALP_predBESSE","TM_FLUXALP_predVALLOIRE")],MARGIN = 1,FUN=max,na.rm=T)
    temp$TM_FLUXALP_predMED=apply(temp[,c("TM_FLUXALP_predMONETIER","TM_FLUXALP_predBESSE","TM_FLUXALP_predVALLOIRE")],MARGIN = 1,FUN=quantile,probs=0.5,na.rm=T)
    
    temp$TM_NIVOGALI_predMEAN=rowMeans(temp[,c("TM_NIVOGALI_predMONETIER","TM_NIVOGALI_predBESSE","TM_NIVOGALI_predVALLOIRE")],na.rm=T)
    temp$TM_NIVOGALI_predMIN=apply(temp[,c("TM_NIVOGALI_predMONETIER","TM_NIVOGALI_predBESSE","TM_NIVOGALI_predVALLOIRE")],MARGIN = 1,FUN=min,na.rm=T)
    temp$TM_NIVOGALI_predMAX=apply(temp[,c("TM_NIVOGALI_predMONETIER","TM_NIVOGALI_predBESSE","TM_NIVOGALI_predVALLOIRE")],MARGIN = 1,FUN=max,na.rm=T)
    temp$TM_NIVOGALI_predMED=apply(temp[,c("TM_NIVOGALI_predMONETIER","TM_NIVOGALI_predBESSE","TM_NIVOGALI_predVALLOIRE")],MARGIN = 1,FUN=quantile,probs=0.5,na.rm=T)
    
    # Linear
    temp$TM_FLUXALP_predMEAN_lm=rowMeans(temp[,c("TM_FLUXALP_predMONETIER_lm","TM_FLUXALP_predBESSE_lm","TM_FLUXALP_predVALLOIRE_lm")],na.rm=T)
    
    
    ## MY OLD STUFF - Calculate linear model without replacing with true value for recent period and average over 3 stations
    # temp$FLUXALP.pred_mean_lm <- temp$FLUXALP.pred
    # lmfit_monet               <- lm(TM_FLUXALP~TM_MONETIER,data = temp)
    # lmfit_bes                 <- lm(TM_FLUXALP~TM_BESSE,data = temp)
    # lmfit_val                 <- lm(TM_FLUXALP~TM_VALLOIRE,data = temp)
    # lmpred_monet              <- predict(lmfit_monet,data.frame(TM_MONETIER=temp$TM_MONETIER))
    # lmpred_bes                <- predict(lmfit_bes,data.frame(TM_BESSE=temp$TM_BESSE))
    # lmpred_val                <- predict(lmfit_val,data.frame(TM_VALLOIRE=temp$TM_VALLOIRE))
    # temp$FLUXALP.pred_mean_lm <- rowMeans(cbind(lmpred_monet,lmpred_bes,lmpred_val),na.rm=T)
    
      # Calculate accumulation of DD -> FXAoverMON.csv (20250124) ----
    ret <- temp[,c("DATE","TM_FLUXALP_predMONETIER")]
    colnames(ret) <- c("DATE","TEMP")
    ret$YEAR <- lubridate::year(ret$DATE)
    ret$WDD <- ifelse(ret$TEMP>0,ret$TEMP,0)
    ret$DOY <- lubridate::yday(ret$DATE)
    
    # complete year 1880 with 0 and compute DD using collapse:: cumsum
    myret <- ret[,c("YEAR","DOY","WDD")]
    #myret <- rbind(data.frame(YEAR=1780,DOY=1:14,WDD=0),myret)
    tmp   <- reshape::cast(myret,YEAR~DOY,value="WDD",add.missing=T)
    tmp   <- t(tmp)
    tmp <- collapse::fcumsum(tmp)
    plot(tmp[,59])
    setwd(WD); write.csv(tmp,"FXAoverMON.DD.csv",row.names=F)
    
        # DO THE SAME on monthly means of MONETIER and HISTALP (20250124) ----
        # using the predicted FXA based on MON !!
    
    setwd("~/CLIMATOmc/DATA/DATA_HISTALP/")
    tmp   <- terra::rast("HISTALP_temperature_1780-2014.nc")
    tmp0  <- data.frame(X_wgs84=c(6.38994,6.41061),Y_wgs84=c(45.05495,45.04128))
    tmp1  <- terra::vect(tmp0,geom=c("X_wgs84","Y_wgs84"),crs="epsg:4326")
    EXT   <- t(terra::extract(tmp,tmp1))[-(1:3),1]
    graphics.off();windows(15,10);plot(EXT,type="l")

    #  extract monthly values
    length(1780:2014)*12
    st <- ymd("1780-01-15")
    en <- ymd("2014-12-15")
    YM <- st %m+% months(seq(0, round(interval(st, en) / months(1)), 1))

    M.HIST <- data.frame(YEAR = lubridate::year(YM),
                         MONTH= lubridate::month(YM),
                         TEMP = EXT
    )
    # delete years after 2010
    M.HIST <- M.HIST[-which(M.HIST$YEAR>2010),]
    
    #  temp$TM_FLUXALP_predMONETIER from previous step !!
    
    df <- temp[,c("MONTH","YEAR","TM_FLUXALP_predMONETIER")]
    M.MON <- aggregate(df$TM_FLUXALP_predMONETIER,list(MONTH=df$MONTH,YEAR=df$YEAR),mean)
    colnames(M.MON) <-c("MONTH","YEAR","TEMP")
    
    DATA <- M.HIST%>%left_join(M.MON,by=c("YEAR","MONTH"),suffix=c(".HIS",".MON"))
    
    plot(DATA$TEMP.MON,type="l")
    
    DATA$TEMP.MON_predHISTALP_rq <- NA
    DATA$TEMP.MON_predHISTALP_lm <- NA

    qrfit <- lmfit <- list(length=12)
    
    for(i in 1:12){
      print(i)
      tmp=DATA[DATA$MONTH==i,]
      qrfit[[i]] <- rq(TEMP.MON~TEMP.HIS,data = tmp)
      
      DATA[DATA$MONTH==i,]$TEMP.MON_predHISTALP_rq <- predict(qrfit[[i]],data.frame(TEMP.HIS=DATA[DATA$MONTH==i,]$TEMP.HIS))
      
      lmfit[[i]] <- lm(TEMP.MON~TEMP.HIS,data = tmp)      
      
      DATA[DATA$MONTH==i,]$TEMP.MON_predHISTALP_lm <- predict(lmfit[[i]],data.frame(TEMP.HIS=DATA[DATA$MONTH==i,]$TEMP.HIS))
    }
    
    
    plot(DATA$TEMP.MON_predHISTALP_lm,DATA$TEMP.MON)
    plot(DATA$TEMP.MON_predHISTALP_rq,DATA$TEMP.MON)
    abline(0,1)
    
    
      # TRANSFORM TO DAILY DATA ---- 
    # (need to use Harmonic or Discrete Fourier Transform)
    
    # removing 4 years !!!
    df <- data.frame(TEMP=DATA[,"TEMP.MON_predHISTALP_rq"],
                     DATE=seq(as.Date('1780-01-15'),as.Date('2010-12-15'),'month')
    )
    
    # order of day
    seq.REF <- seq(as.Date('1780-01-01'),as.Date('2010-12-31'),'day')
    length(seq.REF)
    
    df$ORD <- match(df$DATE,seq.REF)
    
    graphics.off();windows(20,10);plot(df$DATE,df$TEMP,type="l")
    
    x=100+df$TEMP; n=nrow(df)/2; up=30; plot=TRUE; add=FALSE; main=NULL
    # using Fast Fourier Transform
    
    # The direct transformation
    # The first frequency is DC, the rest are duplicated
    
    dff = stats::fft(x)
    # The time
    # t = seq(from = 1, to = length(x))
    t = df$ORD
    
    # Upsampled time
    # nt = seq(from = 1, to = length(x)+1-1/up, by = 1/up)
    nt = df$ORD[1]:df$ORD[nrow(df)]
    
    #New spectrum
    ndff    = array(data = 0, dim = c(length(nt), 1L))
    ndff[1] = dff[1] #Always, it's the DC component
    
    if(n != 0){
      #The positive frequencies always come first
      ndff[2:(n+1)] = dff[2:(n+1)] 
      
      #The negative ones are trickier
      ndff[length(ndff):(length(ndff) - n + 1)] = dff[length(x):(length(x) - n + 1)]
    }
    
    #The inverses
    indff = fft(ndff/nrow(df), inverse = TRUE)
    idff  = fft(dff/nrow(df), inverse = TRUE)
    
    if(plot){
      if(!add){
        plot(x = t, y = x, pch = 16L, xlab = "Time", ylab = "Measurement",
             main = ifelse(is.null(main), paste(n, "harmonics"), main))
        lines(y = Mod(idff), x = t, col = "green")
      }
      lines(y = Mod(indff), x = nt,col="red")
    }
    
    ret = data.frame(ORD = nt, 
                     TEMP = Mod(indff)-100, 
                     DATE = seq.REF[df$ORD[1]:df$ORD[nrow(df)]],
                     YEAR=lubridate::year(seq.REF[df$ORD[1]:df$ORD[nrow(df)]])
                     )
    #  §§§§§§   PAS DE DECALAGE décalage de 15 jours ? ££££
   ret$MONTH <- lubridate::month(ret$DATE)
   
   chk <- aggregate(ret$TEMP,list(MONTH=ret$MONTH,YEAR=ret$YEAR),mean)
   plot(chk$x,DATA$TEMP.MON_predHISTALP_rq)
   abline(0,1)
    
    
    plot(ret$DATE[1:2000],ret$TEMP[1:2000],type='l',col="red")
    points(df$DATE,df$TEMP,pch=21,bg="green")
    abline(h=0)
    

        # FLUXALP RMSE ----
    
    ## Look for for best predictor among MONETIER, BESSE, VALLOIRE, mean, min, max and median of the 3 ; also mean of 3 linear models
    temp_test=subset(temp[!is.na(temp$TM_FLUXALP),],select=-c(TM_NIVOGALI))
    temp_test=temp_test[!rowSums(is.na(temp_test)) > 0,]
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$FLUXALP.pred))# Philippe
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$FLUXALP.pred_mean_lm))# not by month
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predMEAN_lm))# by month
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predMONETIER))
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predVALLOIRE))
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predBESSE))
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predMEAN))
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predMIN))
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predMAX))
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predMED))
    
    ## Similar but removing months after August
    temp_test=temp_test[temp_test$MONTH<=8,]
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$FLUXALP.pred_mean_lm))#not by month
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predMEAN_lm))#by month
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predMONETIER))
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predVALLOIRE))
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predBESSE))
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predMEAN))
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predMIN))
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predMAX))
    print(Metrics::rmse(temp_test$TM_FLUXALP,temp_test$TM_FLUXALP_predMED))
    
    ## Look at best predictor of day with cumulative positive degree-days above THR
    temp_test=subset(temp[!is.na(temp$TM_FLUXALP),],select=-c(TM_NIVOGALI))
    temp_test=temp_test[!rowSums(is.na(temp_test)) > 0,]
    temp_test$DATE=ymd(temp_test$DATE)
    temp_test$DOY=yday(temp_test$DATE)
    
    TM_FLUXALP_cumsum=data.frame(DOY=temp_test$DOY[temp_test$TM_FLUXALP>0],YEAR=temp_test$YEAR[temp_test$TM_FLUXALP>0],val=ave(temp_test$TM_FLUXALP[temp_test$TM_FLUXALP>0], temp_test$YEAR[temp_test$TM_FLUXALP>0], FUN = cumsum))
    TM_FLUXALP_predMEAN_lm_cumsum=data.frame(
      DOY  = temp_test$DOY[temp_test$TM_FLUXALP_predMEAN_lm>0],
      YEAR = temp_test$YEAR[temp_test$TM_FLUXALP_predMEAN_lm>0],
      val  = ave(temp_test$TM_FLUXALP_predMEAN_lm[temp_test$TM_FLUXALP_predMEAN_lm>0],temp_test$YEAR[temp_test$TM_FLUXALP_predMEAN_lm>0], FUN = cumsum)
      )
    
    FLUXALP.pred_mean_lm_cumsum=data.frame(
      DOY=temp_test$DOY[temp_test$FLUXALP.pred_mean_lm>0],
      YEAR=temp_test$YEAR[temp_test$FLUXALP.pred_mean_lm>0],
      val=ave(temp_test$FLUXALP.pred_mean_lm[temp_test$FLUXALP.pred_mean_lm>0], temp_test$YEAR[temp_test$FLUXALP.pred_mean_lm>0], FUN = cumsum))
    
    TM_FLUXALP_predMONETIER_cumsum=data.frame(
      DOY=temp_test$DOY[temp_test$TM_FLUXALP_predMONETIER>0],
      YEAR=temp_test$YEAR[temp_test$TM_FLUXALP_predMONETIER>0],
      val=ave(temp_test$TM_FLUXALP_predMONETIER[temp_test$TM_FLUXALP_predMONETIER>0], temp_test$YEAR[temp_test$TM_FLUXALP_predMONETIER>0], FUN = cumsum))
    
    TM_FLUXALP_predBESSE_cumsum=data.frame(
      DOY=temp_test$DOY[temp_test$TM_FLUXALP_predBESSE>0],
      YEAR=temp_test$YEAR[temp_test$TM_FLUXALP_predBESSE>0],
      val=ave(temp_test$TM_FLUXALP_predBESSE[temp_test$TM_FLUXALP_predBESSE>0], temp_test$YEAR[temp_test$TM_FLUXALP_predBESSE>0], FUN = cumsum))
    
    TM_FLUXALP_predVALLOIRE_cumsum=data.frame(DOY=temp_test$DOY[temp_test$TM_FLUXALP_predVALLOIRE>0],YEAR=temp_test$YEAR[temp_test$TM_FLUXALP_predVALLOIRE>0],val=ave(temp_test$TM_FLUXALP_predVALLOIRE[temp_test$TM_FLUXALP_predVALLOIRE>0], temp_test$YEAR[temp_test$TM_FLUXALP_predVALLOIRE>0], FUN = cumsum))
    TM_FLUXALP_predMEAN_cumsum=data.frame(DOY=temp_test$DOY[temp_test$TM_FLUXALP_predMEAN>0],YEAR=temp_test$YEAR[temp_test$TM_FLUXALP_predMEAN>0],val=ave(temp_test$TM_FLUXALP_predMEAN[temp_test$TM_FLUXALP_predMEAN>0], temp_test$YEAR[temp_test$TM_FLUXALP_predMEAN>0], FUN = cumsum))
    TM_FLUXALP_predMIN_cumsum=data.frame(DOY=temp_test$DOY[temp_test$TM_FLUXALP_predMIN>0],YEAR=temp_test$YEAR[temp_test$TM_FLUXALP_predMIN>0],val=ave(temp_test$TM_FLUXALP_predMIN[temp_test$TM_FLUXALP_predMIN>0], temp_test$YEAR[temp_test$TM_FLUXALP_predMIN>0], FUN = cumsum))
    TM_FLUXALP_predMAX_cumsum=data.frame(DOY=temp_test$DOY[temp_test$TM_FLUXALP_predMAX>0],YEAR=temp_test$YEAR[temp_test$TM_FLUXALP_predMAX>0],val=ave(temp_test$TM_FLUXALP_predMAX[temp_test$TM_FLUXALP_predMAX>0], temp_test$YEAR[temp_test$TM_FLUXALP_predMAX>0], FUN = cumsum))
    TM_FLUXALP_predMED_cumsum=data.frame(DOY=temp_test$DOY[temp_test$TM_FLUXALP_predMED>0],YEAR=temp_test$YEAR[temp_test$TM_FLUXALP_predMED>0],val=ave(temp_test$TM_FLUXALP_predMED[temp_test$TM_FLUXALP_predMED>0], temp_test$YEAR[temp_test$TM_FLUXALP_predMED>0], FUN = cumsum))
    
    
    # evaluate the mean absolute error - in number of days to reach THR
    # clearly the quatile regression with Monetier is the best
    print_rmse_date_thresh=function(df,thresh=1200){
      Metrics::mae(find_first_day_above_thresh(TM_FLUXALP_cumsum,thresh=thresh)[!is.na(find_first_day_above_thresh(TM_FLUXALP_cumsum,thresh=thresh))],find_first_day_above_thresh(df,thresh=thresh)[!is.na(find_first_day_above_thresh(df,thresh=thresh))])
    }
    
    THR=700
    print_rmse_date_thresh(FLUXALP.pred_mean_lm_cumsum,THR)
    print_rmse_date_thresh(TM_FLUXALP_predMEAN_lm_cumsum,THR)
    print_rmse_date_thresh(TM_FLUXALP_predMONETIER_cumsum,THR) # best model
    print_rmse_date_thresh(TM_FLUXALP_predBESSE_cumsum,THR)
    print_rmse_date_thresh(TM_FLUXALP_predVALLOIRE_cumsum,THR)
    print_rmse_date_thresh(TM_FLUXALP_predMEAN_cumsum,THR)
    print_rmse_date_thresh(TM_FLUXALP_predMIN_cumsum,THR)
    print_rmse_date_thresh(TM_FLUXALP_predMAX_cumsum,THR)
    print_rmse_date_thresh(TM_FLUXALP_predMED_cumsum,THR)
    
    THR <- seq(500,2000,50)
    MAE.MONETIER <- rep(NA, length(THR))
    for (i in 1:length(THR)){
      MAE.MONETIER[i] <- print_rmse_date_thresh(TM_FLUXALP_predMONETIER_cumsum,THR[i])
    }
    plot(THR,MAE.MONETIER,ylim=c(0,10))
    
        # NIVOSE RMSE (not the linear by month model) ----
    
    ## Look for best predictor among MONETIER, BESSE, VALLOIRE, mean, min, max and median of the 3 ; also mean of 3 linear models
    temp_test=temp[!is.na(temp$TM_NIVOGALI),]
    temp_test=temp_test[!rowSums(is.na(temp_test)) > 0,]
    print(Metrics::rmse(temp_test$TM_NIVOGALI,temp_test$TM_NIVOGALI_predMONETIER))
    print(Metrics::rmse(temp_test$TM_NIVOGALI,temp_test$TM_NIVOGALI_predVALLOIRE))
    print(Metrics::rmse(temp_test$TM_NIVOGALI,temp_test$TM_NIVOGALI_predBESSE))
    print(Metrics::rmse(temp_test$TM_NIVOGALI,temp_test$TM_NIVOGALI_predMEAN))
    print(Metrics::rmse(temp_test$TM_NIVOGALI,temp_test$TM_NIVOGALI_predMIN))
    print(Metrics::rmse(temp_test$TM_NIVOGALI,temp_test$TM_NIVOGALI_predMAX))
    print(Metrics::rmse(temp_test$TM_NIVOGALI,temp_test$TM_NIVOGALI_predMED))
    
    ## Similar but removing months after August
    temp_test=temp_test[temp_test$MONTH<=8,]
    print(Metrics::rmse(temp_test$TM_NIVOGALI,temp_test$TM_NIVOGALI_predMONETIER))
    print(Metrics::rmse(temp_test$TM_NIVOGALI,temp_test$TM_NIVOGALI_predVALLOIRE))
    print(Metrics::rmse(temp_test$TM_NIVOGALI,temp_test$TM_NIVOGALI_predBESSE))
    print(Metrics::rmse(temp_test$TM_NIVOGALI,temp_test$TM_NIVOGALI_predMEAN))
    print(Metrics::rmse(temp_test$TM_NIVOGALI,temp_test$TM_NIVOGALI_predMIN))
    print(Metrics::rmse(temp_test$TM_NIVOGALI,temp_test$TM_NIVOGALI_predMAX))
    print(Metrics::rmse(temp_test$TM_NIVOGALI,temp_test$TM_NIVOGALI_predMED))
    
    ## Look at best predictor of day with cumulative positive degree-days above 700 and 1200
    temp_test=temp[!is.na(temp$TM_NIVOGALI),]
    temp_test=temp_test[!rowSums(is.na(temp_test)) > 0,]
    temp_test$DATE=ymd(temp_test$DATE)
    temp_test$DOY=yday(temp_test$DATE)
    TM_NIVOGALI_cumsum=data.frame(DOY=temp_test$DOY[temp_test$TM_NIVOGALI>0],YEAR=temp_test$YEAR[temp_test$TM_NIVOGALI>0],val=ave(temp_test$TM_NIVOGALI[temp_test$TM_NIVOGALI>0], temp_test$YEAR[temp_test$TM_NIVOGALI>0], FUN = cumsum))
    TM_NIVOGALI_predMONETIER_cumsum=data.frame(DOY=temp_test$DOY[temp_test$TM_NIVOGALI_predMONETIER>0],YEAR=temp_test$YEAR[temp_test$TM_NIVOGALI_predMONETIER>0],val=ave(temp_test$TM_NIVOGALI_predMONETIER[temp_test$TM_NIVOGALI_predMONETIER>0], temp_test$YEAR[temp_test$TM_NIVOGALI_predMONETIER>0], FUN = cumsum))
    TM_NIVOGALI_predBESSE_cumsum=data.frame(DOY=temp_test$DOY[temp_test$TM_NIVOGALI_predBESSE>0],YEAR=temp_test$YEAR[temp_test$TM_NIVOGALI_predBESSE>0],val=ave(temp_test$TM_NIVOGALI_predBESSE[temp_test$TM_NIVOGALI_predBESSE>0], temp_test$YEAR[temp_test$TM_NIVOGALI_predBESSE>0], FUN = cumsum))
    TM_NIVOGALI_predVALLOIRE_cumsum=data.frame(DOY=temp_test$DOY[temp_test$TM_NIVOGALI_predVALLOIRE>0],YEAR=temp_test$YEAR[temp_test$TM_NIVOGALI_predVALLOIRE>0],val=ave(temp_test$TM_NIVOGALI_predVALLOIRE[temp_test$TM_NIVOGALI_predVALLOIRE>0], temp_test$YEAR[temp_test$TM_NIVOGALI_predVALLOIRE>0], FUN = cumsum))
    TM_NIVOGALI_predMEAN_cumsum=data.frame(DOY=temp_test$DOY[temp_test$TM_NIVOGALI_predMEAN>0],YEAR=temp_test$YEAR[temp_test$TM_NIVOGALI_predMEAN>0],val=ave(temp_test$TM_NIVOGALI_predMEAN[temp_test$TM_NIVOGALI_predMEAN>0], temp_test$YEAR[temp_test$TM_NIVOGALI_predMEAN>0], FUN = cumsum))
    TM_NIVOGALI_predMIN_cumsum=data.frame(DOY=temp_test$DOY[temp_test$TM_NIVOGALI_predMIN>0],YEAR=temp_test$YEAR[temp_test$TM_NIVOGALI_predMIN>0],val=ave(temp_test$TM_NIVOGALI_predMIN[temp_test$TM_NIVOGALI_predMIN>0], temp_test$YEAR[temp_test$TM_NIVOGALI_predMIN>0], FUN = cumsum))
    TM_NIVOGALI_predMAX_cumsum=data.frame(DOY=temp_test$DOY[temp_test$TM_NIVOGALI_predMAX>0],YEAR=temp_test$YEAR[temp_test$TM_NIVOGALI_predMAX>0],val=ave(temp_test$TM_NIVOGALI_predMAX[temp_test$TM_NIVOGALI_predMAX>0], temp_test$YEAR[temp_test$TM_NIVOGALI_predMAX>0], FUN = cumsum))
    TM_NIVOGALI_predMED_cumsum=data.frame(DOY=temp_test$DOY[temp_test$TM_NIVOGALI_predMED>0],YEAR=temp_test$YEAR[temp_test$TM_NIVOGALI_predMED>0],val=ave(temp_test$TM_NIVOGALI_predMED[temp_test$TM_NIVOGALI_predMED>0], temp_test$YEAR[temp_test$TM_NIVOGALI_predMED>0], FUN = cumsum))
    
    
    
    print_rmse_date_thresh=function(df,thresh=1200){
      Metrics::mae(find_first_day_above_thresh(TM_NIVOGALI_cumsum,thresh=thresh)[!is.na(find_first_day_above_thresh(TM_NIVOGALI_cumsum,thresh=thresh))],find_first_day_above_thresh(df,thresh=thresh)[!is.na(find_first_day_above_thresh(df,thresh=thresh))])
    }
    
    THR=400
    print_rmse_date_thresh(TM_NIVOGALI_predMONETIER_cumsum,THR)
    print_rmse_date_thresh(TM_NIVOGALI_predBESSE_cumsum,THR)
    print_rmse_date_thresh(TM_NIVOGALI_predVALLOIRE_cumsum,THR)
    print_rmse_date_thresh(TM_NIVOGALI_predMEAN_cumsum,THR)
    print_rmse_date_thresh(TM_NIVOGALI_predMIN_cumsum,THR)
    print_rmse_date_thresh(TM_NIVOGALI_predMAX_cumsum,THR)
    print_rmse_date_thresh(TM_NIVOGALI_predMED_cumsum,THR)
    
    
        # FINAL outputs ----
    ## Mean of monthly quantile regression (or linear regression) is best predictor for both stations. But monthly quantile regression is better regarding first day exceeding cumulative positive degree day threshold
    
    tmp <-  temp[,c("DATE","TM_MONETIER","TM_BESSE","TM_VALLOIRE","TM_FLUXALP","TM_NIVOGALI","MONTH","YEAR","FLUXALP.pred","TM_FLUXALP_predMONETIER")]
    colnames(tmp) <- c("DATE","TM_MONETIER","TM_BESSE","TM_VALLOIRE","TM_FLUXALP","TM_NIVOGALI","MONTH","YEAR","FLUXALP.pred_Choler","FLUXALP.pred_Reverdy")
    setwd(WD); write.csv(tmp,file="FORCAGE_METEO_reverdy.csv")
    
    graphics.off(); windows(10,10); plot(tmp$FLUXALP.pred_Choler,tmp$FLUXALP.pred_Reverdy,pch=".")
    abline(0,1,col="red")
    
    # II.C. do the same for long time series (OLD NOT UPDATED) ----
  setwd("~/CLIMATOmc/DATA/CLIMATHEQUE/DAILY_MEANS/")
  TMP1 <- read.csv("ANNECY_METEO.csv")
  TMP2 <- read.csv("CRANGEVR_METEO.csv")
  TMP3 <- read.csv("MONETIER_METEO.csv")
  TMP1$TM2   <- 0.5*(TMP1$TX+TMP1$TN)
  TMP1$DATE  <- as.Date(TMP1$DATE)
  TMP1$MONTH <- lubridate::month(TMP1$DATE)
  TMP1$YEAR  <- lubridate::year(TMP1$DATE)
  
  TMP2$TM2   <- 0.5*(TMP2$TX+TMP2$TN)
  TMP2$DATE  <- as.Date(TMP2$DATE)
  TMP2$MONTH <- lubridate::month(TMP2$DATE)
  TMP2$YEAR  <- lubridate::year(TMP2$DATE)
  plot(TMP1$DATE,TMP1$TM2,type="l",col="darkgray")
  lines(TMP2$DATE,TMP2$TM2,col="red")
  
  TMP3$TM2   <- 0.5*(TMP3$TX+TMP3$TN)
  TMP3$DATE  <- as.Date(TMP3$DATE)
  TMP3$MONTH <- lubridate::month(TMP3$DATE)
  TMP3$YEAR  <- lubridate::year(TMP3$DATE)
  
  out   <- split(TMP3 , f = TMP3$YEAR)
  CDD.MONETIER   <- data.frame(t(as.data.frame(lapply(out,FUN=DDf))))[,10:20]
  CDD.MONETIER$YEAR   <- as.numeric(substr(rownames(CDD.MONETIER),2,5))
  
  out   <- split(TMP2 , f = TMP2$YEAR)
  CDD.CRANGEVR   <- data.frame(t(as.data.frame(lapply(out,FUN=DDf))))[,10:20]
  CDD.CRANGEVR$YEAR   <- as.numeric(substr(rownames(CDD.CRANGEVR),2,5))
  
  out   <- split(TMP1 , f = TMP1$YEAR)
  CDD.ANNECY   <- data.frame(t(as.data.frame(lapply(out,FUN=DDf))))[,10:20]
  CDD.ANNECY$YEAR   <- as.numeric(substr(rownames(CDD.ANNECY),2,5))
  
  DATA <- CDD.MONETIER %>%left_join(CDD.CRANGEVR,by="YEAR")
  DATA <- DATA %>% drop_na()
  cor(DATA[,1],DATA[,20])
  
  TMP2  <- lapply(out,FUN=DDf2)
  plot(DATA[,1],DATA[,20])
  abline(c(0,1))
  
  myYEAR <- as.numeric(substr(colnames(TMP),2,5))
  graphics.off();windows(18,5);plot(myYEAR,as.numeric(TMP[20,]),type="p")
  lines(myYEAR,data.table::frollmean(as.numeric(TMP[20,]),n=10,align="center",na.rm=T),col="blue",lwd=2)

    # II.D. Monthly means and anomalies at FLUXALP v ----
    setwd(WD); FM <- read.csv("FORCAGE_METEO_Reverdy.csv")
  
    FM$HYR <- ifelse(FM$MONTH<=6,1,2)
  
  # 1961-1990 monthly means
  TMP  <- FM[which(FM$YEAR>=1961 & FM$YEAR<=1990),c("HYR","FLUXALP.pred_Reverdy")]
  NORm <- aggregate(TMP$FLUXALP.pred_Choler,list(HYR=TMP$HYR),mean)
  
  
  dfm  <- aggregate(FM$FLUXALP.pred_Reverdy,list(YEAR=FM$YEAR,HYR=FM$HYR),mean)
  dfm$ANOm <- NA
  HYR.FLUX <- dfm[dfm$HYR==1,c("YEAR","x")]
  plot(HYR.FLUX,type="b")
  
  for (i in 1:2){
    SEL <- which(dfm$HYR==i)
    dfm$ANOm[SEL] <- dfm$x[SEL]-NORm$x[i]
  }
  
  SEL <- which(dfm$HYR==1)
  plot(dfm$YEAR[SEL],dfm$ANOm[SEL],type="l")
  lines(dfm$YEAR[SEL],data.table::frollmean(dfm$ANOm[SEL],align="center",n=5),type="l",col="red",lwd=2)
  abline(h=0)
  
    # II.E. Accumulation of degree days in HISTALP -> HIST.DD.csv (OLD NOT UPDATED) ----
      # Calculate accumulation of DD using HISraw (20240726) ----
        ret <- HISraw
        # complete year 1880 with 0 and compute DD using collapse:: cumsum
        myret <- ret[,c("YEAR","DOY","WDD")]
        myret <- rbind(data.frame(YEAR=1780,DOY=1:14,WDD=0),myret)
        tmp   <- reshape::cast(myret,YEAR~DOY,value="WDD",add.missing=T)
        tmp   <- t(tmp)
        library(collapse)
        tmp <- collapse::fcumsum(tmp)
        plot(tmp[,59])
        setwd(WD); write.csv(tmp,"HISTfull.DD.csv",row.names=F)

        # break table per year
        #  Split on userid
        out <- split(ret , f = ret$YEAR)
        DDf <- function(x){
          test0 <- x$TEMP
          test0 <- zoo::na.approx(test0,xout=1:length(test0))
          test1 <- ifelse(test0<=0,0,test0)
          test2 <- cumsum(test1)
          
          THR <- seq(100,2000,100)
          RES <- rep(NA,length(THR))
          for (j in 1:length(THR)){
            RES[j] <- which(test2>THR[j])[1]
          }
          return(RES)
        }
        TMP <- t(as.data.frame(lapply(out,FUN=DDf)))
        
        # compute the DD up to 1st August = 214
        ret.l <- split(ret,f=ret$YEAR)
        lapply(ret.l,cumsum(x))
               
        CUM  <- function(x) {
          END <- which(x$DOY==214)
          return(cumsum(x$WDD[1:END])[END])
        }
        df <- data.frame(YEAR  = 1880:2014,
                         DD214 = unlist(lapply(ret.l,CUM))
                        )
        
        setwd(WD); write.csv(df,"HIST.DD.csv",row.names=F)
        
        # preliminary figures
        # simul for WDD=500 (col5) or 1000(col10)
        # scenario 1
        COL = 8   # WDD = 800°C
        THR = 214 # first of August 
        
        # scenario 2
        COL = 5   # WDD = 600°C
        THR = 214 # first of July 
        
        # scenario 2
        COL = 10   # WDD = 1000°C
        THR = 214 # first of July 
        
        graphics.off()
        windows(10,5)
        par(mar=c(4,4,1,0))
        plot(1880:2014,TMP[,COL],pch=21,bg="lightblue",cex=1.2,xlab="",ylab=paste0("Julian Day for WDD=",100*COL,"°C"),xlim=c(1880,2024))
        lines(1880:2014,data.table::frollmean(TMP[,COL],n=9,align="center",na.rm=T),col="blue",lwd=2)
        abline(v=(1880:2014)[is.na(TMP[,COL])],col="black",lwd=4)
        abline(v=(1880:2014)[which(TMP[,COL]>THR)],col="gray",lwd=3)
        abline(h=THR)
        lines(1880:2014,data.table::frollmean(TMP[,COL],n=9,align="center",na.rm=T),col="blue",lwd=2)
        points(1880:2014,TMP[,COL],pch=21,bg="lightblue",cex=1.2)
        # text(2015,THR+1,paste0("DOY=",THR))
        
        setwd(WD); TMPF <- t(read.csv("FLUXALP.DDseq.csv"))
        points(1943:2023,TMPF[,COL],pch=21,bg="lightgreen",cex=0.8)
        lines(1943:2023,data.table::frollmean(TMPF[,COL],n=9,align="center",na.rm=T),col="darkgreen",lwd=2)
        legend("topright",lwd=2,col=c("darkgreen","blue"),legend=c("FluxAlp","HistAlp"))
        # legend("top",lwd=3,col=c("black","gray"),legend=c("No hit",paste0("No hit before DOY=",THR)),bg="white")
        legend("bottomleft",lwd=3,col=c("gray"),legend=c(paste0("No hit before DOY=",THR)),bg="white")
        
      # compare  -----   
    setwd(WD); FM <- read.csv("FORCAGE_METEO_Reverdy.csv")
    M.FLUX  <- aggregate(FM$FLUXALP.pred_Choler,list(YEAR=FM$YEAR,MONTH=FM$MONTH),mean,na.rm=T)
    
    colnames(M.FLUX)[3]<-"TEMP"
    # il manque 1986 !

    
    DATA <- M.HIST %>%left_join(M.FLUX,by=c("YEAR","MONTH"))
    
    # select a period / month
    COR <- rep(NA,12)
    for(i in 1:12) {
      print(i)
      
      
    SEL<-which(DATA$MONTH==i & DATA$YEAR>=1944 & DATA$YEAR<=2010)

    X.FLUX <- DATA$TEMP.x[SEL]
    Y.HIST <- DATA$TEMP.y[SEL]
    COR[i] <- cor(X.FLUX,Y.HIST,"complete.obs")   
    
    }
    
    plot(1:12,COR)
    mean(COR)
    
    
    SEL <- which(DATA$MONTH==5)
    
    graphics.off();windows(20,10);
    plot(DATA$YEAR[SEL],DATA$TEMP.x[SEL],type="b",cex=0.8,pch=21,bg="gray")
    lines(DATA$YEAR[SEL],data.table::frollmean(DATA$TEMP.x[SEL],n=11,align="center"),col="red",lwd=2)
    abline(h=-1)
    abline(v=1980,col="blue")
    abline(v=c(1945,1981),col="darkgreen")
    abline(v=seq(1880,2020,10),lty=2)
    # lines(HYR.FLUX$YEAR,HYR.FLUX$x,col="darkgreen",lwd=2)
    lines(DATA$YEAR[SEL],data.table::frollmean(DATA$TEMP.y[SEL]-0.5,n=11,align="center"),col="darkgreen",lwd=2,type="b",cex=1,pch=21,bg="darkgreen")
    
  
      # 1°grid datasets temperature anomalies ----
    # from the 1901–2000 average
     # NOT TO RUN
      setwd("~/CLIMATOmc/DATA/DATA_HISTALP/")
      TAB <- read.delim("T01-rekonstr_V3_1780-grid-low.dat",header=F)
      # find the center of the closest 1°1° cell
      INI <- which(TAB==" 6    45")
      setwd(WD);write.csv(TAB[(INI+2):(INI+231),],"HISTALP.06.45.csv",row.names=F)
      # time series are Jan-Dez, the 4 seasons, 2 half years and the year
      # the next 4 columns contain the values for 14. spring (March-May), 15. summer (June-August), 16. autumn (September to November) and 17. winter (December to next year's February) and the next two colums contain the values for the 18.summer halfyear (April to September) and 19.winter halfyear (October to next year's March) while the last column contains the 20. yearly value. The values for winter and winter halfyear are assigned to the year the December value belongs to.
      # convert into columns in excel then reimport
      
      setwd(WD);tmp <- read.csv("HISTALP.06.45.csv",sep=";",header=T)
      # plot the summer anomaly 
      plot(tmp$YEAR,tmp$SUM/10,type="l")
      lines(tmp$YEAR,data.table::frollmean(tmp$SUM/10,align="center",n=10),type="l",col="red",lwd=2)
      abline(h=0)
      
      colnames(tmp)<-c("YEAR",substr(101:119,2,3))
      tmp<- tmp[,1:13]
      temp <- reshape2::melt(tmp,id.vars="YEAR",variable.name="MONTH",value.name="TEMP")
      temp <- temp[order(temp$YEAR,temp$MONTH),]
      temp$TEMP <- temp$TEMP/10
      plot(temp$YEAR,temp$TEMP,ylim=c(-7,8),type="l")
      
      # anomalies from Histalp 1961-1990
      setwd("~/CLIMATOmc/DATA/DATA_HISTALP/")
      TMP <- read.csv("T01-SW_rel1961-1990_CRSM.csv",header=T,sep=";")
      # average over the first 6 month
      df <- data.frame(YEAR=TMP$YEAR,ANOm=apply(TMP[,2:7],1,mean))
      plot(df$YEAR,df$ANOm,type="l")
      lines(df$YEAR,data.table::frollmean(df$ANOm,align="center",n=10),col="red",lwd=2,type="l")
      abline(h=0)
      
      # OLD STUFF ----
      
      res = nff(x = y, n = 18L, up = 100L, col = 2L)
      
      plot(cosreg_fit)
  
    # II.F. Accumulation of DD at FLUXALP -> FLUXALP.DD.csv (20250124) ----
      # REMOVE BISSEXTILE DAYS ! and consider NA spline
    setwd(WD);
    ret       <- read.csv("FORCAGE_METEO.csv")
      ret$DATE  <- as.Date(ret$DATE)
      ret$MONTH <- lubridate::month(ret$DATE)
      ret$YEAR  <- lubridate::year(ret$DATE)
      ret$DOY   <- lubridate::yday(ret$DATE)
      ret       <- ret[,c("YEAR","DOY","TM_FLUXALP")]
      colnames(ret)<-c("YEAR","DOY","TEMP")
      ret$WDD   <- ifelse(ret$TEMP>0,ret$TEMP,0)
      myret <- ret[,c("YEAR","DOY","WDD")]
      tmp   <- reshape::cast(myret,YEAR~DOY,value="WDD",add.missing=T)
      tmp   <- t(tmp)
      tmp <- collapse::fcumsum(tmp)
      setwd(WD); write.csv(tmp,"FXA.DD.csv",row.names=F)
  
    # II.G. Accumulation of DD for BIOCLIM/ALL.PTS -> BIOCLIM.DD.csv (20250121) ----
    #  Climate data already aggregated in ALL.PTS. see part I
  
  setwd(DATA_WD); ALL.PTS <- read.csv("ALL.PTS.csv")
  tmp <- as.factor(as.vector(ALL.PTS$COMcor))
  levels(tmp) <-c("2.CF","4.CT","1.EN","7.FP","6.KS","3.PA","5.VM")
  ALL.PTS$COMcor <- as.vector(tmp)
  
  FM <- ALL.PTS[,c(grep("ID",colnames(ALL.PTS)),grep("COMcor",colnames(ALL.PTS)),grep("SMOD",colnames(ALL.PTS)))]
  FM <- FM[,-grep("SMODmean",colnames(FM))]# delete the mean
  colnames(FM)    <- gsub("SMOD","",colnames(FM))
  colnames(FM)[1] <- "PLOT"
  
  tmp <- reshape2::melt(FM,id.vars=1:2,variable.name="YEAR",value.name="SMOD")
  tmp <- tmp %>% drop_na()
  tmp$SITE <- substr(tmp$PLOT,1,3)
  tmp$COM  <- tmp$COMcor # already corrected in part 1
  
  # add degree days estimated at Fluxalp
  tmp$DD <- NA
  setwd(WD); DD  <- read.csv("FLUXALP.DD.csv")
  
  for (i in 1:nrow(tmp)){
    print(i)
    tmpYEAR <- which(colnames(DD)==paste0("X",tmp$YEAR[i]))
    tmp$DD[i] <- DD[tmp$SMOD[i],tmpYEAR]
  }
  hist(tmp$DD)
  
  # select Lautaret-Galibier area : RON - selected already see part 1
  table(tmp$SITE)
  tmp <- tmp[,c("SITE","COM","PLOT","YEAR","SMOD","DD")]
  
  setwd(WD);write.csv(tmp,"BIOCLIM.DD.csv",row.names=F)
  
  
# III. BINARY SNOW COVER MAPS & fSCA estimates from SATELITTE data ----
      # the original files used in Dedieu & al. were not found
      # I downloaded the Theia Snow products for the same dates and other years
      # Catalogue Theia Land - https://catalogue.theia-land.fr/ 
      # Email UGA + Alchemilla@2024
      # In 2015, there were 12 dates from March23 to August30 in Dedieu 
  myYEAR <- 13:23
  
  # III.A PERIOD 2013-2023 using THEIA -> SITE.REC.fSCA.csv & THE.fSCA.csv ----  
    # III.A.1. PREPARE SCAXX.L1 (level 1) - from various RS products ----
    YEAR <- 2023 # loop over year manually
    myREF30 <- my3REF30
    myCONT  <- my3CONT
         # SPOT4TAKE5 (only 2015 for now) ----
  if(YEAR==2015){
    DIR.SPOTndsi <- "D:/DATA_ROCHENOIRE/SPOT4TAKE5_2015/NDSI_10m/"
    setwd(DIR.SPOTndsi)
    
    LF  <- list.files(pattern="NDSI",full=T)[1:8] # only the 8 first files are relevant
    LFs <- list.files(pattern="NDSI",full=F)[1:8]
    
    tmp  <- terra::rast(LF)
    tmp1 <- terra::project(tmp,myREF30)
    m    <- c(-Inf, 0.4, 0, 0.4, +Inf, 1)
    rcl  <- matrix(m, ncol=3, byrow=TRUE)
    SPOT <- terra::classify(tmp1,rcl) 
    
    # tmp10<- terra::project(tmp,myREF10)
    # tmp20 <- terra::classify(tmp10,rcl)
    names(SPOT) <- substr(LFs,11,20)
  }
         # SENTINEL-2 -> SENT ----
  LF   <- list.files("D:/DATA_ROCHENOIRE/SENTINEL_SNOW/",full=T,pattern=paste0("_",YEAR))
  LFs  <- list.files("D:/DATA_ROCHENOIRE/SENTINEL_SNOW/",pattern=paste0("_",YEAR))
  TMP  <- unlist(lapply(strsplit(LFs,"_"),function(x) substr(x[2],1,8))) 
  
  # import GL
  SELGL <- grep("GL",LFs)
  if(length(SELGL>0)){
  tmp       <- terra::rast(LF[grep("GL",LFs)])
  SENTGL    <- terra::project(tmp,myREF30,method="near")
  DOY.SENT  <- lubridate::yday(as.Date(TMP,format="%Y%m%d"))[grep("GL",LFs)]
  DAT.SENT  <- as.Date(TMP,format="%Y%m%d")[grep("GL",LFs)]
  names(SENTGL) <- DAT.SENT
  }else{
    SENTGL <- NULL
  }
  # import GK
  SELGK <- grep("GK",LFs)
  if(length(SELGK>0)){
    tmp      <- terra::rast(LF[grep("GK",LFs)])
    SENTGK   <- terra::project(tmp,myREF30,method="near")
    DOY.SENT <- lubridate::yday(as.Date(TMP,format="%Y%m%d"))[grep("GK",LFs)]
    DAT.SENT <- as.Date(TMP,format="%Y%m%d")[grep("GK",LFs)]
    names(SENTGK) <- DAT.SENT
  }else{
    SENTGK <- NULL
  }
  
  if(!is.null(SENTGL) & !is.null(SENTGK))  SENT <- c(SENTGL,SENTGK)
  if(is.null(SENTGL) & !is.null(SENTGK))   SENT <- SENTGK
  if(!is.null(SENTGL) & is.null(SENTGK))   SENT <- SENTGL
  if(is.null(SENTGL) & is.null(SENTGK))    SENT <- NULL
  
  if(!is.null(SENT)) SENT <- SENT[[order(names(SENT))]]
  
         # LANDSAT -> LAND ----
  LF   <- list.files("D:/DATA_ROCHENOIRE/LANDSAT_SNOW/",full=T,pattern=paste0("XS_",YEAR))
  LFs  <- list.files("D:/DATA_ROCHENOIRE/LANDSAT_SNOW/",pattern=paste0("XS_",YEAR))
  TMP  <- unlist(lapply(strsplit(LFs,"_"),function(x) substr(x[2],1,8)))  
  
  # import GL
  SELGL <- grep("GL",LFs)
  if(length(SELGL>0)){
  tmp <- terra::rast(LF[grep("GL",LFs)])
  LANDGL <- terra::project(tmp,myREF30,method="near")
  DOY.LAND <- lubridate::yday(as.Date(TMP,format="%Y%m%d"))[grep("GL",LFs)]
  DAT.LAND <- as.Date(TMP,format="%Y%m%d")[grep("GL",LFs)]
  names(LANDGL) <- DAT.LAND
  }else{
    LANDGL <- NULL
  }
  
  # import GK
  SELGK <- grep("GK",LFs)
  if(length(SELGK>0)){
  tmp      <- terra::rast(LF[grep("GK",LFs)])
  LANDGK   <- terra::project(tmp,myREF30,method="near")
  DOY.LAND <- lubridate::yday(as.Date(TMP,format="%Y%m%d"))[grep("GK",LFs)]
  DAT.LAND <- as.Date(TMP,format="%Y%m%d")[grep("GK",LFs)]
  names(LANDGK) <- DAT.LAND
  }else{
    LANDGK <- NULL
  }
  
  if(!is.null(LANDGL) & !is.null(LANDGK))  LAND <- c(LANDGL,LANDGK)
  if(is.null(LANDGL) & !is.null(LANDGK))   LAND <- LANDGK
  if(!is.null(LANDGL) & is.null(LANDGK))   LAND <- LANDGL
  if(is.null(LANDGL) & is.null(LANDGK))    LAND <- NULL
  
         # FINAL ASSEMBLAGE ----
  if(YEAR==2015)              LANDSENT <- c(LAND,SPOT)
  if(YEAR<=2017 & YEAR!=2015) LANDSENT <- LAND
  if(YEAR>=2018 & YEAR<2022)  LANDSENT <- c(LAND,SENT)
  if(YEAR>=2022)              LANDSENT <- SENT
  
  LANDSENT <- LANDSENT[[order(names(LANDSENT))]]
  LANDSENT <- terra::subst(LANDSENT,from=c(0,100,205),to= c(0,1,-1))
  LANDSENT[is.na(LANDSENT)]<-(-1)

         # selection of good dates ----
  nlyrs <- dim(LANDSENT)[3]
  for (i in 1:nlyrs) terra::coltab(LANDSENT[[i]]) <- coltb
  if(nlyrs>=16){
  graphics.off()
    windows(20,15);terra::plot(LANDSENT[[1:16]],nc=5,nr=4,fun=function() {terra::plot(myCONT,add=T,border="white");  terra::plot(OBS.PTS,add=T,col="red") })
    windows(20,15);terra::plot(LANDSENT[[17:nlyrs]],nc=5,nr=4,fun=function() {terra::plot(myCONT,add=T,border="white");  terra::plot(OBS.PTS,add=T,col="red") })
  }else{
    graphics.off()
    windows(20,15);terra::plot(LANDSENT,nc=5,nr=4,fun=function() {terra::plot(myCONT,add=T,border="white");  terra::plot(OBS.PTS,add=T,col="red") }) 
  }
         # SAVE SCAXX.L1.tif ----
  ## 2013 OK remise en neige tardive en mai, start on May 20
  LANDSENT <-LANDSENT
  setwd(WD);terra::writeRaster(LANDSENT,"SCA13.L1.tif",overwrite=T) 
  
  ## 2014 OK date du 27 mars à rejeter car trop de nuages
  LANDSENT <-LANDSENT
  setwd(WD);terra::writeRaster(LANDSENT,"SCA14.L1.tif",overwrite=T) 
  
  ## 2015 OK
  LANDSENT <-LANDSENT
  setwd(WD);terra::writeRaster(LANDSENT,"SCA15.L1.tif",overwrite=T) 
  
  ## 2016 OK remise en neige en mars, après fonte très partielle en février mais concernant surtout l'adret de Chaillol, 21 janvier très nuageux
  LANDSENT <-LANDSENT[[-2]]
  setwd(WD);terra::writeRaster(LANDSENT,"SCA16.L1.tif",overwrite=T) 
  
  ## 2017 OK but remise en neige en février, après fonte très partielle en janvier mais concernant surtout l'adret de Chaillol, date produite en GK du 2 juillet inutile, gros gap entre mi-juin et mijuillet ! enlver aussi le 4 septembre car remise en neige au nord du BV 
  LANDSENT <-LANDSENT[[-c(16)]]
  setwd(WD);terra::writeRaster(LANDSENT,"SCA17.L1.tif",overwrite=T)
  
  ## 2018 OK déneigement tardif, belle image le 26 juin
  LANDSENT <-LANDSENT
  setwd(WD);terra::writeRaster(LANDSENT,"SCA18.L1.tif",overwrite=T)
  
  ## 2019 OK fonte partielle en mars et grosse remise en neige en avril ! pas de bonnes images en mai pour Landsat mais OK avec Sentinel
  LANDSENT <-LANDSENT
  setwd(WD);terra::writeRaster(LANDSENT,"SCA19.L1.tif",overwrite=T)
  
  ## 2020 OK delete the first two dates
  LANDSENT <-LANDSENT
  setwd(WD);terra::writeRaster(LANDSENT,"SCA20.L1.tif",overwrite=T)
  
  ## 2021 OK delete the first two dates
  LANDSENT <-LANDSENT[[-19]]
  setwd(WD);terra::writeRaster(LANDSENT,"SCA21.L1.tif",overwrite=T)
  
  ## 2022 OK but few dates in June for Roche Noire haut
  LANDSENT <-LANDSENT
  setwd(WD);terra::writeRaster(LANDSENT,"SCA22.L1.tif",overwrite=T)
  
  ## 2023 OK but remove the last two as a snowfall in september
  LANDSENT <-LANDSENT[[-c(17,18)]]
  setwd(WD);terra::writeRaster(LANDSENT,"SCA23.L1.tif",overwrite=T)
  
         # graphical checking ----
  myY = 15 # pick a year
  SCAtmp.L1 <- terra::rast(paste0("SCA",myY,".L1.tif"))
  SCAtmp.L1[is.na(SCAtmp.L1)]<-(-1)
  for (i in 1:dim(SCAtmp.L1)[3]) terra::coltab(SCAtmp.L1[[i]]) <- coltb
  
  graphics.off();windows(14,16)
  terra::plot(SCAtmp.L1,mar=c(0.1, 0.1, 0.1, 0.1),loc.main="bottomleft",col.main="black",cex.main=1.5,axes=F,alpha=my3ALPHA,maxnl=16,fun=function() {
    terra::plot(my3HILL,add=T,alpha=0.5,col=grey(0:100 / 100),legend=F);
    terra::plot(my3CONT,add=T,border="black"); 
    terra::plot(terra::vect("SAG3lines.shp"),add=T,col="gray");
    terra::plot(OBS.PTS[2],add=T,pch=22,bg="blue")
  })
  
    # III.A.2. PREPARE SCAXX.L2 (level 2) - applying the sequential rules ----
        # loop over years ----
  myYEAR<-13:23
  for(j in 1:length(myYEAR)){
    print(j)
    
    setwd(WD)
    SCA.L1 <- terra::rast(paste0("SCA",myYEAR[j],".L1.tif"))
    SCA.L2 <- SCA.L1
    nlyrs  <- dim(SCA.L2)[3]
    
    # first run
    for (i in 2:nlyrs){
  print(i)
    # RULE1. if snow at t+1 then snow at t
   SEL1 <- which(SCA.L1[[i]][]==1)
   if (length(SEL1)>0) SCA.L2[[i-1]][SEL1] <- 1 
   
   # RULE2. if land at t then land at t-1
   SEL1 <- which(SCA.L1[[i-1]][]==0)
   if (length(SEL1)>0) SCA.L2[[i]][SEL1] <- 0  
  }  
    # second run this could be repeated two times - again on L2
    for (i in 2:nlyrs){
    print(i)
    # RULE1. if snow at t+1 then snow at t
    SEL1 <- which(SCA.L2[[i]][]==1)
    if (length(SEL1)>0) SCA.L2[[i-1]][SEL1] <- 1   
    
    # RULE2. if land at t then land at t-1
    print(i)
    SEL1 <- which(SCA.L2[[i-1]][]==0)
    if (length(SEL1)>0) SCA.L2[[i]][SEL1] <- 0  
  }
     
  setwd(WD);terra::writeRaster(SCA.L2,paste0("SCA",myYEAR[j],".L2.tif"),overwrite=T)     
  
  }
  
        # graphical checking ----
  myY = 15 # pick a year
  SCAtmp.L2 <- terra::rast(paste0("SCA",myY,".L2.tif"))
  SCAtmp.L2[is.na(SCAtmp.L2)]<-(-1)
  for (i in 1:dim(SCAtmp.L2)[3]) terra::coltab(SCAtmp.L2[[i]]) <- coltb
  
  graphics.off();windows(14,16)
  terra::plot(SCAtmp.L2,mar=c(0.1, 0.1, 0.1, 0.1),loc.main="bottomleft",col.main="black",cex.main=1.5,axes=F,alpha=my3ALPHA,maxnl=16,fun=function() {
    terra::plot(my3HILL,add=T,alpha=0.5,col=grey(0:100 / 100),legend=F);
    terra::plot(my3CONT,add=T,border="black"); 
    terra::plot(terra::vect("SAG3lines.shp"),add=T,col="gray");
    terra::plot(OBS.PTS[2],add=T,pch=22,bg="blue")
  })

    # III.A.3. PREPARE SCAXX.L3 (level 3) - gap filling of NA ----
        # loop over years ----
  myYEAR<-13:23
  for(j in 1:length(myYEAR)){
    print(j)
    
    setwd(WD)
    SCA.L2 <- terra::rast(paste0("SCA",myYEAR[j],".L2.tif"))
    SCA.L2[is.na(SCA.L2)]<-(-1)
    SCA.L3 <- SCA.L2
    nlyrs  <- dim(SCA.L2)[3]
  
    for (i in 1:nlyrs){
    
    print(i)
    TMP <- SCA.L2[[i]]
    if(length(which((TMP[]==-1)))>0){
    TMP[TMP[]==-1] <- NA
    plot(TMP,colNA="black")
    table(TMP[],useNA="ifany")
    OBJ <- length(which(is.na(TMP[])))
    while(OBJ>0){
      TMP  <- terra::focal(TMP, w=3, fun="modal",na.rm=T,na.policy="only",silent=F) 
      OBJ  <- length(which(is.na(TMP[])))
    }
    SCA.L3[[i]] <- TMP
    }else{
    SCA.L3[[i]] <- SCA.L2[[i]]
    }
  }
    names(SCA.L3) <- names(SCA.L2)
    setwd(WD);terra::writeRaster(SCA.L3,paste0("SCA",myYEAR[j],".L3.tif"),overwrite=T)     
    
  }
  
  # ERROR at j=11 / 2023
  
  
        # graphical checking ----
  myY = 23 # pick a year
  SCAtmp.L3 <- terra::rast(paste0("SCA",myY,".L3.tif"))
  SCAtmp.L3[is.na(SCAtmp.L3)]<-(-1)
  for (i in 1:dim(SCAtmp.L3)[3]) terra::coltab(SCAtmp.L3[[i]]) <- coltb
  
  graphics.off();windows(14,16)
  terra::plot(SCAtmp.L3,mar=c(0.1, 0.1, 0.1, 0.1),loc.main="bottomleft",col.main="black",cex.main=1.5,axes=F,alpha=my3ALPHA,maxnl=16,fun=function() {
      terra::plot(my3HILL,add=T,alpha=0.5,col=grey(0:100 / 100),legend=F);
      terra::plot(myCONT,add=T,border="black"); 
      terra::plot(terra::vect("SAG3lines.shp"),add=T,col="gray");
      terra::plot(OBS.PTS[2],add=T,pch=22,bg="blue")
    })
  
    # III.A.4. fSCA and DD estimates on L1 -> THE.fSCA.csv (20240718)  ----
  SCA.l <- list(); SCA.l[[1]]<-SCA.l[[2]]<-SCA.l[[3]]<-list()
  names(SCA.l)<- WS
  SCAf <- function(x) length(which(x==1))
  SFAf <- function(x) length(which(x==0))
  NAAf <- function(x) length(which(x==-1))
  setwd(WD); DD  <- read.csv("FLUXALP.DD.csv")
  
  for(j in 1:length(myYEAR)){
    print(j)
    
    setwd(WD)
    tmp <- terra::rast(paste0("SCA",myYEAR[j],".L1.tif"))
    tmp[is.na(tmp)]<-(-1)
    nlyrs  <- dim(tmp)[3]
  
    # extract for each WS
    for (k in 1:3){ 
      EXTtmp <- terra::extract(tmp,my3CONT[k,])[,-1]
      
      df <- cbind(
      apply(EXTtmp,2,SCAf),
      apply(EXTtmp,2,SFAf),
      apply(EXTtmp,2,NAAf)
      )
      
    df <- data.frame(
                   YEAR   = lubridate::year(as.Date(names(tmp))),
                   DOY    = lubridate::yday(as.Date(names(tmp))),
                   snow   = apply(EXTtmp,2,SCAf),
                   land   = apply(EXTtmp,2,SFAf),
                   NAval  = apply(EXTtmp,2,NAAf)
                      )

    # save in a list
    SCA.l[[k]][[j]] <- df
    }

  }
  
  tmp.RNPI <- do.call(rbind,SCA.l[[1]]); tmp.RNPI$WS<-"RNPI"
  tmp.LAUZ <- do.call(rbind,SCA.l[[2]]); tmp.LAUZ$WS<-"LAUZ"
  tmp.MAND <- do.call(rbind,SCA.l[[3]]); tmp.MAND$WS<-"MAND"
  
  all.equal(rownames(tmp.RNPI),rownames(tmp.LAUZ),rownames(tmp.MAND))
  
  # assemble
  tmp.THE <- data.frame(YEAR = tmp.RNPI$YEAR,
                        DOY  = tmp.RNPI$DOY,
                        RNPI = 100*tmp.RNPI$snow/NCELL30.WS3[1],
                        LAUZ = 100*tmp.LAUZ$snow/NCELL30.WS3[2],
                        MAND = 100*tmp.MAND$snow/NCELL30.WS3[3],
    TOT    = round(100*(tmp.MAND$snow+tmp.LAUZ$snow+tmp.RNPI$snow)/sum(NCELL30.WS3),1),
    fNA   = round(100*(tmp.MAND$NAval+tmp.LAUZ$NAval+tmp.RNPI$NAval)/sum(NCELL30.WS3),1)
                        )
  head(tmp.THE)
  dim(tmp.THE)
  hist(tmp.THE$fNA)
  
  # add degree days
  setwd(WD); DD <- read.csv("FLUXALP.DD.csv")
  
  tmp.THE$DD <- NA
  for (i in 1:nrow(tmp.THE)){
    tmpYEAR <- which(colnames(DD)==paste0("X",tmp.THE$YEAR[i]))
    tmp.THE$DD[i] <- DD[tmp.THE$DOY[i],tmpYEAR]
  }
  
  
  # Remove images with NA 
  
  # Select a threshold for fNAval and DOY after 1st of October (DOY 275)
  THR = 5 # meaning 5 % of clouds or NA values in the WS
  SEL <- which(tmp.THE$fNA>=THR) # 75 images
  
  # remove above 40% after the 1st of August that corresponds to transcient late snowfalls
  OUT  <- which(tmp.THE$DOY>=215 & tmp.THE$fSCA>40)
  df <- tmp.THE[-c(SEL,OUT),]
  df <- reshape2::melt(df,id.vars=c(1,2,7,8),variable.name="WS",value.name="fSCA")
  df$SRC <- "THE"
  
  df <- df %>% drop_na()
  
  df <- df[,c("YEAR","DOY","DD","WS","fSCA","fNA","SRC")] 
  
  dim(df)
  table(df$WS)
  
  # sensitivity
  test <- rep(NA,16)
  for(i in 0:15) test[i]<- length(which(tmp.THE$fNAval<=i))
  graphics.off();plot(0:15,100*test/nrow(tmp.THE),ylim=c(0,100))
  
  setwd(WD); write.csv(df,"THE.fSCA.csv",row.names=F)
  setwd(WD); df <- read.csv("THE.fSCA.csv")  

    # III.A.5. PREPARE SCAXX.SMODs using BEYOND MEDIANS (TO BE FIXED) ----
  # the sequence if fine up to the 10th of July - layer 6
  # need to calculate the border to the polygon of snow
  # https://stackoverflow.com/questions/74153263/calculate-distance-of-points-to-polygon-boundary-using-terra-package-in-r
  # need to create polygons of snow-covered cells
        # loop over years (TO BE IMPROVED) ----
  myREF30 <- my3REF30
        for (j in 1:length(myYEAR)){
  
    print(j)
    
    setwd(WD)
    SCA.L3 <- terra::rast(paste0("SCA",myYEAR[j],".L3.tif"))
    SCA.L3[is.na(SCA.L3)]<-(-1)
    nlyrs  <- dim(SCA.L3)[3] 
  
  LAND <- DIST <- DISTtl <- MELT.l <- SMOD.l <- SMODtl.l <- SMODf.l <- list()
  # SMODtl.l[[1]] <- SMOD.l[[1]] <-   SMODf.l <- MELT.l[[1]] <- NULL
  
  # DO I NEED TO CROP HERE ?
    SCA.L3m <- terra::mask(SCA.L3,myCONT)
    SCA.L3m <- SCA.L3
  
  SCA.DOY <- lubridate::yday(as.Date(names(SCA.L3m),format="%Y-%m-%d"))
  
 for (i in 1:(nlyrs-1)){
    print(paste(j,i,sep="-"))
   
    DOY1      <- SCA.DOY[i]
    DOY2      <- SCA.DOY[i+1]
    
    # find the area that melts between t and t+1
    tmp <- myREF30;tmp[]<-NA
    tmp[which(SCA.L3m[[i]][]==1 & SCA.L3m[[i+1]][]==0)]<-1
    tmp1 <- terra::as.polygons(tmp,dissolve=T)
    MELT.l[[i+1]] <-tmp
    
    # distance au nearest snowpatch at t+1 for the area that melted between t and t+1
    tmp         <- terra::mask(SCA.L3m[[i+1]], tmp1, inverse=TRUE)
    DIST[[i]]   <- terra::distance(tmp)*MELT.l[[i+1]]
    # graphics.off(); windows(10,10);terra::plot(DIST[[i]])
    
    # et aussi considérer la distance au plus proche snowfree à t ! -> DISTtl
    tmp <- myREF30;tmp[]<-NA
    tmp[which(SCA.L3m[[i]][]==0)]<-1
    tmp1 <- terra::as.polygons(tmp,dissolve=T)
    test        <- terra::mask(SCA.L3m[[i+1]], tmp1, inverse=FALSE)
    DISTtl[[i]] <- terra::distance(test)*MELT.l[[i+1]]
    # graphics.off(); windows(10,10);terra::plot(DISTtl[[i]])
  
  # Fin the minimum between DISTt and DIStl
    #graphics.off(); windows(10,10);terra::plot(DIST[[i]])
    #graphics.off(); windows(10,10);terra::plot(DISTtl[[i]])
    
    ftp <- sign(DISTtl[[i]]-DIST[[i]])
    plot(ftp,col=c("red","black","green"))
    # if green (+1) distance to snowpatch at t+1 is shorter -> use DIST to compute SMOD
    # if red (-1) distance to snowfree at t is shorter -> use DISTtl to compute SMOD
  
    bin4dist   <- terra::subst(ftp,from=c(-1,0,1),to=c(0,1,1))
    bin4disttl <- terra::subst(ftp,from=c(-1,0,1),to=c(1,0,0))
    
  # A. rescale distance and SMOD to nearest snowpatch at t+1 
  TMP       <- DIST[[i]]
  # issue here a svery long distance ecrase le reste, use IDR range !
  MIN <-   Hmisc::smedian.hilow(TMP,0.9)[1]
  MAX <-   Hmisc::smedian.hilow(TMP,0.9)[3]
  TMP1      <- (TMP - MIN)/ (MAX - MIN)
  TMP1[TMP1[]>1] <- 1 ; TMP1[TMP1[]<0]<-0
  SMOD.l[[i+1]] <- (DOY2-1 + TMP1*(DOY1-DOY2+2))
  #SMOD.l[[i+1]] <- bin4dist*(DOY2-1 + TMP1*(DOY1-DOY2+2))
  # terra::plot(SMOD.l[[i+1]],range=c(DOY1,DOY2))
  
  # B. rescale distance and SMOD to nearest snowfree at t
  TMP       <- DISTtl[[i]]
  # issue here a svery long distance ecrase le reste, use IDR range !
  MIN <-   Hmisc::smedian.hilow(TMP,0.9)[1]
  MAX <-   Hmisc::smedian.hilow(TMP,0.9)[3]
  TMP1      <- (TMP - MIN)/ (MAX - MIN)
  TMP1[TMP1[]>1] <- 1 ; TMP1[TMP1[]<0]<-0
  SMODtl.l[[i+1]] <- (DOY1+1 + TMP1*(DOY2-DOY1+2))
  # SMODtl.l[[i+1]] <- bin4disttl*(DOY1+1 + TMP1*(DOY2-DOY1+2))
  # terra::plot(SMODtl.l[[i+1]],range=c(DOY1,DOY2))
  
  # select the final SMOD based on distance
  SMODf.l[[i+1]] <- terra::app(c(SMOD.l[[i+1]],SMODtl.l[[i+1]]),mean,na.rm=T)
  # SMODf.l[[i+1]] <- terra::app(c(SMOD.l[[i+1]],SMODtl.l[[i+1]]),max,na.rm=T)

 }    
  
  #graphics.off();par(mfrow=c(4,4));for (i in 2:15) hist(DIST[[i]][])  
  #graphics.off();par(mfrow=c(4,4));for (i in 2:16) hist(SMOD.l[[i]][],main=SCA15$DOY[i])
  #graphics.off();par(mfrow=c(4,4));for (i in 3:15) hist(SMODtl.l[[i]][],main=SCA15$DOY[i])
  
  # checking figures
  # MELT <- terra::rast(MELT.l)
  # names(MELT) <-  SCA$DOY
  # graphics.off();windows(18,15)
  # terra::plot(MELT,col="blue",legend=F,mar=c(0.1, 0.1, 0.1, 0.1),loc.main="top",col.main="red",cex.main=1.5,axes=F,maxnl=16,fun=function() {terra::plot(myCONT,add=T,border="black");  terra::plot(OBS.PTS[2],add=T,pch=22,bg="red") })
  # MELT.SCA <- terra::mask(terra::app(MELT,sum,na.rm=T),myCONT)
   
    # distance to front for pixels snowmelting between i and i+1
  SMOD      <- terra::rast(SMODf.l[-1])
  
  # SCA.SMODs <- terra::mask(terra::app(SMOD,max,na.rm=T),myCONT)
  SCA.SMODs <- terra::app(SMOD,max,na.rm=T)
  
  setwd(WD);terra::writeRaster(SCA.SMODs,paste0("SCA",myYEAR[j],".SMODs.tif"),overwrite=T)  
        }  
  
    # graphical checking
  setwd(DATA_WD); terra::plot(terra::rast("SCA22.SMODs.tif"))
  
  # III.B. PERIOD 1984-2015 using DLR and SpotWorldHeritage (20240320) ----
    # data were transferred from Zacharie on 20240319
    # III.B.1. DLR data -> DLR.fSCA.csv ----
    LFs <- list.files("U:/DATA_ROCHENOIRE/SWH_DLR_ZBD/SNOW/DLR/")
    LF  <- list.files("U:/DATA_ROCHENOIRE/SWH_DLR_ZBD/SNOW/DLR/",full=T)
    length(LF)
    
    YEAR = lubridate::year(as.Date(substr(LFs,12,19),format="%Y%m%d"))
    DOY  = lubridate::yday(as.Date(substr(LFs,12,19),format="%Y%m%d"))
    
    # 0 : sol/land
    # 100: neige
    # 205: nuage
    # 255: nodata
      NAval <- matrix(NA,length(LF),3)
      snow  <- matrix(NA,length(LF),3) 
      land  <- matrix(NA,length(LF),3) 
      
    for(i in 1:length(LF)){
      print(i)
      tmp0  <- terra::rast(LF[i])
      tmp <- terra::project(tmp0,my3REF30,method="near")
      # plot(tmp);plot(my3CONT,add=T)
      names(tmp)<-"CLASS"
      
      # count the number of NA, snow and nosnow
        for (k in 1:3){
          NAval[i,k] <- length(which(terra::extract(tmp,my3CONT[k,])$CLASS>=205))
          snow[i,k]  <- length(which(terra::extract(tmp,my3CONT[k,])$CLASS==100))
          land[i,k]  <- length(which(terra::extract(tmp,my3CONT[k,])$CLASS==0))
        }
    }
      
      hist(NAval[,1])
      setwd(WD)
      df      <- data.frame(NAval=NAval[,1],snow=snow[,1],land=land[,1])
      df$YEAR <- YEAR; df$DOY<- DOY
      write.csv(df,"RNPI.DLR.csv")
      
      df$DD <- NA
      for (i in 1:nrow(df)){
        tmpYEAR <- which(colnames(DD)==paste0("X",df$YEAR[i]))
        df$DD[i] <- DD[df$DOY[i],tmpYEAR]
      }
      
      
      df      <- data.frame(NAval=NAval[,2],snow=snow[,2],land=land[,2])
      df$YEAR <- YEAR; df$DOY<- DOY
      write.csv(df,"LAUZ.DLR.csv")
      
      df      <- data.frame(NAval=NAval[,3],snow=snow[,3],land=land[,3])
      df$YEAR <- YEAR; df$DOY<- DOY
      write.csv(df,"MAND.DLR.csv")
      
      
      # add DD and reorganize -> DLR.fSCA.csv
      setwd(WD); DD <- read.csv("FLUXALP.DD.csv")
      RNPI.DLR <- data.frame(read.csv("RNPI.DLR.csv")[,-1],WS="RNPI")
      LAUZ.DLR <- data.frame(read.csv("LAUZ.DLR.csv")[,-1],WS="LAUZ")
      MAND.DLR <- data.frame(read.csv("MAND.DLR.csv")[,-1],WS="MAND") 
      all.equal(nrow(RNPI.DLR),nrow(LAUZ.DLR),nrow(MAND.DLR))
      
      # Check area covered and Add a TOT estimate
      tmp<-RNPI.DLR[,1:3]+LAUZ.DLR[,1:3]+MAND.DLR[,1:3]      
      DIM <- apply(tmp,1,sum)
      TOT.DLR <- data.frame(NAval= tmp[,1],
                            snow = tmp[,2],
                            land = tmp[,3],
                            YEAR = RNPI.DLR$YEAR,
                            DOY  = RNPI.DLR$DOY,
                            WS   = "TOT"
                            )
      SEL <- which(DIM==16580) # checking full area covered      
      tmp.DLR <- rbind(RNPI.DLR[SEL,],LAUZ.DLR[SEL,],MAND.DLR[SEL,],TOT.DLR[SEL,])
      nrow(tmp.DLR)/# 2190 images 

      
      # add degree days and compute fSCA and fNA
      tmp.DLR$DD <- NA
      for (i in 1:nrow(tmp.DLR)){
        tmpYEAR <- which(colnames(DD)==paste0("X",tmp.DLR$YEAR[i]))
        tmp.DLR$DD[i] <- DD[tmp.DLR$DOY[i],tmpYEAR]
      }
      tmp.DLR$fSCA    <- 100*tmp.DLR$snow/(tmp.DLR$snow+tmp.DLR$land+tmp.DLR$NAval)
      tmp.DLR$fNA     <- 100*tmp.DLR$NAval/(tmp.DLR$snow+tmp.DLR$land+tmp.DLR$NAval)
      tmp.DLR$SRC     <- "DLR"
    
      # reorganize and save
      colnames(tmp.DLR)
      DLR <- tmp.DLR[,c("YEAR","DOY","DD","WS","fSCA","fNA","SRC")]
      setwd(WD); write.csv(DLR,"DLR.fSCA.csv",row.names=F)
      
      # NOT AT THIS STAGE !
      # Select a threshold for fNAval and DOY after 1st of October (DOY 275)
      THR = 4 # meaning 5 % of clouds or NA values in the WS
      SEL <- which(tmp.DLR$fNAval<THR)
      test <- rep(NA,16)
      for(i in 0:15) test[i]<- length(which(tmp.DLR$fNAval<=i))
      graphics.off();plot(0:15,100*test/nrow(tmp.DLR),ylim=c(0,30))

    # III.B.2. fSCA and DD estimates -> SWH.fSCA.csv ----
      LFs <- list.files("U:/DATA_ROCHENOIRE/SWH_DLR_ZBD/SNOW/SWH/")
      LF  <- list.files("U:/DATA_ROCHENOIRE/SWH_DLR_ZBD/SNOW/SWH/",full=T)
      length(LF) # including repet !
      
      YEAR = lubridate::year(as.Date(substr(LFs,12,19),format="%Y%m%d"))
      DOY  = lubridate::yday(as.Date(substr(LFs,12,19),format="%Y%m%d"))
      
      # 0 : sol/land
      # 100: neige
      # 205: nuage
      # NA: nodata - beware this is different from DLR tif !
      NAval <- matrix(NA,length(LF),3)
      snow  <- matrix(NA,length(LF),3) 
      land  <- matrix(NA,length(LF),3) 
      
      for(i in 1:length(LF)){
        print(i)
        tmp0  <- terra::rast(LF[i])
        tmp <- terra::project(tmp0,my3REF30,method="near")
        plot(tmp);plot(my3CONT,add=T)
        names(tmp)<-"CLASS"
        
        # count the number of NA, snow and nosnow
        for (k in 1:3){
          EXTtmp <- terra::extract(tmp,my3CONT[k,])$CLASS
          NAval[i,k] <- length(which(EXTtmp>=205 | is.na(EXTtmp)==TRUE))
          snow[i,k]  <- length(which(EXTtmp==100))
          land[i,k]  <- length(which(EXTtmp==0))
        }
      }
      
      df1      <- data.frame(
        YEAR = YEAR, DOY  = DOY, WS="RNPI",
        NAval=NAval[,1], snow=snow[,1],land=land[,1]
        )
      
      df2      <- data.frame(
        YEAR = YEAR, DOY  = DOY, WS="LAUZ",
        NAval=NAval[,2], snow=snow[,2],land=land[,2]
      )
      
      df3      <- data.frame(
        YEAR = YEAR, DOY  = DOY, WS="MAND",
        NAval=NAval[,3], snow=snow[,3],land=land[,3]
      )
      
      tmp.SWH <- rbind(df1,df2,df3)
      tmp.SWH$SRC <- "SWH"
      
      # add DD and reorganize
      setwd(WD); DD <- read.csv("FLUXALP.DD.csv")
      
      tmp.SWH$DD <- NA
      for (i in 1:nrow(tmp.SWH)){
        tmpYEAR <- which(colnames(DD)==paste0("X",tmp.SWH$YEAR[i]))
        tmp.SWH$DD[i] <- DD[tmp.SWH$DOY[i],tmpYEAR]
      }
      
      tmp.SWH$fSCA    <- 100*tmp.SWH$snow/(tmp.SWH$snow+tmp.SWH$land+tmp.SWH$NAval)
      tmp.SWH$fNAval  <- 100*tmp.SWH$NAval/(tmp.SWH$snow+tmp.SWH$land+tmp.SWH$NAval)
      
      # Select a threshold for fNAval and DOY after 1st of October (DOY 275)
      THR = 4 # meaning 5 % of clouds or NA values in the WS
      hist(tmp.SWH$fNAval)
      SEL <- which(tmp.SWH$fNAval<THR)
      test <- rep(NA,16)
      for(i in 0:15) test[i]<- length(which(tmp.SWH$fNAval<=i))
      graphics.off();plot(0:15,100*test/nrow(tmp.SWH),ylim=c(0,30))
      
      SWH <- tmp.SWH[which(tmp.SWH$fNAval<=THR & tmp.SWH$DOY<=275),c("YEAR","DOY","DD","SRC","WS","fSCA")]
      hist(SWH$fSCA)
      
      setwd(WD); write.csv(SWH,"SWH.fSCA.csv",row.names=F)
      
    
  # III.C. PERIOD 1984-2023 using GEE raw NDSI -> GEE.fSCA.csv (20240718)  ----
    # III.C.1. PRELIMINARY TESTS ----
  # Searching data in earthexplorer
  # https://earthexplorer.usgs.gov/
  # login is pcholer ; password is Cardamine2021
  # https://www.usgs.gov/faqs/what-are-band-designations-landsat-satellites
  # using Landsat LT05 level2
  # big issue with badly georeferenced Landsat scenes !! coordinate reference
  # attention: Tier 2 data has not been orthorectified due to issues such as too much cloud coverage (which was the case in this location). For this particular location, there were no Tier 1 level data available for Landsat 5. 
  # si seule source d'images correcte il faudrait les orthorectifier !
  # example for the 30/07/1984
  # NDSI = VIS - SWIR / VIS + SWIR
  # VIS is the reflectance in the green visible spectrum (band 2 in Landsat TM and ETMþ, band 3 in OLI)
  # SWIR denotes the short wave infrared band, being 5 for TM and ETM+ or 6 for OLI
  # 1984 with few available images
  
  # see
  # https://www.usgs.gov/landsat-missions/normalized-difference-snow-index
  # https://eros.usgs.gov/doi-remote-sensing-activities/2022/usgs/landsat-collection-2-normalized-difference-snow-index
  # In Landsat 4-7, NDSI = (Band 2 – Band 5) / (Band 2 + Band 5)
  # In Landsat 8-9, NDSI = (Band 3 – Band 6) / (Band 3 + Band 6)
  # compare with a composite using infrared, near-infrared, and red bands (Bands 5, 4, 3 for Landsat 5 TM)
  
  # test with a T1 product 19840612
  LF <- list.files("D:/DATA_ROCHENOIRE/LANDSAT_USGS/Tier1",full=TRUE,pattern="19840612")
  GREEN <- terra::rast(LF[[grep("SR_B2",LF)]]) # green is band 2
  SWIR  <- terra::rast(LF[[grep("SR_B5",LF)]]) # SWIR is band 5
  RED   <- terra::rast(LF[[grep("SR_B3",LF)]]) # RED is band 3
  NIR   <- terra::rast(LF[[grep("SR_B4",LF)]]) # NIR is band 4
  BLUE  <- terra::rast(LF[[grep("SR_B1",LF)]]) # BLUE is band 1
  
  # create a multiLayered SpatRaster
  ILL <- 25000
  REDadj   <- (RED-ILL)/65000; REDadj[REDadj < 0] = 0; REDadj = REDadj^0.4
  GREENadj <- (GREEN-ILL)/65000; GREENadj[GREENadj < 0] = 0; GREENadj = GREENadj^0.5
  BLUEadj  <- (BLUE-ILL)/65000; BLUEadj[BLUEadj < 0] = 0; BLUEadj=BLUEadj^0.8
  TCtest   <- terra::project(c(REDadj,GREENadj,BLUEadj),my3REF30)
  terra::plotRGB(TCtest, scale=1, zlim=c(0, 1))
  
  test <- terra::project(c(GREEN,SWIR,RED,NIR,BLUE),my3REF30)
  NDSI <- (test[[1]]-test[[2]])/(test[[1]]+test[[2]])
  graphics.off();plot(NDSI);plot(my3CONT,add=T)
  
  # how to rescale to 
  terra::plotRGB(test,r=2,g=4,b=3,scale=65535,smooth=T,stretch="lin");plot(my3CONT,add=T)
  terra::plotRGB(test,r=3,g=4,b=1,scale=65535,stretch="hist",smooth=T);plot(my3CONT,add=T) # weird

  # thresholding with NDSI=0.154
  THR     <- 0.3
  m       <- c(-Inf, THR, 0,  THR, +Inf, 1)
  rclmat  <- matrix(m, ncol=3, byrow=TRUE)
  NDSIbin <- terra::classify(NDSI,rclmat)
  graphics.off();terra::plot(NDSIbin,col=c("darkgreen","lightgray"));plot(my3CONT,add=T)
  
  100*length(which(NDVI[]>=0.3))/NCELL30
  
    # use GEE instead for the same date !
    LF <- list.files("D:/DATA_ROCHENOIRE/GEE/",full=TRUE)
    NDSI <- terra::project(terra::rast(LF[1]),myREF30)
    plot(NDSI);plot(myCONT,add=T)
  
    # III.C.2. fSCA & DD estimates (20240718) ----
    # the first version (Mars 2024) including spectral corrections yields some bizarreries
    LFs <- list.files("D:/DATA_ROCHENOIRE/LANDSAT_USGS/GEE_from_ARTHUR/NDSI_Rochenoire/")
    LF  <- list.files("D:/DATA_ROCHENOIRE/LANDSAT_USGS/GEE_from_ARTHUR/NDSI_Rochenoire/",full=T)
    
    # the second version (July 2024) includes raw NDSI data and a visual selection of cloud mask
    # une centaine d'images cloudfree
    LFs <- list.files("D:/DATA_ROCHENOIRE/LANDSAT_USGS/GEE_from_ARTHUR/NDSI_BRUT_CLOUDFREE/")
    LF  <- list.files("D:/DATA_ROCHENOIRE/LANDSAT_USGS/GEE_from_ARTHUR/NDSI_BRUT_CLOUDFREE/",full=T)
    
    # check availability per year
    barplot(table(as.numeric(substr(LFs,17,20))))
    
    setwd("D:/DATA_ROCHENOIRE/LANDSAT_USGS/GEE_from_ARTHUR/NDSI_WS3_CLOUDFREE/")
  for(i in 1:length(LF)){
      print(i)
      NAME <- gsub("NDSI.tif","NDSI_WS3.tif",LFs[i])
      tmp  <- terra::rast(LF[i])
      tmp1 <- terra::project(tmp,my3REF30)
      terra::writeRaster(tmp1,NAME,overwrite=T)
  }
    
    setwd("D:/DATA_ROCHENOIRE/LANDSAT_USGS/GEE_from_ARTHUR/NDSI_WS3_CLOUDFREE/")
    LF <- list.files(full=T,pattern="19840612")
    graphics.off();windows(10,10);terra::plot(terra::rast(LF));terra::plot(my3CONT,add=T)
    
  # compute fSCA
    LF  <- list.files("D:/DATA_ROCHENOIRE/LANDSAT_USGS/GEE_from_ARTHUR/NDSI_WS3_CLOUDFREE/",full=T)
    LFs  <- list.files("D:/DATA_ROCHENOIRE/LANDSAT_USGS/GEE_from_ARTHUR/NDSI_WS3_CLOUDFREE/",full=F)
    
    SNval    <- as.data.frame(matrix(NA,length(LF),3))
    colnames(SNval)<- values(my3CONT)$Name
    NAval <- fSCA   <- SNval
    NDSI.THR <- 0.1
    
    for (i in 1:length(LF)){
      # print(i)
      tmp <- terra::rast(LF[i])
      for (k in 1:3){
      NAval[i,k] <- length(which(is.na(terra::extract(tmp,my3CONT[k,])$NDSI)==TRUE))
      SNval[i,k]  <- length(which(terra::extract(tmp,my3CONT[k,])$NDSI>=NDSI.THR))
      }
    }
    SNval$TOT <- apply(SNval,1,sum)
    
    for (k in 1:3) fSCA[,k] <- round(100*SNval[,k]/NCELL30.WS3[k],1)
    fSCA$TOT <- round(100*SNval$TOT/sum(NCELL30.WS3),1)
    
    fSCA$YEAR = lubridate::year(as.Date(substr(LFs,17,24),format="%Y%m%d"))
    fSCA$DOY  = lubridate::yday(as.Date(substr(LFs,17,24),format="%Y%m%d"))

    # % of NAval
    fSCA$fNA <- round(100*apply(NAval,1,sum)/sum(NCELL30.WS3),1)
    
    # Compute the thermal time
    setwd(WD); DD <- read.csv("FLUXALP.DD.csv")   
    
    fSCA$DD <- NA
    for (i in 1:nrow(fSCA)){
      tmpYEAR <- which(colnames(DD)==paste0("X",fSCA$YEAR[i]))
      fSCA$DD[i] <- round(DD[fSCA$DOY[i],tmpYEAR],2)
    }
    
    
    # Final selection
    # remove above 40% after the 1st of August that corresponds to transcient late snowfalls
    OUT  <- which(fSCA$DOY>=215 & fSCA$TOT>40)
    fSCA <- fSCA[-OUT,]
    
    graphics.off()
    windows(10,10)
    plot(fSCA$DD,fSCA$MAND,xlim=c(0,1000),pch=21,bg="red",cex=1.5)
    points(fSCA$DD,fSCA$LAUZ,pch=21,bg="lightblue",cex=1.5)
    points(fSCA$DD,fSCA$RNPI,pch=21,bg="green",cex=1.5) 
    
    # export results
    df <- reshape2::melt(fSCA,id.vars=5:8,variable.name="WS",value.name="fSCA")
    df$SRC <- "GEE"
    df <- df[,c("YEAR","DOY","DD","WS","fSCA","fNA","SRC")]
    table(df$WS)
    setwd(WD); write.csv(df,"GEE.fSCA.csv",row.names=F)
    
    setwd(WD);df <- read.csv("GEE.fSCA.csv")
    # this dataset could be used to evaluate the DD-fSCA model ?
    
    # NOT TO RUN USE OF PANCHROMATIC IMAGES (LEAVE IT FOR NOW ?) ----
  # 26-02 ; 23-03; 08-04
  # qui a produit ces images Landsat 8 à 15m
  # should be the 15 meters (panchromatic) product
  # OLI collects data for visible, near infrared, and short wave infrared spectral bands as well as a panchromatic band.
  test <- terra::rast("S:/LECA/PLATEAU-PASTIS/GIS_DATA/Alpes/IMAGERY/L8_Guisane/Snow_bin/L8_apr_18_2013_definiens.img") # 15m resolution ? why
  terra::plot(test)
  
    # II.A. fSCA for available dates 
  # III.D. SMOD and SMOD.fSCA from various sources ----
    # III.D.1. SMOD.fSCA.PC from this study  (20240719) ----  
    # from this study using the 2013:2023 sequence
    # better to do so after a smoothing if the beyond median routine stays as is
  myYEAR  <- 13:23
  myCONT  <- my3CONT
  myPROBS <- seq(0.1,0.9,0.1)
  myPROBS.NAME <- gsub("\\.","",as.character(1-myPROBS)) # beware the 1- !
  SCA.f.l <- list()
  SMOD.fSCA.PC <- data.frame(matrix(NA,length(myYEAR),length(myPROBS)))
  colnames(SMOD.fSCA.PC) <- paste0("fSCA",myPROBS.NAME)
  rownames(SMOD.fSCA.PC) <- 2000+myYEAR
  
  for (i in 1:length(myYEAR)){
    print(i)
    setwd(WD) 
    TMP0        <- terra::mask(terra::rast(paste0("SCA",myYEAR[i],".SMODs.tif")),myCONT)
    TMP1        <- round(terra::focal(TMP0,w=3,fun="mean",na.rm=T))
   
  # issue here as the f05 area varies too much
   # THR <- as.numeric(PARAM[i,4])
   # m    <- c(-Inf, THR, NA, THR, +Inf, 1)
   # rcl  <- matrix(m, ncol=3, byrow=TRUE)
  
  # need to find another routine here with median value
  THR <- quantile(TMP1[],probs=myPROBS,na.rm=T)
  SMOD.fSCA.PC[i,] <- THR
  
  # THR  <- Hmisc::smedian.hilow(TMP1[],0.8,na.rm=T)
  # need to adjust the lower threshold here
    RECLASS <- list()
    for (k in 1:length(myPROBS)){
      rcl           <- matrix(c(-Inf, THR[k], NA, THR[k], +Inf, 1), ncol=3, byrow=TRUE)
      RECLASS[[k]]  <- terra::classify(TMP1,rcl)
      # SCA.f05pol <- terra::as.polygons(SCA.f05==1,dissolve=T) 
    }
  
  SCA.f.l[[i]] <- terra::rast(RECLASS)
  }
  
  # reorder the list by myPROBS
  for (j in 1:length(myPROBS)){
    TMP     <- lapply(SCA.f.l, function(x) x[[j]])
    TMP.sca <- terra::rast(TMP)
    names(TMP.sca) <- 2013:2023
    NAME    <- paste0("SCA.f",myPROBS.NAME[j],".tif")
    setwd(WD);terra::writeRaster(TMP.sca,NAME,overwrite=T)
  }
  
  setwd(WD);write.csv(SMOD.fSCA.PC,"SMOD.fSCA.PC.csv",row.names=F)
  
  # get areas (to update)
  AREA <- matrix(NA,length(myYEAR),length(myPROBS))
  
  # issue here as the area does not vary according to threshold - WHY ?
  for (i in 1:11) {
    AREA[i,1] <- sum(SCA.flow[[i]][],na.rm=T)/NCELL30tot.WS3
    AREA[i,2] <- sum(SCA.fmed[[i]][],na.rm=T)/NCELL30tot.WS3
    AREA[i,3] <- sum(SCA.fhig[[i]][],na.rm=T)/NCELL30tot.WS3
  }
  (AREA)
  

      # graphical checkings ----
  
  graphics.off();windows(8,6)
  par(mfrow=c(1,2))
  terra::plot(terra::rast("SCA.f09.tif"), smooth=F,col="deepskyblue3",mar=c(0.1, 0.1, 0.1, 0.1),legend=F,axes=F,loc.main="bottomleft",cex.main=1.5,
  fun=function() {
    terra::plot(terra::mask(my3HILL,my3CONT),add=T,alpha=0.5,col=grey(0:100 / 100),legend=F);
    terra::plot(terra::vect("SAG3lines.shp"),add=T,col="gray");
    terra::plot(myCONT,add=T,border="black");  
    terra::plot(OBS.PTS[2],add=T,pch=22,bg="red")
    })
              
  LF <- list.files(pattern=".SMODs.tif")
  SMODs <- terra::rast(LF);names(SMODs) <- 2013:2023
  
  graphics.off();windows(8,6)
  par(mfrow=c(1,2))
  terra::plot(terra::mask(SMODs,myCONT), smooth=F,col=SMOD_pal(100),mar=c(0.1, 0.1, 0, 0.1),loc.main="bottomleft",legend=F,axes=F,range=c(40,241),fun=function() {
    terra::plot(terra::mask(HILL,myCONT),add=T,alpha=0.1,col=grey(0:100 / 100),legend=F);
    terra::plot(terra::vect("SAG3lines.shp"),add=T,col="gray");
    terra::plot(myCONT,add=T,border="darkgray");  
    terra::plot(OBS.PTS[1:2],add=T,pch=22,bg="blue")
    })
  
    
    # III.D.2. SMOD.fSCA.AB from Arthur's annual Landsat estimates ----
  LF <- list.files("U:/SMOD_EUALPS/",full=T)
  
  TMP <- list()
  for (i in 1:length(LF)) {
    print(i)
    TMP[[i]] <- terra::project(terra::rast(LF[[i]]),myREF30)
  }
  
  myYEAR <- 1984:2023
  SMOD.ALLA <- terra::mask(terra::rast(TMP),myCONT)
  names(SMOD.ALLA)<- myYEAR
  graphics.off(); windows(20,15);plot(SMOD.ALLA)
  
  SCA.f05.l <- SCA.f01.l <- SCA.f09.l <- list()
  SMOD.fSCA.AB <- data.frame(matrix(NA,length(myYEAR),3))
  colnames(SMOD.fSCA.AB) <- c("fSCA09","fSCA05","fSCA01")
  rownames(SMOD.fSCA.AB) <- 1984:2023
  
  Q01<-0.09;Q05<-0.49;Q09<-0.89
  for (i in 1:length(myYEAR)) {
    print(i)
    TMP1        <- round(SMOD.ALLA[[i]])

  # need to find another routine here with median value
  THR <- quantile(TMP1[],probs=c(Q01,Q05,Q09),na.rm=T)
  SMOD.fSCA.AB[i,] <- THR
  # THR  <- Hmisc::smedian.hilow(TMP1[],0.8,na.rm=T)
  # need to adjust the lower threhold here
  rcl.01  <- matrix(c(-Inf, THR[1], NA, THR[1], +Inf, 1), ncol=3, byrow=TRUE)
  rcl.05  <- matrix(c(-Inf, THR[2], NA, THR[2], +Inf, 1), ncol=3, byrow=TRUE)
  rcl.09  <- matrix(c(-Inf, THR[3], NA, THR[3], +Inf, 1), ncol=3, byrow=TRUE)
  
  SCA.f01.l[[i]]  <- terra::classify(TMP1,rcl.01)
  SCA.f05.l[[i]]  <- terra::classify(TMP1,rcl.05)
  SCA.f09.l[[i]]  <- terra::classify(TMP1,rcl.09)
  }
  names(SCA.f01.l) <-   names(SCA.f05.l)<-   names(SCA.f09.l) <- myYEAR
  
  SCA.f01<- terra::rast(SCA.f01.l); names(SCA.f01)  
  SCA.f05<- terra::rast(SCA.f05.l); names(SCA.f05)  
  SCA.f09<- terra::rast(SCA.f09.l); names(SCA.f09)  
  AREA <- matrix(NA,length(myYEAR),3)
  # issue here as the area does not vary according to thershold - WHY ?
  for (i in 1:length(myYEAR)) {
    AREA[i,1] <- sum(SCA.f01[[i]][],na.rm=T)/NCELL30
    AREA[i,2] <- sum(SCA.f05[[i]][],na.rm=T)/NCELL30
    AREA[i,3] <- sum(SCA.f09[[i]][],na.rm=T)/NCELL30
  }
  (AREA)
  
  setwd(WD);terra::writeRaster(SCA.f05,"SCA.f05.AB.tif",overwrite=T)
  setwd(WD);terra::writeRaster(SCA.f09,"SCA.f09.AB.tif",overwrite=T)
  setwd(WD);terra::writeRaster(SCA.f01,"SCA.f01.AB.tif",overwrite=T)
  setwd(WD);write.csv(SMOD.fSCA.AB,"SMOD.fSCA.AB.csv",row.names=F)
  
    # III.D.3. SMOD.fSCA.SG from Simon's SENTINEL-2 last simulations ----
  # beware: SMOD are estimated from 1st September (DOY=245) of the year before
  # need to substract 120 = 365-245
  DIR.SENT <- "D:/DATA_SNOWBED/SENTINEL/SMOD/RAW"
  setwd(DIR.SENT)
  LF <- list.files(DIR.SENT, pattern="30m.tif",full =T)
  
  FSFD.SENT <- list()
  for (i in 1:length(LF)){
    print(i)
    tmp  <- terra::rast(LF[i])
    tmp1 <- terra::project(tmp,myREF30)
    FSFD.SENT[[i]] <- tmp1-120
  }
  names(FSFD.SENT) <- 2018:2022
  
  FSFD.SENT <- terra::rast(FSFD.SENT)
  setwd(WD);terra::writeRaster(FSFD.SENT,"FSFD.SENT.tif")
  
    # III.D.4. SMOD from SPOT4Take5 as in Dedieu & al. 2016 Remote Sensing ----
  # see manuscript m69.r for details
  # SPOT4TAKE5 at 10m + Landsat 8 at 30m resolution for year 2015
  # all products were resampled at 25m resolution
  # here used as a backbone of the study

  DIR.SPOTsmod <- "D:/DATA_ROCHENOIRE/SPOT4TAKE5_2015/FSFD/"
  setwd(DIR.SPOTsmod)
  tmp <- terra::rast("FSFD_median_2015.tif") # 25m resolution
  terra::crs(tmp)<-"epsg:2154"
  terra::ext(tmp)
  SMOD.DEDIEU <- terra::project(tmp,myREF30)
  
  setwd(WD);terra::writeRaster(SMOD.DEDIEU,"SMOD.DEDIEU.tif",overwrite=T)
  terra::plot(SMOD.DEDIEU)
  terra::plot(my3CONT,add=T)
  terra::plot(OBS.PTS,col="red",add=T)
  
  # Native SPOT4-Take5 image on 5 june 2015 downloaded at 
  # https://tpm-ds.eo.esa.int/smcat/SPOT4-5Take5_ESA/3/France/Ecrins/1/
  # https://tpm-ds.eo.esa.int/smcat/SPOT4-5Take5_ESA/2/18/18/
  
  tmp <- terra::rast("D:/DATA_SPOT/SPOT4TAKE5_2015/SRC/diskA/SPOT5_TAKE5/N1C_TUILE/SPOT5_HRG2_XS_20150605_N1_TUILE_EcrinsFranceD0000B0000/SPOT5_HRG2_XS_20150605_N1_TUILE_EcrinsFranceD0000B0000.TIF")
  tmp1 <- terra::crop(tmp,my3EXT)
  terra::plotRGB(tmp1,r=4,g=3,b=1,scale=822) # same as the Dedieu paper !
  terra::plot(my3CONT,add=T,border="yellow")
  
    # III.D.5. SMOD from Landsat as in Choler & al. NCC ----
  setwd("D:/DATA_SNOWBED/SWALPS/")
  LF <- list.files(pattern="SMOD",full=T)
  
  FSFD.LAND <- list()
  for (i in 1:length(LF)){
    print(i)
    tmp  <- terra::rast(LF[i])
    tmp1 <- terra::project(tmp,myREF30)
    FSFD.LAND[[i]] <- tmp1
  }
  names(FSFD.LAND) <- c("1985-1996","2011-2022")

  FSFD.LAND <- terra::rast(FSFD.LAND)
  setwd(WD);terra::writeRaster(FSFD.LAND,"FSFD.LAND.tif")
  
    # III.D.6. SMOD from ground-level data & S2M re-analyses ----
    # III.D.7. PRELIMINARY FIGURES ----
  setwd(WD)
  nivo <- read.csv("SMOD.NIVO.csv")
  flux <- read.csv("SMOD.FLUX.csv")
  grou <- read.csv("SMOD.SOIL.csv",row.names=1)
  croc <- read.csv("SMOD.CROC.csv",row.names=1)
  
  grou["RON_CT_012",]
  windows(10,10)
  plot(unlist(grou["RON_CT_012",4:10]),nivo[2:8,2],pch=23,cex=3,bg="green",xlim=c(100,160),ylim=c(100,160),xlab="SMOD from soil temp",ylab="SMOD from nivose")
  text(unlist(grou["RON_CT_012",4:10]),nivo[2:8,2],labels=2016:2022)
  abline(0,1)
  abline(10,1,lty=2)
  abline(-10,1,lty=2)
  
  plot(flux[4:10,3],nivo[2:8,2],pch=23,cex=3,bg="green",xlim=c(105,160),ylim=c(105,160))
  text(flux[4:10,3],nivo[2:8,2],labels=2016:2022)
  abline(0,1)
  abline(50,1,lty=2)
  abline(25,1,lty=2)
  
  setwd(WD); FSFD.SPOT <- terra::rast("FSFD.SPOT.tif")
  graphics.off()
  windows(10,10)
  terra::plot(FSFD.SPOT)
  terra::plot(myCONT,add=T)
  
  setwd(WD); FSFD.SENT <- terra::rast("FSFD.SENT.tif")  
  graphics.off()
  windows(10,10)
  terra::plot(FSFD.SENT,3)
  terra::plot(myCONT,add=T,border="red")
  terra::plot(terra::vect("SAGlines.shp"),col="black",add=T)
  
  setwd(WD);myDEM <- terra::rast("myDEM.tif") 
  graphics.off()
  windows(10,10)
  terra::plot(myDEM)
  terra::plot(RNcontL93,add=T)
  terra::plot(PIcontL93[1],add=T)
  points(mySIT$X_L93,mySIT$Y_L93,bg="black",pch=21,cex=1)
  terra::plot(myRWBODY,add=T,lwd=3,col="blue")
  terra::plot(terra::vect("SAG3lines.shp"),add=T)
  
# IV SNOWMELTING model (20240718) ----
  # IV.A. First comparison between Landsat & SPOT ----
  setwd(WD); SMOD.DEDIEU <- terra::mask(terra::rast("SMOD.DEDIEU.tif"),myCONT)
  setwd(WD); SMOD.mod15  <- terra::rast(paste0("SCA15.SMODs.tif"))
  range(SMOD.mod15);range(SMOD.DEDIEU)
  graphics.off(); windows(10,10);par(mfrow=c(2,2))
  plot(SMOD.DEDIEU,range=c(50,210));plot(SMOD.mod15,range=c(50,210),colNA="black")
  hist(SMOD.DEDIEU[],breaks=seq(50,210,5));hist(SMOD.mod15[],breaks=seq(50,210,5))
  
  # Conclusion : selon toute vraisemblance un smooth a été appliqué sur les données Dedieu & al. !
  
  SMOD.mod15sm <- round(terra::focal(SMOD.mod15,w=3,fun="mean",na.rm=T))
  graphics.off(); windows(10,10);
  par(mfrow=c(2,2))
  plot(SMOD.DEDIEU,range=c(60,210));
  plot(SMOD.mod15sm,range=c(60,210))
  
  hist(SMOD.DEDIEU[],breaks=seq(50,210,5));hist(SMOD.mod15sm[],breaks=seq(50,210,5))
  
  
  df     <- terra::spatSample(SMOD.DEDIEU, size=2000, method="random",xy=TRUE)
  colnames(df)[3] <- "DEDIEU" 
  mycell <- terra::cellFromXY(SMOD.DEDIEU,df[,c("x","y")])
  df$mod15 <- terra::extract(SMOD.mod15,mycell)[,1]
  
  graphics.off()
  windows(10,10)
  plot(df$DEDIEU,df$mod15,pch=".",cex=1.5,ylab="Beyond Median",xlab="Dedieu 2016")
  abline(0,1,col="blue",lty=2,lwd=2)
  
  graphics.off()
  windows(10,5)
  hist(df$SPOT,breaks=seq(50,200,1))
  abline(v=DOY.SPOT,col="red",lwd=2)
  abline(v=DOY.LAND,col="blue")
  
  # Conclusion : there is a "stair" effect that calls for further computations than just computing the median between two available dates 
  
  # IV.B. DD-fSCA model ----
    # IV.B.1. UPDATED  procedure using nlstools on TOT only (20250315) ----
    # the selected model is a 3 parameter sigmoid
  setwd(WD); THE <- read.csv("THE.fSCA.csv")   # select training dataset
  # use only THE for training
  df <- THE
  df$fSCA <- df$fSCA/100
  
  df$WS <- as.factor(df$WS)
  levels(df$WS) <- c("3.LAUZ","1.MAND","2.RNPI","4.TOT")
  df$WS <- as.vector(df$WS)
  
  # remove points in September !
  df <- df[which(df$DOY<=245),]
  # remove a few outliers
  OUT <- which(df$fSCA>0.9 & df$DD>1000)
  if(length(OUT)>0) df <- df[-OUT,]
  
  df.l<- split(df,df$WS)
  dim(df.l[[4]]) # 128 images
  
  w=4; graphics.off(); windows(10,10); plot(df.l[[w]]$DD,df.l[[w]]$fSCA)
  
  # optimize parameters
  formula.LL3 <- as.formula(fSCA ~ a/(1 + exp(b*(log(DD)-log(c))))) # LL3
  formula.SIG <- as.formula(fSCA ~ a/(1 + exp(-b*(DD-c))))          # sigmoid
  
  LL3.nls1   <- nls(formula.LL3, start = list(a=1, b = 3, c = 400), data = df.l[[4]])
  SIG.nls1   <- nls(formula.SIG, start = list(a=1, b = -0.1, c = 400), data = df.l[[4]])
  
  nlstools::overview(SIG.nls1) # lower residual sum of squares
  nlstools::overview(LL3.nls1)

  # SELECT a sigmoid with 3 parameters 
  summary(SIG.nls1)
  
  w=4; graphics.off(); windows(10,10); plot(df.l[[w]]$DD,df.l[[w]]$fSCA)
  lines(1:2000,sigmoid(coef(SIG.nls1),1:2000))
  lines(1:2000,LL3(coef(LL3.nls1),1:2000),col="red")
  
  SIG.boot1  <- nlsBoot(SIG.nls1, niter = 1000)
  
  CIrse <- Hmisc::smedian.hilow(SIG.boot1$rse,conf.int=.95)
  ID    <- which(SIG.boot1$rse>CIrse[2] & SIG.boot1$rse<CIrse[3])
  CIPar <- SIG.boot1$coefboot[ID,]
  dim(CIPar) # 950 out of 1000
  
  SIMUL <- list()
  for (i in 1:nrow(CIPar)){
    print(i)
    # SIMUL[[i]] <- 1/(1 + exp(CIPar[i,1]*(log(1:2000)-log(CIPar[i,2]))))
    SIMUL[[i]] <- sigmoid(CIPar[i,],1:2000)
  }
  SIMUL <- do.call(rbind,SIMUL) 
   
  setwd(WD)
  write.csv(SIMUL,"SIG.SIMUL.csv",row.names=F)
  write.csv(SIG.boot1$bootCI,"SIG.PARAM.csv",row.names=F)
  SIG.PARAM <- coef(SIG.nls1)
  
  MED <- apply(SIMUL, 2, median)
  LOW <- apply(SIMUL, 2, function(x) smedian.hilow(x)[2])
  HIG <- apply(SIMUL, 2, function(x) smedian.hilow(x)[3])
  
  w=4; graphics.off(); windows(10,10); plot(df.l[[w]]$DD,df.l[[w]]$fSCA)
  lines(1:2000,MED,col="blue")
  lines(1:2000,sigmoid(SIG.PARAM,1:2000),col="blue")
  lines(1:2000,LOW,col="lightblue",lwd=2)
  lines(1:2000,HIG,col="lightblue",lwd=2)
  
  range(SIG.boot1$coefboot[ID,"b"])
  graphics.off()
  
  windows(5,5)
  ggplot(SIG.boot1$coefboot[ID,], aes(x=b, y=c)) +
    ggdensity::geom_hdr(show.legend=T,probs = c(0.99, 0.95, 0.75, 0.5,0.25))+
    labs(x="b", y="c")+
    scale_x_continuous(lim=range(SIG.boot1$coefboot[ID,"b"]))+
    scale_y_continuous(lim=range(SIG.boot1$coefboot[ID,"c"]))+
    geom_vline(xintercept = SIG.boot1$bootCI[2,], color = "blue",linetype="dotted")+
    geom_hline(yintercept = SIG.boot1$bootCI[3,], color = "blue",linetype="dotted")

  # This representation is generally close to the representation of the 95 percent Beale’s confidence region provided by nlsConfRegions
  
  plot(nlsResiduals(SIG.nls1))
  nlstools::test.nlsResiduals(nlsResiduals(SIG.nls1)) # non normal
  plot(nlsResiduals(LL3.nls1))
  nlstools::test.nlsResiduals(nlsResiduals(LL3.nls1)) # very non normal
  
  # The null-hypothesis of this test is that the population is normally distributed. If the p value is less than the chosen alpha level, then the null hypothesis is rejected and there is evidence that the data tested are not normally distributed.[4]
  
  contmaf <- nlsContourRSS(O2K.nls1)
  plot(contmaf, col = FALSE, nlev = 10,xlim=c(2,4),ylim=c(300,500))
  plot(contmaf$seqPara[,1],contmaf$seqPara[,2])
  # red is the 95 percent Beale's confidence region
  # how to sample parameters in the confidence region ?
  contmaf$seqPara
  boomaf <- nlsBoot(O2K.nls1, niter = 2000)
  jackmaf <- nlstools::nlsJack(O2K.nls1)
  
    # DEPRECATED Train model with 2013-2023 period  ----
    # implement model for the TOT domain - irrespective of years and then calculate error per year ?
  PARAM    <- matrix(NA,4,3)
  fitmodel <- list()
  
  setwd(WD)  
  
  tmp <- read.csv("THE.fSCA.csv")  # select training data set
  
  tmp$fSCA <- tmp$fSCA/100
  TMP   <- split(tmp,f=tmp$WS)
  
  for (i in 1:length(TMP)) {
    print(i)
    x=TMP[[i]]$DD; y=TMP[[i]]$fSCA/100
    # constraint on parameter a=1
    fitmodel[[i]] <- nls(fSCA~1/(1 + exp(-b * (DD-c))), data=TMP[[i]],start=list(b=-.1,c=400))
    
    # bootstraping for non linear model
    # https://cran.r-project.org/web/packages/nlraa/vignettes/Confidence-Bands.html
    
    boot <- nlraa::boot_nls(fitmodel[[i]],fitted) ## This takes a few seconds
    boot.prd <- summary_simulate(t(boot$t),na.rm=T)
    TMP[[i]] <- data.frame(TMP[[i]], 
                         fit = boot.prd[,1], 
                         lwr = boot.prd[,3],
                         upr = boot.prd[,4])
    PARAM[i,] <- c(a=1,coef(fitmodel[[i]]))
  }
  
  graphics.off()
  windows(4,4)
  ggplot(data = TMP[[4]], aes(x = DD, y = fSCA)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "darkgray", alpha = 0.8)+
    geom_line(aes(y = fit),col="black",linewidth=0.8) + 
    geom_point(size=2) + 
    # geom_hline(yintercept=c(0.5,0.2,0.1))+
    scale_x_continuous(limits=c(0,1600),breaks=seq(0,1500,500),expand=c(0.01,0.01))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0.01,0))+
    labs(x="Degree Days", y="fSCA")+
    theme_bw()+
    theme(legend.position = c(0.5,0.5),legend.direction = "horizontal",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))
  
  # better to add an insert for separate watersheds
  
  graphics.off()
  windows(4,4)
  ggplot(data=TMP[[1]],aes(x = DD, y = fSCA)) + 
    geom_ribbon(data=TMP[[1]],aes(ymin = lwr, ymax = upr), fill = "blue", alpha = 0.2)+
    geom_line(data=TMP[[1]],aes(y = fit),col="blue") + 
    geom_ribbon(data=TMP[[2]],aes(ymin = lwr, ymax = upr), fill = "red", alpha = 0.2)+
    geom_line(data=TMP[[2]],aes(y = fit),col="red") + 
    geom_ribbon(data=TMP[[3]],aes(ymin = lwr, ymax = upr), fill = "green", alpha = 0.2)+
    geom_line(data=TMP[[3]],aes(y = fit),col="darkgreen") + 
    # geom_hline(yintercept=c(0.5,0.2,0.1))+
    scale_x_continuous(limits=c(0,1600),breaks=seq(0,1500,500),expand=c(0.01,0.01))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0.01,0))+
    labs(x="Degree Days", y="fSCA")+
    theme_bw()+
    theme(legend.position = c(0.5,0.5),legend.direction = "horizontal",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))
  
    # DEPRECATED linear quantile regression  ----
  # select training dataset
  setwd(WD)
  # DLR <- read.csv("DLR.fSCA.csv")
  # SWH <- read.csv("SWH.fSCA.csv") # use DLR instead
  THE <- read.csv("THE.fSCA.csv") 
  # IGN <- read.csv("IGN.fSCA.csv")
  # df <- rbind(THE,IGN)
  
  # use only THE for training
  df <- THE
  df$fSCA <- df$fSCA/100
  
  df$WS <- as.factor(df$WS)
  levels(df$WS) <- c("3.LAUZ","1.MAND","2.RNPI","4.TOT")
  df$WS <- as.vector(df$WS)
  
  # remove points in September !
  df <- df[which(df$DOY<=245),]
  # remove a few outliers
  OUT <- which(df$fSCA>0.9 & df$DD>1000)
  if(length(OUT)>0) df <- df[-OUT,]
  
  df.l<- split(df,df$WS)
  dim(df.l[[4]]) # 128 images
  
  w=4; graphics.off(); windows(10,10); plot(df.l[[w]]$DD,df.l[[w]]$fSCA)
  
  # Three-parameter log logistic curves (LL3) was fitted to the data. Rees and Long (1993) successfully demonstrated the utility of this model in the description of seed survival processes in several species.

  DDseq <- 1:2000 # remind log(x) in LL3
  
  fitmodel <- fitmodel01 <- fitmodel05 <- fitmodel09 <- list()
  
  for (i in 1:length(df.l)){
    x=df.l[[i]]$DD
    y=df.l[[i]]$fSCA
    
    PARAMini <- list(b=4,c=300)
    # bound a to 1
    fitmodel[[i]]  <- nls(y~1/(1 + exp(b*(log(x)-log(c)))), start=PARAMini)
    
    # using non linear quantile regression
    fitmodel01[[i]] <- nlrq(y~1/(1 + exp(b*(log(x)-log(c)))), start=PARAMini, tau=0.1, trace=TRUE)
    fitmodel05[[i]] <- nlrq(y~1/(1 + exp(b*(log(x)-log(c)))), start=PARAMini, tau=0.5, trace=TRUE)
    fitmodel09[[i]] <- nlrq(y~1/(1 + exp(b*(log(x)-log(c)))), start=PARAMini, tau=0.9, trace=TRUE)
    
    # fitmodel  <- nls(y~a/(1 + exp(b*(log(x)-log(c)))), start=list(a=1,b=1.6,c=200))
  }
  
  
  df.LL3   <- data.frame(DD = DDseq,
                         MAND = LL3(c(1,summary(fitmodel[[1]])$coeff[,1]),DDseq),
                         RNPI = LL3(c(1,summary(fitmodel[[2]])$coeff[,1]),DDseq),
                         LAUZ = LL3(c(1,summary(fitmodel[[3]])$coeff[,1]),DDseq),
                         TOT =  LL3(c(1,summary(fitmodel[[4]])$coeff[,1]),DDseq)
  )
  setwd(WD);write.csv(df.LL3,"df.LL3.csv",row.names=F)
  
  df01.LL3   <- data.frame(DD = DDseq,
                           MAND = LL3(c(1,summary(fitmodel01[[1]])$coeff[,1]),DDseq),
                           RNPI = LL3(c(1,summary(fitmodel01[[2]])$coeff[,1]),DDseq),
                           LAUZ = LL3(c(1,summary(fitmodel01[[3]])$coeff[,1]),DDseq),
                           TOT =  LL3(c(1,summary(fitmodel01[[4]])$coeff[,1]),DDseq)
  )
  setwd(WD);write.csv(df01.LL3,"df01.LL3.csv",row.names=F)
  df05.LL3   <- data.frame(DD = DDseq,
                           MAND = LL3(c(1,summary(fitmodel05[[1]])$coeff[,1]),DDseq),
                           RNPI = LL3(c(1,summary(fitmodel05[[2]])$coeff[,1]),DDseq),
                           LAUZ = LL3(c(1,summary(fitmodel05[[3]])$coeff[,1]),DDseq),
                           TOT =  LL3(c(1,summary(fitmodel05[[4]])$coeff[,1]),DDseq)
  )
  setwd(WD);write.csv(df05.LL3,"df05.LL3.csv",row.names=F)
  df09.LL3   <- data.frame(DD = DDseq,
                           MAND = LL3(c(1,summary(fitmodel09[[1]])$coeff[,1]),DDseq),
                           RNPI = LL3(c(1,summary(fitmodel09[[2]])$coeff[,1]),DDseq),
                           LAUZ = LL3(c(1,summary(fitmodel09[[3]])$coeff[,1]),DDseq),
                           TOT =  LL3(c(1,summary(fitmodel09[[4]])$coeff[,1]),DDseq)
  )
  setwd(WD);write.csv(df09.LL3,"df09.LL3.csv",row.names=F) 
  
  w=4; graphics.off(); windows(10,10); plot(df.l[[w]]$DD,df.l[[w]]$fSCA)
  lines(df.LL3$DD,df.LL3$TOT)
  lines(df.LL3$DD,df05.LL3$TOT,col="blue")
  lines(df.LL3$DD,df01.LL3$TOT,col="lightblue",lwd=2)
  lines(df.LL3$DD,df09.LL3$TOT,col="lightblue",lwd=2)
  
  
  # Simulate a curve for more snowy BV - requiring more heat to melt
  # plot(DDseq,LL3(c(0.97,3.1,383),Xseq),type="l")
  # lines(DDseq,LL3(c(1,3.1,483),Xseq),type="l",col="blue")
  # lines(DDseq,LL3(c(1,3.1,583),Xseq),type="l",col="magenta")
  
    # DEPRECATED ----
  for (i in 1:length(TMP)) {
    print(i)
    x=TMP[[i]]$DD; y=TMP[[i]]$fSCA/100
    # constraint on parameter a=1
    fitmodel[[i]] <- nls(fSCA~1/(1 + exp(-b * (DD-c))), data=TMP[[i]],start=list(b=-.1,c=400))
    
    # bootstraping for non linear model
    # https://cran.r-project.org/web/packages/nlraa/vignettes/Confidence-Bands.html
    
    boot <- nlraa::boot_nls(fitmodel[[i]],fitted) ## This takes a few seconds
    boot.prd <- summary_simulate(t(boot$t),na.rm=T)
    TMP[[i]] <- data.frame(TMP[[i]], 
                           fit = boot.prd[,1], 
                           lwr = boot.prd[,3],
                           upr = boot.prd[,4])
    PARAM[i,] <- c(a=1,coef(fitmodel[[i]]))
  }
  
    # IV.B.2. Test model with GEE & DLR dataset -> NGS FIG. 2c----
  setwd(WD); 
  DLR.fSCA <- read.csv("DLR.fSCA.csv");table(DLR.fSCA$WS) 
  GEE.fSCA <- read.csv("GEE.fSCA.csv");table(GEE.fSCA$WS)
  tmp      <- rbind(DLR.fSCA,GEE.fSCA)  # select testing set
  nrow(tmp)/4 # 2284 images

  # predict first period
  # tmp  <- read.csv("THE.fSCA.csv")
  TMP  <- split(tmp,f=tmp$WS)
    # remove cloudy days for TOT
    myTOT <- TMP[[4]][which(TMP[[4]]$fNA/100<0.1),]
    # remove DOY after 1st of September (DOY 245)
    myTOT <- myTOT[which(myTOT$DOY<215),]
    dim(myTOT) # reste 457 images
    
  PRED <- sigmoid(as.numeric(PARAM[4,]),myTOT$DD)
  MEAS <- myTOT$fSCA/100
  
  df <- data.frame(MEAS=MEAS,PRED=PRED)
  dim(df)
  SMA <- m.eval(df$MEAS,df$PRED) 
  
  EXT <- which(df$MEAS>=0.9 | df$MEAS<=0.1)
  SMAmEXT <- m.eval(df$MEAS[-EXT],df$PRED[-EXT])
  dim(df[-EXT,]) # 56 images remaining
  
  graphics.off()
  windows(4,4)
  ggplot(data = df, aes(x = MEAS, y = PRED)) + 
    geom_point(size=2,col="gray") + 
    # geom_hline(yintercept=c(0.5,0.2,0.1))+
    scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0.01,0))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0.01,0))+
    geom_abline(aes(intercept = 0, slope = 1),linetype="dashed")+
    geom_abline(aes(intercept = SMA[,'elev.SMA'], slope = SMA[,'slope.SMA']),col="red")+
    geom_abline(aes(intercept = SMAmEXT[,'elev.SMA'], slope = SMAmEXT[,'slope.SMA']),col="blue")+
    labs(x="Measured fSCA", y="Predicted fSCA")+
    theme_bw()+
    theme(legend.position = c(0.5,0.5),legend.direction = "horizontal",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))
  
  
  # IV.C. TO USE ? DD-fSCA model per year/site ----
  tmp   <- read.csv("THE.fSCA.csv")
  tmp.l <- split(tmp,tmp$WS)
  WS    <- names(tmp.l)  
  
  # initialize 
  SCA.mod <- list(); SCA.mod[[1]]<-SCA.mod[[2]]<-SCA.mod[[3]]<-SCA.mod[[4]]<-list()
  
  PARAM.l <- list()
  for (i in 1:4){
    PARAM.l[[i]]  <- data.frame(YEAR=2013:2023,matrix(NA,length(myYEAR),3))
    colnames(PARAM.l[[i]]) <- c("YEAR","a","b","c")
  }
  names(PARAM.l) <- WS  

  for (i in 1:length(tmp.l)){
    df <- tmp.l[[i]]
    for(j in 1:length(myYEAR)){
      print(paste(i,j,sep="-"))
      AV <- which(df$YEAR==2000+myYEAR[j])
      # select number of images to model
      if(length(AV)>=3){
        mydf <- df[AV,]
        x=mydf$DOY; y=mydf$fSCA/100
        # constraint on parameter a=1
        fitmodel <- nls(y~1/(1 + exp(-b * (x-c))), start=list(b=-.1,c=135))
        ## get the coefficients using the coef function 
        PARAM.l[[i]][j,2:4] <- c(a=1,coef(fitmodel)) 
        SCA.mod[[i]][[j]]   <- fitmodel
      }
    }
  }
  
  for (i in 1:length(PARAM.l)){
  PARAM.l[[i]]$DELTA <- NA
  # add the time between 0.2 and 0.8 of fSCA
    for (j in 1:length(myYEAR)) {
      tmp         <- invsig(as.numeric(PARAM.l[[i]][j,2:4]),c(0.1,0.9))
      PARAM.l[[i]]$DELTA[j] <- tmp[1]-tmp[2]
    }
  }
  
  # save PARAM
  for (i in 1:length(PARAM.l)){
    setwd(WD); write.csv(PARAM.l[[i]],paste0("PARAM_",names(PARAM.l)[i],".csv"),row.names=F) 
  }
  
  plot(sigmoid(as.numeric(PARAM.l[[1]][1,2:4]),61:270),type="l")
  SEQ<-seq(0.01,0.99,0.01); plot(SEQ,invsig(as.numeric(PARAM.l[[1]][1,2:4]),SEQ),typ="l")
  
    # graphical checkings ----
  setwd(WD); PARAM <- read.csv("PARAM_TOT.csv")
  myYEAR <- 13:23
  graphics.off(); windows(5,5); par(mgp=c(2,0.7,0),mar=c(3,3,0,0))
  XSEQ <- 61:240
  plot(x=NA,y=NA,xlim=c(XSEQ[1],XSEQ[length(XSEQ)]),ylim=c(0,1),type="n",ylab="Fractional Snow Covered Area",xlab="Day Of the Year")
  for (i in 1:length(myYEAR)) {
    print(i)
    lines(XSEQ,sigmoid(as.numeric(PARAM[i,2:4]),XSEQ),type="l",col=rainbow(n=12)[i],lwd=2,cex.main=2)
    text(165,sigmoid(as.numeric(PARAM[i,2:4]),165),labels=myYEAR[i])
  }
  abline(h=1);abline(h=0)
  legend("right",lwd=3,lty=1,col=rainbow(n=12),legend=myYEAR)
  
  graphics.off(); windows(8,8); par(mgp=c(2,0.7,0),mar=c(3,3,0,0))
  plot(x=NA,y=NA,ylim=c(61,270),xlim=c(0,1),type="n",xlab="Fractional Snow Covered Area",ylab="Day Of the Year")
  SEQ<-seq(0.01,0.99,0.01)
  for (i in 1:length(myYEAR)) {
    print(i)
    lines(SEQ,invsig(as.numeric(PARAM[i,2:4]),SEQ),type="l",col=rainbow(n=12)[i],lwd=2,cex.main=2)
    text(0.2,invsig(as.numeric(PARAM[i,2:4]),0.2),labels=myYEAR[i])
  }
  abline(v=c(0.2,0.8))
  
  # years 2016, 2020, 2021 with low slope and long time between low and high fSCA
  graphics.off(); windows(5,5); par(mgp=c(2,0.7,0),mar=c(3,3,0,0))
  plot(PARAM[,"c"],PARAM[,"DELTA"],pch=21,bg=rainbow(n=12),cex=2,xlab="SMOD at fSCA=0.5",ylab="Time between low and high fSCA",ylim=c(30,90))
  text(PARAM[,"c"],PARAM[,"DELTA"]+0.9,labels=myYEAR)
  abline(h=Hmisc::smedian.hilow(PARAM[,"DELTA"],0.75),lty=c(1,2,2))
  
  # compare inflexion point with SMOD.FLUX / SMOD.NIVO
  setwd(WD); SMOD.FLUX <- read.csv("SMOD.FLUX.csv")
  setwd(WD); SMOD.NIVO <- read.csv("SMOD.NIVO.csv")  
  
  graphics.off()
  windows(8,8)
  XLIM <- YLIM <- c(130,180)
  X= PARAM[,"c"]; Y=SMOD.FLUX[,"SMOD_esti"]
  plot(x=X,y=Y,pch=21,cex=2,bg="green",xlim=XLIM,ylim=YLIM,xlab="fSCA=0.5 from Imagery",ylab="SMOD from FLUXALP/NIVOSE")
  text(x=PARAM[,"c"],y=2+SMOD.FLUX[,"SMOD_esti"],labels=myYEAR)
  abline(0,1)
  abline(LM0 <- lm(Y~X))
  abline(LM <- lm(Y[-4]~X[-4]),col="green")
  
  par(new=TRUE)
  plot(x=PARAM[4:10,"c"],y=SMOD.NIVO[2:8,"SMOD_nivose"],pch=21,cex=2,bg="magenta",xlim=XLIM,ylim=YLIM,xlab="",ylab="")
  text(x=PARAM[4:10,"c"],y=2+SMOD.NIVO[2:8,"SMOD_nivose"],labels=myYEAR[4:10])
  abline(0,1)
  abline(-20,1,lty=2)
  
    # SAVE DATA part I -> SNOW_RN.RData ----
  setwd(WD); save.image("SNOW_RN.RData") 
  
# V. ANOMALY OF GREENNESS (20240719) ----
  # V.A. 5M or 10M resolution topography (20250119) ----
  # import DEM & DEM-derived variables at 10M resolution 
  setwd(DATA_WD)
  DEM10 <- terra::aggregate(myDEM,fact=20,"mean")
  DEM05 <- terra::aggregate(myDEM,fact=10,"mean")
  
  tmp1  <- terra::terrain(myDEM,v="slope",neighbors=8,unit="degrees")
  SLO10 <- terra::aggregate(tmp1,fact=20,"mean") # in degrees
  SLO05 <- terra::aggregate(tmp1,fact=10,"mean") # in degrees
  
  DAH10 <- terra::aggregate(myDAH,fact=20,"mean")
  DAH05 <- terra::aggregate(myDAH,fact=10,"mean")
  # SLO10 <- terra::aggregate(mySLOPE,fact=20,"mean") 
  
  # Flat areas using slope in degrees - limit 30/35° ?
  THR      <- c(0,15)
  m        <- c(-Inf, THR[1], NA,THR[1],THR[2], 1, THR[2], +Inf, NA)
  rclmat   <- matrix(m, ncol=3, byrow=TRUE)
  FLAT10.cl1   <- terra::classify(SLO10,rclmat)
  FLAT05.cl1   <- terra::classify(SLO05,rclmat)

  THR      <- c(15,30)
  m        <- c(-Inf, THR[1], NA,THR[1],THR[2], 1, THR[2], +Inf, NA)
  rclmat   <- matrix(m, ncol=3, byrow=TRUE)
  FLAT10.cl2   <- terra::classify(SLO10,rclmat)
  FLAT05.cl2   <- terra::classify(SLO05,rclmat)
  
  THR      <- c(30,+Inf)
  m        <- c(-Inf, THR[1], NA,THR[1],THR[2], 1, THR[2], +Inf, NA)
  rclmat   <- matrix(m, ncol=3, byrow=TRUE)
  FLAT10.cl3   <- terra::classify(SLO10,rclmat)
  FLAT05.cl3   <- terra::classify(SLO05,rclmat)
  
  # High areas
  THR      <- 2300
  m        <- c(-Inf, THR, NA,  THR, +Inf, 1)
  rclmat   <- matrix(m, ncol=3, byrow=TRUE)
  HIGH10   <- terra::classify(DEM10,rclmat)
  HIGH05   <- terra::classify(DEM05,rclmat)
  terra::plot(HIGH05,legend=F)
  setwd(DATA_WD);writeRaster(HIGH05,"HIGH05.tif",overwrite=T)
  
  # V.B. GREENNESS ANOMALY from 10M res SPOT4-Take5 2015 (20240721) ----
  # NOTE: pros 
  # pros : multi-temporal allowing to capture a yearly maximum value, avoiding phenological effects
  # cons : only 10M resolution product
  
  LF  <- list.files("E:/DATA_ROCHENOIRE/SPOT4TAKE5_2015/NDVI_10m/",full=TRUE)
  tmp <- lapply(LF,function(x) terra::crop(terra::rast(x),my3EXT)) #
  tmp <- do.call(c,tmp)
  NDVImax <- terra::app(tmp,max,na.rm=T) # maximum value of NDVImax
  terra::crs(NDVImax) <- "epsg:2154"
  
  graphics.off();
  windows(10,10)
  terra::plot(NDVImax)
  terra::plot(my3CONT,add=T,border="white")
  
  # Assemble products
  # assemble 1M resolution products
  dem  <- terra::crop(terra::mask(DEM10,my3CONT),my3EXT)
  dah  <- terra::crop(terra::mask(DAH10,my3CONT),my3EXT)
  slo  <- terra::crop(terra::mask(SLO10,my3CONT),my3EXT)
  ndvi <- terra::crop(terra::mask(NDVImax,my3CONT),my3EXT)
  
  range(dem[],na.rm=T)
  dem.CUT <- terra::classify(dem,seq(1850,3150,100))
  dem.com <- 100*(19+as.numeric(dem.CUT)) # value is the median class
  plot(dem.com,legend=F)
  
  range(dah[],na.rm=T)
  dah.CUT <- terra::classify(dah,seq(-0.85,0.85,0.1))
  dah.com <- 1+as.numeric(dah.CUT)  # 16 levels from 1 to 17
  plot(dah.com,legend=F)
  
  mycomb <- dem.com+dah.com
  myval  <- sort(na.omit(unique(mycomb[])[,1]))
  length(myval) # 206 values
  
  # step1. extract cell values of topo and NDVImax at 10M res
  # prepare anomaly of NDVI using data.frame df - 2 minutes of computation !
  df <- data.frame(  cell = terra::extract(dem,my3CONT,cells=T)[,3],
                     elev = terra::extract(dem,my3CONT)[,2],
                     dah  = terra::extract(dah,my3CONT)[,2],
                     ndvi = terra::extract(ndvi,my3CONT)[,2]
  )
  dim(df) # 149 188
  df <- df %>% drop_na()
  dim(df) # 149 188
  
  df$elevCUT <- cut(df$elev,seq(1850,3150,100),labels=seq(1900,3100,100))
  df$dahCUT  <- cut(df$dah,seq(-0.85,0.85,0.1),labels=1:17)
  

  # step2. calculate the median NDVImax per class -> mydf
  # calculate the median NDVImax per class
  df.l <- split(df,list(df$elevCUT,df$dahCUT),drop=TRUE,sep=" ")
  length(df.l) # 205 combinations
  NAM  <- strsplit(names(df.l)," ",fixed=T)
  
  # Reminder : For smedian.hilow, conf.int is the coverage probability the outer quantiles should target. When the default, 0.95, is used, the lower and upper quantiles computed are 0.025 and 0.975.
  
  mydf <- data.frame(
    elev = as.vector(as.numeric(unlist(lapply(NAM,function(x) x[[1]])))),
    dah  = as.vector(as.numeric(unlist(lapply(NAM,function(x) x[[2]])))),
    tot  = unlist(lapply(df.l,function(x) nrow(x))),    
    l01  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.8)[2])), # 1st
    l09  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.8)[3])), # 9nt 
    l02  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.6)[2])), # 2nd
    l08  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.6)[3])), # 8th
    l03  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.4)[2])), # 3rd
    l07  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.4)[3])), # 7th
    l04  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.2)[2])), # 4th
    l06  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.2)[3])), # 6th
    l10  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.9999)[3])), # maximum value
    l05  = unlist(lapply(df.l,function(x) median(x$ndvi,na.rm=T))) # median
  )
  
  dim(mydf) # 205 comb x 13 indicators
  
  mydf$comb=mydf$elev+mydf$dah
  mydf <- mydf[order(mydf$comb),]
  
  mydf[which(mydf$comb==2010),]
  head(mydf)
  all.equal(myval,mydf$comb)
  myval[which(is.na(match(myval,mydf$comb)))] # 2401 present dans myval et pas dans mydf$comb
  
  setwd(DATA_WD);write.csv(mydf,"MYDFspot.WS3.csv",row.names=F)

  # step3. replace values of comb with deciles using subst -> ANO.WS3.dec.tif
  # all.equal(myval,mydf$comb)
  
  DEC <- substr(101:110,2,3)
  ano <- NULL
  for (j in 1:length(DEC)){
    print(j)
    selcol <- which(colnames(mydf)==paste0("l",DEC[j]))
    dec    <- terra::subst(mycomb,from=myval,to=mydf[,selcol])
    add    <- sign(ndvi-dec) # set 1 if NDVI is below the specified decile
    m      <- terra::subst(add,from=c(-1,0,1),to=c(j,NA,NA))
    ano    <- c(ano,m)
  }
  
  ANO <- terra::sds(ano)
  ANOspot.WS3 <- terra::app(ANO,min,na.rm=T) # keep the minimum value !
  setwd(DATA_WD); terra::writeRaster(ANOspot.WS3,"ANOspot.WS3.dec.tif",overwrite=T)
  
  windows(10,10);plot(ANOspot.WS3,legend=F,col=GRAY_pal1(10))

  # replace values of comb with deviation to median -> ANOspot.WS3.tif
  
  selcol <- which(colnames(mydf)=="l05")
  dec    <- terra::subst(mycomb,from=myval,to=mydf[,selcol])
  dev    <- (NDVImax-dec) # difference between ndvi and the specified decile
  
  setwd(DATA_WD); terra::writeRaster(dev,"ANOspot.WS3.tif",overwrite=T)
  
  graphics.off(); windows(10,10);plot(dev,legend=T)
  # note il reste des effets de tranche - mieux vaudrait centrer les classes pour chaque pixel
  
  # V.B.EXT. ALTERNATIVE WAY -> "GREEN.ANOMALY.10M.csv" (20250120) ----
  # to avoid border effects -> need to calculate the anomaly and standardized anomaly for each distribution around the focal value

  # starting from df values - for each cell identify similar conditions and extract correspoding values in images
  df$MEAN <- df$MEDIAN <- df$COUNT <- df$SD <- NA
  for (i in 1:nrow(df)){
    print(i)
    ELEVmin <- round(df$elev[i])-50
    ELEVmax <- round(df$elev[i])+50
    m       <- c(-Inf, ELEVmin, 0,  ELEVmin, ELEVmax, 1, ELEVmax, +Inf, 0)
    rclmat  <- matrix(m, ncol=3, byrow=TRUE)
    dem.rec <- terra::classify(dem,rclmat)

    DAHmin <- round(df$dah[i],2)-0.05
    DAHmax <- round(df$dah[i],2)+0.05
    m       <- c(-Inf, DAHmin, 0,  DAHmin, DAHmax, 1, DAHmax, +Inf, 0)
    rclmat  <- matrix(m, ncol=3, byrow=TRUE)
    dah.rec <- terra::classify(dah,rclmat)

    tmp <- dem.rec*dah.rec
    SEL <- which(tmp[]==1)
    df$COUNT[i] <- length(SEL)
    
    if (length(SEL)>0){
      df$COUNT[i] <- length(SEL)
      NDVImax.SEQ <- terra::extract(NDVImax,SEL)$max
      df$MEAN[i]     <- mean(NDVImax.SEQ)
      df$MEDIAN[i]   <- median(NDVImax.SEQ)
      df$SD[i]       <- sd(NDVImax.SEQ)
    }
  }
  
  setwd(DATA_WD);write.csv(df,"GREEN.ANOMALY.10M.csv")
  
  # simplifier le code en regardant les combinaions identiques - avec arrondi à 5/10m et 0.02 dah ?
  
  NDVImax.MEAN <- NDVImax ; NDVImax.MEAN[df$cell] <- df$MEAN
  NDVImax.MEDIAN <- NDVImax ; NDVImax.MEDIAN[df$cell] <- df$MEDIAN
  
  NDVImax.STANO <- NDVImax ; NDVImax.STANO[df$cell] <- (df$ndvi-df$MEAN)/df$SD
  NDVImax.ANO   <- NDVImax ; NDVImax.ANO[df$cell] <- (df$ndvi-df$MEDIAN)
  
  setwd(DATA_WD)
  terra::writeRaster(NDVImax.ANO,"NDVImax.ANOspot.tif",overwrite=T)
  terra::writeRaster(NDVImax.STANO,"NDVImax.STANOspot.tif",overwrite=T)


  # Checking NDVImax distribution and pixel number based on mydf ----
  setwd(DATA_WD)
  mydf<-read.csv("MYDFspot.WS3.csv")
  graphics.off()
  windows(6,6)
  ggplot(data=mydf,aes(x=dah, y=elev,fill=tot)) +
    geom_raster(show.legend=T)+
    scale_fill_viridis_c(name="Number of pixels",option="D",direction=-1,begin=0,end=1)+
    xlab("DAH")+
    ylab("Elevation")+
    theme_bw()+
    theme(legend.position = "bottom",legend.direction = "horizontal",axis.text=element_text(size=12,family=myfont),axis.title=element_text(size=12,family=myfont),strip.text = element_text(size = 12,family=myfont),legend.text=element_text(size=12,family=myfont),legend.title=element_text(size=12,family=myfont,vjust = 1),legend.key.width = unit(1.5, "cm"),legend.key.height = unit(0.4, "cm"),legend.box.margin = margin(-10,-10,-10,-10))
  
  ggplot(data=mydf,aes(x=dah, y=elev,fill=l05)) +
    geom_raster(show.legend=T)+
    scale_fill_viridis_c(name="Median NDVImax",option="D",direction=-1,begin=0,end=1)+
    xlab("DAH")+
    ylab("Elevation")+
    theme_bw()+
    theme(legend.position = "bottom",legend.direction = "horizontal",axis.text=element_text(size=12,family=myfont),axis.title=element_text(size=12,family=myfont),strip.text = element_text(size = 12,family=myfont),legend.text=element_text(size=12,family=myfont),legend.title=element_text(size=12,family=myfont,vjust = 1),legend.key.width = unit(1.5, "cm"),legend.key.height = unit(0.4, "cm"),legend.box.margin = margin(-10,-10,-10,-10))
  
  # V.C. GREENNES ANOMALY from the BD ortho (20240721) ----
  # NOTE : pros : high resolution product; cons : only one date and issues for histogram equalizing
    # step1. harmonized greenness from BD ortho -> NDVI.WS3.tif ----
  # see rochenoire.r section I.E. for production of IRC mosaics by dpt 05 and 73
  # products at 0M20 (2022à and 0M50 (2009) were resampled at 1M
  # better to use IRC and a NDVI estimate
  # there are issues in histogram of values -> need of matching ! (before computaton of NDVI !)

  setwd("D:/DATA_ROCHENOIRE/MY_BD_ORTHO/IRC_2009_0M50/")
  LF  <- list.files(pattern="WS3",full=T)

    # RNPI+MAND (05)
  OP05 <- terra::crop(terra::rast(LF[[1]]),rbind(RNPI.CONT,MAND.CONT)) # only RNPI+MAND
  
  # gap-fill pixels with NA values
  w <- 3
  OP5 <- terra::focal(OP05, w = w, fun = mean, na.policy = "only", na.rm = T)
  filled <- OP05
  to_fill <- any(is.na(values(filled)))
  while(to_fill) {
    w <- w + 2  
    filled <- focal(filled, w = w, fun = mean, na.policy = "only", na.rm = T) 
    to_fill <- any(is.na(values(filled)))
  }
  OP05 <- filled
  OP05 <- terra::mask(OP05,rbind(RNPI.CONT,MAND.CONT)) # only RNPI+MAND
  
  # OP05 <- terra::mask(OP05,rbind(RNPI.CONT,MAND.CONT)) " do not mask before equalizing
  graphics.off()
  windows(10,10)
  terra::plotRGB(OP05,legend=F)
  terra::plot(rbind(RNPI.CONT,MAND.CONT),add=T,border="yellow",lwd=2)
  
  # compute NDVI
  # Bande 1 : infrarouge; Bande 2 : rouge; Bande 3 : vert
  NIRn <- OP05[[1]]/(OP05[[1]]+OP05[[2]]+OP05[[3]])
  REDn <- OP05[[2]]/(OP05[[1]]+OP05[[2]]+OP05[[3]])
  GREn <- OP05[[3]]/(OP05[[1]]+OP05[[2]]+OP05[[3]])
  
  NDVI05 <- (NIRn-REDn)/(NIRn+REDn)
  hist(NDVI05[],seq(-1,1,0.02))
  
  graphics.off()
  windows(10,10)
  SEQ<- seq(-1,0.15,0.025)
  terra::plot(NDVI05,breaks=SEQ,legend=F,col=GRAY_pal1(length(SEQ)),stretch="lin")
  terra::plot(rbind(RNPI.CONT,MAND.CONT),add=T,border="yellow",lwd=2)
  
  terra::plot(NDVI05,legend=F,breaks=SEQ,col=GRAY_pal1(length(SEQ)),stretch="lin")
  terra::plot(rbind(RNPI.CONT,MAND.CONT),add=T,border="black",lwd=2)

    # LAUZETTE D73
  OP73 <- terra::crop(terra::rast(LF[[2]]),LAUZ.CONT) # only RNPI+MAND
  OP73 <- terra::mask(OP73,LAUZ.CONT) # only RNPI+MAND
  graphics.off()
  windows(10,10)
  terra::plotRGB(OP73,legend=F)
  terra::plot(LAUZ.CONT,add=T,border="yellow",lwd=2)
  
   # remove shadow on the IRC image
  GRAY <- (OP73[[1]]+OP73[[2]]+OP73[[3]])/3
  terra::plot(GRAY,col=GRAY_pal1(100))
  hist(GRAY[])
  THR     <- 60
  m       <- c(-Inf, THR, NA,  THR, +Inf, 1)
  rclmat  <- matrix(m, ncol=3, byrow=TRUE)
  GRAY.sh <- terra::classify(GRAY, rclmat)
  graphics.off(); windows(20,10);par(mfrow=c(1,2))
    terra::plot(GRAY,legend=F)
    terra::plot(GRAY.sh,alpha=0.6,col="red")
    
    OP73 <- OP73*GRAY.sh
    terra::plot(OP73,legend=F)

    # NDVI for IRC
    # Bande 1 : infrarouge; Bande 2 : rouge; Bande 3 : vert
    NIRn <- OP73[[1]]/(OP73[[1]]+OP73[[2]]+OP73[[3]])
    REDn <- OP73[[2]]/(OP73[[1]]+OP73[[2]]+OP73[[3]])
    GREn <- OP73[[3]]/(OP73[[1]]+OP73[[2]]+OP73[[3]])
    
    NDVI73 <- (NIRn-REDn)/(NIRn+REDn)
    NDVI73 <- terra::mask(NDVI73,LAUZ.CONT) 
    hist(NDVI73[],seq(-1,1,0.01))
    graphics.off()
    windows(10,10)
    SEQ<- seq(-0.8,0.05,0.025)
    terra::plot(NDVI73,legend=F,breaks=SEQ,col=GRAY_pal1(length(SEQ)),stretch="lin")
    terra::plot(LAUZ.CONT,add=T,border="red",lwd=2)
    
  # histogram equalizing
    
    # original image
    # histogram matching using RStoolbox - à faire sur la m^me portion d'image !
    
    library(RStoolbox)
    NDVI73.hm <- RStoolbox::histMatch(NDVI73,NDVI05,intersectOnly = FALSE,paired=FALSE) # quite long if multi channels
    
    graphics.off(); windows(15,20);par(mfrow=c(3,1))
    hist(NDVI73[],seq(-1,1,0.01))
    hist(NDVI73.hm[],seq(-1,1,0.01))
    hist(NDVI05[],seq(-1,1,0.01))
    
    graphics.off(); windows(20,15);par(mfrow=c(1,3))
    terra::plot(NDVI73)
    terra::plot(NDVI73.hm*GRAY.sh)
    terra::plot(NDVI05)
    
    # mosaic - beware to the extent
    NDVI.WS3 <- terra::mosaic(NDVI73.hm*GRAY.sh,NDVI05)
    hist(NDVI.WS3[],seq(-1,1,0.01))
    terra::ext(NDVI.WS3)
    NDVI.WS3 <- terra::extend(NDVI.WS3,my3EXT)
    
    hist(NDVI.WS3)
    
    graphics.off();windows(10,10)
    SEQ<- seq(-0.8,0.15,0.025)
    terra::plot(NDVI.WS3,legend=F,breaks=SEQ,col=GRAY_pal1(length(SEQ)),stretch="lin")
    terra::plot(my3CONT,border="blue",add=T,lwd=1)
    
    # NOTE : remaining matching issues in the upper part of Lauzette
    
    setwd(DATA_WD); terra::writeRaster(NDVI.WS3,"NDVI.WS3.tif",overwrite=T)
    
      # OLD STUFF ----
    # img05 <- lidaRtRee::raster2Cimg(OP05[[1]])
    # img05.g <- grayscale(img05)
    # grayscale(img05.g) %>% hist(main="Luminance values in boats picture")
    # 
    # img73 <- lidaRtRee::raster2Cimg(NDVI73)
    # 
    # f <- ecdf(img73)
    # f(img73) %>% hist(main="Transformed luminance values")
    # f(img73) %>% as.cimg(dim=dim(img73)) %>% plot(
    #   main="With histogram equalisation")
    # # Hist. equalisation for grayscale
    # hist.eq <- function(im) as.cimg(ecdf(im)(im),dim=dim(im))
    
    # step2. anomaly of greenness -> ANO.WS3.tif (2025019) ----
    # Classify the DEM x DAH combination
    setwd(DATA_WD); NDVI.WS3 <- terra::rast("NDVI.WS3.tif")
    crs(NDVI.WS3) <-"EPSG:2154"
    
    # assemble 1M resolution products
    dem  <- terra::crop(terra::mask(terra::aggregate(myDEM,fact=2),my3CONT),my3EXT)
    dah  <- terra::crop(terra::mask(terra::aggregate(myDAH,fact=2),my3CONT),my3EXT)
    slo  <- terra::crop(terra::mask(terra::aggregate(mySLOPE,fact=2),my3CONT),my3EXT)
    ndvi <- terra::crop(terra::mask(NDVI.WS3,my3CONT),my3EXT)
    
    range(dem[],na.rm=T)
    dem.CUT <- terra::classify(dem,seq(1850,3150,100))
    plot(dem.CUT,legend=F)
    dem.com <- 100*(19+as.numeric(dem.CUT)) # value is the median class
    plot(dem.com,legend=F)
    
    range(dah[],na.rm=T)
    dah.CUT <- terra::classify(dah,seq(-1,1,0.1))
    dah.com <- 1+as.numeric(dah.CUT)  # 20 levels from 1 to 20
    plot(dah.com,legend=F)
    
    mycomb <- dem.com+dah.com
    myval  <- sort(na.omit(unique(mycomb[])[,1]))
    length(myval)
    
    # prepare anomaly of NDVI using data.frame df - 2 minutes of computation !
    df <- data.frame(  cell = terra::extract(dem,my3CONT,cells=T)[,3],
                       elev = terra::extract(dem,my3CONT)[,2],
                       dah  = terra::extract(dah,my3CONT)[,2],
                       ndvi = terra::extract(ndvi,my3CONT)[,2]
    )
    dim(df) # environ 15 millions de cellules
    
    df$elevCUT <- cut(df$elev,seq(1850,3150,100),labels=seq(1900,3100,100))
    df$dahCUT  <- cut(df$dah,seq(-1,1,0.1),labels=1:20)
    
    df <- df %>% drop_na()
    dim(df)
    
    # calculate the median NDVImax per class
    df.l <- split(df,list(df$elevCUT,df$dahCUT),drop=TRUE,sep=" ")
    length(df.l) # 258 combinations
    NAM  <- strsplit(names(df.l)," ",fixed=T)

    # Reminder : For smedian.hilow, conf.int is the coverage probability the outer quantiles should target. When the default, 0.95, is used, the lower and upper quantiles computed are 0.025 and 0.975.
    
    mydf <- data.frame(
      elev = as.vector(as.numeric(unlist(lapply(NAM,function(x) x[[1]])))),
      dah  = as.vector(as.numeric(unlist(lapply(NAM,function(x) x[[2]])))),
      tot  = unlist(lapply(df.l,function(x) nrow(x))),    
      l01  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.8)[2])), # 1st
      l09  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.8)[3])), # 9nt 
      l02  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.6)[2])), # 2nd
      l08  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.6)[3])), # 8th
      l03  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.4)[2])), # 3rd
      l07  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.4)[3])), # 7th
      l04  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.2)[2])), # 4th
      l06  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.2)[3])), # 6th
      l10  = unlist(lapply(df.l,function(x) Hmisc::smedian.hilow(x$ndvi,.9999)[3])), # maximum value
      l05  = unlist(lapply(df.l,function(x) median(x$ndvi,na.rm=T))) # median
   )
    
    dim(mydf) # 258 comb x 13 indicators
    
    mydf$comb=mydf$elev+mydf$dah
    mydf <- mydf[order(mydf$comb),]
    
    mydf[which(mydf$comb==2010),]
    head(mydf)
    
    # replace values of comb with deciles using subst -> ANO.WS3.dec.tif
    all.equal(myval,mydf$comb)
    
    DEC <- substr(101:110,2,3)
    ano <- NULL
    for (j in 1:length(DEC)){
      print(j)
      selcol <- which(colnames(mydf)==paste0("l",DEC[j]))
      dec    <- terra::subst(mycomb,from=myval,to=mydf[,selcol])
      add    <- sign(ndvi-dec) # set 1 if NDVI is below the specified decile
      m      <- terra::subst(add,from=c(-1,0,1),to=c(j,NA,NA))
      ano    <- c(ano,m)
    }
    
    ANO <- terra::sds(ano)
    ano.WS3 <- terra::app(ANO,min,na.rm=T) # keep the minimum value !
    setwd(DATA_WD); terra::writeRaster(ano.WS3,"ANO.WS3.dec.tif",overwrite=T)
    
    graphics.off(); windows(10,10);plot(ANO.WS3,legend=F)
    
    # replace values of comb with deviation to median -> ANO.WS3.tif
    all.equal(myval,mydf$comb)
    
    ano <- NULL
      selcol <- which(colnames(mydf)=="l05")
      dec    <- terra::subst(mycomb,from=myval,to=mydf[,selcol])
      dev    <- (ndvi-dec) # difference between ndvi and the specified decile

    setwd(DATA_WD); terra::writeRaster(dev,"ANO.WS3.tif",overwrite=T)
    
    graphics.off(); windows(10,10);plot(dev,legend=T)

    
  # OLD : Compute anomalies by elev x DAH
    # setwd("D:/DATA_ROCHENOIRE/LIDAR/LIDAR_IGN/")
    # DEM2M <- terra::rast("DEM_2M.tif")
    # DAH2M <- terra::rast("DAH_2M.tif")
    # SLO2M <- terra::rast("SLO_2M.tif")
    # 
    # df <- data.frame(  cell = terra::extract(DEM2M,my3CONT,cells=T)[,3],
    #                    elev = terra::extract(DEM2M,my3CONT)[,2],
    #                    dah  = terra::extract(DAH2M,my3CONT)[,2],
    #                    slo  = terra::extract(SLO2M,my3CONT)[,2]
    # )
    # 
    # dim(df) # 3,730,151 cells
    
  # V.D. NOT TO RUN - compare with excess of red using the refernec SPOT image ----
  tmp  <- terra::rast("E:/DATA_ROCHENOIRE/SPOT/50512591109131049392J0.TIF")
  tmp1 <- terra::crop(tmp,my3CONT)
  terra::plotRGB(tmp1)

  # compute an excess of red
  # SPOT 5 10 m multispectral:
  # Band 1: green (0.50â€“0.59 µm)
  # Band 2: red (0.61â€“0.68 µm)
  # Band 3: near infrared (0.78â€“0.89 µm)
  # Band 4: middle infrared (MIR) (1.58â€“1.75 µm) at 20 m 
  GREn <- tmp1[[1]]/(tmp1[[1]]+tmp1[[2]]+tmp1[[3]])
  REDn <- tmp1[[2]]/(tmp1[[1]]+tmp1[[2]]+tmp1[[3]])
  graphics.off();windows(10,10);terra::plot(2*GREn-REDn)
  # V.E. NOT TO RUN - Excess of green from a 2M res RGB image ----
    # issue with land use and snow patches and unharmonized reflectance
    # Excess/deficit of green for RVB
    setwd("E:/DATA_ROCHENOIRE/MY_BD_ORTHO/RVB_2m/")
    
    RNPI.OP <- terra::rast("RNED_ORTHO_D05_LA93_2m.tif")
    
    # normalize bands (useful ?)
    REDn <- RNPI.OP[[1]]/(RNPI.OP[[1]]+RNPI.OP[[2]]+RNPI.OP[[3]])
    GREn <- RNPI.OP[[2]]/(RNPI.OP[[1]]+RNPI.OP[[2]]+RNPI.OP[[3]])
    BLUn <- RNPI.OP[[3]]/(RNPI.OP[[1]]+RNPI.OP[[2]]+RNPI.OP[[3]])
    
    EXG <- 2*GREn-BLUn-REDn
    hist(EXG[],seq(-1,2,0.01))
    graphics.off()
    windows(10,10)
    terra::plot(EXG,breaks=seq(-0.075,0.25,0.025),stretch="lin",legend=F)
    terra::plot(EXG,breaks=seq(-0.05,0.075,0.025),col=GRAY_pal1(20),stretch="lin",legend=F)
    terra::plot(rbind(RNPI.CONT,MAND.CONT),add=T,border="red",lwd=2)
    
# VI. COMPARING ANOMALY OF GREENNESS AND SNOW COVER MAP ----
    # Range of deviation along snow cover gradient ----
    setwd(DATA_WD); DEV10 <- terra::aggregate(terra::rast("ANO.WS3.tif"),fact=10)
    FLAT10 <- terra::crop(FLAT10,terra::ext(DEV10))
    HIGH10 <- terra::crop(HIGH10,terra::ext(DEV10))    
    
    graphics.off(); windows(10,10);plot(DEV10m <- DEV10*FLAT10*HIGH10,legend=T,colNA="black")
    
    graphics.off(); windows(10,10);plot(DEV05m <- DEV05*FLAT05*HIGH05,legend=T,colNA="black")
    
    terra::plot(my3HILL,add=T,alpha=0.2,col=grey(0:100 / 100),legend=F)
    terra::plot(my3CONT,add=T,border="blue")
    terra::plot(SAG3lines,add=T,col=gray(0.5))
    terra::plot(OBS.PTS[2],col="blue",add=T)
    terra::plot(myROUT,col="black",add=T,lwd=3)
    terra::sbar(1000, type="bar", below="km", label=c(0,0.5,1), cex=.8)
    terra::north(type=2)
    
    
    # 1. to compare with binary 5M res RapidEye images (20250120) ----
    # using BD ortho anomaly of greennness
    setwd(DATA_WD); DEV05ign <- terra::aggregate(terra::rast("ANO.WS3.tif"),fact=5)  
    windows(10,10); plot(DEV05ign)
    
    # using SPOT4take5 anomaly of greenness
    setwd(DATA_WD); tmp <- terra::rast("NDVImax.STANOspot.tif") 
    setwd(DATA_WD); tmp <- terra::rast("NDVImax.ANOspot.tif") 
    ANO <- terra::disagg(tmp,fact=2) # at 5M
    

    # select snow cvover binary map    
    setwd(DATA_WD); REbin <- terra::rast("REbin_20110625.tif")
    setwd(DATA_WD); REbin <- terra::rast("REbin_20110407.tif") # too early    
    
    setwd(DATA_WD); REbin <- erra::rast("REbin_20100707.tif")
    # setwd(DATA_WD); REbin <- terra::rast("REbin_20110529.tif")
    
    terra::crs(REbin) <-"epsg:2154"
    SNOWbin  <- terra::mask(REbin,my3CONT)
    
    ANOc <- terra::crop(ANO,terra::ext(SNOWbin))
    ANOc1 <- ANOc*HIGH05*terra::crop(FLAT05.cl1,terra::ext(SNOWbin))
    ANOc2 <- ANOc*HIGH05*terra::crop(FLAT05.cl2,terra::ext(SNOWbin))
    ANOc3 <- ANOc*HIGH05*terra::crop(FLAT05.cl3,terra::ext(SNOWbin))
    
    graphics.off()
    windows(3,6)
    boxplot(ANOc2[]~SNOWbin[],outline=F)
    abline(h=0)  
    # to be done by slope classes
    
    # do not use this / split anomaly per class of slope
    # DEV05  <- terra::crop(DEV05ign, terra::ext(REbin))
    # FLAT05 <- terra::crop(FLAT05,terra::ext(REbin))
    # HIGH05 <- terra::crop(HIGH05,terra::ext(REbin))
    # DEV05m <- DEV05*FLAT05*HIGH05    
    # REY05m <- terra::mask(REbin*FLAT05*HIGH05,my3CONT)
    
    # 2. to compare with SMOD 2015 at 10M ----
    setwd(DATA_WD); SMOD10 <- terra::project(terra::rast("SCA13.SMODs.tif"),DEV10)
    # compute deciles and reckassify
    SEQ <- seq(0.05,0.95,0.05)
    DEC <- as.numeric(global(SMOD10m <- SMOD10*FLAT10*HIGH10,fun=quantile,probs=SEQ,na.rm=T))
    length(DEC)
    
    df <- as.matrix(cbind(A=c(-Inf,DEC),
                B=c(DEC,+Inf),
                C=1:(1+length(DEC))))

    TMP    <- terra::classify(SMOD10m,df) 
    
    # for figure 1
    boxplot(DEV10m[]~TMP[])
    abline(h=0)
      
    
    length(SMOD10Rm[])  # 280665
    myDEC <- Hmisc::smedian.hilow(x$ndvi,.4)
    
    
    RAN <- range(SMOD10R[],na.rm=T)
    hist(SMOD10R[],breaks=seq(RAN[1],RAN[2],5))
    
    SMOD10R[SMOD10R>220] <- 220
    SMOD10R[SMOD10R<100] <- 100
    hist(SMOD10R[],breaks=seq(95,240,5),xlim=c(95,240))
    
    graphics.off(); windows(10,10);plot(SMOD10Rm<-SMOD10R*FLAT10*HIGH10,legend=T,colNA="black")
    terra::plot(my3HILL,add=T,alpha=0.2,col=grey(0:100 / 100),legend=F)
    terra::plot(my3CONT,add=T,border="blue")
    terra::plot(SAG3lines,add=T,col=gray(0.5))
    terra::plot(OBS.PTS[2],col="blue",add=T)
    terra::plot(myROUT,col="black",add=T,lwd=3)
    terra::sbar(1000, type="bar", below="km", label=c(0,0.5,1), cex=.8)
    terra::north(type=2)
    

    
    # see the post https://jakubnowosad.com/posts/2024-10-13-spatcomp-bp1/
    # and also https://isprs-archives.copernicus.org/articles/XLVIII-4-W12-2024/127/2024/
    

    # RASTER A - greenness anaomaly a 10M   
    Atmp <- ANO02.10M*FLAT.10M*HIGH.10M
    # smooth A
    A <- terra::focal(Atmp,w=3,fun="modal")
    
    terra::plot(A,col=c("white","red"),legend=F)
    terra::plot(my3HILL,add=T,alpha=0.4,col=grey(0:100 / 100),legend=F)
    terra::plot(my3CONT,add=T,border="blue")
    terra::plot(SAG3lines,add=T,col=gray(0.5))
    terra::plot(OBS.PTS[2],col="blue",add=T)
    terra::plot(myROUT,col="black",add=T,lwd=3)
    terra::sbar(1000, type="bar", below="km", label=c(0,0.5,1), cex=.8)
    terra::north(type=2)
    
    
    # RASTER B - proba de neige a10M
    PRSNO02.10M <- terra::project(PRSNO.02,ANO02.10M)
    B <- PRSNO02.10M*FLAT.10M*HIGH.10M
    terra::plot(B,col=c("white","red"),legend=F)
    terra::plot(my3HILL,add=T,alpha=0.4,col=grey(0:100 / 100),legend=F)
    terra::plot(my3CONT,add=T,border="blue")
    terra::plot(terra::vect("SAG3lines.shp"),add=T,col=gray(0.5))
    terra::plot(OBS.PTS[2],col="blue",add=T)
    terra::plot(myROUT,col="black",add=T,lwd=3)
    terra::sbar(1000, type="bar", below="km", label=c(0,0.5,1), cex=.8)
    terra::north(type=2)
    
    # compare A & B
    # the sabre package (Nowosad and Stepinski, 2018) returns two maps and three values: V-measure, homogeneity, and completeness
    library(sabre)
    mapcurves_calc(A, B)
    vmeasure_calc(A, B)
    
    # The homogeneity measures how well regions from the first map fit inside of regions from the second map, the completeness measures how well regions from the second map fit inside of regions from the first map, and the V-measure is the weighted harmonic mean of homogeneity and completeness. All of them range from 0 to 1, where larger values indicate better spatial agreement.
    
    # diffeR (Pontius Jr. and Santacruz, 2023) offers a set of metrics to compare continuous and categorical rasters at original and multiple resolutions, such as the overall difference and mean absolute deviation (MAD).
    library(diffeR)
    mad <- MAD(A, B, strata = NULL, eval = "multiple")
    MADscatterplot(A, B, strata = NULL)
     
    # landscapemetrics package (Hesselbarth et al., 2019) offers more than 130 metrics to describe the composition and configuration of spatial patterns
    library(landscapemetrics)
    
    landscape <- terra::rast(landscapemetrics::landscape)
    
    show_patches(A)
    
    show_lsm(A, what = "lsm_p_area", directions = 8)
    show_lsm(B, what = "lsm_p_area", directions = 4)
    show_lsm(A, what = "lsm_p_shape", class = c(1, 2), label_lsm = TRUE)
    show_lsm(A, what = "lsm_p_circle", class = 3, labels = TRUE)
    
    # The categorical rasters can be compared using spatialEco (Evans and Murphy, 2023) that provides a function raster.change() with various statistics calculated in a moving window, including two-tailed paired t-test and cross-entropy loss function
    library(spatialEco)
    AB <- raster.change(A, B, s = 3,stat="cross-entropy")
    # the larger the value, the more different the two rasters are./ quite long
    graphics.off();windows(10,10);terra::plot(AB)
    
# VI. GEOREFERENCEMENT DES PHOTOS AERIENNES ----
  # https://www.univ-st-etienne.fr/wikimastersig/doku.php/fonctions:integration:georeferencement:images:arcgis
# VII. MAIN FIGURES (20250120) ----
  # NGS FIGURE 1. ANOMALIES OF GREENNESS and SMOD / Landscape scale ----  
    # NGS FIG. 1a Study Area with 10m res. SPOT image (20250314) ----
  setwd(DATA_WD) 
  tmp <- terra::rast("50512591109131049392J0.TIF")
  tmp1 <- terra::crop(tmp,my3CONT)
  
  largeCONT <- terra::ext(964500-800, 969460+2000, 6443020-500, 6448690+500)
  tmp1 <- terra::crop(tmp,largeCONT)
  
  # Create graticule object with sf
  grat <- sf::st_graticule(lon=seq(5.025,7.025, 0.05),
                           lat = seq(44.025,47.025, 0.05),margin=0.2) %>%
    vect() %>%
    project("EPSG:2154") 
  
  graphics.off();windows(5,4.8)
  terra::plotRGB(tmp1,alpha=0.9)
  plot(my3CONT,add=T,border="black",lwd=1)
  points(OBS.PTS[2],bg="black",pch=22,cex=1)
  terra::plot(myROUT,col=gray(0.6),add=T,lwd=3)
  terra::plot(grat,add=TRUE,col="yellow",lty=2,lwd=0.5)
  terra::north(xy=c(963900,6448100),type=2,col="black")

  terra::plot(myHYDRO,add=T,col="skyblue3")
  terra::sbar(xy=c(963900,6442800),d=2000, type="bar", below="km", label=c(0,0.5,1), cex=1.2,lwd=2,col="black")  
  
    # NGS FIG. 1b Standardized anomaly based on SPOT4Take5 (20250314) ----
  NDVImax.STANO <- terra::mask(terra::rast("NDVImax.STANOspot.tif"),my3CONT)
  BORNE <- terra::global(NDVImax.STANO,min,na.rm=T)
  BORNE <- terra::global(NDVImax.STANO,quantile,probs=c(0.05,0.95),na.rm=T)
  test2 <- terra::stretch(NDVImax.STANO, minq=0.05, maxq=0.95)
  test2 <- terra::stretch(NDVImax.STANO,histeq=T)
  
  graphics.off();windows(4.3,4.8)
  terra::plot(test2,col=DIVcbs_VGpal(20),axes=T,clip=T,
              mar=c(0,0,0,0),
              # pax=list(side=1:2),
              plg=list( # parameters for drawing legend
                title = "", 
                title.adj=0, # cannot left justify !
                xjust=0,
                loc="bottom",
                horiz=TRUE,
                ext = c(964600, 967000, 6443400, 6443600),
                title.cex = 1, # Legend title size
                cex = 1 # Legend text size
              ))
  terra::plot(my3HILL,add=T,alpha=0.2,col=grey(0:100 / 100),legend=F)
  terra::plot(my3CONT,add=T,border="black")
  terra::plot(SAG3lines,add=T,col=gray(0.5))
  points(OBS.PTS[2],bg="black",pch=22,cex=1)
  terra::plot(myROUT,col="gray",add=T,lwd=3)
  # terra::sbar(1000, type="bar", below="km", label=c(0,0.5,1), cex=.8,col="white")
  # terra::north(type=2)
  
    # NGS FIG. 1c Example of location of late lying snowfields (20250314) ----
    setwd(DATA_WD); REbin <- terra::rast("REbin_20100707.tif")
    terra::crs(REbin) <-"epsg:2154"
    REY05m <- terra::mask(REbin,my3CONT)
    
    graphics.off();windows(4.3,4.8)
    terra::plot(REY05m,axes=T,clip=T,col=c(NOSNOWcol,SNOWcol),
                mar=c(0,0,0,0),
                plg=list( # parameters for drawing legend
                  title = "", 
                  title.adj=0, # cannot left justify !
                  horiz=TRUE,
                  # ext = c(964600, 967000, 6443400, 6443600), # does not work for discrete
                  title.cex = 1, # Legend title size
                  cex = 1 # Legend text size
                ))
    terra::add_legend(x=964600,y=6444000,legend=c("snow free","snow covered"),pt.bg=c(NOSNOWcol,SNOWcol),pch=22,bty="n",y.intersp=0.6,pt.cex=2)
    terra::plot(my3HILL,add=T,alpha=0.2,col=grey(0:100 / 100))
    terra::plot(my3CONT,add=T,border="black")
    terra::plot(SAG3lines,add=T,col=gray(0.5))
    points(OBS.PTS[2],bg="black",pch=22,cex=2)
    terra::plot(myROUT,col=gray(0.6),add=T,lwd=3)
    # terra::sbar(1000, type="bar", below="km", label=c(0,0.5,1), cex=.8,col="black")
    # terra::north(type=2)
  
    # NGS FIG. 1d Distribution of anomaly per snow cover (20250314) ----
    # using SPOT4take5 anomaly of greenness
    # nee to run part V.A.
    setwd(DATA_WD) 
    tmp    <- terra::rast("NDVImax.STANOspot.tif") 
    ANO    <- terra::disagg(tmp,fact=2) # disaggregate at 5M
    MYEXT <- ext(ANO)
    
    # select snow cover binary map of Fig. 1c    
    REbin <- terra::rast("REbin_20100707.tif")
    # setwd(DATA_WD); REbin <- terra::rast("REbin_20110625.tif")
    terra::crs(REbin) <-"epsg:2154"
    REbin <- terra::extend(terra::project(terra::mask(REbin,my3CONT),ANO),MYEXT)
    
    HIGH05 <- terra::crop(terra::rast("HIGH05.tif"),MYEXT)
    SNOWbin  <- HIGH05*REbin
    
    ANOc  <- terra::crop(ANO,terra::ext(SNOWbin))
    ANOc1 <- ANOc*HIGH05*terra::crop(FLAT05.cl1,terra::ext(SNOWbin))
    ANOc2 <- ANOc*HIGH05*terra::crop(FLAT05.cl2,terra::ext(SNOWbin))
    ANOc3 <- ANOc*HIGH05*terra::crop(FLAT05.cl3,terra::ext(SNOWbin))
    
    tmp1 <- na.omit(terra::extract(ANOc1,which(SNOWbin[]==1)))
    df1 <- data.frame(ANO=tmp1,SNOW="3.snow")
    tmp2 <- na.omit(terra::extract(ANOc1,which(SNOWbin[]==0)))
    df2 <- data.frame(ANO=tmp2,SNOW="2.nosnow")
    tmp3 <- na.omit(terra::extract(ANOc1,which(SNOWbin[]>=0)))
    df3 <- data.frame(ANO=tmp3,SNOW="1.all")
    
    df <- rbind(df1,df2,df3)
    colnames(df) <- c("ANO","SNOW")
    df$SLOP <- "1.LOW"
    df$SNOW <-as.factor(df$SNOW)
    
    tmp1 <- na.omit(terra::extract(ANOc2,which(SNOWbin[]==1)))
    df1 <- data.frame(ANO=tmp1,SNOW="3.snow")
    tmp2 <- na.omit(terra::extract(ANOc2,which(SNOWbin[]==0)))
    df2 <- data.frame(ANO=tmp2,SNOW="2.nosnow")
    tmp3 <- na.omit(terra::extract(ANOc2,which(SNOWbin[]>=0)))
    df3 <- data.frame(ANO=tmp3,SNOW="1.all")
    
    dm <- rbind(df1,df2,df3)
    colnames(dm) <- c("ANO","SNOW")
    dm$SLOP <- "2.MID"
    dm$SNOW <-as.factor(dm$SNOW)
    
    tmp1 <- na.omit(terra::extract(ANOc3,which(SNOWbin[]==1)))
    df1 <- data.frame(ANO=tmp1,SNOW="3.snow")
    tmp2 <- na.omit(terra::extract(ANOc3,which(SNOWbin[]==0)))
    df2 <- data.frame(ANO=tmp2,SNOW="2.nosnow")
    tmp3 <- na.omit(terra::extract(ANOc3,which(SNOWbin[]>=0)))
    df3 <- data.frame(ANO=tmp3,SNOW="1.all")
    
    dg <- rbind(df1,df2,df3)
    colnames(dg) <- c("ANO","SNOW")
    dg$SLOP <- "3.HIG"
    dg$SNOW <-as.factor(dg$SNOW)
    
    DF <- rbind(df,dm,dg)
    
    # add a test
    head(DF)
    summary(lm(ANO~SNOW,data=DF[which(DF$SLOP=="1.LOW"),]))
    summary(lm(ANO~SNOW,data=DF[which(DF$SLOP=="2.MID"),]))
    summary(lm(ANO~SNOW,data=DF[which(DF$SLOP=="3.HIG"),]))
    
    COLOR <- c("gray",NOSNOWcol,SNOWcol)
    graphics.off()
    windows(6,5)
    ggplot(data=DF, aes(x=SNOW, y=ANO, fill=SNOW)) + 
      geom_violin(alpha=0.4,show.legend=F)+
      # geom_point(size=0.5,position = position_jitterdodge(seed = 1, dodge.width = 0.9),show.legend=F)+
      ylab("Standardized Anomaly of Greenness") + xlab("")+
      scale_x_discrete(labels=c("all","no\nsnow", "snow"))+
      # geom_boxplot(width=0.3,notch=T)+
      stat_summary(fun.data = f, geom="boxplot", 
                   position=position_dodge(1),
                   width=0.2,linewidth=0.7,show.legend=F)+
      geom_hline(yintercept=c(0),linetype="solid",color="black")+
      # stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.2 )+
      # stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red")+
      scale_fill_manual(values=COLOR)+
      facet_wrap2(vars(SLOP),
                  labeller=as_labeller(c("1.LOW"="< 15°","2.MID"="[15°-30°]","3.HIG"="> 30°")))+
      scale_y_continuous(expand=c(0.02,0.01),limits=c(-4,6),breaks=seq(-4,6,2))+
      theme_bw()+
      theme(legend.position = "bottom",legend.direction = "horizontal",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))
    
  # NGS FIGURE 2. ANOMALIES OF GREENNESS, SMOD & SOM / Plot scale ----
    # NGS FIG. 2a Image transect EN->CF ----
    # NGS FIG. 2b Anomaly of greenness for BIOCLIM plots (20250314) ----
   setwd(DATA_WD)
   ALL.PTS        <- read.csv("ALL.PTS.csv") 
   NDVImax.STANO <- terra::mask(terra::rast("NDVImax.STANOspot.tif"),my3CONT)
   ALL.PTS$ANO   <- terra::extract(NDVImax.STANO,ALL.PTS[,c("X_L93","Y_L93")])[,2]

   tmp <- as.factor(as.vector(ALL.PTS$COMcor))
   levels(tmp) <-c("2.CF","4.CT","1.EN","7.FP","6.KS","3.PA","5.VM")
   ALL.PTS$COMcor <- as.vector(tmp)
   
   # issue with one CF / RON_CF_148 located in a very small area
   ALL.PTS[which(ALL.PTS$COMcor=="2.CF"),c("ID","ANO")]
   
   # ggplot 
   graphics.off(); windows(3,4)
   ggplot(data=ALL.PTS, aes(x=COMcor, y=ANO, fill=COMcor)) + 
     geom_violin(alpha=0.4,show.legend=F)+
     geom_point(position = position_jitter(seed = 1, width = 0.2),shape=16,size = 1.5,show.legend=F)+
     ylab("Anomaly of greenness") + xlab("Habitat")+
     scale_x_discrete(labels=rev(c("FP", "KS", "VM", "CT","PA","CF","EN")))+
     # geom_boxplot(width=0.3,notch=T)+
     stat_summary(fun.data = f, geom="boxplot", 
                  position=position_dodge(1),
                  width=0.25,linewidth=0.7,show.legend=F)+
     geom_point(size=1,position = position_jitterdodge(seed = 1, dodge.width = 0.9),show.legend=F)+
     geom_hline(yintercept=c(0),linetype="solid",color="black")+
     # stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.2 )+
     # stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red")+
     scale_fill_manual(values=HABcol)+
     scale_y_continuous(expand=c(0.02,0.01),limits=c(-0.9,2.1),breaks=seq(-1,2.5,0.5))+
     theme_bw()+
     theme(legend.position = "bottom",legend.direction = "horizontal",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))
   
    # NGS FIG. 2c SMOD & SOD for ALL plots (20250314) ----
      # Prepare data
   setwd(DATA_WD); ALL.PTS <- read.csv("ALL.PTS.csv")
   
   tmp <- as.factor(as.vector(ALL.PTS$COMcor))
   levels(tmp) <-c("2.CF","4.CT","1.EN","7.FP","6.KS","3.PA","5.VM")
   ALL.PTS$COMcor <- as.vector(tmp)
   
   ALL.SMOD <- ALL.PTS[,c(grep("COMcor",colnames(ALL.PTS)),grep("SMOD",colnames(ALL.PTS)))]
   ALL.SMOD <- ALL.SMOD[,-grep("SMODmean",colnames(ALL.SMOD))]# delete the mean
   tmp <- reshape2::melt(ALL.SMOD,id.vars="COMcor",variable.name="YEAR",value.name="SMOD")
   df.smod <- tmp[!is.na(tmp$SMOD),]
   
   ALL.SOD <- ALL.PTS[,c(grep("COMcor",colnames(ALL.PTS)),grep("SOD",colnames(ALL.PTS)))]
   ALL.SOD <- ALL.SOD[,-grep("SODmean",colnames(ALL.SOD))]# delete the mean
   tmp <- reshape2::melt(ALL.SOD,id.vars="COMcor",variable.name="YEAR",value.name="SOD")
   df.sod <- tmp[!is.na(tmp$SOD),]
   
      # ggplot 
   graphics.off();windows(3,4)
   ggplot() + 
    # ggplot(data=df.smod, aes(x=COMcor, y=SMOD, fill=COMcor)) + 
     geom_violin(data=df.smod, aes(x=COMcor, y=SMOD, fill=COMcor),alpha=0.4,show.legend=F)+
     geom_point(data=df.smod, aes(x=COMcor, y=SMOD, fill=COMcor),position = position_jitter(seed = 1, width = 0.2),shape=18,size = .5,show.legend=F)+
     ylab("SMOD & SOD (DOY)") + xlab("Habitat")+
     scale_x_discrete(labels=rev(c("FP", "KS", "VM", "CT","PA","CF","EN")))+
     # geom_boxplot(width=0.3,notch=T)+
     stat_summary(data=df.smod, aes(x=COMcor, y=SMOD, fill=COMcor),
                  fun.data = f, geom="boxplot", 
                  position=position_dodge(1),
                  width=0.15,linewidth=0.7,show.legend=F)+

     geom_violin(data=df.sod, aes(x=COMcor, y=SOD, fill=COMcor),alpha=0.4,show.legend=F)+
     geom_point(data=df.sod, aes(x=COMcor, y=SOD, fill=COMcor),position = position_jitter(seed = 1, width = 0.2),shape=18,size = .5,show.legend=F)+
     # geom_boxplot(width=0.3,notch=T)+
     stat_summary(data=df.sod, aes(x=COMcor, y=SOD, fill=COMcor),
                  fun.data = f, geom="boxplot", 
                  position=position_dodge(1),
                  width=0.15,linewidth=0.7,show.legend=F)+
     # geom_hline(yintercept=c(600,1000),linetype="dashed",color="black")+
     # add 227 = 15-08
     geom_hline(yintercept=c(227),linetype="solid",color="black")+
     scale_fill_manual(values=HABcol)+
     scale_y_continuous(expand=c(0.02,0.01),limits=c(50,350),breaks=seq(50,3650,25))+
     theme_bw()+
     theme(legend.position = "bottom",legend.direction = "horizontal",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))

   
    # NGS FIG. 2d biplot of %C x CN ratio (TO UPDATE) ----
   
   setwd(DATA_WD); ALL.PTS <- read.csv("ALL.PTS.csv")
   tmp <- ALL.PTS[,c("COMcor","X.C","CNratio")]
   tmp$YUKI <- ifelse(tmp$COMcor=="EN","EN","non EN")
   df <- tmp %>% drop_na()
   
   tmp0 <- as.factor(as.vector(df$COMcor))
   levels(tmp0) <-c("2.CF","4.CT","1.EN","7.FP","6.KS","3.PA","5.VM")
   df$COMcor <- as.vector(tmp0)
   
   dim(df) # 38 soils
   
   df1 <- df[,1:3]
   df2 <- df[,2:4]
   
   graphics.off()
   windows(5,5)
   ggplot()+
     geom_point(data=df2,aes(x=X.C,y=CNratio),size = 2,show.legend=FALSE)+
     stat_ellipse(data=df2,aes(x=X.C,y=CNratio,fill=YUKI),color="gray",
                  show.legend = FALSE,
                  geom="polygon", 
                  alpha = 0.1,
                  level = 0.9)+
     scale_fill_manual(values=c("gray","gray"),name="EN vs. non EN")+

     guides(fill = "none" )+
     
     geom_point(data=df1,aes(x=X.C,y=CNratio,color=COMcor),size = 3,show.legend=T)+
   
     scale_color_discrete(name="Habitat")+
   ylab("C:N ratio") + xlab("% C")+
     guides(color = guide_legend(override.aes = list(size = 3) ) )+
     scale_alpha(guide = 'none')+
     theme_bw()+
     theme(legend.position = c(0.9,0.3),legend.direction = "vertical",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))
   
   
  # NGS FIG. 3 YUKIGATA  ----
    # NGS FIG. 3a-b Yukigata pictures (20250314) ----
    # NGS FIG. 3c-e Probability maps of snow covered at various fSCA (20250314) ----
    # focus on LSM and show fSCA05,fSCA02 and fSCA01
  setwd(DATA_WD)
  TMP <- terra::rast("SCA.f01.tif") # pick a fSCA level. see part II.D.1
  TMP <- terra::rast("SCA.f02.tif") # pick a fSCA level. see part II.D.1
  TMP <- terra::rast("SCA.f05.tif") # select the fSCA
  
  for (i in 1:11) TMP[[i]]<- terra::mask(terra::subst(TMP[[i]],NA,0),my3CONT)
  PRSNO    <-  terra::app(TMP,mean)
  PRSNOcat <-  PRSNO; PRSNOcat[PRSNO>=0.9]<-1;PRSNOcat[PRSNO<0.9]<-2; PRSNOcat[PRSNO<=0.1]<-3
  TEST <- terra::project(PRSNOcat,PRSNO)
 
 # main plot is a map
 myCOL <- rev(c('#fee090','#abd9e9','#4575b4'))
 graphics.off();windows(4.3,4.8)
 terra::plot(PRSNOcat,add=F,alpha=0.7,col=myCOL,legend=F,axes=F,mar=c(0,0,0,0))
 terra::plot(my3HILL,add=T,alpha=0.4,col=grey(0:100 / 100),legend=F)
 terra::plot(my3CONT,add=T,border="black",lwd=1)
 points(OBS.PTS[2],bg="black",pch=22,cex=1)
 terra::plot(myROUT,col=gray(0.6),add=T,lwd=3)
 terra::plot(SAG3lines,add=T,col=gray(0.5))
 #terra::sbar(1000, type="bar", below="km", label=c(0,0.5,1), cex=.8)
 #terra::north(type=2)
 
 BP <- as.matrix(table(values(PRSNOcat))/sum(table(values(PRSNOcat))))
 graphics.off();windows(1,2.5);par(mar=c(0.5,0.5,0.5,0.5),mgp=c(1,0.5,0));barplot(BP,horiz=F,beside=F,col=alpha(myCOL,0.7),las=2,yaxt="n")
 myY <- c(0.5*BP[1,1],BP[1,1]+0.5*BP[2,1],(BP[2,1]+BP[1,1])+0.5*(1-(BP[2,1]+BP[1,1])))
 for(i in 2:3) text(0.7,myY[i],round(BP[i,],2),adj=0.5,cex=1.5)
 
 # # add a ggpot pie chart  
 # df <- data.frame(
 #   group  = c("no snow","int.","snow"),
 #   values = c(length(which(PRSNO[]<=0.1)),
 #              length(which(PRSNO[]>0.1 & PRSNO[]<0.9)),
 #              length(which(PRSNO[]>=0.9))))
 # df <- df %>% 
 #   arrange(desc(group)) %>%
 #   mutate(prop = values / sum(df$values) *100) %>%
 #   mutate(ypos = cumsum(prop)- 0.5*prop )
 # 
 # graphics.off();windows(3,3)
 # ggplot(data=df, aes(x="", y=prop, fill=group)) +
 #   geom_bar(stat="identity", width=1,show.legend=F,alpha=0.6) +
 #   coord_polar("y", start=0)+
 #   scale_fill_manual(values=SMOD_pal(3)[c(2,1,3)])+
 #   # geom_text(aes(y = ypos, label = round(prop,1)), color = "black", size=6) +
 #   theme_void()+
 #   theme(plot.margin=margin(0,0,0,0))
 
  # NGS FIG. 4. THERMAL MODEL of fSCA (20250315) ----
 setwd(WD) 
 DLR.fSCA <- read.csv("DLR.fSCA.csv");table(DLR.fSCA$WS) 
 GEE.fSCA <- read.csv("GEE.fSCA.csv");table(GEE.fSCA$WS)
 tmp      <- rbind(DLR.fSCA,GEE.fSCA)  # select testing set
 nrow(tmp)/4 # 2284 images
 
 # tmp  <- read.csv("THE.fSCA.csv")
 TMP  <- split(tmp,f=tmp$WS)
 # remove cloudy days for TOT
 myTOT <- TMP[[4]][which(TMP[[4]]$fNA/100<0.1),]
 # remove DOY after 1st of September (DOY 245)
 myTOT <- myTOT[which(myTOT$DOY<215),]
 dim(myTOT) # reste 291+80 = 371 images
 
 myTOT$MEAS <- myTOT$fSCA/100
 plot(myTOT$DD,myTOT$MEAS,pch=".",cex=3)
 
 # import PARAMETERS & PREDICT
 setwd(WD); SIG.PARAM <- read.csv("SIG.PARAM.csv")
 myTOT$PRED <- sigmoid(SIG.PARAM[,1],myTOT$DD)  
 
  THE.fSCA     <- read.csv("THE.fSCA.csv")
  THE_TOT      <- THE.fSCA[which(THE.fSCA$WS=="TOT"),]
  THE_TOT$PRED <- sigmoid(SIG.PARAM[,1],THE_TOT$DD)
  THE_TOT$MEAS <- THE_TOT$fSCA/100
  
  # assemble THEIA and old (DLR,GEE) period
  dfl <- rbind(myTOT,THE_TOT)
  plot(dfl$DD,dfl$MEAS,pch=".",cex=3)
  
  setwd(WD); BIO.DD <- read.csv("BIOCLIM.DD.csv")
  graphics.off();windows(10,10)
  boxplot(DD~COM,data=BIO.DD)
  abline(h=c(600,1000))
  tmp    <- aggregate(BIO.DD$DD,list(BIO.DD$COM),quantile,probs=c(0.1,0.5,0.9))
  
  df.seg <- data.frame(COM=tmp$Group.1,
                       x1=tmp$x[,1],
                       x2=tmp$x[,3],
                       y1=sigmoid(SIG.PARAM[,1],tmp$x[,2]),
                       y2=sigmoid(SIG.PARAM[,1],tmp$x[,2])
  )
  
  # prepare ribbon of model
  setwd(WD); SIMUL <- read.csv("SIG.SIMUL.csv")
  MED <- apply(SIMUL, 2, median)
  LOW <- apply(SIMUL, 2, function(x) smedian.hilow(x)[2])
  HIG <- apply(SIMUL, 2, function(x) smedian.hilow(x)[3])
  BOOT <- data.frame(DD=1:2000,lwr=LOW,upr=HIG,fSCA=MED)
  
  
  graphics.off()
  windows(6,6)
  ggplot(data = BOOT, aes(x = DD, y = fSCA)) + 
    geom_ribbon(data=BOOT, aes(ymin = lwr, ymax = upr), fill = "gray", alpha = 0.8)+
    geom_line() + 
    geom_point(data=dfl, aes(x=DD,y=MEAS),size=1,col="gray") +
    geom_segment(data=df.seg,aes(x=x1,xend=x2,y=y1,yend=y2,col=COM),alpha=0.9,show.legend=T,linewidth=3)+
    # scale_color_manual(HABcol)+
    scale_color_discrete(name="Habitat",type=HABcol,labels = c("EN", "CF", "PA","CT","VM","KS","FP"))+
    # geom_hline(yintercept=c(0.5,0.2,0.1))+
    scale_x_continuous(limits=c(0,1600),breaks=seq(0,1500,200),expand=c(0.01,0.01))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0.01,0))+
    labs(x="Degree Days", y="fSCA")+
    #geom_vline(xintercept=c(800,900),linetype="dotted",color="black")+
    theme_bw()+
    theme(legend.position = c(0.85,0.75),legend.direction = "vertical",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),legend.title=element_text(size=14,,family="MS Reference Sans Serif"),legend.text=element_text(size=14,,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))
  
  # NGS FIG. 5 LONG-TERM TRENDS OF DD (20250124) ----
    # NGS FIG. 5a Long-term changes in GDD ----
  DOYc <- 227 # mid August
  setwd(WD); 
  tmp    <- read.csv("HISwdd.csv")
  HIS212 <- data.frame(YEAR   = as.numeric(gsub("X","",colnames(tmp))),
                       DD     = unlist(as.vector(tmp[DOYc,]))
                        )
  
  tmp <- read.csv("MONwdd.csv")
  MON212 <- data.frame(YEAR  = as.numeric(gsub("X","",colnames(tmp))),
                       DD    = as.vector(unlist(as.vector(tmp[DOYc,])) )
                       )
  
  tmp <- read.csv("FXAwdd.csv")
  FXA212 <- data.frame(YEAR  = as.numeric(gsub("X","",colnames(tmp))),
                       DD    = as.vector(unlist(as.vector(tmp[DOYc,])) )
                        )
  
  # assemble data
  INI212    <- data.frame(YEAR=1780:2025)
  DATA      <- INI212%>%left_join(HIS212,by="YEAR")
  DATA      <- DATA%>%left_join(MON212,by="YEAR")
  DATA      <- DATA%>%left_join(FXA212,by="YEAR")
  colnames(DATA)[2:4] <- c("DD.HIS","DD.MON","DD.FXA")
  
  BUG <- mean(DATA$DD.HIS-DATA$DD.MON,na.rm=T) # -49 DD ! where do they come from !
  
  # add a gaussian filter 
  DATA$FIL.MON    <- smoother::smth.gaussian(DATA$DD.MON,8)
  DATA$FIL.HIS    <- smoother::smth.gaussian(DATA$DD.HIS,8)
  DATA$FIL.FXA    <- smoother::smth.gaussian(DATA$DD.FXA,3)
  
  # ggplot
  graphics.off()
  windows(6,8)
  ggplot(data=DATA,aes(x=YEAR,y=DD.MON))+
    geom_point(color="gray",size=2)+     
    geom_point(aes(x=YEAR,y=DD.FXA),color="blue",size=2)+ 
    geom_point(aes(x=YEAR,y=DD.HIS),color="black",size=1)+
    geom_line(aes(x=YEAR,y=FIL.MON),linewidth=2,alpha=2,color="gray")+
    geom_line(aes(x=YEAR,y=FIL.HIS),linewidth=2,alpha=2,color="black")+
    geom_line(aes(x=YEAR,y=FIL.FXA),linewidth=1,alpha=2,color="blue")+
    ylab("Warming Degree Days (from 01-01 to 15-08)") + xlab("")+
    geom_hline(yintercept=1000,linetype="solid",color="gray",lwd=0.5)+
    scale_x_continuous(expand=c(0.02,0.01),limits=c(1780,2025),breaks=seq(1780,2020,20))+
    scale_y_continuous(expand=c(0.02,0.01),limits=c(700,1650),breaks=seq(700,1650,100))+
    theme_bw()+
    theme(legend.position = "bottom",legend.direction = "horizontal",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))+
    coord_flip()

    # NGS FIG. 5b DOY for reaching specific DD ----
    DOYc <- 227 
    setwd(WD); 
    BIO.DD <- read.csv("BIOCLIM.DD.csv")
    MED.HABITAT    <- aggregate(BIO.DD$DD,list(BIO.DD$COM),quantile,probs=c(0.5))
 
    setwd(WD); fh    <- read.csv("HISwdd.csv")
    
   # DOY when DD equals the median of snowbeds habitat
   
   CRIT <- function(x) {length(which(x==(-1)))}
   dfh <- data.frame(YEAR   = as.numeric(gsub("X","",colnames(fh))),
                     DOY.EN = unlist(lapply(sign(fh-round(MED.HABITAT$x)[1]), CRIT)),
                     DOY.CF = unlist(lapply(sign(fh-round(MED.HABITAT$x)[2]), CRIT)),
                     DOY.PA = unlist(lapply(sign(fh-round(MED.HABITAT$x)[3]), CRIT)),
                     DOY800 = unlist(lapply(sign(fh-800), CRIT)),
                     DOY900 = unlist(lapply(sign(fh-900), CRIT)),
                     DOY.CT = unlist(lapply(sign(fh-round(MED.HABITAT$x)[4]), CRIT))
   )
   
   # how to process years which do not hit 365
   dfh$DOY.EN <- ifelse(dfh$DOY.EN==365,NA,dfh$DOY.EN)
   dfh$YEARf   <- as.factor(dfh$YEAR)
   
   # comment centrer les couleurs 0 ?
   # dfh$DIFF <- dfh$DOY.EN-214
   # tmp <- cut(dfh$DIFF, quantile(dfh$DIFF, breaks=seq(0,1,0.2),na.rm=T))
   # DECAL=20
   # levels(tmp) <- brewer.pal(11, "RdBu")[2+c(2,4,7,8,9)]
   # dfh$COL214  <- as.character(as.vector(tmp))
   # 
   # CUR.EN <- mean(dfh[116:135,"DOY.EN"])
   # dfh[121:135,"YEAR"]
   
   
   # DEALING with missing values
   # This function applies a filter to a time series with potentially missing data 
   
   filter_with_NA(dfh$DOY.EN)
   
   graphics.off(); windows(6,8)
   ggplot(data=dfh,linewidth=2)+
     geom_point(aes(x=YEAR,y=DOY.EN),alpha=0.8,size=1)+
     geom_segment(aes(x=YEAR,xend=YEAR,y=DOYc,yend=DOY.EN),color=HABcol[1],linewidth=1,alpha=0.5,show.legend=F)+
     geom_line(aes(x=YEAR,y=filter_with_NA(DOY.EN)),linewidth=2,alpha=0.8,color=HABcol[1])+
     geom_line(aes(x=YEAR,y=smoother::smth.gaussian(DOY.CF,8)),linewidth=2,alpha=0.8,color=HABcol[2])+
     geom_line(aes(x=YEAR,y=smoother::smth.gaussian(DOY800,8)),linewidth=1,alpha=0.8,color="black",linetype="dashed")+
     geom_line(aes(x=YEAR,y=smoother::smth.gaussian(DOY900,8)),linewidth=1,alpha=0.8,color="black",linetype="dotdash")+
     geom_line(aes(x=YEAR,y=smoother::smth.gaussian(DOY.PA,8)),linewidth=2,alpha=0.8,color=HABcol[3])+
     scale_x_continuous(expand=c(0.02,0.01),limits=c(1780,2025),breaks=seq(1780,2020,20))+
     scale_y_continuous(expand=c(0.02,0.01),limits=c(150,255),breaks=seq(160,260,20))+
     geom_hline(yintercept=DOYc,linetype="solid",color="black")+
     # geom_hline(yintercept=mean(dfh[121:135,"DOY.EN"]),linetype="dashed",color="black")+
     ylab("SMOD") + xlab("")+
     theme_bw()+
     theme(legend.position = "bottom",legend.direction = "horizontal",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))+
     coord_flip()
   
     # tidyquant::geom_ma(n = 10,linetype = "solid",size=2,color = "black")+
     # tidyquant::geom_bbands(data=df,aes(x=YEAR,y=DD214,high=high,low=low,close=close),n = 10,sd=1)+

     # geom_segment(data= df.seg,aes(x=x1,xend=x2,y=y1,yend=y2),linewidth=3,alpha=0.5,col="green")+       
     # geom_point(data=dfa,aes(x=YEAR,y=DD214),size=2,col="darkgray",alpha=0.6)+
     # geom_line(data=dfa,aes(y=zoo::rollmean(DD214, 10, na.pad=TRUE,align="right")),col="darkgray",linewidth=2.5,alpha=0.8)
     # tidyquant::geom_ma(data=dfa,n = 10,linetype = "solid",color = "blue")+
     # tidyquant::geom_bbands(data=dfa,aes(x=YEAR,y=DD214,high=high,low=low,close=close),n = 10,sd=1)+

   
   # TRY a warming stripes using R climate chart representation -----
   col_strip <- brewer.pal(11, "RdBu")
   brewer.pal.info
   
   graphics.off()
   windows(15,1)
   ggplot(df, aes(x = YEAR, y = 1, fill = DD214-800))+
     geom_tile(show.legend=F)+
     ylab("")+xlab("")+
     scale_y_continuous(expand = c(0, 0))+
     scale_x_continuous(expand = c(0, 0))+
     scale_fill_gradientn(colors = rev(col_strip))+
     guides(fill = guide_colorbar(barwidth = 1))+
     theme_bw()
   
   
   # OLD STUFF ----
   # DD on the First of August (DOY. 214)
   dfa        <- data.frame(YEAR=1943:2023,DD214=as.numeric(fa[214,]))
   dfa$low    <- dfa$DD214-2; dfa$high=dfa$DD214+2
   dfa$close  <- dfa$DD214
   
   
   setwd(WD); df<- read.csv("HIST.DD.csv")
   df$high    <- df$DD214+2; df$low <- df$DD214-2
   df$close   <- df$DD214
   df$DIFF800 <- df$DD214-800
   df$DD800   <- 800
   df$YEARf   <- as.factor(df$YEAR)
   
   tmp <- cut(df$DIFF800, quantile(df$DIFF800, seq(0,1,0.1)))
   levels(tmp) <- rev(brewer.pal(10, "RdBu"))
   df$COL800  <- as.character(as.vector(tmp))
   
   df.seg<- data.frame(x1=c(1889,1945,1995),
                       x2=rep(2023,3),
                       y1=seq(600,1000,200),
                       y2=seq(600,1000,200)
   )
   
  # compute the Accumulation of degree days at FLUXALP
  setwd(WD); FM <- read.csv("FORCAGE_METEO.csv")
  FM$DATE  <- as.Date(FM$DATE)
  FM$MONTH <- lubridate::month(FM$DATE)
  FM$YEAR  <- lubridate::year(FM$DATE)
  FM <- FM[,c("DATE","YEAR","MONTH","FLUXALP.pred")]
  colnames(FM)<-c("DATE","YEAR","MONTH","TM2")
  
  plot(FM$DATE[which(FM$YEAR==1986)],FM$TM2[which(FM$YEAR==1986)],type="l")
  
  # break table per year
  #  Split on userid
  out <- split(FM , f = FM$YEAR)
  # remove years with no data and start at 1945
  out <- out[16:97]
  
  DDf <- function(x){
    test0 <- x$TM2
    test0 <- na.approx(test0,xout=1:length(test0))
    test1 <- ifelse(test0<=0,0,test0)
    test2 <- cumsum(test1)
    
    THR <- seq(100,2000,100)
    RES <- rep(NA,length(THR))
    for (j in 1:length(THR)){
      RES[j] <- which(test2>THR[j])[1]
    }
    return(RES)
  }
  
  TMP <- as.data.frame(lapply(out,FUN=DDf))

  # using segments  time between 1500 and 2000 °C of accumulation
  # using segments  time between 1500 and 2000 °C of accumulation
  i=6;j=10;XLIM <- c(160,265); DXLIM <- c(-8,13)
  graphics.off();windows(5,10);par(mar=c(3,3,1,1),mgp=c(1.2,0.4,0))
  plot(x=0,y=0,type="n",xlim=XLIM,ylim=-c(2024,1943),ylab="",xlab="DOY",las=1,yaxt="n")
  abline(h=-seq(1951,2021,10),lty=2)
  segments(x0=as.numeric(TMP[i,81:1]),y0=-(2023:1943),x1=as.numeric(TMP[j,81:1]),lwd=3,col="gray")
  points(as.numeric(TMP[j,81:1]),-(2023:1943),pch=21,bg="lightblue",cex=1.5)
  points(as.numeric(TMP[i,81:1]),-(2023:1943),pch=21,bg="pink2",cex=1.5)
  lines(data.table::frollmean(as.numeric(TMP[j,81:1]),n=3,align="center"),-(2023:1943),type="l",col="blue",lwd=3) 
  lines(data.table::frollmean(as.numeric(TMP[i,81:1]),n=3,align="center"),-(2023:1943),type="l",col="red",lwd=3) 
  abline(v=seq(180,260,40),lty=2)
  abline(v=225,lwd=2)
  points(IGN.fSCA$DOY,-HIST.fSCA$YEAR,pch=22,bg="black",cex=1.5)
  axis(2,ylim=-c(2023,1943),labels=seq(2023,1943,-2),at=-seq(2023,1943,-2),las=1,tck=0.05)
  axis(side = 2, ylim=-c(2023,1943),at=-seq(2022,1944,-2), labels = FALSE, tck=0.02) 
  
  graphics.off();windows(3.5,10);par(mar=c(3,3,1,1),mgp=c(2.1,0.4,0))
  SHIFT      <- as.numeric(TMP[j,81:1])-as.numeric(TMP[i,81:1])
  NORM.SHIFT <- SHIFT-mean(SHIFT[15:35])
  A<-barplot(NORM.SHIFT,names=1943:2023,horiz=T,las=1,yaxt="n",xlim=DXLIM,col="lightgray",xlab="Anomaly of the elapsed time\n between 600°C and 1000°C")
  abline(v=0)
  lines(data.table::frollmean(NORM.SHIFT,n=3,align="center"),A[,1],lwd=3)
  axis(2,at=A[seq(1,81,2),1],labels=seq(2023,1943,-2),las=1,tck=0.05)
  axis(2, at=A[seq(2,80,2),1], labels = FALSE, tck=0.02) 
  
# VIII. SUPPLEMENTARY MATERIAL (20250123) ----
  # SUPP FIG. 1 HYPSOMETRY (20250314) ----
     # SUPP. FIG. 1a TOPOGRAPHY ----

  graphics.off(); windows(10,12); 
  terra::plot(myDEM,axes=F,alpha=0.8,col=SRTM_pal2(100),
              mar=c(5,0,0,0),
              plg=list( # parameters for drawing legend
                title = "", 
                title.adj=0, # cannot left justify !
                xjust=0,
                loc="bottom",
                horiz=TRUE,
                ext = c(965700, 968000, 6442800, 6442900),
                title.cex = 1, # Legend title size
                cex = 1,
                bty="o",
                bg="white"# Legend text size
              ))
  terra::plot(myHILL,add=T,alpha=0.5,col=grey(0:100 / 100),legend=F)
  terra::plot(my3CONT,add=T,border="black",lwd=2)
  terra::plot(SAG3lines,add=T,col=gray(0.5))
  terra::sbar(1000, xy="bottomleft", type="bar", below="km", label=c(0,0.5,1), cex=.8)
  terra::north(xy="topleft",type=2)
  points(OBS.PTS[2],bg="black",pch=22,cex=1)
  terra::plot(myROUT,col=gray(0.6),add=T,lwd=3)
  
     # SUPP. FIG. 1b Elevation-DAH for the 3 WS ----
    # compute on 2m res to reduce computation time
  setwd(DATA_WD)
  dem2m <- terra::rast("DEM_10M.tif")
  dah2m <- terra::rast("DAH_10M.tif")
  
  WS <- c("2.Roche Noire", "3.Lauzette", "1.Mandette")
  DEM.l <- list()
  for (i in 1:3){
    print(i)
    tmpCONT <- my3CONT[i,]
    tmp1 <- terra::mask(dem2m,tmpCONT); 
    tmp2 <- terra::mask(dah2m,tmpCONT)
    df <- data.frame(ELEV=tmp1[],DAH=tmp2[])
    dfs <- df %>% drop_na()
    colnames(dfs) <- c("ELEV","DAH")
    dfs$WS <- WS[i]
    DEM.l[[i]] <- dfs
  }
  
  df <- do.call(rbind,DEM.l)
  
  graphics.off()
  windows(5,9)
  ggplot(df, aes(x=DAH, y=ELEV)) +
    ggdensity::geom_hdr(show.legend=T,probs = c(0.99, 0.95, 0.75, 0.5,0.25))+
    labs(x="DAH", y="Elevation")+
    geom_vline(xintercept=0)+
    scale_y_continuous(limits=c(1800,3100),breaks=seq(1800,3100,200))+
    scale_x_continuous(limits=c(-0.75,0.75),breaks=seq(-0.5,0.5,0.5))+
    facet_wrap(~WS,ncol= 1,labeller=as_labeller(c("1.Mandette"="Mandette","2.Roche Noire"="Roche Noire","3.Lauzette"="Lauzette")))+
    theme_bw()+
    theme(legend.position = "bottom", legend.background = element_rect(fill = "white", colour = NA),axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"),legend.box.background = element_rect(color = "black"),legend.box.margin = margin(t = 2, l = 2),legend.direction="horizontal")
  
  # SUPP FIG. 2  CONSISTENCY indices of estimated fSCA05 maps (20250314) ----
    # SUPP. FIG. 2a for fSCA01 ----
  SCA.f05 <- terra::rast("SCA.f01.tif")
  for (i in 1:11) SCA.f05[[i]]<- terra::mask(terra::subst(SCA.f05[[i]],NA,0),my3CONT)
  SCA.comp <- SCA.f05
  SCA.comp <- SCA.f05
  OA <- data.frame(matrix(NA,11,11))
  rownames(OA)<-colnames(OA)<-2013:2023
  KC <- OA
  for (j in 1:11){
    print(j)
    SCA.comp <- SCA.f05
    for (i in j:11){
      TP = length(which(SCA.f05[[i]][]==1 & SCA.f05[[j]][]==1))
      TN = length(which(SCA.f05[[i]][]==0 & SCA.f05[[j]][]==0))
      FP = length(which(SCA.f05[[i]][]==0 & SCA.f05[[j]][]==1))
      FN = length(which(SCA.f05[[i]][]==1 & SCA.f05[[j]][]==0))
      O  = (TP+TN)/N
      R  = (((TP+FP)*(TP+FN) )+((FN+TN)*(FP+TN) ) )/(N*N)
      
      SCA.comp[[i]] <- abs(SCA.f05[[i]]-SCA.f05[[j]]) # map of overall accuracy
      OA[i,j]       <- O
      KC[i,j]       <- (O-R)/(1-R)
    }
  }
  rownames(OA)<-colnames(OA)<-2013:2023
  
  OA$YEAR <- rownames(OA)
  
  longData<-reshape2::melt(OA,id.vars="YEAR")
  longData<-longData[longData$value<=1,]
  colnames(longData)<- c("YEAR1","YEAR2","SIMIL")
  longData$SIMIL <-as.numeric(longData$SIMIL)
  DOWN <- longData %>% drop_na()
  
  KC <- as.data.frame(t(KC))
  KC$YEAR <- rownames(KC)
  longData<-reshape2::melt(KC,id.vars="YEAR")
  longData<-longData[longData$value<=1,]
  colnames(longData)<- c("YEAR1","YEAR2","SIMIL")
  longData$SIMIL <-as.numeric(longData$SIMIL)
  UP <- longData %>% drop_na()
  
  DATA <- rbind(UP,DOWN)
  DATA <- DATA[DATA$SIMIL<1,]
  
  graphics.off(); windows(6,6)
  ggplot(data=DATA,aes(x = YEAR1, y = YEAR2)) + 
    geom_raster(aes(fill=SIMIL),show.legend=F)+
    xlab("")+ylab("")+
    geom_text(aes(label = round(SIMIL,2)),size=4)+
    scale_fill_gradient(low='#e7d4e8', high='#7fbf7b')+
    theme_bw() + 
    theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
          axis.text.y=element_text(size=10),
          plot.title=element_text(size=11))
  
    # SUPP. FIG. 2b for fSCA05 ----
  # use of two indicators of accuracy
  # Overall accuracy  OA = TP + TN / N
  # Kappa coefficient KC = (O) − (R) / 1 − (R) 
 # where TP=True Positive; TN=True negative
  SCA.f05 <- terra::rast("SCA.f05.tif")
  for (i in 1:11) SCA.f05[[i]]<- terra::mask(terra::subst(SCA.f05[[i]],NA,0),my3CONT)
  plot(SCA.f05[[1]],legend=F)
  N <- length(which(!is.na(values(SCA.f05[[1]]))==TRUE)) # 17 003 non NA cells
  SCA.comp <- SCA.f05
  OA <- data.frame(matrix(NA,11,11))
  rownames(OA)<-colnames(OA)<-2013:2023
  KC <- OA
  for (j in 1:11){
    print(j)
    SCA.comp <- SCA.f05
  for (i in j:11){
    TP = length(which(SCA.f05[[i]][]==1 & SCA.f05[[j]][]==1))
    TN = length(which(SCA.f05[[i]][]==0 & SCA.f05[[j]][]==0))
    FP = length(which(SCA.f05[[i]][]==0 & SCA.f05[[j]][]==1))
    FN = length(which(SCA.f05[[i]][]==1 & SCA.f05[[j]][]==0))
    O  = (TP+TN)/N
    R  = (((TP+FP)*(TP+FN) )+((FN+TN)*(FP+TN) ) )/(N*N)
    
    SCA.comp[[i]] <- abs(SCA.f05[[i]]-SCA.f05[[j]]) # map of overall accuracy
    OA[i,j]       <- O                # overall accuracy
    KC[i,j]       <- (O-R)/(1-R)      # Kappa coefficient
  }
  }
  rownames(OA)<-colnames(OA)<-2013:2023
  
  OA$YEAR <- rownames(OA)
  
  longData<-reshape2::melt(OA,id.vars="YEAR")
  longData<-longData[longData$value<=1,]
  colnames(longData)<- c("YEAR1","YEAR2","SIMIL")
  longData$SIMIL <-as.numeric(longData$SIMIL)
  DOWN <- longData %>% drop_na()
  
  KC <- as.data.frame(t(KC))
  KC$YEAR <- rownames(KC)
  longData<-reshape2::melt(KC,id.vars="YEAR")
  longData<-longData[longData$value<=1,]
  colnames(longData)<- c("YEAR1","YEAR2","SIMIL")
  longData$SIMIL <-as.numeric(longData$SIMIL)
  UP <- longData %>% drop_na()
  
  DATA <- rbind(UP,DOWN)
  DATA <- DATA[DATA$SIMIL<1,]
  
  HABcol <- c('#762a83','#af8dc3','#e7d4e8','#f7f7f7','#d9f0d3','#7fbf7b','#1b7837')
  graphics.off(); windows(6,6)
  ggplot(data=DATA,aes(x = YEAR1, y = YEAR2)) + 
    geom_raster(aes(fill=SIMIL),show.legend=F)+
    xlab("")+ylab("")+
    geom_text(aes(label = round(SIMIL,2)),size=4)+
    scale_fill_gradient(low='#e7d4e8', high='#7fbf7b')+
    theme_bw() + 
    theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=10),
                       plot.title=element_text(size=11))
  
    # SUPP. FIG. 2c for fSCA09 ----
  SCA.f05 <- terra::rast("SCA.f09.tif")
  for (i in 1:11) SCA.f05[[i]]<- terra::mask(terra::subst(SCA.f05[[i]],NA,0),my3CONT)
  N <- length(which(!is.na(values(SCA.f05[[1]]))==TRUE)) # 17 003 non NA cells
  SCA.comp <- SCA.f05
  OA <- data.frame(matrix(NA,11,11))
  rownames(OA)<-colnames(OA)<-2013:2023
  KC <- OA
  for (j in 1:11){
    print(j)
    SCA.comp <- SCA.f05
    for (i in j:11){
      TP = length(which(SCA.f05[[i]][]==1 & SCA.f05[[j]][]==1))
      TN = length(which(SCA.f05[[i]][]==0 & SCA.f05[[j]][]==0))
      FP = length(which(SCA.f05[[i]][]==0 & SCA.f05[[j]][]==1))
      FN = length(which(SCA.f05[[i]][]==1 & SCA.f05[[j]][]==0))
      O  = (TP+TN)/N
      R  = (((TP+FP)*(TP+FN) )+((FN+TN)*(FP+TN) ) )/(N*N)
      
      SCA.comp[[i]] <- abs(SCA.f05[[i]]-SCA.f05[[j]]) # map of overall accuracy
      OA[i,j]       <- O
      KC[i,j]       <- (O-R)/(1-R)
    }
  }
  rownames(OA)<-colnames(OA)<-2013:2023
  
  OA$YEAR <- rownames(OA)
  
  longData<-reshape2::melt(OA,id.vars="YEAR")
  longData<-longData[longData$value<=1,]
  colnames(longData)<- c("YEAR1","YEAR2","SIMIL")
  longData$SIMIL <-as.numeric(longData$SIMIL)
  DOWN <- longData %>% drop_na()
  
  KC <- as.data.frame(t(KC))
  KC$YEAR <- rownames(KC)
  longData<-reshape2::melt(KC,id.vars="YEAR")
  longData<-longData[longData$value<=1,]
  colnames(longData)<- c("YEAR1","YEAR2","SIMIL")
  longData$SIMIL <-as.numeric(longData$SIMIL)
  UP <- longData %>% drop_na()
  
  DATA <- rbind(UP,DOWN)
  DATA <- DATA[DATA$SIMIL<1,]
  
  graphics.off(); windows(6,6)
  ggplot(data=DATA,aes(x = YEAR1, y = YEAR2)) + 
    geom_raster(aes(fill=SIMIL),show.legend=F)+
    xlab("")+ylab("")+
    geom_text(aes(label = round(SIMIL,2)),size=4)+
    scale_fill_gradient(low='#e7d4e8', high='#7fbf7b')+
    theme_bw() + 
    theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
          axis.text.y=element_text(size=10),
          plot.title=element_text(size=11))

  # SUPP FIG. 3 EVALUATION of RS SMOD (20250123) ----
    # SUPP. FIG. 3a COMPARE GROUND and RS SMOD at FLUXALP and NIVOSE (20250314) ----
  setwd(DATA_WD); LF <- list.files(pattern=".SMODs.tif")
  SMODs <- terra::rast(LF);names(SMODs) <- 2013:2023
  DATA.RS <- t(terra::extract(terra::focal(SMODs,w=3,"mean"),OBS.PTS))[-1,]
  # DATA.RS <- t(terra::extract(SMODs,OBS.PTS))[-1,]
  
  df <- data.frame(YEAR=2013:2023,DATA.RS)
  colnames(df)[2:3]<- c("NIVOSE","FLUXALP")
  dfs <- reshape2::melt(df,id.vars=1,measure.vars=c(2,3))
  colnames(dfs)<-c("YEAR","SITE","RS")
  dfs$GROUND <- NA
  
  setwd(WD); DATA1 <- read.csv("SMOD.FLUX.csv")[,c("YEAR","SMOD_ndvi")]
  setwd(WD); DATA2 <- read.csv("SMOD.NIVO.csv")[,c("YEAR","SMOD_nivose")]
  
  dfs$GROUND[12:21] <- DATA1$SMOD_ndvi        # 2013:2022
  dfs$GROUND[4:10]  <- DATA2$SMOD_nivose[2:8] # 2016:2022
  
  dfs <- dfs %>% drop_na()
  
  # run the smatr test
  mySMA  <- smatr::sma(RS ~ GROUND, data=dfs)
  summary(mySMA)
  SMAfit <- data.frame(t(mySMA$coef[[1]]))
  SMAfit$GROUP <- 1:3
  
  library(ggpmisc)
  graphics.off();windows(5,5);
  ggplot(data=dfs,aes(x=GROUND,y=RS))+
  geom_point(size=2)+
    xlab("SMOD from weather Station")+ ylab("SMOD from satellite imagery")+  
    geom_abline(slope=1,intercept=0,linetype="dashed")+
    stat_ma_line() +
    stat_ma_eq(method="SMA",use_label(c("eq","R2")),rr.digits=2,parse=T)+
    # geom_point(data=dfs,aes(x=GROUND,y=RS,color=SITE),size=3)+
    # geom_abline(data=SMAfit[2,],aes(slope=slope,intercept=elevation,group=GROUP),linewidth=1)+
    scale_x_continuous(limits=c(120,160),expand=c(0,0))+
    scale_y_continuous(limits=c(105,175),expand=c(0,0))+
    theme_bw()+
    theme(legend.position = c(0.7,.06), legend.background = element_rect(fill = "white", colour = NA),axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 11,family="MS Reference Sans Serif"),legend.box.background = element_rect(color = "black"),legend.box.margin = margin(t = 1, l = 1),legend.direction="horizontal",plot.margin = margin(6,10,6,6))
  
  
    # SUPP. FIG. 3b COMPARE GROUND and RS at SOIL TEMP sites (20250314) ----
  # better to keep the original value of 
  setwd(DATA_WD); LF <- list.files(pattern=".SMODs.tif")
  SMODs <- terra::rast(LF);names(SMODs) <- 2013:2023
  DATA.RS <- t(terra::extract(terra::focal(SMODs,w=3,"mean"),SOIL.PTS))[-1,]
   #DATA.RS <- t(terra::extract(SMODs,SOIL.PTS))[-1,]
  colnames(DATA.RS)<- SOIL.PTS$ID
  # remove the plots out myCONT
  DATA.RS <- DATA.RS[,-which(is.na(DATA.RS[1,]))]
  boxplot(t(DATA.RS))
  
  setwd(WD); tmp     <- t(read.csv("SMOD.SOIL.csv",row.names=1))
  rownames(tmp) <- sub("X","",  rownames(tmp))
  
  SEL1    <- colnames(tmp)%in%colnames(DATA.RS)
  DATA.SOIL <- tmp[,SEL1]
  
  SEL     <- colnames(DATA.RS)%in%colnames(tmp)
  DATA.RS <- DATA.RS[,SEL] # 65 sites remaining
  DATA.RS <- DATA.RS[,order(colnames(DATA.RS))]
  
  COM <- substr(colnames(DATA.RS),5,6)
  mySEL <- which(COM=="VM" | COM=="PA" | COM=="CF")
  all.equal(colnames(DATA.SOIL),colnames(DATA.RS))
  
  table(c(!is.na(DATA.SOIL))) # 404 sites x year data  

  cor(c(DATA.RS),c(DATA.SOIL),use="complete.obs") # 0.68
  
  dfs= data.frame(reshape2::melt(DATA.RS),reshape2::melt(DATA.SOIL))[,c(1:3,6)]
  colnames(dfs) <- c("YEAR","SITE","RS","GROUND")
  dfs$COM <- substr(dfs$SITE,5,6)
  
  dfs <- dfs%>%drop_na()
  
  mySEL <- which(dfs$COM=="VM" | dfs$COM=="PA" | dfs$COM=="CF" | dfs$COM=="SH")
  mySEL <- which(dfs$COM=="VM" | dfs$COM=="PA" | dfs$COM=="CF" | dfs$COM=="TR")
  dfsel <- dfs[mySEL,]
  cor(dfsel[,"RS"],dfsel[,"GROUND"]) # 0.85
  
  # smatr test
  (mySMA  <- smatr::sma(RS ~ GROUND, data=dfsel))
  summary(mySMA)
  SMAfit <- data.frame(t(mySMA$coef[[1]]))
  
  library(ggpmisc)
  graphics.off();windows(5,5);
  ggplot(data=dfsel,aes(x=GROUND,y=RS))+
    geom_point(size=2)+
    xlab("SMOD from Soil Temperature")+ ylab("SMOD from Satellite imagery")+  
    geom_abline(slope=1,intercept=0,linetype="dashed")+
    stat_ma_line()+
    stat_ma_eq(method="SMA",use_label(c("eq","R2")),rr.digits=2,parse=T)+
    # geom_point(data=dfsel,aes(x=GROUND,y=RS,color=COM),size=3)+
    theme_bw()+
    theme(legend.position = c(0.66,.06), legend.background = element_rect(fill = "white", colour = NA),axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 11,family="MS Reference Sans Serif"),legend.box.background = element_rect(color = "black"),legend.box.margin = margin(t = 1, l = 1),legend.direction="horizontal",plot.margin = margin(6,10,6,6))
  
  
  # PB <- abs(DATA.RS-DATA.SOIL)
  # boxplot(PB,las=2);abline(h=c(5,10),lty=2)
      
  # SUPP FIG. 4 ANNUAL VARIATION OF SMOD ----
    # SUPP. FIG. 4a SMOD MAPS (20250123) -----
  setwd(DATA_WD); LF <- list.files(pattern=".SMODs.tif") 
  SMODs <- terra::rast(LF);names(SMODs) <- 2013:2023
  for(i in 1:length(myYEAR)) SMODs[[i]][which(SMODs[[i]][]<=41)]<-41
  
  graphics.off();windows(7,6)
  terra::plot(terra::mask(SMODs,my3CONT), smooth=F,stretch="hist",col=rev(DIVcbs_VGpal(100)),mar=c(0.1, 0.1, 0, 0.1),loc.main="bottomleft",cex.main=1.5,legend=F,axes=F,range=c(40,270),fun=function() {
    terra::plot(terra::mask(my3HILL,my3CONT),add=T,alpha=0.1,col=grey(0:100 / 100),legend=F);
    terra::plot(SAG3lines,add=T,col=gray(0.5));
    terra::plot(my3CONT,add=T,border="black");  
    points(OBS.PTS[2],bg="black",pch=22,cex=1)
  })
  
    # SUPP. FIG. 4b VIOLIN PLOTS of SMOD for 2013-2023 (20210123) ----
    # DO NOT WORK ANY LONGER on 20250120 - FIXED with bypassing a tidy step !!!!
  setwd(DATA_WD); LF <- list.files(pattern=".SMODs.tif")
  SMODs <- terra::rast(LF);names(SMODs) <- 2013:2023
  tmp0  <- terra::mask(SMODs,my3CONT)
  tmp1  <- round(terra::focal(tmp0,w=3,fun="mean",na.rm=T))
  df    <- as.data.frame(tmp1[])
  df    <- df %>% drop_na()

  df.gsb <- reshape2::melt(df)
  colnames(df.gsb) <- c("Year","SMOD")
  df.gsb$Year <- as.factor(df.gsb$Year)
  
  # simple volin - without filling the volin 
  ggplot(data=df.gsb, aes(x=Year, y=SMOD)) + 
    geom_violin(trim=F)
  
  # more complicated violin - with filling
  # This is all you need for the fill:   
  p<- ggplot(data=df.gsb, aes(x=Year, y=SMOD,fill=SMOD)) + 
    geom_violin(trim=F)
  mywidth <- .5 # bit of trial and error

  vl_fill <- data.frame(ggplot_build(p)$data) %>%
    mutate(xnew = x- mywidth*violinwidth, xend = x+ mywidth*violinwidth) 
  
  vl_poly <- vl_fill[,c("xnew", "xend", "y", "group")]
  vl_poly <- pivot_longer(vl_poly,cols=-c(y, group), names_to = "oldx", values_to = "x")
  vl_poly <- arrange(vl_poly,y)
  vl_pol <-  split(vl_poly,"oldx")
  
  vl_poly <- 
  vl_pol %>%
  map(., function(x) {
    if(all(x$oldx == "xnew")) x <- arrange(x, desc(y))
    x
  }) %>%
    bind_rows()
  
    
  # vl_poly <-
  #   vl_fill %>%
  #   select(xnew, xend, y, group) %>%
  #   pivot_longer(-c(y, group), names_to = "oldx", values_to = "x") %>%
  #   arrange(y) %>%
  #   split(., .$oldx) %>% # ISSUE HERE
  #   map(., function(x) {
  #     if(all(x$oldx == "xnew")) x <- arrange(x, desc(y))
  #     x
  #   }) %>%
  #   bind_rows()
  
  
  df.gsb$Year <- as.numeric(as.vector(df.gsb$Year))-2012
  table(df.gsb$Year)
  
  graphics.off()
  windows(10,5)
  ggplot() +
    geom_polygon(data = vl_poly, aes(x, y, group = group), 
                 color= "gray", linewidth = 1, fill = NA,show.legend=F) +  
    ylab("SMOD")+ xlab("")+   
    scale_color_gradientn(colours = rev(DIVcbs_VGpal(50)))+
    geom_segment(data = vl_fill, aes(x = xnew, xend = xend, y = y, yend = y, color = y),show.legend=F)+
    scale_y_continuous(breaks=seq(50,250,50),limits=c(30,260))+
    scale_x_continuous(breaks=1:11,labels=2013:2023,limits=c(0.5,11.5),expand = c(0,0))+
    stat_summary(data=df.gsb,aes(x=Year, y=SMOD,group=Year),fun.data = f, 
                 geom="boxplot", 
                 position=position_dodge(1),
                 width=0.1,linewidth=0.5)+
    theme_bw()+
    theme(axis.text=element_text(size=14,family="MS Reference Sans Serif"),axis.title=element_text(size=14,family="MS Reference Sans Serif"),strip.text = element_text(size = 14,family="MS Reference Sans Serif"),plot.margin = margin(1,1,2,2))
  
  # SUPP FIG. 5 fSCA-DD validation for the period 1984-2013 (20250123) ----
    # SUPP. FIG. 5a fSCA-DOY for TOT - training period (20250315) ----
  # fSCA data
  setwd(WD); tmp   <- read.csv("THE.fSCA.csv")
  tmp.l <- split(tmp,tmp$WS)
  df    <- tmp.l[[4]] 
  df$YEAR <- as.factor(df$YEAR)
  
  # parameters of model
  setwd(WD); PARAM <- read.csv("PARAM_TOT.csv")
  XSEQ <- 61:240
  df.mod <- NULL
  for (i in 1:nrow(PARAM)){
    print(i)
    mydf <- data.frame(YEAR=PARAM$YEAR[i],
               DOY = XSEQ,
               MOD = sigmoid(as.numeric(PARAM[i,2:4]),XSEQ)
               )
    df.mod <- rbind(df.mod,mydf)
  }
  df.mod$YEAR <- as.factor(df.mod$YEAR)
  
  graphics.off()
  windows(4,4)
  ggplot(data = df, aes(x = DOY, y = fSCA/100, col=YEAR)) + 
    geom_point(size=2) + 
    geom_line(data=df.mod,aes(x = DOY, y = MOD, col=YEAR))+
    # geom_hline(yintercept=c(0.5,0.2,0.1))+
    scale_x_continuous(limits=c(0,260),breaks=seq(0,260,50),expand=c(0.01,0))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0.01,0))+
    labs(x="DOY", y="fSCA")+
    theme_bw()+
    guides(col = guide_legend(ncol = 2))+
    theme(legend.position = c(0.25,0.3),legend.direction = "vertical",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))
  
    # SUPP. FIG. 5b fSCA-DD for TOT (20250315) ----
  setwd(WD); THE <- read.csv("THE.fSCA.csv")   # select training dataset
  df <- THE
  df$fSCA <- df$fSCA/100
  df$WS <- as.factor(df$WS)
  levels(df$WS) <- c("3.LAUZ","1.MAND","2.RNPI","4.TOT")
  df$WS <- as.vector(df$WS)
  
  # remove points in September !
  df <- df[which(df$DOY<=245),]
  # remove a few outliers
  OUT <- which(df$fSCA>0.9 & df$DD>1000)
  if(length(OUT)>0) df <- df[-OUT,]
  
  df.l <- split(df,df$WS)
  TMP  <- df.l[[4]]
  dim(TMP) # 128 images



  (myEVAL.TOT  <- m.eval(TMP$fSCA,TMP$PRED)) # to be reported on figure
  
  (SysErr   <- 100*myEVAL.TOT[,"MSEs"]/myEVAL.TOT[,"MSE"])
  (UnSysErr <- 100*myEVAL.TOT[,"MSEu"]/myEVAL.TOT[,"MSE"])
  
  # resampling of parameters and average
  # confidence interval on parameters
  setwd(WD); SIG.PARAM <- read.csv("SIG.PARAM.csv")
  TMP$PRED <- sigmoid(SIG.PARAM[,1],TMP$DD)  
  
  setwd(WD); SIMUL <- read.csv("SIG.SIMUL.csv")
  MED <- apply(SIMUL, 2, median)
  LOW <- apply(SIMUL, 2, function(x) smedian.hilow(x)[2])
  HIG <- apply(SIMUL, 2, function(x) smedian.hilow(x)[3])
  BOOT <- data.frame(DD=1:2000,lwr=LOW,upr=HIG,fSCA=MED)
  
  # ggplot
  graphics.off()
  windows(4,4)
  ggplot(data = BOOT, aes(x = DD, y = fSCA)) + 
    geom_ribbon(data=BOOT, aes(ymin = lwr, ymax = upr), fill = "gray", alpha = 0.8)+
    geom_line() + 
    scale_x_continuous(limits=c(0,1600),breaks=seq(0,1500,500),expand=c(0.01,0.01))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0.01,0))+
    geom_point(data=TMP,aes(x=DD,y=fSCA))+
    labs(x="Degree Days", y="fSCA")+
    theme_bw()+
    theme(legend.position = c(0.85,0.4),legend.direction = "vertical",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))
  
    # SUPP. FIG. 5c MODEL EVALUATION for the period 1984-2013 (20250315) ----
  setwd(WD) 
  DLR.fSCA <- read.csv("DLR.fSCA.csv");table(DLR.fSCA$WS) 
  GEE.fSCA <- read.csv("GEE.fSCA.csv");table(GEE.fSCA$WS)
  tmp      <- rbind(DLR.fSCA,GEE.fSCA)  # select testing set
  nrow(tmp)/4 # 2284 images
  
  # tmp  <- read.csv("THE.fSCA.csv")
  TMP  <- split(tmp,f=tmp$WS)
  # remove cloudy days for TOT
  myTOT <- TMP[[4]][which(TMP[[4]]$fNA/100<0.1),]
  # remove DOY after 1st of September (DOY 245)
  myTOT <- myTOT[which(myTOT$DOY<215),]
  dim(myTOT) # reste 291+80 = 371 images
  
  myTOT$MEAS <- myTOT$fSCA/100
  plot(myTOT$DD,myTOT$MEAS,pch=".",cex=3)
    
  # import PARAMETERS & PREDICT
  setwd(WD); SIG.PARAM <- read.csv("SIG.PARAM.csv")
  myTOT$PRED <- sigmoid(SIG.PARAM[,1],myTOT$DD)  
  
  # evaluation of the model with or without extreme balues
  (myEVAL   <- m.eval(myTOT$MEAS,myTOT$PRED))
  (SysErr   <- 100*SMA[,"MSEs"]/SMA[,"MSE"])
  (UnSysErr <- 100*SMA[,"MSEu"]/SMA[,"MSE"])
  
  EXT <- which(myTOT$MEAS>=0.95 | myTOT$MEAS<=0.15)
  SMAmEXT <- m.eval(myTOT$MEAS[-EXT],myTOT$PRED[-EXT])
  dim(myTOT[-EXT,]) # 56 images remaining
  
  graphics.off()
  windows(4,4)
  ggplot(data = myTOT, aes(x = MEAS, y = PRED)) + 
    geom_point(size=1.5,col="black") + 
    # geom_hline(yintercept=c(0.5,0.2,0.1))+
    scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0.01,0))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0.01,0))+
    geom_abline(aes(intercept = 0, slope = 1),linetype="dashed")+
    geom_abline(aes(intercept = SMA[,'elev.SMA'], slope = SMA[,'slope.SMA']),col="red")+
    # geom_abline(aes(intercept = SMAmEXT[,'elev.SMA'], slope = SMAmEXT[,'slope.SMA']),col="blue")+
    geom_abline(slope=1,intercept=0,linetype="dashed")+
    ggpmisc::stat_ma_line()+
    ggpmisc::stat_ma_eq(method="SMA",ggpmisc::use_label(c("eq","R2")),rr.digits=2,parse=T)+
    # geom_point(data=dfsel,aes(x=GROUND,y=RS,color=COM),size=3)+
    labs(x="Measured fSCA", y="Predicted fSCA")+
    theme_bw()+
    theme(legend.position = c(0.5,0.5),legend.direction = "horizontal",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))
  
   
    # SUPP. FIG. 5d Residuals per year (20250315) ----
  myTOT$RES <- myTOT$PRED-myTOT$MEAS
  df <- myTOT[which(myTOT$YEAR>=1985),]

  graphics.off();  windows(4,4)
  ggplot(df, aes(x=YEAR, y=RES,group=YEAR))+
    geom_boxplot(fill=c(rep("gray",27),rep("black",11)))+
    geom_hline(yintercept=0,linetype="solid",color="gray")+
    scale_y_continuous(limits=c(-.3,.3),expand=c(0,0.01))+
    scale_x_continuous(breaks=seq(1990,2020,10),expand=c(0.01,0.03))+
    xlab("")+ylab("Residual")+
    theme_bw()+
    theme(legend.position = c(0.5,0.5),legend.direction = "horizontal",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))
  
  # SUPP. FIG. 6 SNOW/NO SNOW distribution using HISTORICAL AERIAL PHOTOGRAPHS ----
  # 1980-08-24
  tmp <- terra::rast("D:/DATA_ROCHENOIRE/BD_HISTORIQUE/OR_DATA/GRAY/IGNF_PVA_1-0__1980-08-24_orthorect.tif")
  crs(tmp) <-"EPSG:2154"
  tmp1 <- terra::crop(tmp,my3CONT)
  terra::plot(tmp1,col=grey(0:100 / 100),legend=F,axes=F)
  plot(my3CONT,add=T,border="yellow")
  
  setwd(WD);read.csv("HIST.fSCA.csv")
  # First example on cold year 24 August 1980 with fSCA at 6%
  setwd(WD); tmp <- terra::rast("1980-08-24_snowbin_1M.tif")
  terra::plot(tmp,col=c("wheat","darkblue"),legend=F,axes=F,colNA="white")
  plot(my3CONT,add=T,border="black")

  # 1960-08-22
  tmp <- terra::rast("E:/DATA_ROCHENOIRE/BD_HISTORIQUE/OR_DATA/STUFF_ARTHUR/ROCHE NOIRE/images/1960-08-22/IGNF_PVA_1-0__1960-08-22_orthorect.tif")
  crs(tmp) <-"EPSG:2154"
  tmp0 <- terra::mask(tmp,my3CONT)
  
  plot(tmp0)
  tmp1 <- terra::crop(tmp,my3CONT)
  terra::plot(tmp1,col=grey(0:100 / 100),legend=F,axes=F)
  plot(my3CONT,add=T,border="yellow")
  
  setwd(WD); tmp <- terra::rast("1960-08-22_snowbin_1M.tif")
  terra::plot(tmp,col="darkblue",legend=F,axes=F,colNA="lightgray")
  plot(myCONT,add=T,border="black")
  
  # 1994-08-15
  tmp <- terra::rast("E:/DATA_ROCHENOIRE/BD_HISTORIQUE/OR_DATA/STUFF_ARTHUR/ROCHE NOIRE/images/1994-08-15/IGNF_PVA_1-0__1994-08-15_orthorect.tif")
  crs(tmp) <-"EPSG:2154"
  tmp1 <- terra::crop(tmp,my3CONT)
  terra::plot(tmp1,col=grey(0:80 / 100),legend=F,axes=F)
  plot(my3CONT,add=T,border="yellow")
  
  # 1978-09-14
  tmp <- terra::rast("E:/DATA_ROCHENOIRE/BD_HISTORIQUE/OR_DATA/STUFF_ARTHUR/ROCHE NOIRE/images/1978-09-14/IGNF_PVA_1-0__1978-09-14_orthorect.tif")
  crs(tmp) <-"EPSG:2154"
  tmp1 <- terra::crop(tmp,my3CONT)
  terra::plot(tmp1,col=grey(0:100 / 100),legend=F,axes=F)
  plot(my3CONT,add=T,border="yellow")
  
  # 2010-07-07
  LF <- list.files("S:/LECA/PLATEAU-PASTIS/GIS_DATA/Alpes/IMAGERY/RAPIDEYE/_05",pattern=".tif",full=T)
  TMP <- terra::crop(terra::rast(LF[13]),my3CONT)
  TMP1 <- terra::subset(TMP,3:1)
  MAX  <- max(terra::global(TMP1,max))
  TMP2 <- 255*terra::subset(TMP1,3:1)/MAX # need to rescale
  terra::plotRGB(TMP2)
  plot(my3CONT,add=T,border="yellow")

  # unused
  # extract the hypsometry of points
  setwd("U:/DATA_ROCHENOIRE/LIDAR/LIDAR_IGN/")
  dem <- terra::rast("DEM_1M.tif")
  slo <- terra::rast("SLO_1M.tif")
  
  sno <- dem[tmp==1][,1];range(sno)
  # sno <- sno[-which(sno<2300)]
  df1 <- data.frame(ELEV=sno);df1$SNOWY="YES" 
  df1$SLOPE <- slo[tmp==1][,1]
  
  nosno <- dem[tmp==0][,1];range(nosno)
  df2 <- data.frame(ELEV=nosno);df2$SNOWY="NO"    
  df2$SLOPE <- slo[tmp==0][,1]
  
  df <- rbind(df1,df2)
  dff <- reshape2::melt(df,id.vars="SNOWY",variable.name="TOPO")
  
  hist(dff$value)
  
  # ggplot for slope/elevation
  ggplot(dff, aes(x=value, fill=SNOWY)) +
    geom_density(alpha=0.6)+
    scale_fill_manual(values=c("beige","darkblue"),name="Cover",labels=c("no snow","snow"))+
    theme_bw()+
    facet_wrap(vars(TOPO),ncol=2,scales = 'free')+
    theme(legend.position = c(0.5,0.9),legend.direction = "vertical",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),legend.box.background = element_rect(color = "black"),legend.box.margin = margin(t = 1, l = 1),strip.text = element_text(size = 12,family="MS Reference Sans Serif"),strip.background=element_rect(colour="black",fill=c("gray","yellow")))
  
  # using ggdensity - quite long - consider subsampling
  dfs <- df[sample(1:nrow(df),100000),]
  ggplot(dfs, aes(ELEV, SLOPE, fill = SNOWY)) +
    ggdensity::geom_hdr(show.legend=T)+
    scale_fill_manual(values=c("yellow","darkblue"),name="Cover",labels=c("no snow","snow"))+
    xlab("Elevation")+
    theme_bw()+
    theme(legend.position ="top",legend.direction = "vertical",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),legend.box.background = element_rect(color = "black"),legend.box.margin = margin(t = 1, l = 1),strip.text = element_text(size = 12,family="MS Reference Sans Serif"),strip.background=element_rect(colour="black",fill=c("gray","yellow")))
  

  # SUPP. FIG. 8 BUBLE PLOT OF LATE SNOWMELT SITES (20250315) ----
  
   # aggregate to cells of 30x30 km = 900 km2 (fact=1000)
   SCALE  <- 30 # in squared in km
   pix_tot <- (SCALE*1000)^2/30^2  # nbre total de pixels de 30m dans SCALE 
   
   myEXT  <- terra::ext(3900000,4800000,2200000, 2800000)
   INIT   <- terra::rast(myEXT,crs="epsg:3035",res=30) # 20000 x 30000

  # NOT TO RUN 
  # LATE SMOD from Landsat average 2004-2023 with SMOD >195 prepared by Arthur in February
   setwd(DATA_WD)
   # tmp    <- terra::rast("SMOD_2004_2023_EUALPS_sup195.tif") # at 30m
   tmp    <- terra::rast("SMOD_2004_2023_smod195to245_EUALPS.tif") # at 30m
   LLSF   <- terra::project(tmp,INIT,method="near")
   setwd(DATA_WD); writeRaster(LLSF,"LLSF_suppmat.tif")
  
  # DEM above 2000m
  dem      <- terra::project(terra::rast("U:/DATA_SNOWBED/EUALPS/DEM_30M_EUALPS.tif"),INIT)
  bin.dem  <- dem; bin.dem[dem<2000]<-0; bin.dem[dem>=2000]<-1
  setwd(DATA_WD); writeRaster(bin.dem,"DEM2k_suppmat.tif")

  # CORINE land cover 2018 class 34 is glacier and perpetual snow
  clc <- terra::project(terra::rast("U:/DATA_LANDCOVER/CORINE_LAND_COVER/CLC_Alpes_2018.tif"),INIT,method="near")
  clc34 <- clc; clc34[clc==34]<-1;clc34[clc!=34]<-0
  setwd(DATA_WD); writeRaster(clc34,"CLC34_suppmat.tif")
  
  # CGLS cover of ice & permanent snow
  ice <- terra::project(terra::rast("U:/DATA_SNOWBED/EUALPS/ICE_CGLS_LC100M_30M_EUALPS.tif"),INIT,method="near")
  setwd(DATA_WD); writeRaster(ice,"CGLSice_suppmat.tif")
 # end of NOT TO RUN
  
  # Randolph Glacier inventory
  RGI <- crop(terra::project(terra::vect("D:/DATA_YUKIGATA/RGI/11_rgi60_CentralEurope.shp"),"epsg:3035"),myEXT)
  RGIrast <- rasterize(RGI,INIT)
  setwd(DATA_WD); writeRaster(RGIrast,"RGIice_suppmat.tif")
  # graphics.off(); windows(15,10);plot(RGIrast)
  
  # COMPUTE AREAS
  # nombre de pixels de 30m dans l'emprise de SCALE
  pix_lsm  <- terra::aggregate(LLSF,fact=SCALE*1000/30,fun="sum",na.rm=T) # 20 x 30
  lsm_area <- 900*pix_lsm/10^6 # area in km2
  
  pix_alp  <- terra::aggregate(bin.dem,fact=SCALE*1000/30,fun="sum",na.rm=T)
  alp_area <- 900*pix_alp/10^6 # en km2
  names(alp_area) <- "alp_area"

  pix_ilc  <- terra::aggregate(clc34,fact=SCALE*1000/30,fun="sum",na.rm=T)
  ilc_area <- 900*pix_ilc/10^6
  
  pix_ice  <- terra::aggregate(ice,fact=SCALE*1000/30,fun="sum",na.rm=T)
  ice_area <- 900*pix_ice/10^6

  pix_rgi  <- terra::aggregate(RGIrast,fact=SCALE*1000/30,fun="sum",na.rm=T)
  rgi_area <- 900*pix_rgi/10^6
  

# assemble into DATAbub
  DATAbub <- data.frame(terra::xyFromCell(lsm_area,which(lsm_area[]>0)))
  tmp  <- terra::vect(DATAbub[,1:2],geom=c("x","y"),crs="epsg:3035")
  tmp1 <- terra::project(tmp,"epsg:4326")
  tmp2 <- terra::crds(tmp1)
  DATAbub <- cbind(tmp2,DATAbub)
  colnames(DATAbub) <- c("x_wgs84","y_wgs84","x","y")

  DATAbub$alp_area <- terra::extract(alp_area,which(lsm_area[]>0))[,1]
  DATAbub$lsm_area <- terra::extract(lsm_area,which(lsm_area[]>0))[,1]
  DATAbub$ilc_area <- terra::extract(ilc_area,which(lsm_area[]>0))[,1]
  DATAbub$ice_area <- terra::extract(ice_area,which(lsm_area[]>0))[,1]
  DATAbub$rgi_area <- terra::extract(rgi_area,which(lsm_area[]>0))[,1]
  DATAbub$ice_area[is.na(DATAbub$ice_area)]<-0 
  DATAbub$rgi_area[is.na(DATAbub$rgi_area)]<-0 

# fraction de lsm dans l'alpin non englacé & save
  DATAbub$lsm_frac <- 100*(DATAbub$lsm_area)/(DATAbub$alp_area)
  DATAbub$ice_frac <- 100*(DATAbub$ilc_area)/(DATAbub$alp_area)

  SEL <- is.finite(DATAbub$lsm_frac)
  DATAbub <- DATAbub[SEL,]
  head(DATAbub)
  setwd(WD); write.csv(DATAbub,"DATAbub.csv",row.names=F)


# area covered in km2
# The highest axial sector still hosts most of the 3500 Alpine glaciers that cover a total area of about 370 km2.
sum(DATAbub$alp_area) # 35455 km2
sum(DATAbub$lsm_area) #  1748 km2

# compare with CGLS ice
sum(DATAbub$ice_area) # 2822 km2

# compare with CORINE ice
sum(DATAbub$ilc_area) # 1970 km2

# compare with RGI ice
sum(DATAbub$rgi_area) # 2087 km2


# prepare dem
tmp.dem <- terra::crop(terra::aggregate(dem,fact=30),myEXT)
names(tmp.dem)<- "elevation"

# add graticule
GRAT <- sf::st_graticule(lon=seq(6,15, 3),
                         lat = seq(44,50, 3),
                         ndiscr = 5000,
                         margin = 0.00000000001)%>%
  st_transform("epsg:3035")
GRAT$crs <- "epsg:3035"


graphics.off(); windows(15,8)
ggplot() +
  geom_raster(data=tmp.dem,aes(x=x,y=y,fill=elevation),show.legend=F)+
  scale_fill_gradient(low = "white", high = gray(0.3))+
  geom_point(data=DATAbub, aes(x=x, y=y, size=lsm_area),colour='#762a83',alpha=0.6)+
  scale_x_continuous(limits = c(3900000,4800000),expand = c(0, 0))+
  scale_y_continuous(limits = c(2200000, 2800000),expand = c(0, 0))+
  geom_sf(data = GRAT, size = 0.3, color = "#424242", alpha = 0.6)+
  theme_bw()+
  xlab("")+ylab("")+
  scale_size_area(max_size=10,name="Area of late-lying snow fields (km2)")+
  theme(legend.position = c(0.5,.05),legend.direction = "horizontal",legend.title=element_text(size=12,family="MS Reference Sans Serif"),axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))

graphics.off(); windows(15,8)
ggplot() +
  geom_raster(data=tmp.dem,aes(x=x,y=y,fill=elevation),show.legend=F)+
  scale_fill_gradient(low = "white", high = gray(0.3))+
  geom_point(data=DATAbub, aes(x=x, y=y, size=rgi_area),colour='#762a83',alpha=0.6)+
  scale_x_continuous(limits = c(3900000,4800000),expand = c(0, 0))+
  scale_y_continuous(limits = c(2200000, 2800000),expand = c(0, 0))+
  geom_sf(data = GRAT, size = 0.3, color = "#424242", alpha = 0.6)+
  theme_bw()+
  xlab("")+ylab("")+
  scale_size_area(max_size=10,name="Area of ice and permanent snow (km2)")+
  theme(legend.position = c(0.5,.05),legend.direction = "horizontal",legend.title=element_text(size=12,family="MS Reference Sans Serif"),axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))


setwd(WD); ggsave("FIGsup8.svg")

   
  # SUPP FIG. 7  GAP - FILLING of SMOD - EXAMPLE of image processing for 2015 (20250122) ----
    # SUPP FIG. 7a Using level 2 products ----
      myY = 15 # pick a year
      SCAtmp.L2 <- terra::rast(paste0("SCA",myY,".L2.tif"))
      SCAtmp.L2[is.na(SCAtmp.L2)]<-(-1)
      for (i in 1:dim(SCAtmp.L2)[3]) terra::coltab(SCAtmp.L2[[i]]) <- coltb
  
    graphics.off();windows(14,16)
    terra::plot(SCAtmp.L2,mar=c(0.1, 0.1, 0.1, 0.1),loc.main="bottomleft",col.main="black",cex.main=1.5,axes=F,alpha=my3ALPHA,maxnl=16,fun=function() {
    terra::plot(my3HILL,add=T,alpha=0.5,col=grey(0:100 / 100),legend=F);
    terra::plot(my3CONT,add=T,border="black"); 
    terra::plot(terra::vect("SAG3lines.shp"),add=T,col="gray");
    terra::plot(OBS.PTS[2],add=T,pch=22,bg="blue")
  })
  
  
    # SUPP FIG. 7b Using level 3 products ----
    myY = 15 # pick a year
    SCAtmp.L3 <- terra::rast(paste0("SCA",myY,".L3.tif"))
    SCAtmp.L3[is.na(SCAtmp.L3)]<-(-1)
    for (i in 1:dim(SCAtmp.L3)[3]) terra::coltab(SCAtmp.L3[[i]]) <- coltb
    
    graphics.off();windows(14,16)
    terra::plot(SCAtmp.L3,mar=c(0.1, 0.1, 0.1, 0.1),loc.main="bottomleft",col.main="black",cex.main=1.5,axes=F,alpha=my3ALPHA,maxnl=16,fun=function() {
      terra::plot(my3HILL,add=T,alpha=0.5,col=grey(0:100 / 100),legend=F);
      terra::plot(my3CONT,add=T,border="black"); 
      terra::plot(terra::vect("SAG3lines.shp"),add=T,col="gray");
      terra::plot(OBS.PTS[2],add=T,pch=22,bg="blue")
    })
    
  # SUPP FIG. 9 ESTIMATED MAPS of fSCA0.5 (20250123) ----
  setwd(DATA_WD); SCA.f05 <- terra::rast("SCA.f05.tif")
  for (i in 1:11) SCA.f05[[i]]<- terra::mask(terra::subst(SCA.f05[[i]],NA,0),my3CONT)
  
  graphics.off();windows(7,6)
  terra::plot(SCA.f05, smooth=F,col=c("lightgreen","deepskyblue3"),mar=c(0.1, 0.1, 0.1, 0.1),legend=F,axes=F,loc.main="bottomleft",cex.main=1.5,
              fun=function() {
                terra::plot(terra::mask(my3HILL,my3CONT),add=T,alpha=0.5,col=grey(0:100 / 100),legend=F);
                terra::plot(terra::vect("SAG3lines.shp"),add=T,col="gray");
                terra::plot(my3CONT,add=T,border="black");  
                terra::plot(OBS.PTS[2],add=T,pch=22,bg="red")
              })
  
  # SUPP FIG. 11 BIPLOT fSCA x SMOD (2025012) ----
    setwd(WD); TMP <- read.csv("SMOD.fSCA.PC.csv")
    LABEL <- paste0("fSCA",substr(101:109,2,3))
   graphics.off();windows(6,5);par(mar=c(4,4,1,1),mgp=c(2,0.6,0));
   boxplot(TMP,xlab="fSCA",ylab="SMOD",names=LABEL,las=1)
   
   graphics.off();windows(2,5);
   df <- reshape2::melt(TMP)
   ggplot(data=df,aes(x=variable,y=value))+
   geom_boxplot(outlier.colour="black", outlier.shape=16,
                outlier.size=2, notch=FALSE)+
     xlab("")+ylab("SMOD")+
     scale_y_continuous(breaks=seq(100,200,50),limits=c(100,220))+
     scale_x_discrete(labels=LABEL)+
     theme_bw()+
     theme(axis.text=element_text(size=14,family="MS Reference Sans Serif"),axis.title=element_text(size=14,family="MS Reference Sans Serif"),strip.text = element_text(size = 14,family="MS Reference Sans Serif"),plot.margin = margin(0,0,0.2,0.2))

  # SUPP FIG. 12 boxplot of SMOD x fSCA  ----  
  setwd(WD); TMP <- read.csv("SMOD.fSCA.PC.csv")
  TMP$YEAR <- 2013:2023
  TMP$DELTA <- TMP[,3]-TMP[,1]
  setwd(WD); DD  <- read.csv("FLUXALP.DD.csv")
  
  TMP.DD <- TMP
  for (i in 1:nrow(TMP)){
    tmpYEAR <- which(colnames(DD)==paste0("X",TMP$YEAR[i]))
    TMP.DD[i,1] <- DD[TMP[i,1],tmpYEAR]
    TMP.DD[i,2] <- DD[TMP[i,2],tmpYEAR]
    TMP.DD[i,3] <- DD[TMP[i,3],tmpYEAR]
  }
  TMP.DD$DELTA <- TMP.DD[,3]-TMP[,1]
  
  graphics.off();windows(4,4);
  df <- reshape2::melt(TMP.DD[,1:3])
  ggplot(data=df,aes(x=variable,y=value))+
    geom_boxplot(outlier.colour="black", outlier.shape=16,
                 outlier.size=2, notch=FALSE)+
    xlab("fSCA")+ylab("SMOD")+
    scale_y_continuous(breaks=seq(0,1200,200),limits=c(0,1200))+
    scale_x_discrete(labels=c("0.9","0.5","0.1"))+
    theme_bw()+
    theme(axis.text=element_text(size=14,family="MS Reference Sans Serif"),axis.title=element_text(size=14,family="MS Reference Sans Serif"),strip.text = element_text(size = 14,family="MS Reference Sans Serif"),plot.margin = margin(0,10,0.2,0.2))
  
  
  graphics.off();windows(4,6);par(mfrow=c(2,1),mar=c(4,3,0,1),mgp=c(2,0.7,0))
  plot(TMP$DELTA,TMP.DD$DELTA,pch=21,bg="gray",cex=1.5,ylim=c(100,900),ylab="Difference in DD",xlab="sigma SMOD")
  abline(h=Hmisc::smedian.hilow(TMP.DD$DELTA,0.8),lty=c(1,2,2))
  
  # same stuff with Arthur's data
  setwd(WD); TMP <- read.csv("SMOD.fSCA.AB.csv"); TMP$YEAR <- 1984:2023
  TMP$DELTA <- TMP[,3]-TMP[,1]
  setwd(WD); DD  <- read.csv("FLUXALP.DD.csv")
  
  TMP.DD <- TMP
  for (i in 1:nrow(TMP)){
    tmpYEAR <- which(colnames(DD)==paste0("X",TMP$YEAR[i]))
    TMP.DD[i,1] <- DD[TMP[i,1],tmpYEAR]
    TMP.DD[i,2] <- DD[TMP[i,2],tmpYEAR]
    TMP.DD[i,3] <- DD[TMP[i,3],tmpYEAR]
  }
  TMP.DD$DELTA <- TMP.DD[,3]-TMP[,1]
  
  graphics.off();windows(4,6);par(mfrow=c(2,1),mar=c(4,3,0,1),mgp=c(2,0.7,0))
  boxplot(TMP.DD[30:40,1:3],ylab="Degree Days",ylim=c(0,900))
  plot(TMP$DELTA[30:40],TMP.DD$DELTA[30:40],pch=21,bg="gray",cex=1.5,ylim=c(100,900),ylab="Difference in DD",xlab="sigma SMOD")
  abline(h=Hmisc::smedian.hilow(TMP.DD$DELTA,0.8),lty=c(1,2,2))
  
  
  
  
  # SUPP FIG. 13 Anomaly of grenness per class of fSCA (TO UPDATE) ----
     # WHERE ARE THE TIFs !!
     # matching a sequence of fSCA and greenness anomaly
   setwd(DATA_WD)
   DEC <- substr(101:109,2,3)
   df <- NULL
   for (i in 1:length(DEC)){
     print(i)
     FILE <- paste0("SCA.f",DEC[i],".tif")
     setwd(WD); SCA <- terra::rast(FILE) # see part II.D.1
     for (j in 1:11) SCA[[j]]<- terra::mask(terra::subst(SCA[[j]],NA,0),my3CONT)
     PRSNO <- terra::app(SCA,mean)
     dat   <- ifelse(PRSNO[][,1]<0.5,0,1)
     df    <- cbind(df,dat)
   }
   colnames(df)<-paste0("f",DEC)
   df <- data.frame(df)
   
   setwd(WD); ANO <- terra::rast("ANO30.WS3.tif")
   df <- data.frame(cbind(df,ANO=ANO[][,1]))
   
   df <- na.omit(df)
   dim(df)
   
   mydf <- reshape2::melt(df,measure.vars=1:9,value.name="PROB",variable.name="fSCA")
   mydf$PROB <- as.factor(mydf$PROB)
   
   boxplot(ANO~PROB+fSCA,data=mydf)
   
    # ggplot
     graphics.off(); windows(6,5)
     ggplot(data=mydf, aes(x=fSCA, y=ANO,fill=PROB)) + 
     # geom_violin(alpha=0.4,show.legend=F)+
     ylab("Anomaly of greenness") + xlab("fSCA")+
     # scale_x_discrete(labels=rev(c("FP", "KS", "VM", "CT","PA","CF","EN")))+
     # geom_boxplot(width=0.3,notch=T)+
     stat_summary(fun.data = f, geom="boxplot", 
                  position=position_dodge(0.5),
                  width=0.6,linewidth=0.7,show.legend=F,alpha=0.6)+
     # stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.2 )+
     # stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red")+
     # scale_fill_manual(values=colran4)+
     # scale_y_continuous(expand=c(0.02,0.01),limits=c(0,1400),breaks=seq(0,1500,200))+
     theme_bw()+
     theme(legend.position = "bottom",legend.direction = "horizontal",axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"))

   
# V. UNUSED SUPP FIGURES ----
  # SUPP. FIG. 4b DEGREE DAYS - fSCA model for 2013-2023 for WS3  ----
    # Prepare longitudinal data  ----
  setwd(WD)
  DLR <- read.csv("DLR.fSCA.csv")
  # SWH <- read.csv("SWH.fSCA.csv") # use DLR instead
  THE <- read.csv("THE.fSCA.csv") 
  IGN <- read.csv("IGN.fSCA.csv")
  
  df <- rbind(THE,IGN)
  
  # only THE
  df <- THE
  df$fSCA <- df$fSCA/100

 df$WS <- as.factor(df$WS)
 levels(df$WS) <- c("3.LAUZ","1.MAND","2.RNPI","4.TOT")
 df$WS <- as.vector(df$WS)
 
 # remove points in September !
 df <- df[which(df$DOY<=245),]
 # remove a few outliers
 OUT <- which(df$fSCA>0.9 & df$DD>1000)
 if(length(OUT)>0) df <- df[-OUT,]
 
 df.l<- split(df,df$WS)
 
 w=4; graphics.off(); windows(10,10); plot(df.l[[w]]$DD,df.l[[w]]$fSCA)
 
    # Adjust log logistic curves using non linear quantile regression per WS ----
 # Three-parameter log logistic curves was fitted to the data. Rees and Long (1993) successfully demonstrated the utility of this model in the description of seed survival processes in several species.
  LL3 <- function(params,x){
      params[1]/(1+exp(params[2]*(log(x)-log(params[3]))))
    }
    
  DDseq <- 0:2000
    
  fitmodel <- fitmodel01 <- fitmodel05 <- fitmodel09 <- list()
  
  for (i in 1:length(WS)){
      x=df.l[[i]]$DD
      y=df.l[[i]]$fSCA
      
      PARAMini <- list(b=4,c=300)
      # bound a to 1
      fitmodel[[i]]  <- nls(y~1/(1 + exp(b*(log(x)-log(c)))), start=PARAMini)
      
      # using non linear quantile regression
fitmodel01[[i]] <- nlrq(y~1/(1 + exp(b*(log(x)-log(c)))), start=PARAMini, tau=0.1, trace=TRUE)
fitmodel05[[i]] <- nlrq(y~1/(1 + exp(b*(log(x)-log(c)))), start=PARAMini, tau=0.5, trace=TRUE)
fitmodel09[[i]] <- nlrq(y~1/(1 + exp(b*(log(x)-log(c)))), start=PARAMini, tau=0.9, trace=TRUE)
      
      # fitmodel  <- nls(y~a/(1 + exp(b*(log(x)-log(c)))), start=list(a=1,b=1.6,c=200))
    }
    
    # Simulate a curve for more snowy BV - requiring more heat to melt
    # plot(DDseq,LL3(c(0.97,3.1,383),Xseq),type="l")
    # lines(DDseq,LL3(c(1,3.1,483),Xseq),type="l",col="blue")
    # lines(DDseq,LL3(c(1,3.1,583),Xseq),type="l",col="magenta")
    
    df.LL3   <- data.frame(DD = DDseq,
                           MAND = LL3(c(1,summary(fitmodel[[1]])$coeff[,1]),DDseq),
                           RNPI = LL3(c(1,summary(fitmodel[[2]])$coeff[,1]),DDseq),
                           LAUZ = LL3(c(1,summary(fitmodel[[3]])$coeff[,1]),DDseq)
    )
    setwd(WD);write.csv(df.LL3,"df.LL3.csv",row.names=F)
    
    df01.LL3   <- data.frame(DD = DDseq,
                           MAND = LL3(c(1,summary(fitmodel01[[1]])$coeff[,1]),DDseq),
                           RNPI = LL3(c(1,summary(fitmodel01[[2]])$coeff[,1]),DDseq),
                           LAUZ = LL3(c(1,summary(fitmodel01[[3]])$coeff[,1]),DDseq)
    )
    df05.LL3   <- data.frame(DD = DDseq,
                             MAND = LL3(c(1,summary(fitmodel05[[1]])$coeff[,1]),DDseq),
                             RNPI = LL3(c(1,summary(fitmodel05[[2]])$coeff[,1]),DDseq),
                             LAUZ = LL3(c(1,summary(fitmodel05[[3]])$coeff[,1]),DDseq)
    )
    df09.LL3   <- data.frame(DD = DDseq,
                             MAND = LL3(c(1,summary(fitmodel09[[1]])$coeff[,1]),DDseq),
                             RNPI = LL3(c(1,summary(fitmodel09[[2]])$coeff[,1]),DDseq),
                             LAUZ = LL3(c(1,summary(fitmodel09[[3]])$coeff[,1]),DDseq)
    )
    
    # Evaluate model performance using m.eval ----
    EVAL <- FIT <- list()
    for (i in 1:length(WS)){
      OBS  <- df.l[[i]]$fSCA
      PRED <- LL3(c(1,summary(fitmodel[[i]])$coeff[,1]),df.l[[i]]$DD)
      EVAL[[i]] <- m.eval(OBS,PRED)
      EVAL[[i]] <- cbind(EVAL[[i]],t(summary(fitmodel[[i]])$coeff[,1]))
      FIT[[i]]  <- cbind(df.l[[i]],PRED=PRED)
      FIT[[i]]$RES <- PRED-OBS
    }
    EVAL.LL3 <- as.data.frame(do.call(rbind,EVAL))
    EVAL.LL3$WS <- WS 
    names(FIT)<- WS
    
    FIT.LL3 <- as.data.frame(do.call(rbind,FIT))
    boxplot(FIT.LL3$RES~FIT.LL3$YEAR)
    boxplot(FIT.LL3$RES~FIT.LL3$SRC,ylim=c(-0.1,0.15))  
    abline(h=0)
    
    plot(FIT.LL3$fSCA,FIT.LL3$PRED)

    setwd(WD);write.csv(EVAL.LL3,"EVAL.LL3.csv",row.names=F)
    setwd(WD);write.csv(FIT.LL3,"FIT.LL3.csv",row.names=F)
    
    # Main ggplot ----
    df.LL <- reshape2::melt(df.LL3,id.vars="DD",variable.name="WS",value.name="fSCA")
    df.LL$WS <- as.factor(df.LL$WS)
    levels(df.LL$WS)<- c("1.MAND","2.RNPI","3.LAUZ")
    df.LL$WS <- as.vector(df.LL$WS)
    
    df01.LL <- reshape2::melt(df01.LL3,id.vars="DD",variable.name="WS",value.name="fSCA")
    df01.LL$WS <- as.factor(df01.LL$WS)
    levels(df01.LL$WS)<- c("1.MAND","2.RNPI","3.LAUZ")
    df01.LL$WS <- as.vector(df01.LL$WS)
    
    df05.LL <- reshape2::melt(df05.LL3,id.vars="DD",variable.name="WS",value.name="fSCA")
    df05.LL$WS <- as.factor(df05.LL$WS)
    levels(df05.LL$WS)<- c("1.MAND","2.RNPI","3.LAUZ")
    df05.LL$WS <- as.vector(df05.LL$WS)
    
    df09.LL <- reshape2::melt(df09.LL3,id.vars="DD",variable.name="WS",value.name="fSCA")
    df09.LL$WS <- as.factor(df09.LL$WS)
    levels(df09.LL$WS)<- c("1.MAND","2.RNPI","3.LAUZ")
    df09.LL$WS <- as.vector(df09.LL$WS)
    
    # where params[1] is fitted initial germination (percentage), params[2] is the rate of viability loss in the rapidly declining section of the curve and t is the accumulated time in the controlled aging test (CAT);p50 is the time taken in aging days for viability to be reduced to 50%. 
  graphics.off();windows(10,5)
  ggplot()+
  geom_point(data=df,aes(x=DD,y=fSCA,color=WS),size=1,show.legend=F)+
    geom_line(data=df05.LL,aes(x=DD,y=fSCA,color=WS),linewidth=1.2,show.legend=F)+
    geom_line(data=df01.LL,aes(x=DD,y=fSCA,color=WS),linewidth=1.2,alpha=0.5,show.legend=F)+
    geom_line(data=df09.LL,aes(x=DD,y=fSCA,color=WS),linewidth=1.2,alpha=0.5,show.legend=F)+
    scale_color_manual(values=c(COL.WS,"black"))+
    geom_vline(xintercept=c(500,1000,1500),color="gray",linetype="dashed")+
    geom_hline(yintercept=c(0.1,0.5,0.9),linetype="dashed",color="gray")+
    xlab("Degree Days")+ylab("fSCA")+
    scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1),expand=c(0,0.01))+
    scale_x_continuous(limits=c(0,1600),breaks=seq(0,1500,500),expand=c(0,0.01))+
    facet_wrap2(~WS,strip=STRIP4,labeller=as_labeller(c("1.MAND"="Mandette","2.RNPI"="Roche Noire","3.LAUZ"="Lauzette","4.TOT"="Total area")))+
    theme_bw()+
    theme(legend.title=element_blank(),legend.box.background=element_rect(fill=NA),legend.key = element_rect(fill = NA),legend.text = element_text(size = 15),legend.position = c(0.5,0.5),axis.text=element_text(size=14,family="MS Reference Sans Serif"),axis.title=element_text(size=14,family="MS Reference Sans Serif"),strip.text = element_text(size = 14,family="MS Reference Sans Serif"),plot.margin = margin(0,20,10,10))
  
  
    # OLD STUFF - Other trials ----
   # an idea would be to set up a nonparametric Quantile Regression
   # see part 8 of https://cran.r-project.org/web/packages/quantreg/vignettes/rq.pdf
  
  # first need to performs the logistic transformation in Galarza et.al.(2020) (see references) for estimating quantiles for a bounded response. Once the response is transformed, it uses the best.lqr function.
  
  # quantile regression. we have a bounded response between 0 and 1
  # issue here would be better to use a sigmoid curve !
  # cf Robust Logistic Linear Quantile Regression 
  
  # test the nonlinear regression analysis of curves and model parameters using the drc package in R (Version 2.10, http://www.R-project.org). 
  
  # library(quantreg)    # quantile regression
  # df2  <- df[which(df$fSCA<0.95 & df$fSCA>0.05),]
  
  # test data
  x=df[which(df$WS=="3.LAUZ"),"DD"]
  y=df[which(df$WS=="3.LAUZ"),"fSCA"] 
  plot(x,y)
  LL3 <- function(x,a,b) y~1/(1 + exp(b*(log(x)-log(c))))
  Dat.nls  <- nls(y~1/(1 + exp(b*(log(x)-log(c)))), start=list(b=1.6,c=200))
  lines(1:2000, predict(Dat.nls, newdata=list(x=1:2000)), col="blue")
  # then fit the median using nlrq
  # the formul should ne in nls format
  Dat.nlrq.95 <- nlrq(y~1/(1 + exp(b*(log(x)-log(c)))), start=list(b=1.6,c=200), tau=0.9, trace=TRUE)
  Dat.nlrq.05 <- nlrq(y~1/(1 + exp(b*(log(x)-log(c)))), start=list(b=1.6,c=200), tau=0.1, trace=TRUE)
  lines(1:2000, predict(Dat.nlrq.95, newdata=list(x=1:2000)), col=2)
  lines(1:2000, predict(Dat.nlrq.05, newdata=list(x=1:2000)), col=2)
  Dat.nlrq.95

   # SUPP FIG. 2 Example of binary snow cover map from a SPOT4 Take 5 image ----
  # Native SPOT4-Take5 image on 5 june 2015 was downloaded at 
  # https://tpm-ds.eo.esa.int/smcat/SPOT4-5Take5_ESA/3/France/Ecrins/1/
  # https://tpm-ds.eo.esa.int/smcat/SPOT4-5Take5_ESA/2/18/18/
  tmp <- terra::rast("U:/DATA_SPOT/SPOT4TAKE5_2015/SRC/diskA/SPOT5_TAKE5/N1C_TUILE/SPOT5_HRG2_XS_20150605_N1_TUILE_EcrinsFranceD0000B0000/SPOT5_HRG2_XS_20150605_N1_TUILE_EcrinsFranceD0000B0000.TIF")
  tmp1 <- terra::crop(tmp,my3EXT)
  terra::plotRGB(tmp1,r=4,g=3,b=1,scale=822,stretch="hist") # same as the Dedieu paper !
  terra::plot(my3CONT,add=T,border="yellow")
  
  # pour le NDSI (panel c)
  NDSI <- (tmp1[[1]]-tmp1[[4]])/(tmp1[[1]]+tmp1[[4]])
  THR     <- 0.2 # instead of 0.4 in Dedieu & al. !
  m       <- c(-Inf, THR, 0,  THR, +Inf, 1)
  rclmat  <- matrix(m, ncol=3, byrow=TRUE)
  BSCM    <- terra::classify(NDSI,rclmat) 
  
  graphics.off(); windows(10,10); 
  terra::plot(BSCM,axes=F,col=c("white","white","steelblue2"),legend=F) # same as the Dedieu paper !
  terra::plot(my3CONT,add=T,border="black")
  
   # SUPP FIG. 3 FLUXALP data of SNOW COVER / DEGREE DAYS (20250123)  ----
  setwd(WD)
  DATA1     <- read.csv("FLUXALP.daily.csv")
  DATA2     <- read.csv("NIVOSE.daily.csv")
  SMOD.FLUX <- read.csv("SMOD.FLUX.csv")
  
  SEL1 <- na.omit(match(DATA2$DATE,DATA1$DATE))
  SEL2 <- na.omit(match(DATA1$DATE,DATA2$DATE))
  
  DATA1$SNOW_NIVOSE <- NA
  DATA1$TM_NIVOSE <- NA
  
  DATA1$SNOW_NIVOSE[SEL1] <- DATA2$SNOW[SEL2]/100 # attention depth in cm
  DATA1$TM_NIVOSE[SEL1]   <- DATA2$TM[SEL2]
  
  DATA <- DATA1
  hist(DATA2$SNOW)
  
  DATA$DATE <- as.Date(DATA$DATE,format="%Y-%m-%d")
  DATA$SNOW <- ifelse(DATA$NDVI>0,0,DATA$SNOW)
  SEL <- which(lubridate::year(as.Date(DATA$DATE))>=2013)
  DATA$SNOWsm <- data.table::frollmean(DATA$SNOW,n=5,fill=NA,align="center",na.rm=T)
  
  # add the cumulative degree days
  setwd(WD); met <- read.csv("FORCAGE_METEO.csv")[,c("DATE","FLUXALP.pred")]
  met$DATE <- as.Date(met$DATE)
  DATA <- DATA %>%left_join(met,by="DATE")
  
  DATA$DD <- ifelse(DATA$FLUXALP.pred<=0,0,DATA$FLUXALP.pred)
  DATA$YEAR <- lubridate::year(DATA$DATE)
  DATA$DOY <- lubridate::yday(DATA$DATE)
  YEARvec <- 2012:2023
  for (i in 2:nrow(DATA)){
    if(DATA$DOY[i]>DATA$DOY[i-1]) {
      DATA$DD[i]<- DATA$DD[i]+DATA$DD[i-1][]
    }else{
      DATA$DD[i]<-0
    }
  }
  DATA$DD[which(DATA$YEAR==2012)]<-NA
  
  # find the DD at FLUXALP.SMOD
  setwd(WD)
  SMOD.FLUX <- read.csv("SMOD.FLUX.csv")
  DD        <- read.csv("FLUXALP.DD.csv")
  SMOD.FLUX$DD <- NA
  
  
  for (i in 1:nrow(SMOD.FLUX)){
    tmpYEAR <- which(colnames(DD)==paste0("X",SMOD.FLUX$YEAR[i]))
    SMOD.FLUX$DD[i] <- DD[SMOD.FLUX$SMOD_esti[i],tmpYEAR]
  }
  SMOD.FLUX$DATE <- as.Date(paste0(SMOD.FLUX$YEAR,SMOD.FLUX$SMOD_esti),format="%Y%j")
  
  # set vector of dates
  DATEvec <- seq(from = as.Date("2012-10-01"), to = as.Date("2023-10-01"), by = "days")
  SMODvec <- as.Date(paste0(SMOD.FLUX$YEAR,SMOD.FLUX$SMOD_esti),format="%Y%j")
  COEFF   <- 400

  # mind the argument na.rm=T !!
  graphics.off(); windows(20,8)
  ggplot(data=DATA,aes(x=DATE,y=SNOWsm))+
    geom_area(fill=SMOD_pal(10)[10],na.rm=T)+
    geom_line(aes(x=DATE,y=DD/COEFF),na.rm=T,linewidth=1.5,color=SMOD_pal(10)[1])+
    geom_point(data=SMOD.FLUX,aes(x=DATE,y=DD/COEFF),size=3)+
    scale_y_continuous(limits=c(0,2.5),sec.axis = sec_axis(~.*COEFF, name="Cumulative Degree Days (°C)"),expand=c(0,0))+
    ylab("Snow Depth (m)")+xlab("")+
    geom_hline(yintercept=0)+
    geom_hline(yintercept=Hmisc::smedian.hilow(SMOD.FLUX$DD)/COEFF,linetype=c("solid","dashed","dashed"))+
    geom_vline(xintercept=SMODvec,linetype=2)+
    scale_x_date(date_breaks="1 years",limits=c(min(DATEvec),max(DATEvec)),date_labels = "%Y",expand=c(0,0))+
    theme_bw()+
    theme(axis.text=element_text(size=14,family="MS Reference Sans Serif"),axis.title=element_text(size=14,family="MS Reference Sans Serif"),strip.text = element_text(size = 14,family="MS Reference Sans Serif"),plot.margin = margin(5,5,5,5))
  
  # compute a more sophisticated model of SMOD ?
  
  # OLD FIG. 1c HYPSOMETRY ---- 
  setwd(WD); myDEM <- terra::rast("myDEM.tif"); myHILL <- terra::rast("myHILL.tif")
  myDAH <- terra::rast("myDAH.tif")
  tmp1 <- terra::mask(myDEM,myCONT); tmp2 <- terra::mask(myDAH,myCONT)
  
  tmp10 <- terra::aggregate(tmp1,fact=10)
  tmp20 <- terra::aggregate(tmp2,fact=10)
  df <- data.frame(ELEV=tmp10[],DAH=tmp20[])
  dfs <- df[sample(1:nrow(df),20000),]
  dfs <- dfs %>% drop_na()
  colnames(dfs) <- c("ELEV","DAH")
  
  graphics.off()
  windows(5,5)
  ggplot(dfs, aes(DAH, ELEV)) +
  ggdensity::geom_hdr(show.legend=T,probs = c(0.99, 0.95, 0.75, 0.5,0.25))+
  labs(x="DAH", y="Elevation")+
  geom_vline(xintercept=0)+
      theme_bw()+
    theme(legend.position = "bottom", legend.background = element_rect(fill = "white", colour = NA),axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 11,family="MS Reference Sans Serif"),legend.box.background = element_rect(color = "black"),legend.box.margin = margin(t = 1, l = 1),legend.direction="horizontal")
  
  
  # DEPRECATED Mapping of Greenness negative anomalies (DEPRECATED - see Fig. 1b) ----
   # late melting sites from SMOD.mod15sm (0,1)
   setwd(DATA_WD)
   SMOD.mod15  <- terra::rast(paste0("SCA15.SMODs.tif"))
   SMOD.mod15s <- round(terra::focal(SMOD.mod15,w=3,fun="mean",na.rm=T))
   SMOD.mod15d <- terra::disagg(SMOD.mod15s,fact=30,method="bilinear")
   THR     <- 140
   m       <- c(-Inf, THR, 0,  THR, +Inf, 1)
   rclmat  <- matrix(m, ncol=3, byrow=TRUE)
   LSM     <- terra::classify(SMOD.mod15d,rclmat) 
   all.equal(my3EXT,terra::ext(LSM))
   terra::plot(LSM,legend=F)
   setwd(WD); terra::writeRaster(LSM,"LSM.WS3.tif",overwrite=T)
   
   # flat areas (0,1)
   setwd("U:/DATA_ROCHENOIRE/LIDAR/LIDAR_IGN/")
   slo <- terra::mask(terra::rast("SLO_1M.tif"),my3CONT)
   slo <- terra::crop(slo,my3EXT)
   terra::plot(slo)
   range(slo[],na.rm=TRUE)
   setwd(WD); terra::writeRaster(slo,"FLAT.WS3.tif",overwrite=T)
   
   THR     <- 0.6
   m       <- c(-Inf, THR, 1,  THR, +Inf, 0)
   rclmat  <- matrix(m, ncol=3, byrow=TRUE)
   FLAT    <- terra::classify(slo, rclmat)
   terra::plot(FLAT,legend=F)
   terra::plot(FLAT*LSM,legend=F)
   
   # Anomaly of vegetation grenness
   setwd(DATA_WD);ano.WS3 <- terra::rast("ANO.WS3.tif")

   # apply a smoothing parameter and assemble
   # smo.WS3 <- terra::focal(ano.WS3,w=3,fun="modal",na.rm=T)
   # ANO <- terra::crop(smo.WS3,my3EXT)
   ANO <- terra::crop(ano.WS3,my3EXT)
   graphics.off();windows(10,10);terra::plot(ANO)
   all.equal(my3EXT,terra::ext(ANO))
   
   ANO30 <- terra::aggregate(ANO,fact=30)
   setwd(WD); terra::writeRaster(ANO30,"ANO30.WS3.tif",overwrite=T)
   
   # Greenness negative anomalies
   THR     <- 0
   m       <- c(-Inf,THR, 1,THR,+Inf,0)
   rclmat  <- matrix(m, ncol=3, byrow=TRUE)
   # x <- terra::subst(myrast, NA, 999)
   # x[x[]<=3] <- 1; x[x[]>3]<- NA
   ANO.neg <- terra::classify(ANO,rclmat) 
   terra::plot(ANO.neg,legend=F,col=c("lightgray","blue"))    
   
   # assemble
   myrast  <- FLAT*LSM*ANO
   terra::plot(myrast,legend=F) 
   
   myrast.neg  <- FLAT*LSM*ANO.neg
   
   # final plot
   graphics.off();windows(8,8)
   terra::plot(myrast.neg,legend=F,axes=F,col=c("gray70","tan4")) 
   plot(my3CONT,add=T,border="white",lwd=2)
   terra::plot(OBS.PTS[2],col="black",add=T)
   terra::sbar(xy=c(970100,6442800),d=1000, type="bar", below="km", label=c(0,0.5,1), cex=.9,col="yellow")
   # terra::north(xy=c(971000,6448200),type=2,col="yellow")
   terra::plot(myROUT,add=T,lwd=2,col="yellow")
   terra::plot(myHYDRO,add=T,col="skyblue3")

  # SUPP Example of a DOY - fSCA model for 2015 (20250123) ----
  myYEAR <- 13:23
  WS     <- c("RNPI","LAUZ","MAND")
  
  j=3 # 2015
  tmpPARAM <- matrix(NA,3,3)
  tmpfSCA  <- list()
  for(i in 1:3){
    setwd(WD)
    tmpPARAM[i,]   <- as.numeric(read.csv(paste0("PARAM_",WS[i],".csv"))[j,2:4])
    tmpfSCA[[i]]   <- read.csv(paste0("REC.fSCA",myYEAR[j],"_",WS[i],".csv"))
  }
  
  XSEQ <- 50:250
  df <- data.frame(DOY = XSEQ,
                   fSCA.RNPI = sigmoid(as.numeric(tmpPARAM[1,1:3]),XSEQ),
                   fSCA.LAUZ = sigmoid(as.numeric(tmpPARAM[2,1:3]),XSEQ),
                   fSCA.MAND = sigmoid(as.numeric(tmpPARAM[3,1:3]),XSEQ)
                           )
  
  dfs <- reshape2::melt(df,id.vars="DOY")
  
  for (i in 1:3) tmpfSCA[[i]]$WS <- paste0("fSCA.",WS[i])
  
  fSCA.thr <- matrix(NA,3,3)
  for (i in 1:3) {
    fSCA.thr[i,1]<- XSEQ[1]+which.min(abs(df[,i+1]-0.1))-1
    fSCA.thr[i,2]<- XSEQ[1]+which.min(abs(df[,i+1]-0.5))-1
    fSCA.thr[i,3]<- XSEQ[1]+which.min(abs(df[,i+1]-0.9))-1
  }
  
  df.fSCA <- do.call(rbind,tmpfSCA)
  
  graphics.off();windows(4,5)
  ggplot()+
    geom_line(data=dfs,aes(x=DOY,y=value,color=variable),linewidth=1.5)+
    geom_point(data=df.fSCA,aes(x=DOY,y=fSCA,color=WS),size=2)+
    geom_hline(yintercept=0.5,lty=2)+
    geom_segment(aes(x=fSCA.thr[1,3],xend=fSCA.thr[1,1],y=0.49,yend=0.49),linewidth=1.1,color=SMOD_pal(10)[c(7)])+
    geom_segment(aes(x=fSCA.thr[2,3],xend=fSCA.thr[2,1],y=0.5,yend=0.5),linewidth=1.1,color=SMOD_pal(10)[c(9)])+
    geom_segment(aes(x=fSCA.thr[3,3],xend=fSCA.thr[3,1],y=0.51,yend=0.51),linewidth=1.1,color=SMOD_pal(10)[c(3)])+
    scale_color_manual(values=SMOD_pal(10)[c(7,9,3)],labels=c("Roche Noire","Lauzette","Mandette"))+
    xlab("DOY")+ylab("fSCA")+
    scale_x_continuous(limits=c(50,225),expand=c(0,0))+
    scale_y_continuous(limits=c(0,1),expand=c(0.01,0.05))+
    theme_bw()+
    theme(legend.position = c(.24,.12), legend.background = element_rect(fill = "white", colour = NA),legend.title = element_blank(),axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 12,family="MS Reference Sans Serif"),legend.box.background = element_rect(color = "black"),legend.box.margin = margin(t = 1, l = 1),legend.direction="vertical")
  
  # SUPP TRENDS OF DOY for given fSCA and sigma SMOD for the three watersheds ----
    # define the accumulation of degree days for decades of reference 
  setwd(WD); DD <- read.csv("FLUXALP.DD.csv")
  setwd(WD); df.LL3 <- read.csv("df.LL3.csv")  
  
  YEARvec <- 1943:2023
  # Extract the average DD for moving windows MW of X YEAR
  MW = 10; YEARmed <- seq(1943+MW/2,2023-MW/2,1)
  df.DD <- NULL
  for (i in 1943:(2023-MW)){
    myPER <- i:(i+10)
    SEL   <- YEARvec%in%myPER
    TMP   <- apply(DD[,SEL],1,mean)
    df.DD <- cbind(df.DD,TMP)
  }
  colnames(df.DD)<-YEARmed
  df.DDlong <- reshape2::melt(df.DD)
  colnames(df.DDlong)<-c("DOY","YEAR","DD")
  
  # Compute the corresponding fSCA for a given DD-fSCA relationship
  df.fSCA.RNPI <- NULL
  for (i in 1:ncol(df.DD)){
    print(i)
    TMP     <- df.LL3[match(round(df.DD[,i]),df.LL3$DD),"RNPI"]
    df.fSCA.RNPI <- cbind(df.fSCA.RNPI,TMP)
  }

  df.fSCA.LAUZ <- NULL
  for (i in 1:ncol(df.DD)){
    TMP     <- df.LL3[match(round(df.DD[,i]),df.LL3$DD),"LAUZ"]
    df.fSCA.LAUZ <- cbind(df.fSCA.LAUZ,TMP)
  }
  
  df.fSCA.MAND <- NULL
  for (i in 1:ncol(df.DD)){
    TMP     <- df.LL3[match(round(df.DD[,i]),df.LL3$DD),"MAND"]
    df.fSCA.MAND <- cbind(df.fSCA.MAND,TMP)
  }
  
  # Estimate the DOY for critical values of fSCA
  fSCA.cv <- c(0.5,0.05,0.1,0.9,0.95)
  
  # for RNPI
  fSCA.THR.RNPI <- matrix(NA,ncol(df.fSCA.RNPI),5)
  for (i in 1:ncol(df.fSCA.RNPI)){
    for (j in 1:length(fSCA.cv)){
      fSCA.THR.RNPI[i,j]<- which.min(abs(df.fSCA.RNPI[,i]-fSCA.cv[j]))
    }
  }
  
  # put NA if superior to first of October -> no snowmelt
  fSCA.THR.RNPI[which(fSCA.THR.RNPI>270)] <- NA
  # transform into data.frame
  fSCA.THR.RNPI<- data.frame(fSCA.THR.RNPI)
  colnames(fSCA.THR.RNPI) <- paste0("fSCA",fSCA.cv)
  fSCA.THR.RNPI$YEAR <- YEARmed
  fSCA.THR.RNPI$DIFF    <- -fSCA.THR.RNPI$fSCA0.9+fSCA.THR.RNPI$fSCA0.1
  fSCA.THR.RNPI$DIFFext <- -fSCA.THR.RNPI$fSCA0.95+fSCA.THR.RNPI$fSCA0.05  
  
  # for LAUZ
  fSCA.THR.LAUZ <- matrix(NA,ncol(df.fSCA.LAUZ),5)
  for (i in 1:ncol(df.fSCA.LAUZ)){
    for (j in 1:length(fSCA.cv)){
      fSCA.THR.LAUZ[i,j]<- which.min(abs(df.fSCA.LAUZ[,i]-fSCA.cv[j]))
    }
  }
  
  # put NA if superior to first of October -> no snowmelt
  fSCA.THR.LAUZ[which(fSCA.THR.LAUZ>270)] <- NA
  # transform into data.frame
  fSCA.THR.LAUZ<- data.frame(fSCA.THR.LAUZ)
  colnames(fSCA.THR.LAUZ) <- paste0("fSCA",fSCA.cv)
  fSCA.THR.LAUZ$YEAR <- YEARmed
  fSCA.THR.LAUZ$DIFF    <- -fSCA.THR.LAUZ$fSCA0.9+fSCA.THR.LAUZ$fSCA0.1
  fSCA.THR.LAUZ$DIFFext <- -fSCA.THR.LAUZ$fSCA0.95+fSCA.THR.LAUZ$fSCA0.05
  
  # MAND
  fSCA.THR.MAND <- matrix(NA,ncol(df.fSCA.MAND),5)
  for (i in 1:ncol(df.fSCA.MAND)){
    for (j in 1:length(fSCA.cv)){
      fSCA.THR.MAND[i,j]<- which.min(abs(df.fSCA.MAND[,i]-fSCA.cv[j]))
    }
  }
  
  # put NA if superior to first of October -> no snowmelt
  fSCA.THR.MAND[which(fSCA.THR.MAND>270)] <- NA
  # transform into data.frame
  fSCA.THR.MAND<- data.frame(fSCA.THR.MAND)
  colnames(fSCA.THR.MAND) <- paste0("fSCA",fSCA.cv)
  fSCA.THR.MAND$YEAR <- YEARmed
  fSCA.THR.MAND$DIFF    <- -fSCA.THR.MAND$fSCA0.9+fSCA.THR.MAND$fSCA0.1
  fSCA.THR.MAND$DIFFext <- -fSCA.THR.MAND$fSCA0.95+fSCA.THR.MAND$fSCA0.05

  
  # ggplot
  DATA1 <- reshape2::melt(fSCA.THR.RNPI[,1:6],id.vars="YEAR")
  DATA1$WS <-"2.RNPI"
  DATA2 <- reshape2::melt(fSCA.THR.LAUZ[,1:6],id.vars="YEAR")
  DATA2$WS <-"3.LAUZ"
  DATA3 <- reshape2::melt(fSCA.THR.MAND[,1:6],id.vars="YEAR")
  DATA3$WS <-"1.MAND"
  DATA <- rbind(DATA1,DATA2,DATA3)
  
  DATA1 <- reshape2::melt(fSCA.THR.RNPI[,6:8],id.vars="YEAR")
  DATA1$WS <-"2.RNPI"
  DATA2 <- reshape2::melt(fSCA.THR.LAUZ[,6:8],id.vars="YEAR")
  DATA2$WS <-"3.LAUZ"
  DATA3 <- reshape2::melt(fSCA.THR.MAND[,6:8],id.vars="YEAR")
  DATA3$WS <-"1.MAND"
  DIFF <- rbind(DATA1,DATA2,DATA3)
  
  TMP1<-DATA[DATA$variable=="fSCA0.1",];colnames(TMP1)[3]<-"MAX"
  TMP2<-DATA[DATA$variable=="fSCA0.9",];colnames(TMP2)[3]<-"MIN"
  SIGMA <- cbind(TMP1[,c("YEAR","MAX","WS")],TMP2[,c("MIN")])
  colnames(SIGMA) <- c("YEAR","MAX","WS","MIN")
  
  TMP1<-DATA[DATA$variable=="fSCA0.05",];colnames(TMP1)[3]<-"MAX"
  TMP2<-DATA[DATA$variable=="fSCA0.95",];colnames(TMP2)[3]<-"MIN"
  SIGMAext <- cbind(TMP1[,c("YEAR","MAX","WS")],TMP2[,c("MIN")])
  colnames(SIGMAext) <- c("YEAR","MAX","WS","MIN")
  
  graphics.off();windows(8,4)
  ggplot()+
    geom_ribbon(data=SIGMAext,aes(x=YEAR,ymin=MIN,ymax=MAX), fill=gray(0.9), alpha=0.8) +
    geom_ribbon(data=SIGMA,aes(x=YEAR,ymin=MIN,ymax=MAX), fill=gray(0.6), alpha=0.8) +
    geom_line(data=DATA,aes(x=YEAR,y=value,linetype=variable,linewidth=variable))+
    scale_linetype_manual(values=c("solid","dashed","dashed","dashed","dashed","dashed"))+
    scale_linewidth_manual(values=c(1.5,.5,.5,.1,.1))+    
    scale_y_continuous(limits=c(85,275),expand=c(0.,0))+
    scale_x_continuous(breaks=seq(1950,2010,20),expand=c(0.05,0.05))+
    ylab("DOY at fSCA")+xlab("")+
    geom_hline(yintercept=c(225,270),linewidth=rep(c(1,1.5),3))+
    facet_wrap2(~WS, ncol=3,strip=STRIP,labeller=as_labeller(c("1.MAND"="Mandette","2.RNPI"="Roche Noire","3.LAUZ"="Lauzette")))+
    theme_bw()+
    theme(legend.title=element_blank(),legend.direction = "horizontal",legend.box.background=element_rect(fill=NA),legend.key = element_rect(fill = NA),legend.text = element_text(size = 10),legend.position = c(0.5,0.08),axis.text=element_text(size=14,family="MS Reference Sans Serif"),axis.title=element_text(size=14,family="MS Reference Sans Serif"),strip.text = element_text(size = 14,family="MS Reference Sans Serif"),plot.margin = margin(0,20,10,10))
  
  graphics.off();windows(8,3)
  ggplot()+
    geom_line(data=DIFF,aes(x=YEAR,y=value,linetype=variable,linewidth=variable))+
    scale_linetype_manual(values=c("solid","dashed"),labels=c("from fSCA0.9 to fSCA0.1","from fSCA0.95 to fSCA0.05"))+
    scale_linewidth_manual(values=c(1.5,.5),labels=c("from fSCA0.9 to fSCA0.1","from fSCA0.95 to fSCA0.05"))+    
    scale_y_continuous(limits=c(35,105),expand=c(0.,0))+
    scale_x_continuous(breaks=seq(1950,2010,20),expand=c(0.05,0.05))+
    ylab("Sigma SMOD")+xlab("")+
    facet_wrap2(~WS, ncol=3,strip=STRIP,labeller=as_labeller(c("1.MAND"="Mandette","2.RNPI"="Roche Noire","3.LAUZ"="Lauzette")))+
    theme_bw()+
    theme(legend.title=element_blank(),legend.direction = "horizontal",legend.box.background=element_rect(fill=NA),legend.key = element_rect(fill = NA),legend.text = element_text(size = 10),legend.position = c(0.65,0.13),axis.text=element_text(size=14,family="MS Reference Sans Serif"),axis.title=element_text(size=14,family="MS Reference Sans Serif"),strip.text = element_text(size = 14,family="MS Reference Sans Serif"),plot.margin = margin(0,20,10,10))
  

  # SUPP LogLOGISTIC curves for the 2013-2023 years ----
  setwd(WD); PARAM <- read.csv("PARAM.csv")
  myYEAR <- 13:23
  XSEQ   <- 50:250
  
  df.l <- list()
  for (i in 1:length(myYEAR)) df.l[[i]] <- data.frame(
    YEAR = 2000+rep(myYEAR[i],length(XSEQ)),
    DOY  = XSEQ,
    FIT  = sigmoid(as.numeric(PARAM[i,2:4]),XSEQ)
  )
  
  df  <- do.call(rbind,df.l)
  dfr <- do.call(cbind,df.l[c(1,10)])[,c(1,2,3,6)]
  colnames(dfr)[3:4] <- c("lower","upper")

graphics.off();windows(5,5)
ggplot()+
  geom_ribbon(data=dfr,aes(x=DOY, ymax=lower,ymin=upper), fill="gray", alpha=.5)+
  geom_line(data=dfr,aes(x=DOY,y=lower),linewidth=2,col="blue")+
  geom_line(data=dfr,aes(x=DOY,y=upper),linewidth=2,col="darkgreen")+
  geom_line(data=df,aes(x=DOY,y=FIT,group=YEAR))+
  scale_y_continuous(breaks=seq(0,1,0.25),expand = c(0.01,0.01))+
  scale_x_continuous(expand = c(0,0))+
  xlab("DOY")+ylab("fSCA")+
  theme_bw()+
  theme(legend.position = c(0.66,.06), legend.background = element_rect(fill = "white", colour = NA),axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 11,family="MS Reference Sans Serif"),legend.box.background = element_rect(color = "black"),legend.box.margin = margin(t = 1, l = 1),legend.direction="horizontal",plot.margin = margin(6,10,6,6))

  # SUPP BIPLOT fSCA05 x delta SMOD ----
  setwd(WD); PARAM <- read.csv("PARAM.csv")
  myYEAR <- 13:23

  graphics.off();windows(5,5)
  ggplot(data=PARAM,aes(x=c,y=DELTA))+
    geom_point(size=0.2)+
    geom_text(aes(label = YEAR),size=5)+
    xlab("fSCA0.5")+ylab("Delta SMOD")+
    geom_hline(yintercept=Hmisc::smedian.hilow(PARAM[,"DELTA"],0.75),linetype=c("solid","dashed","dashed"))+
        theme_bw()+
    theme(legend.position = c(0.66,.06), legend.background = element_rect(fill = "white", colour = NA),axis.text=element_text(size=12,family="MS Reference Sans Serif"),axis.title=element_text(size=12,family="MS Reference Sans Serif"),strip.text = element_text(size = 11,family="MS Reference Sans Serif"),legend.box.background = element_rect(color = "black"),legend.box.margin = margin(t = 1, l = 1),legend.direction="horizontal",plot.margin = margin(6,10,6,6))
  
  
  # SUPP NAIVE MODEL OF S2M SMOD (TO UPDATE ) ----
  setwd(WD); SMOD.CROC <- read.csv("SMOD.CROC.csv")
  
  SMOD.CROC$ANO_SnowDepth_JFM <- SMOD.CROC$SnowDepth_JFM - mean(SMOD.CROC$SnowDepth_JFM[33:52])
  SMOD.CROC$LPF_SnowDepth <- data.table::frollmean(SMOD.CROC$ANO_SnowDepth_JFM,n=5,align="center")  
  
  SMOD.CROC$ANO_Tmean_AMJ <- SMOD.CROC$Tmean_AMJ  - mean(SMOD.CROC$Tmean_AMJ[33:52])
  SMOD.CROC$LPF_Tmean <- data.table::frollmean(SMOD.CROC$ANO_Tmean_AMJ,n=5,align="center")
  
  SMOD.CROC$ANO_SMOD <- SMOD.CROC$SMOD  - mean(SMOD.CROC$SMOD[33:52])
  SMOD.CROC$LPF_SMOD <- data.table::frollmean(SMOD.CROC$ANO_SMOD,n=5,align="center")
  
  SMOD.CROC$RS <- "NO" ; SMOD.CROC$RS[55:65] <-"YES"
  
  df1 <- reshape2::melt(SMOD.CROC,id.vars=c(2,8),measure.vars=c(7,10,14))
  df2 <- reshape2::melt(SMOD.CROC,id.vars=c(2,8),measure.vars=c(11:13))
  
  df <- cbind(df1,df2$value)
  colnames(df)[4:5] <- c("ANO","LPF")
  
  graphics.off()
  windows(8,10)
  ggplot(data=df, aes(x=Year,y=ANO,fill=RS))+
    geom_col(position = "dodge",show.legend=F)+
    geom_line(aes(x=Year,y=LPF),linewidth=1.2)+
    labs(x="", y="Anomaly")+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=seq(1960,2020,10),col=gray(0.8),linewidth=0.1)+
    # scale_y_continuous(breaks=seq(50,250,50),limits=c(30,260))+
    scale_x_continuous(breaks=seq(1960,2020,10),labels=seq(1960,2020,10),expand = c(0.01,0.01))+
    facet_wrap(vars(variable),ncol=1,scales = "free_y",labeller=as_labeller(c("ANO_SnowDepth_JFM"="Mean Snow Depth (January-March) (m)","ANO_Tmean_AMJ"="Mean Air Temperature (April-June) (°C)","ANO_SMOD"="SMOD")))+
    theme_bw()+
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
    theme(axis.text=element_text(size=14,family="MS Reference Sans Serif"),axis.title=element_text(size=14,family="MS Reference Sans Serif"),strip.text = element_text(size = 14,family="MS Reference Sans Serif"),plot.margin = margin(4,4,4,4))
  
  
  graphics.off()
  windows(6,5)
  par(mar=c(3,4,0.5,4))
  ggplot(data=SMOD.CROC, aes(x=Year,y=ANO_Tmean_AMJ,fill=RS))+
    geom_col(position = "dodge",show.legend=F)+
    geom_line(aes(x=Year,y=LPF_Tmean),linewidth=1.2,col=gray(0.3))+
    labs(x="", y="Anomaly of Mean Air Temperature\n(April-June) (°C)")+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=seq(1960,2020,10),col=gray(0.8),linewidth=0.1)+
    # scale_y_continuous(breaks=seq(50,250,50),limits=c(30,260))+
    scale_x_continuous(breaks=seq(1960,2020,10),labels=seq(1960,2020,10),expand = c(0.01,0.01))+
    theme_bw()+
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
    theme(axis.text=element_text(size=14,family="MS Reference Sans Serif"),axis.title=element_text(size=14,family="MS Reference Sans Serif"),strip.text = element_text(size = 14,family="MS Reference Sans Serif"),plot.margin = margin(4,4,4,4))
  
  
  # old style figures
  plot(SMOD.CROC$Year,SMOD.CROC$SnowDepth_JFM,type="h",col="darkblue",lwd=5,ylab="JFM Mean Snow Depth (m)",xlab="year")
  abline(h=Hmisc::smedian.hilow(SMOD.CROC$SnowDepth_JFM,0.75),lty=c(1,2,2),col="blue")
  
  par(new=TRUE)
  plot(SMOD.CROC$Year,SMOD.CROC$SMOD,pch=21,bg="lightgreen",type="b",axes=F,xlab="",ylab="")
  axis(side = 4, at = pretty(range(SMOD.CROC$SMOD)))      # Add second axis
  mtext("SMOD", side = 4, line = 3)   
  abline(h=Hmisc::smedian.hilow(SMOD.CROC$SMOD,0.75),lty=c(1,2,2),col="darkgreen")  
  
  graphics.off()
  windows(10,10)
  plot(SMOD.CROC$SnowDepth_JFM,SMOD.CROC$SMOD,pch=21,cex=0.1,bg="lightgray",xlab="JFM Mean Snow Depth (m)",ylab="Average SMOD (1950-2250)")
  abline(LM <- lm(SMOD~SnowDepth_JFM,data=SMOD.CROC))
  text(SMOD.CROC$SnowDepth_JFM,SMOD.CROC$SMOD,labels=1959:2023,cex=1.1)
  
  plot(1959:2023,residuals(LM),type="h",xlab="",ylab="Residuals")
  lines(1959:2023,data.table::frollmean(residuals(LM), n=5,align="center"),lwd=2)
  abline(h=0)
  
  SMOD.CROC[order(SMOD.CROC$SMOD),]
  
  LM <- lm(SMOD~SnowDepth_JFM+Tmean_AMJ,data=SMOD.CROC)
  plot(1959:2023,residuals(LM),type="h",xlab="",ylab="Residuals")
  lines(1959:2023,data.table::frollmean(residuals(LM), n=5,align="center"),lwd=2)
  abline(h=0)
  
  plot(SMOD.CROC$SMOD,fitted(LM),xlim=c(100,170),ylim=c(100,170),pch=21,bg="gray")
  abline(0,1)
  abline(-10,1);abline(10,1)
  text(SMOD.CROC$SMOD,fitted(LM),1959:2023)
  # year typology
  graphics.off()
  windows(10,3)
  par(mar=c(4,3,0,0))
  # A <- boxplot(SMOD.CROC$SMOD,horizontal = T)
  boxplot(SMOD.CROC$SMOD,horizontal = T,add=F,col="gray",boxwex=1,xlab="SMOD (DOY)",ylim=c(100,165))
  boxplot(SMOD.CROC$SMOD[33:62],horizontal = T,add=T,col="white",boxwex=0.5) 
  snowyY <- SMOD.CROC$Year%in%c(2013,2018,2019,2021)
  dryY   <- SMOD.CROC$Year%in%c(2014,2015,2020,2022)
  mesicY <- SMOD.CROC$Year%in%c(2017,2016,2023) 
  vsnowyY <- SMOD.CROC$Year%in%c(1984,1991,1995)  # very snowy years in the Landsat period
  vdryY <- SMOD.CROC$Year%in%c(2007,2011)  # very snowy years in the Landsat period
  
  abline(v=SMOD.CROC[snowyY,"SMOD"],col="blue",lwd=2)
  abline(v=SMOD.CROC[mesicY,"SMOD"],col="green",lwd=2)    
  abline(v=SMOD.CROC[dryY,"SMOD"],col="red",lwd=4)
  abline(v=SMOD.CROC[vsnowyY,"SMOD"],lwd=2,lty=3,col="blue")  
  abline(v=SMOD.CROC[vdryY,"SMOD"],lwd=2,lty=3,col="red")  
  
  text(SMOD.CROC[55:65,"SMOD"]-0.5,1,labels=2013:2023,srt=90)
  text(SMOD.CROC[vsnowyY,"SMOD"]-0.5,1,labels=c(1984,1991,1995),srt=90)  
  text(SMOD.CROC[vdryY,"SMOD"]-0.5,1,labels=c(2007,2011),srt=90)  

################################################################################
# FROM PREVIOUS THOUGHTS ----
# snow - nosnow binary classification (20181130) ----
  # Objective : snow - nosnow binary classification of Persistent Snow Spring using late-season RGB aerial photographs
  # the aerial photograph campaign of 11/08/2013 (snowy year) in the Hautes-Alpes seems particularly suitable for this
  # A review of existing approaches is made in Fedorov, R., A. Camerada, P. Fraternali, and M. Tagliasacchi. 2016. Estimating Snow Cover From Publicly Available Images. Ieee Transactions on Multimedia 18:1187-1200.
  # seems that a thresholding technique based on the blue channel would be the most efficient way
  
  # Pilot study
  # UPDATE : la BD ortho est sur 
  DIR_ORTHO_IGN <- "K:/macroeco/GIS_DATA/Alpes/IGN/BD-ORTHO/D005/1_DONNEES_LIVRAISON_2010-10-00066/BDO_RVB_0M50_ECW_LAMB93_D05-ED03"  
  
  # En local
  MY_BDORTHO <- "C:/Users/cholerp/Documents/SIG/ROCHE_NOIRE/MY_BD_ORTHO/"
  LF         <- list.files(MY_BDORTHO,full=T)
  INFILE     <- LF[1]
  
  # not to run
  shell(cmd=paste('gdalinfo',INFILE,sep=' '))  
  
  
  BRIGHT <- stackApply(TEST, indices=c(1,1,1), fun=sum)
  BLUE   <- raster(TEST,3)
  # canal bleu normalis? par la brightness
  BLUEN <- BLUE/BRIGHT 
  plot(BLUEN,zlim=c(0.2,0.4))
  
  graphics.off();windows(10,10);dh <- hist(BLUEN[],breaks=seq(0,1,0.001))
  
  plot(b1)
  plot(raster(TEST,3))
  windows(10,10); plotRGB(TEST)
  
  graphics.off();windows(10,10);dh <- hist(raster(TEST,3)[],breaks=seq(0,255,1))
  ins  <- dh[["density"]]
  ss  <- which(diff(sign(diff(ins)))==2)+1
  abline(v=dh[["mids"]][ss],col="red")
  
  THR     <- 178
  m       <- c(-Inf, THR, 0,  THR, +Inf, 1)
  rclmat  <- matrix(m, ncol=3, byrow=TRUE)
  exg3    <- reclassify(raster(TEST,3), rclmat)
  windows(10,10); plot(exg3)
  windows(10,10); plot(exg3*BRIGHT)
  
  THR     <- 550
  m       <- c(-Inf, THR, 0,  THR, +Inf, 1)
  rclmat  <- matrix(m, ncol=3, byrow=TRUE)
  exg1    <- reclassify(BRIGHT, rclmat)
  windows(10,10); plot(exg1)   
  
  
  library(dplyr)
  mutate(ins, local.minima = if_else(lag(x) > x & lead(x) > x, TRUE, FALSE))
  
  nbins <- length(ins)
  ss <- which(rank(ins)%in%seq(from=nbins-2,to=nbins))
  dh[["mids"]][ss]
  abline(v=dh[["mids"]][ss])
  # the blue component corresponding to snow covered areas is always = 127 (accrding to Salvatori 2011)
  
  
  LF      <- list.files(DIR_ORTHO_IGN,pattern=".ecw")
  LF_FULL <- list.files(DIR_ORTHO_IGN,pattern=".ecw",full.names=T)
  
  TMP     <- match(DALLES_RNED,LF)
  

  # FIGURE 1 : SEASONAL NDVI SPOT x FLORISTIC CLUSTER ----
  graphics.off(); windows(25,15); 
  par(mfrow=c(2,3),mar=c(2,3,0,0),mgp=c(1.8,0.5,0))
  boxplot(FSFD~CLUST[,NC],col=COL)
  mtext(side=3,"FirstSnowFreeDay",line=-1.5) 
  boxplot(NDVI04~CLUST[,NC],ylim=c(0,1),col=COL);abline(h=c(0.2,0.4,0.6),lty=2)
  mtext(side=3,"SPOT NDVI April,11",line=-1.5)
  boxplot(NDVI05~CLUST[,NC],ylim=c(0,1),col=COL);abline(h=c(0.2,0.4,0.6),lty=2)
  mtext(side=3,"SPOT NDVI May, 11",line=-1.5)
  boxplot(NDVI06~CLUST[,NC],ylim=c(0,1),col=COL);abline(h=c(0.2,0.4,0.6),lty=2)
  mtext(side=3,"SPOT NDVI June, 5",line=-1.5)
  boxplot(NDVI07~CLUST[,NC],ylim=c(0,1),col=COL);abline(h=c(0.2,0.4,0.6),lty=2)
  mtext(side=3,"SPOT NDVI End of July",line=-1.5)
  boxplot(NDVI08~CLUST[,NC],ylim=c(0,1),col=COL);abline(h=c(0.2,0.4,0.6),lty=2)
  mtext(side=3,"SPOT NDVI End of August",line=-1.5)
  
  # FIGURE 2 : LEGEND of FIGURE 1 ----
  graphics.off(); windows(10,10); 
  plot(1,1,type="n",axes=F,xlab="",ylab="")
  legend("center",pt.bg=COL,pt.cex=4,pch=22,bty="n",legend=LEGEND,title="Vegetation Types",x.intersp=3,y.intersp=1.5)
  
  # FIGURE 2. NUMBER OF RELEVES SOPHIE & CBNA
  NC   <- 12
  REP  <- matrix(NA,2,NC); IDclu <- list()
  for (i in 1:NC){
    IDclu[[i]] <- which(CLUST[,NC]==i)
    N          <- length(IDclu[[i]])
    Nsop       <- length(which(substr(rownames(CLUST)[IDclu[[i]]],1,6)=="SOPHIE"))
    REP[2,i]   <- Nsop
    REP[1,i]   <- N-Nsop
  }
  
  graphics.off();windows(10,10);barplot(REP,names=1:NC,col=c("black","gray"),ylab="Number of relev?s")
  legend("top",pch=22,pt.bg=c("black","gray"),bty="n",legend=c("CBNA","SOPHIE"))
  
  # FIGURE 3 : FSFD x NDVI SPOT JULY focusing on snowbeds ----
  graphics.off();windows(15,10)
  par(mfrow=c(1,2),mar=c(2,3,0,0),mgp=c(1.8,0.5,0))
  plot(NDVI06,NDVI07,pch=21,bg="lightgray",cex=0.6)
  points(NDVI06[c(IDclu[[6]],IDclu[[7]])],NDVI07[c(IDclu[[6]],IDclu[[7]])],pch=21,bg="blue",cex=1.5)
  plot(NDVI06,NDVI07,pch=21,bg="lightgray",cex=0.8)
  points(NDVI06[IDclu[[12]]],NDVI07[IDclu[[12]]],pch=21,bg="lightblue",cex=1.5)
  
  graphics.off();windows(15,10)
  par(mfrow=c(1,2),mar=c(2,3,3,0),mgp=c(1.8,0.5,0))
  hist(NDVI06[c(IDclu[[6]],IDclu[[7]])],ylim=c(0,150),breaks=seq(-0.1,1,0.025),col='blue',main="Nival ouvert (cluster 6 & 7)")
  hist(NDVI06[IDclu[[12]]],ylim=c(0,150),breaks=seq(-0.1,1,0.025),col='lightblue',main="Nival dense (cluster 12)")
  
  # There are a number of sites (39) for which the NDVI06 values are already high in June
  # issue in cluster 12 defined too broadly
  ID.CLU12.ANO <- which(NDVI06>=0.2 & CLUST[,NC]==12) # 97 valeurs si 0.2 !
  
  # FIGURE 4 : SEASONAL NDSI x FLORISTIC CLUSTERS ----
  graphics.off(); windows(25,15); 
  par(mfrow=c(2,3),mar=c(2,3,0,0),mgp=c(1.8,0.5,0))
  boxplot(FSFD~CLUST[,NC],col=COL)
  mtext(side=3,"FirstSnowFreeDay",line=-1.5) 
  boxplot(NDSI04~CLUST[,NC],ylim=c(-1,1),col=COL);abline(h=c(0),lty=2)
  mtext(side=3,"SPOT NDSI April",line=-1.5)
  boxplot(NDSI05~CLUST[,NC],ylim=c(-1,1),col=COL);abline(h=c(0),lty=2)
  mtext(side=3,"SPOT NDSI Mid-May",line=-1.5)
  boxplot(NDSI06~CLUST[,NC],ylim=c(-1,1),col=COL);abline(h=c(0),lty=2)
  mtext(side=3,"SPOT NDSI Mid June",line=-1.5)
  boxplot(NDSI07~CLUST[,NC],ylim=c(-1,1),col=COL);abline(h=c(0),lty=2)
  mtext(side=3,"SPOT NDSI End of July",line=-1.5)
  boxplot(NDSI08~CLUST[,NC],ylim=c(-1,1),col=COL);abline(h=c(0),lty=2)
  mtext(side=3,"SPOT NDSI End of August",line=-1.5)
  
  
  # FIGURE 5 : PREDICTED MAPS for grassy snowbeds over SPOT and RON ----
  # FIGURE 5A over the SPOT4TAKE5 extent
  graphics.off();windows(14,4);par(mfrow=c(1,3),mar=c(0,0,0,0))
  plot(COND1,col=c("gray","blue"),axes=F,legend=F); 
  plot(RON.EXT,add=T); plot(PNE,add=T); plot(PNV,add=T)
  mtext(side=3,line=-1,text="Cond1 : NDVI June,5 < 0.2")
  plot(COND2,col=c("gray","blue"),axes=F,legend=F)
  plot(RON.EXT,add=T); plot(PNE,add=T); plot(PNV,add=T)
  mtext(side=3,line=-1,text="Cond2 : NDVI Mid-August > 0.5")
  plot(SNOWBET,col=c("gray","blue"),axes=F,legend=F)
  plot(RON.EXT,add=T); plot(PNE,add=T); plot(PNV,add=T)
  mtext(side=3,line=-1,text="Snowbet : Cond2 x Cond2")
  
  # FIGURE 5B over the RON extent  
  graphics.off();windows(14,4);par(mfrow=c(1,3),mar=c(0,0,0,0))
  plot(crop(COND1,RON.EXT),col=c("gray","blue"),axes=F,legend=F); 
  plot(ARA.EXT.400,add=T); plot(PNE,add=T); plot(PNV,add=T)
  mtext(side=3,line=-1.5,text="Cond1 : NDVI June,5 < 0.2")
  plot(crop(COND2,RON.EXT),col=c("gray","blue"),axes=F,legend=F); 
  plot(ARA.EXT.400,add=T); plot(PNE,add=T); plot(PNV,add=T)
  mtext(side=3,line=-1.5,text="Cond2 : NDVI Mid-August > 0.5")
  plot(crop(SNOWBET,RON.EXT),col=c("gray","blue"),axes=F,legend=F); 
  plot(ARA.EXT.400,add=T); plot(PNE,add=T); plot(PNV,add=T)
  mtext(side=3,line=-1.5,text="Snowbet : Cond2 x Cond2")
  
  # Condition #3 the utility of using NDSI is not demonstrated
  
  # TABLE  1 : MATCH SNOWBED PLOTS and MAPS ----
  
  SNOWBET <- raster("SNOWBET_005_05.tif")
  TMP     <- extract(SNOWBET,SIT_TOT[,3:4])
  TAB     <- table(CLUST[,NC],TMP)
  round(100*TAB[,2]/(TAB[,2]+TAB[,1]))
  
  TMP <- extract(SNOWBET,SIT_TOT[-ID.CLU12.ANO,3:4])
  TAB <- table(CLUST[-ID.CLU12.ANO,NC],TMP)
  round(100*TAB[,2]/(TAB[,2]+TAB[,1]))
  
  # bcp de faux n?gatifs avec CLU12 original, much better if removing ID.CLU12.ANO
  
  SNOWBET.AGG <- list(); SNOWBET.AGG[[1]] <- SNOWBET
  for (i in 2:25) {
    print(i)
    SNOWBET.AGG[[i]]  <- aggregate(SNOWBET, fact=i, fun=max)
  }
  
  SEQ <- seq(10,250,10)
  BET <- matrix(NA,length(SEQ),12)
  for (i in 1:length(SEQ)) {
    print(i)
    TMP     <- extract(SNOWBET.AGG[[i]],SIT_TOT[,3:4])
    TAB     <- table(CLUST[,NC],TMP)
    BET[i,] <- round(100*TAB[,2]/(TAB[,2]+TAB[,1]))
  }
  apply(BET,2,mean)
  
  BET12 <- matrix(NA,length(SEQ),5)
  table(PARTSNOW[[1]][,5])
  for (i in 1:length(SEQ)) {
    print(i)
    TMP      <- extract(SNOWBET.AGG[[i]],SIT_TOT[ID.CLUSNO,3:4])
    TAB       <- table(PARTSNOW[[1]][,5],TMP)
    BET12[i,] <- round(100*TAB[,2]/(TAB[,2]+TAB[,1]))
  }
  
  SEQ <- seq(10,250,10)
  BET2 <- matrix(NA,length(SEQ),12)
  for (i in 1:length(SEQ)) {
    TMP     <- extract(SNOWBET.AGG[[i]],SIT_TOT[-ID.CLU12.ANO,3:4])
    TAB     <- table(CLUST[-ID.CLU12.ANO,NC],TMP)
    BET2[i,] <- round(100*TAB[,2]/(TAB[,2]+TAB[,1]))
  }
  apply(BET2,2,mean)
  
  plot(BET[,12],BET2[,12]);abline(0,1)
  
  graphics.off();windows(10,10)
  plot(SEQ,BET[,2],type="b",ylim=c(0,80),pch=21,bg="lightblue",xlab="Resolution (m)",ylab="Percent of correct matches",cex=3)
  # points(SEQ,BET[,12],type="b",ylim=c(0,80),pch=21,bg="lightblue",xlab="Resolution (m)",ylab="Percent of correct matches",cex=1.5)
  points(SEQ,BET[,11],type="b",pch=21,bg="magenta")
  points(SEQ,BET[,10],type="b",pch=21,bg="brown")
  points(SEQ,apply(BET[,6:7],1,mean),type="b",pch=21,cex=2,bg="darkblue")
  points(SEQ,apply(BET[,2:3],1,mean),type="b",pch=21,bg="orange")
  points(SEQ,apply(BET[,c(1,5,8)],1,mean),type="b",pch=21,bg="darkgreen")
  points(SEQ,BET[,4],type="b",pch=21,bg="green",cex=2)
  abline(h=50,lty=2)
  abline(v=c(10,100,250)) # resolution of SPOT, LANDSAT & MODIS
  
  graphics.off();windows(10,10)
  plot(SEQ,apply(BET12,1,mean),type="b",ylim=c(0,100),pch=21,bg="lightblue",xlab="Resolution (m)",ylab="Percent of correct matches",cex=3)
  for (i in 1:5) points(SEQ,BET12[,i],type="b",pch=23,bg=rainbow(5)[i],cex=1)
  abline(h=c(50,90),lty=2)
  abline(v=c(10,100,250)) # resolution of SPOT, LANDSAT & MODIS
  
  # Extract for the different clusters of CLU12
  tmp <- extract(SNOWBET.AGG[[10]],SIT_TOT[ID.CLUSNO,3:4])
  table(tmp)  
  tmp <- extract(SNOWBET.AGG[[10]],SIT_TOT[ID.CLU121,3:4])
  table(tmp)
  tmp <- extract(SNOWBET.AGG[[10]],SIT_TOT[ID.CLU122,3:4])
  table(tmp)
  tmp <- extract(SNOWBET.AGG[[10]],SIT_TOT[ID.CLU123,3:4])
  table(tmp)
  tmp <- extract(SNOWBET.AGG[[10]],SIT_TOT[ID.CLU124,3:4]) # p?le nival
  table(tmp)
  
  # FIGURE 6 ----
  SNOWBET100 <- aggregate(SNOWBET,10,sum)
  plot(SNOWBET100)
  plot(PNE,add=T);plot(PNV,add=T)
  hist(SNOWBET100[SNOWBET100[]>=10],xlab="% de secteur PAMN dans la maille 100m")
  
  ALL   <- xyFromCell(SNOWBET100,which(SNOWBET100[]!=0))
  
  ID    <- which(SNOWBET100[]>10)
  SNO   <- xyFromCell(SNOWBET100,ID)
  
  TESU<-180*ptransform(SNO,LAMB93@projargs,WGS84@projargs)/pi
  
  PTsp<-SpatialPointsDataFrame( TESU[,1:2], 
                                data.frame(COVER=extract(SNOWBET100,ID)),
                                proj4string = WGS84
  )
  
  writeOGR(PTsp, "SNO100_10.kml","SNO100_10",driver="KML")
  
  
  
  save.image("SNOWBED_RS.RData")
  
# From SPOT to MODIS ----
  EXT  <- extent(830000,1080000,6210000,6610000) # 250 km x 400 km
  M250 <- raster(EXT, resolution=250, crs=LAMB93, vals=1)
  
  SPOT250 <- crop(M250,extent(SNOWBET)) # 324 rows and 364 col
  
  SNOWBET250 <- 100*raster::aggregate(SNOWBET,fact=25, fun=sum,na.rm=T)/625
  PERSEQ <- seq(20,90,10); MOD.SNOPER <-list()
  for(i in 1:8) MOD.SNOPER[[i]] <- which(SNOWBET250[]>=PERSEQ[i])
  unlist(lapply(MOD.SNOPER, function(x) length(x)))
  
  COO           <-  xyFromCell(SNOWBET250,MOD.SNOPER[[4]]) # above 50% - 81 pixels
  SNOWBET_ID250 <- cellFromXY(M250,COO)
  PTS           <- (1/pi)*180*ptransform(COO,LAMB93,WGS84@projargs)[,1:2]
  DF            <- data.frame(ID250 = SNOWBET_ID250, 
                              COVSNOW = SNOWBET250[SNOWBET250>=50], 
                              PTS)
  PTsp<-SpatialPointsDataFrame( DF[,3:4], 
                                data.frame(A=DF[,1:2]),
                                proj4string = WGS84
  )
  
  writeOGR(PTsp, "SNOWBET.kml","SNOWBET",driver="KML")
  graphics.off(); plot(PTS)
  
  # Extract MODIS yearly metrics for SNOWBET sites
  SRC.DIR         <- "C:/Users/cholerp/Documents/SIG/MOD09Q1.006.h18v04/NDVI_RES250/"
  SNOWBET_NDVI250 <-  gfsm(DF[,1],SRC.DIR,RETURN="ALL")
  
  # Compute inter-annual variation for DIAGNOSTIC ANALYSIS
  METRICSsno <- SNOWBET_NDVI250$METRICS
  Nsno       <- length(SNOWBET_ID250)
  ONSETsno   <- matrix(unlist(lapply(METRICSsno, function(x) x[,"ONSET"])),Nsno,18,byrow=T)
  OFFSETsno  <- matrix(unlist(lapply(METRICSsno, function(x) x[,"OFFSET"])),Nsno,18,byrow=T) 
  OFFSET2sno <- matrix(unlist(lapply(METRICSsno, function(x) x[,"OFFSET2"])),Nsno,18,byrow=T) 
  INT02sno   <- matrix(unlist(lapply(METRICSsno, function(x) x[,"NDVIint02"])),Nsno,18,byrow=T)
  INT04sno   <- matrix(unlist(lapply(METRICSsno, function(x) x[,"NDVIint04"])),Nsno,18,byrow=T)
  PEAKsno    <- matrix(unlist(lapply(METRICSsno, function(x) x[,"TNDVImax"])),Nsno,18,byrow=T)  
  PLAT90sno  <- matrix(unlist(lapply(METRICSsno, function(x) x[,"NDVIplateau90"])),Nsno,18,byrow=T)
  PLAT80sno  <- matrix(unlist(lapply(METRICSsno, function(x) x[,"NDVIplateau80"])),Nsno,18,byrow=T)
  MINsno     <- matrix(unlist(lapply(METRICSsno, function(x) x[,"NDVImin"])),Nsno,18,byrow=T) 
  MAXsno     <- matrix(unlist(lapply(METRICSsno, function(x) x[,"NDVImax"])),Nsno,18,byrow=T) 
  MAXMsno    <- matrix(unlist(lapply(METRICSsno, function(x) x[,"NDVImaxm"])),Nsno,18,byrow=T) 
  GRDURsno   <- PEAKsno - ONSETsno        # Duration of growing period
  SEN30sno   <- matrix(unlist(lapply(METRICSsno, function(x) x[,"SEN30"])),Nsno,18,byrow=T) 
  SEN45sno   <- matrix(unlist(lapply(METRICSsno, function(x) x[,"SEN45"])),Nsno,18,byrow=T) 
  REG30sno   <- matrix(unlist(lapply(METRICSsno, function(x) x[,"REG30"])),Nsno,18,byrow=T)
  REG45sno   <- matrix(unlist(lapply(METRICSsno, function(x) x[,"REG45"])),Nsno,18,byrow=T)
  
  # Criteria #1 : NDVImaxm <  0.5 ? supprimer (DEL1) - > NO PIXELS, FINE
  graphics.off();windows(15,20);par(mar=c(3,4,0.5,0.5),mfrow=c(2,1))
  tmp <- apply(MAXsno,1,max); range(tmp)
  DEL1sno <- which(tmp<0.5 | tmp>1.1, arr.ind=T) # fine
  
  # Criteria #2 : amplitude of response too low (DEL2) -> NO PIXELS, FINE
  graphics.off();windows(15,20);par(mar=c(3,4,0.5,0.5),mfrow=c(2,1))
  tmp <- apply(MAXsno-MINsno,1,mean)
  hist(tmp,breaks=seq(0,1.3,0.01),col="darkgreen",main="",xlab="Ampltitude du verdissement")
  LOWAMP <- which(tmp<0.4)
  
  # Criteria #3 : average ONSETsno TOO EARLY -> 11 PIXELS
  ONSETm <- apply(ONSETsno,1,mean)
  hist(ONSETm)
  TOOEAR <- which(ONSETm<=150)
  
  # Criteria #4 : GRDUR misestimated (several years with GRDURsno>100)
  hist(GRDURsno,breaks=c(0,400,5))
  tmp <- which(GRDURsno>=100,arr.ind=T)
  apply(GRDURsno,1,mean)
  TOOLON <- as.numeric(names(which(table(tmp[,1])>1)))
  
  # Merging all criteria -> DELSNOW with 11 pixels, reste 70 pixels
  DELsno <- unique(c(LOWAMP,TOOEAR,TOOLON))
  
  (Nsnof <- Nsno - length(DELsno))
  
  # FIGURE 6 : MODIS phenology for snowbed pixels ----
  graphics.off()
  windows(15,10)
  par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(4,6,3,0))
  YLIM=c(0,120)
  hist(PEAKsno[-DELsno,],ylim=YLIM, breaks=seq(120,340,2),add=F,col="gray",main="",xlab="Jour Julien") # TNDVImax
  hist(OFFSETsno[-DELsno,],breaks=seq(120,340,2),col=colors()[38],add=T,main="")  
  hist(ONSETsno[-DELsno,], breaks=seq(120,340,2),add=T,xlim=c(0,370),ylim=c(0,50000),col=colors()[86],main="")
  abline(v=365)
  abline(v=median(ONSETsno[-DELsno,]),col="darkgreen",lwd=2) # 171 - 12 juin - conforme ? SPOT !
  # abline(v=median(PEAKsno[-DELsno,]),col="black",lwd=2) # 208 
  abline(v=median(OFFSETsno[-DELsno,]),col="red",lwd=2) # 294 
  
  legend("top", pch=22,cex=1,pt.cex=2,pt.bg=c(colors()[86],"gray",colors()[38]), c("Date de d?but de saison","Date du pic de v?g?tation","Date de fin de saison")) 
  
  # FIGURE 7 : INTERANNUAL VARIABILITY OF MODIS PHENOLOGY ----
  YEAR <- 2000:2017
  SCA<-0.45
  graphics.off()
  windows(SCA*20, SCA*13)
  par(mar=c(5,5,1,1),cex.lab=1.5,mgp=c(3.5,1,0))
  YLIM <- c(140,320)
  
  boxplot(ONSETsno[-DELsno,],ylim=YLIM,col=colors()[86],outline=F,names=YEAR,las=2)
  abline(h=median(apply(ONSETsno[-DELsno,],2,median,na.rm=T),na.rm=T),col=c('darkgreen'),lwd=3) 
  
  boxplot(PEAKsno[-DELsno,],add=T,col="gray",xaxt="n",yaxt="n",outline=F)
  abline(h=median(apply(PEAKsno[-DELsno,],2,median,na.rm=T),na.rm=T),col=c('darkgray'),lwd=3) 
  
  boxplot(OFFSETsno[-DELsno,],add=T,col=colors()[38],xaxt="n",yaxt="n",outline=F,ylab="Jour de l'ann?e (Jour julien)",xlab="Ann?e")
  abline(h=median(apply(OFFSETsno[-DELsno,],2,median,na.rm=T),na.rm=T),col=c('brown'),lwd=3)
  
  # FIGURE 8 : ONSET x TNDVImax - GROWTH PHASE ----
  graphics.off()
  windows(SCA*20, SCA*20)
  plot(ONSETsno[-DELsno,],PEAKsno[-DELsno,],ylab='TNDVImax', xlab='Onset',bg='black',pch=21,cex=0.6,ylim=c(150,260))
  abline(0,1,lty=2)
  abline(22,1,lty=2,col="blue")
  abline(LM <- lm(c(PEAKsno[-DELsno,])~c(ONSETsno[-DELsno,])))
  mtext(side=3,line=-1,text="y = 67.8 + 0.82x, R2 = 0.39")
  
  # Using of Teil-Sen median slope instead of lm does not give stronger results
  library(mblm)
  
  SLO <- INT <- SLOts <- INTts <- RSQ <- RSQts <- rep(NA,nrow(ONSETsno))
  for (i in 1:nrow(ONSETsno)){
    print(i)
    Y <- PEAKsno[i,]; X <- ONSETsno[i,]
    LRL <- mblm(Y ~ X)
    INTts[i] <- coef(LRL)[1]; SLOts[i] <- coef(LRL)[2]; RSQts[i] <- anova(LRL)$Sum[1]/sum(anova(LRL)$Sum)
    LM <- lm(Y~X)
    INT[i] <- coef(LM)[1]; SLO[i] <- coef(LM)[2]; RSQ[i] <- summary(LM)$r.squared
  }
  
  graphics.off()
  windows(20,15)
  par(mfcol=c(3,2),oma=c(2,2,2,2),mar=c(3,6,3,0))
  
  hist(SLO[-DELsno],breaks=seq(0,2,0.1),xlim=c(0,2)) # a majority of sites with slope<1 (=rattrappage ph?no)
  abline(v=median(SLO[-DELsno],na.rm=T))
  abline(v=1,lwd=2,col='red')
  
  hist(INT[-DELsno],breaks=seq(-100,300,10));abline(v=median(INT),lwd=2) # intercept median = 52j
  hist(INT[which(RSQ>=0.4)],breaks=seq(-100,300,10),add=T,col="green");abline(v=median(INT[which(RSQ>=0.4)]),lwd=2) # intercept median = 52j
  
  hist(RSQ[-DELsno],breaks=seq(0,1,0.1),xlim=c(0,1))
  
  hist(SLOts[-DELsno],breaks=seq(0,2,0.1),xlim=c(0,2)) # a majority of sites with slope<1 (=rattrappage ph?no)
  abline(v=median(SLOts[-DELsno],na.rm=T))
  abline(v=1,lwd=2,col='red')
  
  hist(INTts[-DELsno],breaks=seq(-100,300,10));abline(v=median(INT),lwd=2) # intercept median = 52j
  hist(INTts[which(RSQ>=0.4)],breaks=seq(-100,300,10),add=T,col="green");abline(v=median(INT[which(RSQ>=0.4)]),lwd=2) # intercept median = 52j
  
  
  # FIGURE 9 : NDVImax and NDVIint02 ---- 
  
  YEAR <- 2000:2017
  SCA<-0.45
  graphics.off()
  windows(SCA*20, SCA*20)
  par(mar=c(4,3,1,1),cex.lab=1.5,mgp=c(3.5,1,0),mfcol=c(3,1))
  
  boxplot(MAXMsno[-DELsno,],ylim=c(0.4,0.9),add=F,col="gray",outline=F,names=YEAR,las=2)
  abline(h=median(apply(MAXMsno[-DELsno,],2,median,na.rm=T),na.rm=T),col=c('black'),lwd=2) 
  
  boxplot(INT02sno[-DELsno,],ylim=c(30,100),col=colors()[86],outline=F,names=YEAR,las=2)
  abline(h=median(apply(INT02sno[-DELsno,],2,median,na.rm=T),na.rm=T),col=c('black'),lwd=2)
  
  boxplot(PLAT90sno[-DELsno,],add=F,ylim=c(10,100),col="gray",outline=F,names=YEAR,las=2)
  abline(h=median(apply(PLAT90sno[-DELsno,],2,median,na.rm=T),na.rm=T),col=c('black'),lwd=2) 
  
  
  
  
  
# II. FLORISTIC DATA ----
  
  # in brief
  CLUST  <- read.csv("20180225_Dom_Flo_Clust.csv",header=T,row.names=1)
  # ACHTUNG ! need to reorder PARTa after the randomization
  CLUST  <- CLUST[order(rownames(CLUST)),]
  CLUST$SRC <- unlist(lapply(strsplit(rownames(CLUST),"_"), function(x) x[1]))
  ID.CLUSNO <- which(CLUST[,12]==2)
  SIT_TOT  <- read.csv("20180225_CBNA_SOPHIE_SIT.csv")
  
  # II.A. DATASET extracted for SOPHIE - SPOT4TAKE5 (20180225) ----
  # see 
  # M69.r     I.A.
  # alpages.r V.A. CBNA DB for SPOT4TAK5 data set (20180209)
  
  setwd("C:/Users/cholerp/Documents/PROJETS/SNOWBED/DATA_ANALYSIS/")  
  
  SIT_TOT <- read.csv("20180225_SNOWBED_SIT.csv") # row names not relevant in the data set
  REL_TOT <- read.csv("20180225_SNOWBED_REL.csv",row.names=1) # Attention ? colnames " " -> "."
  SPL_TOT <- read.csv("20180225_SNOWBED_SPL.csv",stringsAsFactors = F)[,-1] # Attention ? colnames " " -> "."
  
  FLO(REL_TOT,1) ; (dim(REL_TOT.db)) # 3105 rel x 487 species 
  REL_TOTf <- REL_TOT[-REL_TOT.i[[4]],]
  SIT_TOTf <- SIT_TOT[-REL_TOT.i[[4]],]
  SPL_TOTf <- SPL_TOT[-REL_TOT.i[[3]],]
  
  # NOT TO RUN
  write.csv(REL_TOTf,"20180225_CBNA_SOPHIE_REL.csv",row.names=T)
  write.csv(SIT_TOTf,"20180225_CBNA_SOPHIE_SIT.csv",row.names=F)
  write.csv(SPL_TOTf,"20180225_CBNA_SOPHIE_SPL.csv",row.names=T)    
  
  # Checking
  sort(REL_TOT[500,],decreasing=T)[1:10] # Checking  is fine
  
  # Estimating distances between relev?s
  SIT_TOTf$REL.NEAR <-   SIT_TOTf$DIS.NEAR <- SIT_TOTf$DAT.NEAR <- NA
  for (i in 1:nrow(SIT_TOTf)){
    XI  <- SIT_TOTf$X[i]
    YI  <- SIT_TOTf$Y[i]
    D   <- sqrt((XI-SIT_TOTf$X)^2+(YI-SIT_TOTf$Y)^2)
    tmp <- order(D)
    SIT_TOTf$DIS.NEAR[i]  <- D[tmp[2]]
    SIT_TOTf$REL.NEAR[i]  <- as.vector(SIT_TOTf[tmp[2],"Id_station"])
    SIT_TOTf$DAT.NEAR[i]  <- SIT_TOTf[tmp[2],"Date"]
  }
  
  hist(SIT_TOTf$DIS.NEAR,breaks=seq(0,5000,5),xlim=c(0,500))
  CUT <- cut(SIT_TOTf$DIS.NEAR,breaks=seq(0,50,5))
  plot(CUT)
  
  write.csv(SIT_TOTf,"20180225_CBNA_SOPHIE_SIT_NEAR.csv",row.names=F)
  
  # II.B. FLORISTIC CLUSTERS (20180225) ----
  nclu <- 15 # number of clusters
  NR   <- nrow(REL_TOT.db)
  
  # any influence of relev? order in the process ?
  TRY  <- REL_TOT.db[sample(1:NR,NR),]
  PART <- part(TRY,nclu)
  
  TRY   <- REL_TOT.da[sample(1:NR,NR),]
  PARTa <- part(TRY,nclu,typedist="jaccard")
  
  # Writing outputs
  CURDAY <- gsub("-","",Sys.Date())
  write.csv(PARTa[[1]]      , paste(CURDAY,"Dom_Flo_Clust.csv",sep="_"))
  write.csv(PARTa[[3]][[12]], paste(CURDAY,"Dom_Flo_Part_12.csv",sep="_")) # " " becomes .
  write.csv(PARTa[[3]][[15]], paste(CURDAY,"Dom_Flo_Part_15.csv",sep="_"))
  write.csv(PARTa[[4]][[12]], paste(CURDAY,"Indic_Spec_12.csv",sep="_")) # " " becomes .
  write.csv(PARTa[[4]][[15]], paste(CURDAY,"Indic_Spec_15.csv",sep="_")) # " " becomes .
  write.csv(PARTa[[5]][[12]], paste(CURDAY,"SignAssoc_Spec_12.csv",sep="_")) 
  write.csv(PARTa[[5]][[15]], paste(CURDAY,"SignAssoc_Spec_15.csv",sep="_"))
  
  # II.C. SNOWBED CLUSTERS (20180225) ----
  
  INDIC     <- read.csv("20180225_Indic_Spec_12.csv")
  SP.NC12   <- read.csv("20180225_Dom_Flo_Part_12.csv")
  CLUST     <- read.csv("20180225_Dom_Flo_Clust.csv",header=T,row.names=1)
  # ACHTUNG ! need to reorder PARTa after the randomization
  CLUST     <- CLUST[order(rownames(CLUST)),]
  CLUST$SRC <- unlist(lapply(strsplit(rownames(CLUST),"_"), function(x) x[1]))
  
  REL_TOT  <- read.csv("20180225_CBNA_SOPHIE_REL.csv")
  SIT_TOT  <- read.csv("20180225_CBNA_SOPHIE_SIT_NEAR.csv")
  SPL_TOT  <- read.csv("20180225_CBNA_SOPHIE_SPL.csv")[,-1]
  
  INDIC[300,];SPL_TOT[300,]
  INDIC[400,];SPL_TOT[400,]
  
  INDIC$"CD_REF" <- SPL_TOT$"CD_REF"
  
  # Check presence of snowbed species within clusters
  ID.SNOSP  <- c(88515,120057,81140,123176,100515,89986)
  SNOSP     <- SPL_TOT[match(ID.SNOSP,SPL_TOT$"CD_REF"),] # list of snowbed specialists
  INDIC.SNO <- INDIC[match(ID.SNOSP,SPL_TOT$"CD_REF"),]
  apply(INDIC.SNO[,2:13],1,which.max)
  
  # Richness of snowbed relev?s
  NC <- 12
  graphics.off()
  hist(REL_TOT.i[[2]][which(CLUST[,NC]==2)],col="lightblue",breaks=seq(0,50,2))  # rather high 
  # hist(REL_TOT.i[[2]][which(CLUST[,NC]==6)],col="blue",add=T,breaks=seq(0,60,2))  # rather high comp?red to what I know from Aravo and Roche Noire
  # hist(REL_TOT.i[[2]][which(CLUST[,NC]==7)],col="darkblue",add=T,breaks=seq(0,60,2))  
  
  # LOOK FOR RESURVEYS per CLUSTER
  SIT_TOT$CLUST <- CLUST[,NC]
  SIT_TOT$CLUST_NEAR <- SIT_TOT$CLUST[match(SIT_TOT$REL.NEAR,SIT_TOT$Id_station)]
  
  write.csv(SIT_TOT,"20180225_CBNA_SOPHIE_SIT_NEAR.csv",row.names=F)
  
  # Identification des Relev?s appari?s
  PAIR          <- SIT_TOT[which(SIT_TOT$DIS.NEAR<=30),]
  hist(abs(PAIR$Date-PAIR$DAT.NEAR))
  PAIR$DIFF.DAT <- abs(PAIR$Date-PAIR$DAT.NEAR)
  graphics.off();windows(10,10);plot(PAIR$DIS.NEAR,PAIR$DIFF.DAT,pch=21,bg="gray")
  
  table(PAIR$CLUST,PAIR$CLUST_NEAR)
  
  length(which(PAIR$DIFF.DAT>=20)) # 57 relev?s si 10 ans
  
  # Distribution of relev?s
  graphics.off();windows(7.5,10);
  # plot(EXTspot,add=F,col="red")
  plot(SIT_TOT[,3:4],pch="+")
  points(SIT_TOT[ID.CLU12,3:4],pch=21,bg="lightblue",cex=1.5)
  plot(PNE,add=T);plot(PNV,add=T)
  # plot(EXTspot,add=T,col="red")
  
  save.image("SNOWBED_RS.RData")
  
  # CONCLUSION on RUN20180214 (3106 rel x 488 sp., abundance data)  
  # 1. Three clusters (6,7 and 12) with nival conditions when 12 prescribed groups. Cluster 12 with grassy snowbeds all snowbed specialists in the 10 most frequent species and there are other cluster with stony, open nival systems
  # 2. The number of snowbed clusters is constant between 11 and 15 clusters with around 150-200 relev?s. 12-15 clusters is a good compromise to get a high number of relev?s  while keeping a strong snowbed identity  
  # 3. NDVI phenology will suggest that the snowbed cluster is too broadly defined see part IV
  
  
  
  # II.D. DIVERSITY patterns in snowbed releves (20180221) ----
  
  REL_TOT  <- read.csv("20180225_CBNA_SOPHIE_REL.csv",row.names = 1) 
  REL_SNO  <- REL_TOT[which(CLUST[,NC]==2),] # 245 relev?s
  
  FLO(REL_SNO,0.5)
  dim(REL_SNO.db) # 244 x 179 species
  dim(REL_SNO.da) # 244 x 179 species / 82 species if 5%
  
  # Partitionning snowbed relev?s
  nsnow     <- 5 # number of clusters
  PARTSNOW  <- part(REL_SNO.da,nsnow,NSPEC=10,typedist="jaccard")
  # PARTSNOW  <- part(REL_SNO.db,nsnow,NSPEC=20)
  table(PARTSNOW[[1]][,5])
  dim(PARTSNOW[[5]][[5]])
  SP <- rownames(PARTSNOW[[5]][[5]])
  SP[which(PARTSNOW[[5]][[5]][,6]==3)]
  write.csv(PARTSNOW[[5]][[5]],"PARTSNOW5.csv")
  
  
  ID.CLU121<- match(names(which(PARTSNOW[[1]][,5]==1)),rownames(REL_TOT))
  ID.CLU122<- match(names(which(PARTSNOW[[1]][,5]==2)),rownames(REL_TOT))
  ID.CLU123<- match(names(which(PARTSNOW[[1]][,5]==3)),rownames(REL_TOT))
  ID.CLU124<- match(names(which(PARTSNOW[[1]][,5]==4)),rownames(REL_TOT))
  ID.CLU125<- match(names(which(PARTSNOW[[1]][,5]==5)),rownames(REL_TOT))
  
  # Extract for the different clusters of CLU12
  tmp <- extract(SNOWBET.AGG[[25]],SIT_TOT[ID.CLU12,3:4])
  table(tmp)  
  tmp <- extract(SNOWBET.AGG[[25]],SIT_TOT[ID.CLU121,3:4])
  table(tmp)
  tmp <- extract(SNOWBET.AGG[[25]],SIT_TOT[ID.CLU122,3:4])
  table(tmp)
  tmp <- extract(SNOWBET.AGG[[25]],SIT_TOT[ID.CLU123,3:4])
  table(tmp)
  tmp <- extract(SNOWBET.AGG[[25]],SIT_TOT[ID.CLU124,3:4]) # p?le nival
  table(tmp)
  tmp <- extract(SNOWBET.AGG[[25]],SIT_TOT[ID.CLU125,3:4]) # p?le nival
  table(tmp)
  
  
  # Non Metric Muldidimensionnal Scaling (vegan) - approche par groupe d'exp?ces
  library(vegan)
  mdloto<-metaMDS(vegdist(t(REL_SNO.da), method='bray'), k=4)
  
  X1=1;X2=2
  graphics.off();windows(10,10);plot(mdloto$points[,X1],mdloto$points[,X2],pch=".",cex=4)
  
  MATCH1<-na.omit(match("CAREX.FOETIDA.ALL...1785",rownames(mdloto$points)))
  points(mdloto$points[MATCH1,X1],mdloto$points[MATCH1,X2],pch='.',cex=10,col='blue')
  MATCH2<-na.omit(match("SALIX.HERBACEA.L...1753",rownames(mdloto$points)))
  points(mdloto$points[MATCH2,X1],mdloto$points[MATCH2,X2],pch='.',cex=10,col='blue')
  MATCH3<-na.omit(match("ALCHEMILLA.PENTAPHYLLEA.L...1753",rownames(mdloto$points)))
  points(mdloto$points[MATCH3,X1],mdloto$points[MATCH3,X2],pch='.',cex=10,col='green')
  MATCH4<-na.omit(match("SIBBALDIA.PROCUMBENS.L...1753",rownames(mdloto$points)))
  points(mdloto$points[MATCH4,X1],mdloto$points[MATCH4,X2],pch='.',cex=10,col='red')
  MATCH5<-na.omit(match("GNAPHALIUM.SUPINUM.L...1768",rownames(mdloto$points)))
  points(mdloto$points[MATCH5,X1],mdloto$points[MATCH5,X2],pch='.',cex=10,col='magenta')
  
  # axe 1 esp?ces est un axe thermique - ?boulis pierrier
  sort(mdloto$points[,2])[1:20]
  sort(mdloto$points[,2],decreasing=T)[1:20]
  
  # Non Metric Muldidimensionnal Scaling (vegan) - approche par groupe de relev?s
  mdloto<-metaMDS(vegdist(REL_SNO.da, method='bray'), k=4)
  # mdloto<-metaMDS(vegdist(REL_SNO.db, method='jaccard'), k=4)
  
  X1=1;X2=2
  graphics.off();windows(10,10);plot(mdloto$points[,X1],mdloto$points[,X2],pch=".",cex=6,col=PARTSNOW[[2]][,5],xlab="NMDS 1", ylab="NMDS 2")
  # text(mdloto$points[,X1],mdloto$points[,X2],labels=PARTSNOW[[1]][,5])
  PARTSNOW[[3]][[4]]
  
  SNO100 <- extract(SNOWBET.AGG[[10]],SIT_TOT[ID.CLU12,3:4])
  
  plot(mdloto$points[,X1],mdloto$points[,X2],pch=".",cex=0.1,col=PARTSNOW[[2]][,5],xlab="NMDS 1", ylab="NMDS 2")
  text(mdloto$points[,X1],mdloto$points[,X2],labels=SNO100)
  
  
  # Conclusion en regardant les extremes de la distribution dans le plan esp?ce
  # tendance m?sophile acidiphile avec C. curvula, Euphrasia minima, Carex semp. etc
  # tendance neutrophile frais nitrophile avec Trifolium badium/thalii, Ranunculus kuepferi
  # tendance ?boulis avec Cryptogramma, Doronicum, Epilobium
  
  # CONCLUSION : 
  
  save.image("SNOWBED_RS.RData")
  
  # II.E. STAGE MARGAUX (not updated) ----
  # (beware of the transformation of initial csv file)
  
  rm(list=ls())
  FDATA <- read.table("especes_presence_absence_PC.csv",header=T,sep=";",row.names=1)
  
  # rename species to get shorter names -> SPEC
  SPEC <-unlist(lapply(strsplit(colnames(FDATA),split="_"),function(x) paste(substr(x[1],1,4),substr(x[2],1,4),sep="_")))
  
  graphics.off()
  par(mar=c(5,5,1,0))
  barplot(sort(apply(FDATA,2,sum),decreasing=T)[1:49]/129,las=2,cex.names=0.8,ylab="% of occurrence",ylim=c(0,0.5)) 
  # 49 species with at least two occurences
  
  # CBNA DB Importation des donn?es floristiques combe ? neige
  DIR.FLO <- "K:/fyse/ALPAGES_SENTINELLES/COMBES_A_NEIGE_CBNA_RAW_DATA/"
  HAB.files <- list.files(DIR.FLO,pattern="_HAB.txt",full=T)
  
  count.fields(HAB[1],sep="\t",quote="")
  scan(HAB[1],what="integer",sep="\t")
  
  HAB <- NULL
  for (i in 1:7){
    print(i)
    TMP <- read.table(HAB.files[i],header=T,sep="\t",quote="")
    HAB <- rbind(HAB,TMP)
  }
  
  dim(HAB)
  
  PREC <- which(as.numeric(HAB$"lmprecision.l931")<=10)
  
  NIV <- extract(NDSI,HAB[PREC,30:31])
  NDV <- extract(NDVI,HAB[PREC,30:31])
  
  plot(HAB[,30],HAB[,31])
  
  graphics.off()
  hist(NDSI[],freq=F,ylim=c(0,5),breaks=seq(-1,1,0.02))
  hist(NIV,freq=F,add=T,col='blue',breaks=seq(-1,1,0.02))
  
  
  graphics.off()
  hist(NDVI[],freq=F,ylim=c(0,3),breaks=seq(-1,1,0.02))
  hist(NDV,freq=F,add=T,col='blue',breaks=seq(-1,1,0.02))
  
  NIVele <- extract(DEM_30,HAB[,30:31])
  
  # Ordination using NMDS - Non Metric Muldidimensionnal Scaling
  mdloto    <- metaMDS(FDATA, distance ='bray', k=2) # to get species scores 
  
  # Display  of sites (black) and species (red)
  graphics.off()
  par(mar=c(4,4,1,1))
  plot(mdloto,display=c("sites","species"),pch="")
  text(mdloto$points,rownames(FDATA),cex=0.6)
  text(mdloto$species,SPEC,col='red',cex=0.6)
  
  # I.C. Clustering of releves & indicator species analysis
  nclu <- 4 # number of clusters
  PART <- part(FDATA,nclu)
  
  # assocaitaed species to cluster 1
  row.names(ISAtmp1)[which(ISAtmp1[,1]>=0.2)] # poa supina
  row.names(ISAtmp1)[which(ISAtmp1[,2]>=0.2)] # nardion
  row.names(ISAtmp1)[which(ISAtmp1[,3]>=0.2)] # alop alpi
  row.names(ISAtmp1)[which(ISAtmp1[,4]>=0.2)] # salix retu reti
  
  
  
# VI. MODIS values on CBNA grid (not updated) ----
  # III.A. SELECTED SITES ----
  
  # III.B. CBNA grid cells at 100m ----
  
  CELLS <- read.csv("centroids.csv")[,2:3]
  CELLS <- read.csv("centroids_last.csv")[,4:5]
  
  
  # reference grid CBNA - see MODIS_ARAVO for details
  
  # les donn?es MODIS updated sont sur le disque externe WD6
  SRC.DIR <- "H:/MODIS.ODYSSEE/MOD09Q1.005.ODYSSEE/NDVI_RES100/"
  FILES   <- list.files(SRC.DIR,full=T)
  test    <- raster(FILES[1])
  
  PTS     <- cellFromXY(test,CELLS)
  
  TEST    <- gfsm(PTS=PTS,SRC.DIR=SRC.DIR)
  
  write.csv(cbind(PTS,CELLS,t(TEST$TSNOWmelt)),"TSNOWmelt.csv")
  write.csv(cbind(PTS,CELLS,t(TEST$NDVImax)),"NDVImax.csv")
  
  
# VII. LONG-TERM GREENING (LANDSAT & MODIS) & SNOWCOVER DURATION ----
  # VII.1 LANDSAT data (LRR.L) ----
  # combined analysis using Landsat & MODIS
  # the browning trend for very long snow-cover trend needs to be investigated. pb of partial pixel snowcoverage within original Landsat scenes ?
  CAN.LRR <- cbind(extract(raster("LRR_P1P3.tif"),CAN_COO),
                        extract(raster("LRR_P1P2.tif"),CAN_COO),
                        extract(raster("LRR_P2P3.tif"),CAN_COO)
                        )
  
  boxplot(CAN.LRR);abline(h=0)
  t
  LRR_P1P3 <- projectRaster(raster("LRR_P1P3.tif"),REF25)*MASKfsfd*MASKaoi
  LRR_P1P2 <- projectRaster(raster("LRR_P1P2.tif"),REF25)*MASKfsfd*MASKaoi
  LRR_P2P3 <- projectRaster(raster("LRR_P2P3.tif"),REF25)*MASKfsfd*MASKaoi
  
  LRR.L    <- stack(LRR_P1P2,LRR_P2P3,LRR_P1P3) 
  
  SMELT    <- FSFDaoi*MASKfsfd*MASKaoi # snowmelting date
  
  x <-  SMELT[]
  y <-  raster(LRR.L,3)[]
  MIS <- which(is.na(x) | is.na(y)) # remove all NA values
  x <- x[-MIS];  y <- y[-MIS]
  FSFDseq <- seq(90,190,10)
  
  xc <- cut(x,FSFDseq)
  boxplot(y~xc,las=2,outline=F,ylim=c(-0.8,0.8),ylab="LRR",xlab="FSFD")
  abline(h=0)
  
  
  # detailed analysis per period

  GREM <- matrix(NA,length(FSFDseq)-1,3)
  for(j in 1:3){
    x <-  SMELT[]
    y <-  raster(LRR,j)[]
    MIS <- which(is.na(x) | is.na(y)) # remove all NA values
    x <- x[-MIS];  y <- y[-MIS]
    
    for(i in 1:length(FSFDseq)-1){
      PSS     <- which(x>FSFDseq[i] & x<=FSFDseq[i+1])    
      GREM[i,j] <- median(y[PSS])
    }
  }
  
  plot(x, y, pch=".", col="grey55",xlab="FSFD", ylab="LRR",ylim=c(-1.5,1.5),xlim=c(90,200))
  points(seq(95,195,10),GREM[,1],pch=21,bg="black",type="b")
  points(seq(95,195,10),GREM[,2],pch=21,bg="red",type="b")
  points(seq(95,195,10),GREM[,3],pch=21,bg="blue",type="b")
  abline(h=0,lwd=1, lty="dashed")
  
  f1 <- kde2d(x, y, n = 200, lims = c(120,200, range(y)))
  contour(f1, add=T,  drawlabels=F)
  abline(h=0,lwd=1, lty="dashed")
  
  # identify hotspot of greening
  BREAKS <- seq(-2,2,0.01)
  data   <- na.omit(raster(LRR.L,3)[])
  H      <- hist(data,freq=F,breaks=seq(-2,2,0.01))
  
  QUA90 <- BREAKS[which(cumsum(H$density)>90)[1]] # 0.32
  abline(v=QUA90)
  
  GAU  <- fitdistr(data, "cauchy")[1]
  
  plot(function(x) dcauchy(x, location=GAU$estimate[1], scale=GAU$estimate[2]),-2,2,add=T,col="red")
  plot(function(x) dcauchy(x, location=GAU$estimate[1], scale=0.1),-2,2,add=T,col="red")
  
  QUAN <- qnorm(c(0.1,0.9), mean=GAU$estimate[1], sd=GAU$estimate[2], lower.tail = TRUE, log.p = FALSE)
  
  abline(v=QUAN,col="red")
  
  THR      <- QUA90
  m        <- c(-Inf, THR, 0,  THR, +Inf, 1)
  rclmat   <- matrix(m, ncol=3, byrow=TRUE)
  LRR90    <- reclassify(raster(LRR.L,3), rclmat)
  windows(10,10); plot(LRR90)

   save.image("SNOWBED_RS.RData") 
  
  # VII.2 MODIS data (mblm & LRR) ----
  
   YEAR        <- 2000:2018
   SRC.DIR     <- "K:/LECA/biogeochem/MODIS/MOD09Q1.006.ALPS/NDVI_RES250/"
   # SRC.DIR   <- "G:/REMOTE_SENSING/MOD09Q1.006.ALPS/NDVI_RES250/"
   LFmod09Q1   <- list.files(SRC.DIR,full=T) # 841 files
   
   # For CAN sites
   CAN_ID250   <- cellFromXY(raster(LFmod09Q1[100]),CAN_COO) # extract ID
   CAN_NDVI250 <- gfsm(PAR_ID250,SRC.DIR,RETURN="ALL")
   
   plot(CAN_NDVI250[[3]][,2]) # checking
   
   VAR <- lapply(CAN_NDVI250[[4]],function(x) colnames(x))[[1]]
   
   setwd("~/PROJETS/SNOWBED/PHENO/")
   
   for (i in 1:length(VAR)) {
     TMP  <- matrix(unlist(lapply(CAN_NDVI250[[4]], function(x) x[,VAR[i]])),length(CAN_ID250),19,byrow=T)
     colnames(TMP) <- 2000:2018
     rownames(TMP) <- CAN[,1]
     write.csv(TMP,paste(VAR[i],".csv",sep=""))   
   }
   
   setwd("~/PROJETS/SNOWBED/PHENO/") 
   LF     <- list.files(full=T)
   LFs    <- list.files(full=F)
   NAME   <- gsub(".csv","",LFs)
   MEDIAN_PIX <- CV_PIX <- matrix(NA,length(CAN_ID250),length(LFs))
   MEDIAN_YR  <- CV_YR <- matrix(NA,length(YEAR),length(LFs))
   colnames(CV_PIX) <-colnames(MEDIAN_PIX)<- NAME
   rownames(CV_PIX) <-rownames(MEDIAN_PIX)<- read.csv(LFs[1])[,1]
   colnames(CV_YR)  <-colnames(MEDIAN_YR) <- NAME
   rownames(CV_YR)  <-rownames(MEDIAN_YR) <- YEAR
   
   for (i in 1:length(LFs)){
     print(i)
     TMP           <- as.matrix(read.csv(LF[i],row.names = 1))
     MEDIAN_PIX[,i] <- apply(TMP,1,median,na.rm=T)
     CV_PIX[,i]     <- round(100*apply(TMP,1,sd,na.rm=T)/apply(TMP,1,mean,na.rm=T))
     MEDIAN_YR[,i] <- apply(TMP,2,median,na.rm=T)
     CV_YR[,i]     <- round(100*apply(TMP,2,sd,na.rm=T)/apply(TMP,2,mean,na.rm=T))
   }

  
  setwd("~/PROJETS/SNOWBED/PHENO/") 
  
  myPEAK   <- as.matrix(read.csv("TNDVImax.csv",row.names = 1))
  myOFFSET <- as.matrix(read.csv("OFFSET.csv",row.names = 1))
  myONSET  <- as.matrix(read.csv("ONSET.csv",row.names = 1))
  myNDVImax   <- read.csv("NDVImax.csv",row.names = 1)
  myNDVIint   <- read.csv("NDVIint02.csv",row.names = 1)
  myPLAT      <- read.csv("NDVIplateau80.csv",row.names = 1)
  
    # A. HISTOGRAMS of ONSET,OFFSET & PEAK ----
  graphics.off()
  windows(15,10)
  par(mfrow=c(2,1),oma=c(2,2,2,2),mar=c(3,6,3,0))
  YLIM=c(0,250); XLIM=c(0,365)
  hist(myPEAK,ylim=YLIM,xlim=XLIM, breaks=seq(1,468,2),add=F,col="gray",main="",xlab="Jour Julien") # TNDVImax
  hist(myOFFSET,breaks=seq(1,468,2),col=colors()[38],add=T,main="")  
  hist(myONSET, breaks=seq(1,468,2),add=T,xlim=c(0,370),ylim=c(0,50000),col=colors()[86],main="")
  abline(v=365)
  
  legend("topleft", pch=22,cex=1,pt.cex=2,pt.bg=c(colors()[86],"gray",colors()[38]), c("Onset","Peak","Offset")) 
  
  hist(myPEAK-myONSET,breaks=seq(-400,350,2),col="green",xlim=c(0,300),main="") # lenght of growth phase -> unimodal
  hist(myOFFSET-myPEAK,breaks=seq(-400,450,2),add=T,col="darkgreen") # lenght of growth phase -> unimodal
  
  legend("topright", pch=22,cex=1,pt.cex=2,pt.bg=c("green","darkgreen"), c("Length of growth phase","Length of senescence phase"))
  
    # B. BOXPLOTS of key dates ----
  
  SCA<-0.45
  graphics.off()
  windows(SCA*22, SCA*15)
  par(mar=c(5,5,1,1),cex.lab=1.5,mgp=c(3.5,1,0))
  YLIM <- c(40,400)
  
  boxplot(myONSET,ylim=YLIM,col=colors()[86],outline=F,names=YEAR,las=2)
  abline(h=median(apply(myONSET,2,median,na.rm=T),na.rm=T),col=c('darkgreen'),lwd=3) 
  
  boxplot(myPEAK,add=T,col="gray",xaxt="n",yaxt="n",outline=F)
  abline(h=median(apply(myPEAK,2,median,na.rm=T),na.rm=T),col=c('darkgray'),lwd=3) 
  
  boxplot(myOFFSET,add=T,col=colors()[38],xaxt="n",yaxt="n",outline=F,ylab="Jour de l'ann?e (Jour julien)",xlab="Ann?e")
  abline(h=median(apply(myOFFSET,2,median,na.rm=T),na.rm=T),col=c('brown'),lwd=3)
  
  legend("topleft", pch=22,cex=1,pt.cex=2,pt.bg=c(colors()[86],"gray",colors()[38]), c("Onset","Peak","Offset"))
  
  
    # C. BOXPLOTS of length of phases ----
  
  SCA<-0.45
  graphics.off()
  windows(SCA*22, SCA*15)
  par(mar=c(5,5,1,1),cex.lab=1.5,mgp=c(3.5,1,0))
  YLIM <- c(0,400)
  
  boxplot(myPEAK-myONSET,ylim=YLIM,col=colors()[86],outline=F,names=YEAR,las=2,ylab="Number of days")
  abline(h=median(apply(myPEAK-myONSET,2,median,na.rm=T),na.rm=T),col=c('darkgreen'),lwd=3)
  boxplot(-myONSET+myOFFSET,ylim=YLIM,col="gray",outline=F,names=YEAR,las=2,add=T)
  abline(h=median(apply(-myONSET+myOFFSET,2,median,na.rm=T),na.rm=T),col=c('darkgreen'),lwd=3)
  
  legend("topleft", pch=22,cex=1,pt.cex=2,pt.bg=c(colors()[86],"gray"), c("Length of growth phase","Length of senescence phase"))
  
    # D. BIVARIATE RELATIONSHIP ONSET & PEAK ----
  
  SLO <- RSQ <- rep(NA,nrow(myONSET))
  for (i in 1:nrow(myONSET)){
    print(i)
    Y <- myPEAK[i,]; X <- myONSET[i,]
    LM <- lm(Y~X)
    SLO[i] <- coef(LM)[2]
    RSQ[i] <- summary(LM)$r.squared
  }
  
  graphics.off()
  windows(15,8)
  par(mfrow=c(1,2),oma=c(0,0,0,0),mar=c(4,4,1,1))
  
  plot(as.matrix(myONSET),as.matrix(myPEAK),ylab='TNDVImax',xlab='Onset',pch=".",cex=0.8,ylim=c(0,350))
  abline(0,1)
  points(MEDIAN_PIX[,"ONSET"],MEDIAN_PIX[,"TNDVImax"],pch=".",cex=1,col="green")
  points(MEDIAN_YR[,"ONSET"],MEDIAN_YR[,"TNDVImax"],pch=21,cex=1,bg="darkgreen")
  legend("bottomright", pch=21,cex=0.8,pt.cex=2,pt.bg=c("gray","green","darkgreen"), c("All pixel-year combination","Median per pixel","Median per year"))
  
  hist(SLO,breaks=seq(-15,15,0.02),xlim=c(-3,3),main="",xlab="Slope of TNDVImax~Onset") # a majority of sites with slope<1 (=rattrappage ph?no)
  abline(v=median(SLO,na.rm=T),col="red")
  abline(v=1,lwd=1,lty=2)
  
  
    # E. NDVImax, NDVIint & PLATEAU length ----
  
  SCA<-0.45
  graphics.off()
  windows(SCA*25, SCA*20)
  par(mar=c(4,5,0,1),cex.lab=1.5,mgp=c(3.5,1,0),mfcol=c(3,1))
  
  boxplot(myNDVImax,ylim=c(0.4,1),ylab="NDVImax",add=F,col="gray",outline=F,names=YEAR,las=2)
  abline(h=median(apply(myNDVImax,2,median,na.rm=T),na.rm=T),col=c('darkgray'),lwd=3) 
  boxplot(myNDVIint,ylab="integrated NDVI",ylim=c(0,250),col=colors()[86],outline=F,names=YEAR,las=2)
  abline(h=median(apply(myNDVIint,2,median,na.rm=T),na.rm=T),col=c('darkgreen'),lwd=3)
  boxplot(myPLAT,add=F,ylim=c(0,200),col="gray",outline=F,ylab="Length of plateau",names=YEAR,las=2)
  abline(h=median(apply(myPLAT,2,median,na.rm=T),na.rm=T),col=c('darkgray'),lwd=3) 
  
  # hist(RSQ,breaks=seq(0,1,0.02),xlim=c(0,1))
  
    # F. TRENDS ON NDVImax and NDVIint ----
  RTA  <- matrix(NA,length(CAN_ID250),2)
  for (i in 1:length(CAN_ID250)){
    print(i)
    Y <- as.numeric(myNDVImax[i,])
    Z <- YEAR
    MY <- which(is.na(Y)==TRUE)
    if (length(MY)==0){
      # LM <- lm(Y~X)
      # LM <- mblm(Y~X,repeated=FALSE)
      # LR <- lm(residuals(LM)~Z)
      # Y2 <- as.vector(residuals(LM))
      LR       <- mblm(Y~Z,repeated=FALSE)
      RTA[i,1] <- summary(LR)$coefficients[2,1]
      RTA[i,2] <- summary(LR)$coefficients[2,4]
    }
  }
  
  rownames(RTA) <- CAN_ID250
  RAN    <- round(range(1000*RTA[,1]))
  BREAKS <- seq(RAN[1]-1,RAN[2]+1,0.25)

  NS.g <- which(RTA[,2]>0.05 & RTA[,1]>0)
  S.g  <- which(RTA[,2]<=0.05 & RTA[,1]>0)
  SS.g  <- which(RTA[,2]<=0.01 & RTA[,1]>0)
  SSS.g  <- which(RTA[,2]<=0.001 & RTA[,1]>0)
  NS.b <- which(RTA[,2]>0.05 & RTA[,1]<0)
  S.b  <- which(RTA[,2]<=0.05 & RTA[,1]<0)
  SS.b  <- which(RTA[,2]<=0.01 & RTA[,1]<0)
  SSS.b  <- which(RTA[,2]<=0.001 & RTA[,1]<0)
  
  graphics.off()
  windows(SCA*25, SCA*20)
  par(mar=c(5,4,0,1),cex.lab=1.5,mgp=c(2.5,1,0),mfcol=c(1,1))
  hist(1000*RTA[,1],breaks=BREAKS,xlab="slope NDVImax~Year",main="") # tr?s nette tendance au verdissemen
  hist(1000*RTA[S.g,1],breaks=BREAKS,col="lightgreen",add=T)
  hist(1000*RTA[SS.g,1],breaks=BREAKS,col="green",add=T)
  hist(1000*RTA[SSS.g,1],breaks=BREAKS,col="darkgreen",add=T) 
  hist(1000*RTA[S.b,1],breaks=BREAKS,col="orange",add=T)
  hist(1000*RTA[SS.b,1],breaks=BREAKS,col="red",add=T)
  hist(1000*RTA[SSS.b,1],breaks=BREAKS,col="darkred",add=T) 
  legend("topright", pch=22,cex=1,pt.cex=2,pt.bg=c("white","lightgreen","green","darkgreen"), c("ns","P<0.05","P<0.01","P<0.001"),title="Greening trends")
  legend("topleft", pch=22,cex=1,pt.cex=2,pt.bg=c("white",colors()[57],"red","darkred"), c("ns","P<0.05","P<0.01","P<0.001"),title="Browning trends")
  
  
    # MAPS OF GREENING / BROWNING ----
  CEX <- 2
  graphics.off(); windows(10,10);par(mar=c(4,4,1,1))
  plot(CAN_COO,pch=22,bg="lightgray",col="lightgray",cex=CEX,xlab="X_etrs89",ylab="Y_etrs89")
  
  points(CAN_COO[S.g,],pch=22,bg="lightgreen",col="lightgreen",cex=CEX)
  points(CAN_COO[SS.g,],pch=22,bg="green",col="green",cex=CEX) 
  points(CAN_COO[SSS.g,],pch=22,bg="darkgreen",col="darkgreen",cex=CEX)
  
  points(CAN_COO[S.b,],pch=22,bg=colors()[57],col=colors()[57],cex=CEX)
  points(CAN_COO[SS.b,],pch=22,bg="red",col="red",cex=CEX) 
  points(CAN_COO[SSS.b,],pch=22,bg="darkred",col="darkred",cex=CEX)
  
  
  # extract ID of CAN sites
  DEM250    <- raster("C:/Users/cholerp/Documents/PROJETS/ALPAGES_SENTINELLES/DATA_ANALYSIS/DEM250.tif")
  COO250 <- coordinates(DEM250)
  
  # Identify picxel ID wihin AOIext
  SNOB_ID250 <- which(!is.na(extract(AOIras,COO250))) # 9554 cells
  
  plot(AOIras)
  points(COO250[SNOB_ID250,],pch="+")
  
  SNOB_COO <- COO250[SNOB_ID250,]
  SNOB_SP  <- SpatialPoints(SNOB_COO); SNOB_SP@proj4string <- LAMB93
  
  # call gsfm function
  SRC.DIR <- "C:/Users/cholerp/Documents/SIG/MOD09Q1.006.h18v04/NDVI_RES250/"
  SNOB_NDVI250   <- gfsm(SNOB_ID250,SRC.DIR,RETURN="ALL")
  
  
  # d?pouillement des sorties de
  METRICS <- SNOB_NDVI250$METRICS
  N       <- length(SNOB_ID250) # 9554 cells
  ONSET   <- matrix(unlist(lapply(METRICS, function(x) x[,"ONSET"])),N,18,byrow=T)
  OFFSET <- matrix(unlist(lapply(METRICS, function(x) x[,"OFFSET"])),N,18,byrow=T) 
  INT02  <- matrix(unlist(lapply(METRICS, function(x) x[,"NDVIint02"])),N,18,byrow=T)
  INT04  <- matrix(unlist(lapply(METRICS, function(x) x[,"NDVIint04"])),N,18,byrow=T)
  PEAK   <- matrix(unlist(lapply(METRICS, function(x) x[,"TNDVImax"])),N,18,byrow=T)  
  PLAT90 <- matrix(unlist(lapply(METRICS, function(x) x[,"NDVIplateau90"])),N,18,byrow=T)
  PLAT80 <- matrix(unlist(lapply(METRICS, function(x) x[,"NDVIplateau80"])),N,18,byrow=T)
  MIN    <- matrix(unlist(lapply(METRICS, function(x) x[,"NDVImin"])),N,18,byrow=T) 
  MAX    <- matrix(unlist(lapply(METRICS, function(x) x[,"NDVImax"])),N,18,byrow=T) 
  MAXM   <- matrix(unlist(lapply(METRICS, function(x) x[,"NDVImaxm"])),N,18,byrow=T) 
  GRDUR  <- PEAK - ONSET        # Duration of growing period
  SEN30  <- matrix(unlist(lapply(METRICS, function(x) x[,"SEN30"])),N,18,byrow=T) 
  SEN45  <- matrix(unlist(lapply(METRICS, function(x) x[,"SEN45"])),N,18,byrow=T) 
  REG30  <- matrix(unlist(lapply(METRICS, function(x) x[,"REG30"])),N,18,byrow=T)
  REG45  <- matrix(unlist(lapply(METRICS, function(x) x[,"REG45"])),N,18,byrow=T)
  
  # NDVI trend using a calculation of LRR as for Landsat 
  MAXINI <- apply(MAX[,1:3],1,max)   # max over 2000:2001:2002
  MAXEND <- apply(MAX[,14:16],1,max) # max over 2013:2014:2015
  LRR.M  <- log(MAXEND/MAXINI)
  
  LRRMrast <- DEM250; LRRMrast[]<-NA; LRRMrast[SNOB_ID250]<-LRR.M; LRRMrast <- crop(LRRMrast,extent(REF25))
  plot(LRRMrast)

  
  # NDVI trend using mblm
  RTA  <- matrix(NA,N,2)
  YEAR <- 2000:2017
  for (i in 1:N){
    print(i)
    Y <- MAX[i,]
    Z <- YEAR
    MY <- which(is.na(Y)==TRUE)
    if (length(MY)==0){
      LR       <- mblm(Y~Z,repeated=FALSE)
      RTA[i,1] <- summary(LR)$coefficients[2,1]
      RTA[i,2] <- summary(LR)$coefficients[2,4]
    }
  }
  
  rownames(RTA) <- SNOB_ID250
  1000*range(RTA[,1],na.rm=T)
  
  RTArast <- DEM250; RTArast[]<-NA; RTArast[SNOB_ID250]<-1000*RTA[,1]; RTArast <- crop(RTArast,extent(REF25))
  plot(RTArast)

  graphics.off()
  windows(10,10)
  BREAKS <- seq(-7,14,0.25)
  hist(1000*RTA[,1],freq=T,breaks=BREAKS,xlab="mean change of NDVImax per year (x1000)",main="") # tr?s nette tendance au verdissemen
  hist(1000*RTA[which(RTA[,2]<0.05),1],freq=T,breaks=BREAKS,col="lightgreen",add=T)
  hist(1000*RTA[which(RTA[,2]<0.01),1],freq=T,breaks=BREAKS,col="green",add=T)
  hist(1000*RTA[which(RTA[,2]<0.001),1],freq=T,breaks=BREAKS,col="darkgreen",add=T) 
  
  # Identify hotspot of greening using LRR.M
  BREAKS   <- seq(-1,2,0.01)
  H        <- hist(LRR.M,freq=F,breaks=BREAKS)
  LRR.M90  <- BREAKS[which(cumsum(H$density)>90)[1]] # 0.32
  abline(v=LRR.M90)
  
  THR       <- LRR.M90
  m         <- c(-Inf, THR, 0,  THR, +Inf, 1)
  rclmat    <- matrix(m, ncol=3, byrow=TRUE)
  LRRM90    <- reclassify(LRRMrast, rclmat)
  windows(10,10); plot(LRRM90)
  
  
  # Identify hotspot of greening using MBLM.M  
  BREAKS <- seq(-7,14,0.25)
  H <- hist(1000*RTA[,1],freq=F,breaks=BREAKS)
  
  mean(1000*RTA[,1])
  GAU <- fitdistr(1000*RTA[,1], "normal")
  
  plot(function(x) dnorm(x, mean=GAU$estimate[1], sd=GAU$estimate[2]),-5,10,add=T,col="red")
  
  QUAN <- qnorm(c(0.2,0.8), mean=GAU$estimate[1], sd=GAU$estimate[2], lower.tail = TRUE, log.p = FALSE)
  
  abline(v=QUAN,col="red")
  
  # use the 10 or 20% top quantile to show responsiveness of snowbeds ?

  THR      <- QUAN[2]
  m        <- c(-Inf, THR, 0,  THR, +Inf, 1)
  rclmat   <- matrix(m, ncol=3, byrow=TRUE)
  MBLM90    <- reclassify(RTArast, rclmat)
  windows(10,10); plot(MBLM90)
  
  
  save.image("SNOWBED_RS.RData")
    
  # III.3.Correspondence between LRR.L and MODIS ----
  y <-  raster::aggregate(raster(LRR.L,3),fact=10,fun=mean)
  
  plot(RTAf[],y[],pch=".")
  abline(lm(y[]~RTAf[]))
  
  
  save.image("SNOWBED_RS.RData")
  
  # III.4. Link MOD/LANDSAT and ORTHO ----
  # First test using the reference tile of figure 1C (exg3 at 2m)
  # BD_ORTHO provides a way to locate Persistent Summer Snow Patches
  
  EXG25 <- crop(projectRaster(100*aggregate(exg3,fact=12,fun=sum,na.rm=T)/144,REF25),extent(exg3))
  graphics.off();windows(10,10);plot(EXG25)
  
  GRE25 <- crop(raster(LRR.L,3),extent(exg3))
  
  EXG25C <- cut(EXG25[],seq(0,50,5))
  boxplot(GRE25[]~EXG25C,outline=F)
  abline(h=0)
  
  EXGPOL <- as(extent(EXG250),"SpatialPolygons")
  EXGPOL@proj4string <- LAMB93

  EXG250_ID <- which(!is.na(over(SNOB_SP,EXGPOL)))
  
  EXG250    <- aggregate(EXG25,10,mean,na.rm=T)
  plot(EXG250)
  points(SNOB_COO)
  points(SNOB_COO[EXG250_ID,],pch=21,bg="red")
  
  X  <- extract(EXG250,SNOB_COO[EXG250_ID,])
  Y  <- RTA[EXG250_ID,1]
  XC <- cut(X,seq(0,20,2))
  boxplot(Y~XC)
  
  save.image("SNOWBED_RS.RData")

# IV. FIGURES ----
  # FIG 1A ----   

  # graphics.off() 
  # windows(20,17.8) 
  # plot(DEM25,col=DEM_pal1(100),zlim=c(0,4000),alpha=0.3,add=F)
  # plot(RWBODYl93,add=T,col="blue")
  # plot(DEM25*POI,col=DEM_pal1(100),add=T)
  # plot(hill,add=T,col=GRAY_pal1(100),alpha=0.4)
  
  b_poly <- as(raster::extent(bbox(DEM25)), "SpatialPolygons")
  RWBODYl93clip <- gIntersection(RWBODYl93, b_poly, byid = TRUE)
  
  
  # legend.mar 	Width in characters of legend margin that has the axis. Default is 5.1 for a vertical legend and 3.1 for a horizontal legend.
  
  XLIM <- c(950000,990000)
  graphics.off() 
  windows(ncol(hill)/100,nrow(hill)/100) 
  # windows(20,17.5) 
  par(mar=c(0,0,0,0),oma=c(0,0,0,0),usr=c(950000,990000,6430000,6465000),pty="m")
  plot(hill,col=GRAY_pal1(100),alpha=1,add=F, legend=F,smallplot= c(1,1,0,0))
  
  plot(DEM25,col=GRAY_pal1(100),zlim=c(0,4000),alpha=0.5,add=T, legend=F,xlim=XLIM)
  abline(v=seq(950000,990000,5000),h=seq(6430000,6465000,5000)) 
  plot(RWBODYl93clip,add=T,col="blue",legend=F,xlim=XLIM,legend.mar=0)

  plot(DEM25*POI,col=DEM_pal1(18),alpha=0.7,add=T,legend.mar=0,xlim=XLIM,zlim=c(1700,3500),legend=F,smallplot= c(.85,.9,0,1),graphics.reset=F)
  
  
  # tout ce qui vien apr?s ne se superpose plus correctement si pas de controle sur smallplot
  abline(v=seq(950000,990000,5000),h=seq(6430000,6465000,5000)) 
  plot(RWBODYl93clip,add=T,col="blue",legend=F,legend.mar=0,smallplot= c(.85,.9,0,1))
  # plot(RWBODYl93,add=T,col="white") # ne se superpose pas ! La solutoin passe par smallplot avec des coordonn?es pour le faire disparaitre
  
  plot(DEM25*POI,col=DEM_pal1(38), add=T,legend.mar=0,xlim=XLIM,zlim=c(1700,3500),legend.only=T,horizontal=F,smallplot= c(.92,.94,0.5,0.9))
 
  
  # plot(DEM25*MASKfsfd*MASKaoi,zlim=c(0,4000),col=DEM_pal1(100),add=T,legend=F)
  # plot(MASSIFS,add=T) 
  # plot(D05,add=T)
  # plot(D73,add=T) 

  # FIG 1B ----
  H <- hist((DEM25*POI)[],breaks=seq(1700,3500,100),col=DEM_pal1(18),xlab="Elevation (m)",main="")
  NAME <- c("1700m-\n1800m","1800m-\n1900m","1900m-\n2000m","2000m-\n2100m","2100m-\n2200m","2200m-\n2300m","2300m-\n2400m","2400m-\n2500m","2500m-\n2600m","2600m-\n2700m","2700m-\n2800m","2800m-\n2900m","2900m-\n3000m","3000m-\n3100m","3100m-\n3200m","3200m-\n3300m","3300m-\n3400m","3400m-\n3500m")
  windows(10,5); par(mar=c(8,5,1,1),mgp=c(4,1,0))
  barplot(H$counts*25*25/10000,col=DEM_pal1(19),xlab="Elevation (m)",ylab="Surface (ha)",ylim=c(0,5000),names.arg = NAME,las=2)
  
  # FIG 1B
  # COL <- adjustcolor(OSOcol2, alpha.f = 0.3) 
  # graphics.off() 
  # windows(20,17.8) 
  # plot(SNOBoso2,col=COL,legend=F) # do not use alpha directlty in plot function
  # plot(SNOBoso2*MASKfsfd*MASKaoi,add=T,legend=T,col=OSOcol2)
  # abline(v=seq(950000,990000,5000),h=seq(6430000,6465000,5000))  
  
  # graphics.off() 
  # windows(20,17.8) 
  # plot(FSFDl2e,col=DEM_pal1(100))
  # abline(v=seq(900000,940000,5000),h=seq(2000000,2030000,5000))
  
  
  # FIG 1C (example with tile 0980-6450 located in Thabor)
  ID03            <- grep("05-2003-0980-6450-LA93",LFtif03)  
  ID09            <- grep("05-2009-0980-6450-LA93",LFtif09)
  ID99            <- grep("05-1999-0980-6450-LA93",LFtif99)
  
  graphics.off()
  windows(30,10)
  par(mfrow=c(1,3),oma=c(0,0,0,0),mar=c(0,0,0,0))
  plotRGB(TEST03 <- stack(LFtif03[ID03]))
  plotRGB(TEST09 <- stack(LFtif09[ID09]))
  plotRGB(TEST99 <- stack(LFtif99[ID99]))
  
  BRIGHT <- stackApply(TEST99, indices=c(1,1,1), fun=sum)
  BLUE   <- raster(TEST99,3)
  
  graphics.off();windows(10,10)
  dh    <- hist(raster(TEST99,3)[],breaks=seq(0,255,1))
  ins   <- dh[["density"]]
  ss    <- which(diff(sign(diff(ins)))==2)+1
  abline(v=dh[["mids"]][ss],col="red")
  
  THR     <- 163
  m       <- c(-Inf, THR, 0,  THR, +Inf, 1)
  rclmat  <-matrix(m, ncol=3, byrow=TRUE)
  exg3    <- reclassify(raster(TEST99,3), rclmat)
  windows(10,10); plot(exg3)
  
  # canal bleu normalis? par la brightness
  
  plot(BLUEN,zlim=c(0.2,0.4))
  
  # FIG 2 ---- 
  SNO <- raster("SNOW_OCCU_07.tif")+raster("SNOW_OCCU_06.tif")
  SNOt <- SNO*POI 
  
  # Snow occurrence from a combination of sentinel-2 and landsat
  XLIM <- c(950000,990000)
  graphics.off() 
  windows(ncol(hill)/100,nrow(hill)/100) 
  # windows(20,17.5) 
  par(mar=c(0,0,0,0),oma=c(0,0,0,0),usr=c(950000,990000,6430000,6465000),pty="m")
  plot(hill,col=GRAY_pal1(100),alpha=1,add=F, legend=F,smallplot= c(1,1,0,0))
  
  plot(DEM25,col=GRAY_pal1(100),zlim=c(0,4000),alpha=0.5,add=T, legend=F,xlim=XLIM)

  
  plot(SNOt,col=rev(THEpal(31)),alpha=0.6,add=T,legend.mar=0,xlim=XLIM,legend=F,smallplot= c(1,1,0,0),zlim=c(1,61))

  abline(v=seq(950000,990000,5000),h=seq(6430000,6465000,5000)) 
  plot(RWBODYl93clip,add=T,col="blue",legend=F,xlim=XLIM,legend.mar=0)
  
  # FSFD map from sentinel-2 over the area of interest and FSFD-NDVImax relationship

  graphics.off() 
  windows(20,17.8) 
  plot(FSFDt,col=GRAY_pal1(100),alpha=1,add=F)
  
  LLS <- which(FSFDt[]>=160 & FSFDt[]<=210) # 102271 points
  
  hist(FSFDt[LLS])
  
  NDVIt1 <- NDVImaxspotaoi1*POI
  NDVIt2 <- NDVImaxspotaoi2*POI
  
  plot(NDVIt1[LLS],NDVIt2[LLS],pch=".")
  
  hist(NDVIt2[LLS],breaks=seq(-0.3,0.75,0.01))
  
  
  SEQ <- c(seq(160,200,5),210)
  CUT <- cut(FSFDt[LLS],breaks=SEQ)
  HIS <- hist(FSFDt[LLS],breaks=SEQ)

  QA <- matrix(NA,9,9) 
  for (i in 1:9) {
    Q<- 0.1*i
    QUA <- function(x) quantile(x,Q,na.rm=T)
    QA[,i] <- aggregate(NDVIt2[LLS],list(Q=CUT),"QUA")$x
  }
  
  graphics.off();windows(10,10)
  plot(FSFDt[LLS],NDVIt2[LLS],pch=".",col="gray",cex=2,ylim=c(0,0.6))
  for(i in c(5,8,9)) lines(HIS$mids,QA[,i],type="b",pch=21,bg="black")
  
  
  # adding marginal histograms
  scatterhist = function(x, y, xlab="", ylab=""){
    zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
    layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
    xhist = hist(x, plot=FALSE,breaks=seq(160,210,5))
    yhist = hist(y, plot=FALSE,breaks=seq(0,0.65,0.01))
    top = max(c(xhist$counts, yhist$counts))
    par(mar=c(3,3,1,1))
    plot(x,y,pch=".",col="gray",cex=2,xlim=c(160,210),ylim=c(0,0.65))
    G=21
    for(i in c(5,8,9)) {
      lines(HIS$mids,QA[,i],type="b",pch=G,cex=1.5,bg="black")
      G=G+1
    }
    legend("topright",pch=c(21,22,23),legend=c("Quantile 0.5","Quantile 0.8","Quantile 0.9"),pt.bg="black",cex=1.5)
    par(mar=c(0,3,1,1))
    barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
    par(mar=c(3,0,1,1))
    barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
    par(oma=c(3,3,0,0))
    mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
          at=0.3)
    mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
          at=0.3)
  }
  
  ds <- data.frame(FSFD=FSFDt[LLS], NDVI=NDVIt2[LLS])
  ds <- ds[-which(ds[,2]<0 | ds[,2]>0.650),]
  with(ds, scatterhist(FSFD, NDVI, xlab="FSFD (Julian Day)", ylab="NDVImax"))
  
# density contour line plot using ggplot
  df <- data.frame(FSFD=FSFDt[LLS], NDVI=NDVIt[LLS])
  ggplot(df,aes(x=FSFD,y=NDVI)) + 
    geom_point(colour="blue", alpha=0.2) + 
    geom_density2d(colour="black")
  
  graphics.off(); windows(10,10)
  ggplot(data=df,aes(FSFD,NDVI)) + 
    geom_point(size=0.1,alpha=0.4, colour="black") +
    scale_fill_continuous(low="green",high="red") +
    guides(alpha="none")+
    stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') 

  
  # FIG 3 ----
  # Bivariate relationship between snowcover duration and long-term greening
  # better to have a synthetic boxplot analysis
  
  x <- SNOt[]
  y <- (POI*raster(LRR.L,3))[]
  graphics.off();plot(x,y,pch='.')
  
  graphics.off();plot(SNOt[],FSFDt[],pch=".")
  
  # optimisatino sur DAH avec 210? !
  
  df  <- na.omit(data.frame(snot=SNOt[],fsfd=FSFDt[],dem=DEMt[],dah=DAHt[],she=SHELTt[],tpi=TPIt[]))
  
  DAH         <- cos(210*pi/180-aspect)*atan(slope)
  DAHt        <- DAH*POI
  
  dfrn  <- na.omit(data.frame(fsfd=crop(FSFDt,RNse)[],dem=crop(DEMt,RNse)[],dah=crop(DAHt,RNse)[],she=crop(SHELTt,RNse)[],tpiL=crop(TPIs25t,RNse)[],tpiS=crop(TPIs5t,RNse)[]))
  
  library(lme4)
  model <- lmer(fsfd ~ tpiL + dah +(1+tpiL | elevc) + (1+dah | elevc), data=dfrn)
 
   plot(dfrn$fsfd,fitted(model),pch='.')
  abline(0,1)
  AIC(model)
  
  # use category for elevation
  graphics.off();hist(dfrn$dem)
  dfrn$elevc <- cut(dfrn$dem,breaks=seq(1800,3200,100)) # 1 valeur au-dessus
  dfrn<- na.omit(dfrn)
  
  ELEVC <- levels(dfrn$elevc);nreg <- length(ELEVC)
  filter <- dfrn$elevc==ELEVC[i]
  
  regdata <- NULL 
  for (i in 1:nreg) { 
    filter <- dfrn$elevc==ELEVC[i]
    y <- dfrn$fsfd[filter] 
    X <- cbind(1,      # intercept placeholder 
               dfrn$tpiL[filter], 
               dfrn$dah[filter]) 
    regdata[[i]] <- list(y=y, X=X) 
  } 
  
  Data <- list(regdata=regdata) 
  Mcmc <- list(R=10000) 
  
  library(bayesm)
  out <- bayesm::rhierLinearModel(Data=Data,Mcmc=Mcmc)

  
  
  
  (BETA <- apply(out$betadraw, c(1,2), mean))
  dfrn$fsfdPB <- NA
  
  for (i in 1:nreg) { 
    filter <- dfrn$elevc==ELEVC[i]
    dfrn$fsfdPB[filter]  <- BETA[i,1]+BETA[i,2]*dfrn$tpiL[filter]+BETA[i,3]*dfrn$dah[filter]
  }
  plot(dfrn$fsfd,dfrn$fsfdPB,pch=".")
  abline(0,1)
  
  
  
  
  plot(out$betadraw)
  
  summary(LM1 <- lm(fsfd~dem,data=dfrn)) # 54 % de variance expliqu?e
  dfrn$R1  <- residuals(LM1)

  
  summary(LM1var <- lm(R1~tpiL+dah,data=dfrn)) # 42% de variance expliqu?e
  plot(dfrn$fsfd,predict.lm(LM1,data=dfrn),pch=".")
  plot(dfrn$R1,predict.lm(LM1var,dfrn),pch=".")
  
  dfrn$fsfdP <- predict.lm(LM1,dfrn)+predict.lm(LM1var,dfrn)
  plot(dfrn$fsfd,dfrn$fsfdP,pch=".")
  abline(lm(fsfdP~fsfd,data=dfrn),col="red")
  abline(0,1)
  summary(lm(fsfdP~fsfd,data=dfrn))

  summary(LM <- lm(fsfd~dem+tpiS+dah,data=dfrn)) # 74 % de variance expliqu?e
  plot(dfrn$fsfd,predict.lm(LM,dfrn),pch=".")
  abline(0,1)
  
  windows(10,8); plot(crop(DAH,RNse))
  plot(crop(TPI25,RNse))
  graphics.off(); 
  windows(10,8); plot(crop(DAHt,RNse),col=SNOWpal(50))
  windows(10,8); plot(crop(SHELTt,RNse),col=rev(SNOWpal(50)))

  graphics.off(); windows(10,8); plot(crop(TPIt,RNse),col=SNOWpal(50))

  windows(10,8); plot(crop(FSFDt,RNse),col=rev(GRAY_pal1(50)))
  windows(10,8); plot(crop(POI,RNse),col=SNOWpal(50))
  windows(10,8); plot(crop(SNOBoso2,RNse),col=heat.colors(8))
  

  graphics.off() 
  windows(10,10) 
  plot(dfrn$dah,dfrn$she,pch=".")
  plot(dfrn$dah,dfrn$tpi,pch=".")
  plot(DEMt[],FSFDt[],pch=".")
  plot(DEMt[LLS],FSFDt[LLS],pch=".")
  plot(DAHt[LLS],FSFDt[LLS],pch=".")
  plot(TPIt[LLS],FSFDt[LLS],pch=".")
  plot(SHEt[LLS],FSFDt[LLS],pch=".")
  

  
  w <- 0.5
  dfrn$var <- (1-w)*dfrn$tpi+w*dfrn$dah
  windows(10,8); plot(crop(TPIt,RNse)+crop(DAHt,RNse),col=SNOWpal(50))  
  
  
  summary(LM1var <- lm(R1~tpi+dah,data=dfrn)) # 35% de variance expliqu?e
  summary(LM1dah <- lm(R1~dah,data=dfrn))     # 26% de variance expliqu?e
  
    LM1she <- lm(R1~she,data=dfrn)
    LM1dah <- lm(R1~dah*she,data=dfrn) # 26% de varaiance expliqu?e
  LM1tpi <- lm(R1~tpi,data=dfrn) 
  graphics.off();windows(10,10);plot(dfrn$var,R1,pch=".");abline(v=0,col="red");abline(h=0,col="red")
  graphics.off();windows(10,10);plot(dfrn$dah,R1,pch=".");abline(v=0,col="red");abline(h=0,col="red")
  
    plot(df$tpi,R1,pch=".",xlim=c(-1,1));  abline(v=0,col="red");abline(h=0,col="red")
  plot(df$she,R1,pch=".");  abline(v=0,col="red");abline(h=0,col="red")
  w=0.5; form <- (1-w)*df$tpi+w*df$dah
  plot(form,R1,pch=".",xlim=c(-1,1));  abline(v=0,col="red");abline(h=0,col="red")
  
  LM2 <- lm(fsfd~dem+tpi+dah,data=df)
  plot(df$fsfd,predict.lm(LM2,df),pch=".")
  abline(0,1,col="red")
  
  LM1 <- lm(FSFDt[LLS]~DEMt[LLS])
  R1  <- residuals(LM1)
  LM2 <- lm(FSFDt[LLS]~DEMt[LLS]+DAHt[LLS])
  LM2 <- lm(FSFDt[LLS]~DEMt[LLS]+DAHt[LLS]+TPIt[LLS])
  summary(LM2)
  anova(LM1,LM2)
  
  plot(DAHt[LLS],R1,pch=".")
  w=0.7;
  plot(w*TPIt[LLS]+(1-w)*DAHt[LLS],R1,pch=".")
  abline(v=0,col="red");abline(h=0,col="red")
  

# OTHER PARTS NOT INCLUDED YET 


# OTHER STUFF ----
  # IV.E. IMPACT OF HEATWAVES ON SNOWBEDS ----
  
  windows(10,10);plot(X <- CANIC80_END,Y<- apply(OFFSETsno[-DELsno,],2,mean)[-c(17,18)] - apply(PEAKsno[-DELsno,],2,mean)[-c(17,18)],pch="" )
  text(X,Y,labels=2000:2015)
  
  
  # IV.C. CLASSIFICATION from TIME-SERIES -----
  
  TS <- matrix(NA,nrow(SIT_TOT),11)
  rownames(TS) <- SIT_TOT$Id_station
  for(i in 1:11){
    print(i)
    TS[,i] <- extract(raster(LF[i]),SIT_TOT[,3:4])
  }
  
  TSI <- matrix(NA,nrow(SIT_TOT),11)
  rownames(TS) <- SIT_TOT$Id_station
  for(i in 1:11){
    print(i)
    TSI[,i] <- extract(raster(LFS[i]),SIT_TOT[,3:4])
  }
  
  
  require(zoo);require(signal);require(gtools)
  TSf <- TS ; TSIf <- TS
  NAL <- rep(NA,nrow(SIT_TOT))
  THR <- 3
  for(i in 1:nrow(SIT_TOT)) {
    print(i)
    NAL[i] <- length(which(is.na(TS[i,])==TRUE))
    if (NAL[i]<=THR) TSf[i,]  <- na.spline(TS[i,])
    NAL[i] <- length(which(is.na(TS[i,])==TRUE))
    if (NAL[i]<=THR) TSIf[i,] <- na.spline(TSI[i,])
  }
  TSf   <- TSf[-which(NAL>THR),]
  TSIf  <- TSIf[-which(NAL>THR),]
  DEL   <- which(NAL>THR) # 2046 si THR=2
  (length(DEL))
  
  PARTf <- PARTa[[1]][-DEL,NC]
  
  TSCLUST <- tsclust(cbind(TSf,TSIf),k=12,centroid="pam")
  table(TSCLUST@cluster)
  
  table(PARTf[which(TSCLUST@cluster==11)]) # spot 10 <->  veget 5 & 11
  table(PARTf[which(TSCLUST@cluster==3)])  # spot 3  <->  veget 5
  table(PARTf[which(TSCLUST@cluster==4)])  # spot 4  <->  veget 5
  table(PARTf[which(TSCLUST@cluster==1)])  # spot 1  <->  veget 9, 1
  
  graphics.off(); windows(15,15); par(mfrow=c(2,2),mar=c(1,3,0,0))
  boxplot(NDVI05[-DEL]~TSCLUST@cluster,ylim=c(0,1));abline(h=c(0.2,0.4,0.6),lty=2)
  boxplot(NDVI06[-DEL]~TSCLUST@cluster,ylim=c(0,1));abline(h=c(0.2,0.4,0.6),lty=2)
  boxplot(NDVI07[-DEL]~TSCLUST@cluster,ylim=c(0,1));abline(h=c(0.2,0.4,0.6),lty=2)
  boxplot(NDVI08[-DEL]~TSCLUST@cluster,ylim=c(0,1));abline(h=c(0.2,0.4,0.6),lty=2)
  
  # k-means clustering of time series
  
  # Issue : classification tree from longitudinal analyses - how to proceed
  # Regression trees from longitudinal and multiresponse data
  # Classification of repeated measurements data using tree-based ensemble methods
  # Previous algorithms for constructing regression tree models for longitudinal and multiresponse data have mostly followed the CART approach
  # CART : Computational difficulties in estimating the covariance matrices limit the method to data observed at equally-spaced time points
  # CART-based RE-EM tree
  
  library(REEMtree)
  help(package="REEMtree")
  
  # https://cran.r-project.org/web/views/FunctionalData.html
  
  # Trees for modelling longitudinal data by means of random effects is offered by package REEMtree
  
  data(simpleREEMdata)
  REEMresult<-REEMtree(Y~D+t+X, data=simpleREEMdata, random=~1|ID)
  print(REEMresult)
  
  # Clustering of time series
  # Time series clustering is implemented in TSclust, dtwclust, BNPTSclust and pdc. 
  
  library(Funclustering)
  # funclust : This function run clustering algorithm for multivariate functional data.
  
  # fdakma: Functional Data Analysis: K-Mean Alignment
  library(fdakma)
  
  data(kma.data)
  x <- kma.data$x # abscissas
  y0 <- kma.data$y0 # evaluations of original functions
  y1 <- kma.data$y1 # evaluations of original function first derivatives
  
  library(dtwclust)
  
  
  tsclust(series, k = 4L, distance = "L2", centroid = "pam",seed = 3247, trace = TRUE, control = partitional_control(nrep = 10L))
  
  
# I. BIOCLIMATE ----

  # I.A. ARAVO DB ----
  MOD_ARA <- read.csv('C:/Users/cholerp/Documents/PROJETS/ARAVO/DATA_ANALYSIS/MODIS_ARA.csv')

  # Ground validation of Aravo snowmelt
  # Was a proposed supplementary figure for PPES 2018, deleted after review
  ARA.SNOW <- read.csv("ARA_SNOW.txt",head=T,sep="\t")
  windows(10,10); plot(ARA.SNOW[1:57,4],FSFD[,2]);abline(0,1) # fine
   FSFDsf <- 0.5*(ARA.SNOW[,5]+ARA.SNOW[,4])

  PRED <- FSFDsf; 
  PRED <- ARA.SNOW[,5]
  range(PRED,na.rm=T)
  graphics.off()
  windows(8,4)
  par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.8,0.5,0))

  plot(ARA.SNOW[,3],PRED,pch=21,bg='green',ylim=c(150,200),xlim=c(150,200),xlab='Observed First Snow Free Day',cex=0,ylab='Predicted First Snow Free Day')
  text(ARA.SNOW[,3],PRED,labels=ARA.SNOW$"ANNEE",cex=0.7)
  abline(LM<-lm(PRED~ARA.SNOW[,3]))
  summary(LM)
  abline(0,1,lty=2)
  mtext(side=3,adj=0.01,"(A)",line=0)


  line.cis(PRED,ARA.SNOW[,3])
  line.cis(ARA.SNOW[,2],ARA.SNOW[,5])
  line.cis(ARA.SNOW[,2],ARA.SNOW[,3])

# Compared distribution of predicetd values (recent and long period)

  PAS=10
  hist(PRED,breaks=seq(130,220,PAS),main="",ylab="Number of years",xlab="Predicted First Snow Free Day")
  hist(PRED[42:57],breaks=seq(130,220,PAS),add=T,col='gray')
  mtext(side=3,adj=0.01,"(B)",line=0)
  
  # I.B. FORCING SAFRAN at 2400m ----
  setwd("~/PROJETS/SNOWBED/DATA_ANALYSIS")
  SFALL <- read.csv("snow_beg_date.csv",sep=";")  # First day of snow >10 cm since August =JD 214
  SFREE <- read.csv("snow_free_date.csv",sep=";") # First snowfree day
# Attention au d?calage d'ann&e
  SFALL_SCL  <- c(apply(SFALL[42:57,5:6],1,mean),NA,NA)+214
  names(SFALL_SCL)<-1999:2016

  SFREE_SCL  <- c(apply(SFREE[41:57,5:6],1,mean),NA)
  names(SFREE_SCL)<-1999:2016

  GSL_SCL    <- SFALL_SCL-SFREE_SCL  
  
  # to be compared with
  RAW <- read.csv("C:/Users/cholerp/Documents/CLIMATO/DATA/CROCUS/PRO_3/oisans_2400_2700_FL_0_SNOW.csv")
  RAW <- read.csv("C:/Users/cholerp/Documents/CLIMATO/DATA/CROCUS/PRO_3/oisans_2700_3000_FL_0_SNOW.csv")
   YEAR <- as.numeric(substr(RAW[,1],1,4))
 FSFD.f <- function(x) which(x<=0.05)[1]
 FSFA.f <- function(x) {X <- which(x<=0.1); X[length(X)]} # seuil de 10 cm
 FSFD <- aggregate(RAW[,2],list(Y=YEAR),"FSFD.f")[-1,]
 windows(10,10); plot(FSFD[,1],FSFD[,2],type='l',lwd=2,col='blue');abline(v=1964,lty=2)
 # compare with Fig. 4 in Dedieu & al.
 FSFA <- aggregate(RAW[,2],list(Y=YEAR),"FSFA.f")[-1,]
 
 ARA.SNOW <- read.csv("ARA_SNOW.txt",head=T,sep="\t")
 windows(10,10); plot(ARA.SNOW[1:57,4],FSFD[,2]);abline(0,1) # fine

# Estimate of SAFRAN_CROCUS GSL - 2400-2700m - p?riode r?cente
  graphics.off()
  plot(SFALL[41:57,1],SFREE[41:57,5],type="b",ylim=c(100,350))

  plot(1999:2016,SFREE_SCL,type="b",ylim=c(100,350))
  abline(v=c(2004,2013,2015),col="red",lty=2)
  lines(1999:2016,SFALL_SCL,type='b')

  plot(1999:2016,GSL_SCL,type="b")
  abline(h=mean(GSL_SCL,na.rm=T))

  # I.C. HEATWAVES in OISANS (2400-2700m) - FLAT ----
    # Heat Waves (see alpages.r)
  load("C:/Users/cholerp/Documents/PROJETS/ALPAGES_SENTINELLES/DATA_ANALYSIS/HEATWAVE_90_ALPS.RData")
  names(HW)
  
  # CONDION
  DATA <- HW[[127]] # 125 is Oisans 1800_2100
  DATA <- HW[[125]] # 125 is Oisans 1800_2100
  
  YEAR<-2000:2017
  
  TMP <- NULL
  for (i in 1:length(YEAR)){
    X <- DATA[[i]]
    if(!is.na(X[1])){
      X <- as.data.frame(X)
      X$"YEAR" <- names(DATA)[i]
      X$"END"  <- X$"START"+X$"LENGTH"-1
      TMP <- rbind(TMP,X)
    }
  }
  TMP90   <- TMP
  CANIC90 <- TMP[which(TMP[,"LENGTH"]>=5 & TMP[,"START"]>=120 & TMP[,"START"]<=350),]
  
  load("C:/Users/cholerp/Documents/PROJETS/ALPAGES_SENTINELLES/DATA_ANALYSIS/HEATWAVE_80_ALPS.RData")
  names(HW)
  DATA <- HW[[125]]
  
  TMP <- NULL
  for (i in 1:length(YEAR)){
    X <- DATA[[i]]
    if(!is.na(X[1])){
      X <- as.data.frame(X)
      X$"YEAR" <- names(DATA)[i]
      X$"END"  <- X$"START"+X$"LENGTH"-1
      TMP <- rbind(TMP,X)
    }
  }
  TMP80   <- TMP
  CANIC80 <- TMP[which(TMP[,"LENGTH"]>=5 & TMP[,"START"]>=120 & TMP[,"START"]<=350),]
  
  # I.D. Synthetic figure of HEAT WAVES events over 2000:2017 (20180216) ----
  YEAR<-2000:2017
  graphics.off()
  windows(20,20);
  plot(CANIC90$START,as.numeric(CANIC90$YEAR),type="n",ylim=c(2000,2017),xlim=c(120,350),xlab="Time (Julian Days)",ylab="Year")
  mtext(side=3,line=1,text="Heat Waves Oisans 1800-2100 m",cex=2)
  legend("topright",legend=c("Seuil 80","Seuil 90"),col=c("gray","black"),lty=1,lwd=4)
  segments(x0=CANIC80$START,y0=as.numeric(CANIC80$YEAR),x1=CANIC80$END,y1=as.numeric(CANIC80$YEAR),lwd=6,col="darkgray")  
  segments(x0=CANIC90$START,y0=as.numeric(CANIC90$YEAR),x1=CANIC90$END,y1=as.numeric(CANIC90$YEAR),lwd=6)  
  abline(h=2000:2017,lty=2)
  
  # Add FSFD
  points(FSFD[42:57,2],2000:2015,pch='|',col='darkgreen',lwd=4,cex=1.5) # CROCUS
  points(FSFA[42:56,2],2000:2014,pch='|',col='blue',lwd=4,cex=1.5)
  
  segments(x0=FSFD[42:57,2],y0=2000:2014,x1=FSFA[42:56,2],y1=2000:2014,lwd=1)
  
  points(ARA.SNOW[42:57,3],2000:2015,pch='|',col='green',lwd=4,cex=1.5) # in situ snowmelt
  # points(ARA.SNOW[42:57,2],2000:2015,pch='|',col='gray',lwd=4,cex=1.5) # MODIS 
  points(MOD_ARA[1:16,"ONSET"],2000:2015,pch='|',col='magenta',lwd=4,cex=1.5) # MODIS onset
  points(MOD_ARA[1:16,"TNDVImax"],2000:2015,pch='|',col='red',lwd=4,cex=1.5) # MODIS peak 
  points(MOD_ARA[1:16,"OFFSET"],2000:2015,pch='|',col='magenta',lwd=4,cex=1.5) # MODIS offset 
  
  CANIC80_END <- rep(0,16)
  CANIC80_BEG <- rep(0,16)
  YEAR <- 2000:2015
  for (i in 1:length(YEAR)){
    EV <- MOD_ARA[which(MOD_ARA$X==YEAR[i]),]
    ID <- which(TMP80$YEAR==YEAR[i] & TMP80$START>=EV$TNDVImax & TMP80$END<=EV$OFFSET)
    if(length(ID)>0) {CANIC80_END[i] <- sum(TMP80[ID,"MAGN"])}
    ID <- which(TMP80$YEAR==YEAR[i] & TMP80$START>=EV$ONSET & TMP80$END<=EV$TNDVImax) 
    if(length(ID)>0) {CANIC80_BEG[i] <- sum(TMP80[ID,"MAGN"])}
  }
  
  CANIC90_END <- rep(0,16)
  YEAR <- 2000:2015
  for (i in 1:length(YEAR)){
    EV <- MOD_ARA[which(MOD_ARA$X==YEAR[i]),]
    ID <- which(TMP90$YEAR==YEAR[i] & TMP90$START>=EV$TNDVImax & TMP90$END<=EV$OFFSET)
    if(length(ID)>0) {CANIC90_END[i] <- sum(TMP90[ID,"MAGN"])}
  }
  
  
  graphics.off();windows(10,10);plot(X <- CANIC80_END,Y<- MOD_ARA[1:16,"OFFSET"]-MOD_ARA[1:16,"TNDVImax"]);text(X,Y,labels=YEAR)
  windows(10,10);plot(X <- CANIC80_END,Y<- MOD_ARA[1:16,"SEN45"])
  text(X,Y,labels=YEAR)
  
  windows(10,10);plot(X <- CANIC80_BEG,Y<- MOD_ARA[1:16,"TNDVImax"]-MOD_ARA[1:16,"ONSET"])
  text(X,Y,labels=YEAR)
  windows(10,10);plot(X <- CANIC80_BEG,Y<- MOD_ARA[1:16,"NDVImax"])
  text(X,Y,labels=YEAR)
  
  plot(X <- CANIC80_END,Y<- MOD_ARA[1:16,"OFFSET"]-MOD_ARA[1:16,"TNDVImax"])
  text(X,Y,labels=YEAR)
  
  # plot(ARA.SNOW[,2],ARA.SNOW[,3],pch=21,bg='green',ylab='Observations',cex=0.5,xlab='Estimations MODIS')
  # abline(LM<-lm(ARA.SNOW[,3]~ARA.SNOW[,2]))
  # text(ARA.SNOW[,2],ARA.SNOW[,3],labels=ARA.SNOW$"ANNEE",cex=0.7)
  # abline(0,1,lty=2)
  

# PART 1. OLD STUFF from OLD SCRIPT VERSION ----
  # OLD STUF ----
  graphics.off()
  windows(10,10)
  p <- radial.plot(lengths=META$ELEV_LAST,radial.pos=2*pi*META$ASPECT_MNT10m/360,start=pi/2,clockwise=F,rp.type="s",point.col="gray",point.symbols=16)
radial.plot.labels(META$ELEV_LAST,2*pi*META$ASPECT_MNT10m/360,labels=META$IDs,start=pi/2,clockwise=F,cex=0.7)
  # OLD STUFF Comparing First Snow Free Day (FSFD) ----
  TMPs <- read.csv("~/CLIMATO/DATA_ANALYSIS/BIOCLIM_SOIL_TEMP_S_FSFD.csv",row.names=1)
  TMPh <- read.csv("~/CLIMATO/DATA_ANALYSIS/BIOCLIM_SOIL_TEMP_H_FSFD.csv",row.names=1)

  SITc  <- intersect(rownames(TMPs),rownames(TMPh)) # 69 sites communs
  FSFDs <- TMPs[match(SITc,rownames(TMPs)),]
  FSFDh <- TMPh[match(SITc,rownames(TMPh)),-c(1,21)] # no simul for 1999 and 2019

  com.tmp <- as.factor(substr(row.names(FSFDs),5,6))
  CP      <- colors()
  levels(com.tmp)
  levels(com.tmp) <- c(be=CP[33],cc=CP[33],cf=CP[68],cr=CP[38],ct=CP[86],en=CP[128],es=CP[17],fg=CP[116],fh=CP[86],fl=CP[33],fp=CP[255],hs=CP[120],kd=CP[37],ks=CP[41],ld=CP[24],pa=CP[8],rk=CP[8],sc=CP[142],sh=CP[72],sr=CP[94],tr=CP[28],vm=CP[91])
  com.col <- as.vector(as.character(com.tmp))
  
  
  plot(unlist(FSFDh),unlist(FSFDs))
  abline(0,1)
  
  DIFF <- (FSFDh)-(FSFDs)
  graphics.off();windows(20,10); par(mar=c(6,5,1,1)); boxplot(t(DIFF),cex.axis=0.8,las=2,col=com.top,ylab="FSFDobs-FSFDmod")
  abline(h=seq(-50,50,10),lty=2)
  abline(h=0,lwd=2)
  abline(v=1:69,lty=3)
  
  # OLD STUFF Comparing Growing Degree Days (snowmelt to 31/07) (GDD) ----  
  TMPs <- read.csv("~/CLIMATO/DATA_ANALYSIS/BIOCLIM_SOIL_TEMP_S_GDD0731.csv",row.names=1)
  TMPh <- read.csv("~/CLIMATO/DATA_ANALYSIS/BIOCLIM_SOIL_TEMP_H_GDD0731.csv",row.names=1)
  
  SITc  <- intersect(rownames(TMPs),rownames(TMPh)) # 69 sites communs
  GDDs <- TMPs[match(SITc,rownames(TMPs)),]
  GDDh <- TMPh[match(SITc,rownames(TMPh)),-c(1,21)] # no simul for 1999 and 2019
  
  graphics.off(); windows(20,10); plot(unlist(GDDh),unlist(GDDs))
  abline(0,1)
  
  DIFF <- (GDDh)-(GDDs)
  graphics.off();windows(20,10); par(mar=c(6,5,1,1)); boxplot(t(DIFF),cex.axis=0.8,las=2,col=com.col,ylab="GDDobs-GDDmod")
  abline(h=seq(-500,500,100),lty=2)
  abline(h=0,lwd=2)
  abline(v=1:69,lty=3)
  
  # OLD STUFF Comparing Freezing Degree Days (FDD) ----  
  TMPs <- read.csv("~/CLIMATO/DATA_ANALYSIS/BIOCLIM_SOIL_TEMP_S_FDDjfm.csv",row.names=1)
  TMPh <- read.csv("~/CLIMATO/DATA_ANALYSIS/BIOCLIM_SOIL_TEMP_H_FDDjfm.csv",row.names=1)
  
  SITc  <- intersect(rownames(TMPs),rownames(TMPh)) # 69 sites communs
  FDDs  <- TMPs[match(SITc,rownames(TMPs)),]
  FDDh  <- TMPh[match(SITc,rownames(TMPh)),-c(1,21)] # no simul for 1999 and 2019
  
  graphics.off(); windows(20,10); plot(unlist(GDDh),unlist(GDDs))
  abline(0,1)
  
  DIFF <- (FDDh)-(FDDs)
  NONA <- !is.na(apply(DIFF,1,mean,na.rm=T))
  com.tmp <- 
  graphics.off();windows(20,10); par(mar=c(6,5,1,1)); boxplot(t(DIFF),cex.axis=0.8,las=2,col=com.col,ylab="FDDobs-FDDmod")
  abline(h=seq(-500,500,100),lty=2)
  abline(h=0,lwd=2)
  abline(v=1:69,lty=3)
  
  

  
  # OLD STUFF PLOTS ----
  
    # Area Of Interest (drawn on GoogleEarth) -> AOIext
  tmp             <- readOGR("THAB.kml")
  AOI             <- spTransform(tmp, CRSobj=LAMB93) # fine
  AOIexttmp       <- extent(AOI)
  AOI@proj4string <- LAMB93
  
  AOIl2e          <- spTransform(tmp, CRSobj=LAMBIIe) # fine
  extent(AOIl2e)
  
  # #adjust extent to ORTHO grid - en multple de 5000 m
  AOIext    <- extent(950000,990000,6430000,6465000)
  AOIextl2e <- extent(900000,945000,2000000,2035000)
  
  # Suivi territoire - CBNA N. Fort -> CAN
  CAN <- read.csv("C:/Users/cholerp/Documents/PROJETS/SNOWBED/DATA/CBNA_SUIVI_HABITAT/sites.csv")
  CAN_COO<- CAN[,4:5] # LAMB93
  graphics.off();windows(10,10);plot(CAN_COO,pch="+",lwd=2,col="red")
  plot(AOIexttmp,add=T)

  # Roche Noire for checking
  RNse  <- raster::extent(292200, 298700, 4990000, 4994200)
  RN10m<-raster(RNse, 
                nrows=(RNse@ymax-RNse@ymin)/10, 
                ncols=(RNse@xmax-RNse@xmin)/10, 
                crs=UTM32)
  RN10_LAMB93 <- projectExtent(RN10m,LAMB93)
  RNse  <- extent(RN10_LAMB93)
  
    
  # OLD STUFF CONTOURS ADMIN & HYDRO ----
  # Departments
  DPT <-readOGR("C:/Users/cholerp/Documents/SIG/DPTS/shp//dep_france_dom_region.shp"); DPT@proj4string<-LAMBIIe
  DPT <- spTransform(DPT, CRSobj=LAMB93) # fine
  D05 <-DPT[which(DPT@data[,8]=='HAUTES-ALPES'),] # extract Is?re
  D73 <-DPT[which(DPT@data[,8]=='SAVOIE'),] # extract Is?re
  
  # Parcs nationaux
  PN<-readOGR("C:\\Users\\cholerp\\Documents\\SIG\\INPN\\pn2010.shp")
  PN@proj4string<-LAMB93
  
  PNE   <- PN[grep('Ecrins',PN@data$NOM),]  # coeur et zone d'adh?sion
  PNV   <- PN[grep('Vanoise',PN@data$NOM),] # coeur et zone d'adh?sion
  
  # secteur Lautaret
  LAUT <- extent(c(938000,998000,6428000,6478000))
  
  # OLD STUFF DEM for POI ----
  DAHt  <- DAH*POI
  TPIt  <- TPI25*POI
  DEMt  <- DEM25*POI
  SHEt  <- SHELT25*POI
  
    # Aggregate at 10 m
  graphics.off(); windows(10,5); plot(raster("TPI2M.tif"))
  
  # OLD DEM version of SINTEGRA at 10m - not sourced ! ----
  RAW<-read.table("C:\\Users\\cholerp\\Documents\\SIG\\DEM\\LAUTARET10m\\DEM10m.xyz")
  RNdem10<-rasterFromXYZ(RAW,crs=LAMBIIe, digits=5)
  # image(dem10,axes=T)
  # fait penser ? un d?calage syst?matique selon l'axe des x de 1 pixel 50m
  
  # reprojection en utm32
  TMP <-projectExtent(RNdem10,UTM32)      # initialize a raster with boudary in utm32
  dem10.utm32<- projectRaster(RNdem10,TMP) # quite fast
  dem10      <- crop(dem10.utm32,RNse)
  dem10      <- extend(dem10,RNse) 
  # expand  replaced by extend to avoid a name conflict with Matrix package
  dem10      <- crop(dem10,RNse)
  #plot(dem10,axes=T)
  #points(mask0[,1],mask0[,2],type='l',lwd=2)
  
  # interpolate sur une grille de 10 x 10
  RNdem10s   <- raster::resample(dem10,RN10m,method='bilinear')
  writeRaster(RNdem10s,'RNdem10s.tif',format="GTiff",overwrite=T)
  
  # OLD STUFF deux alternatives pour calculer un max ----
  NDVImaxspot1           <- stackApply(stack(LF),1,fun=max)
  NDVImaxspot1@crs       <- LAMB93
  NDVImaxspotaoi1        <- projectRaster(NDVImaxspot1,REF25)
  NDVImaxspotl2e1        <- projectRaster(NDVImaxspot1,REF25l2e)  
  writeRaster(NDVImaxspotaoi1,"NDVImaxspotaoi1.tif",format="GTiff",overwrite=T)  
  
  NDVImaxspot2           <- stackApply(stack(LF[7:9]),1,fun=mean,na.rm=T)
  NDVImaxspot2@crs       <- LAMB93
  NDVImaxspotaoi2        <- projectRaster(NDVImaxspot2,REF25)
  NDVImaxspotl2e2        <- projectRaster(NDVImaxspot2,REF25l2e)
  writeRaster(NDVImaxspotaoi2,"NDVImaxspotaoi2.tif",format="GTiff",overwrite=T)

  FSFD      <- raster("C:/Users/cholerp/Documents/SIG/SPOT/SPOT4TAKE5_2015/FSFD/FSFD_median_2015.tif")
  FSFD@crs  <- LAMB93
  FSFDaoi   <- projectRaster(FSFD,REF25)
  FSFDl2e   <- projectRaster(FSFD,REF25l2e)
  plot(FSFDaoi)
  
  MASKfsfd <- FSFDaoi; MASKfsfd[!is.na(FSFDaoi[])]<-1
  plot(MASKfsfd)
  (EXT.FSFD <- extent(FSFD)) # 25m resolution
  graphics.off(); windows(20,20); plot(FSFD)
  
  # Final raster mask -> AOIras
  AOIras <- MASKaoi*MASKfsfd
  plot(AOIras)
  
  
  
  FSFDt <- FSFDaoi*POI 
  graphics.off(); windows(20,20); plot(FSFDt)
  graphics.off(); windows(10,8) ; plot(crop(FSFDt,RNse),col=SNOWpal(50))
  
  # Massifs SAFRAN
  tmp <- readOGR("C:/Users/cholerp/Documents/SIG/CODAGE_SAFRAN/massifs_WGS84.shp")
  MASSIFS <- spTransform(tmp, CRSobj=LAMB93) # fine
  
  # OLD STUFF NVDImax from SPOT4TAKE5 ----
  DIR.SPOT.NDVI <- "C:/Users/cholerp/Documents/SIG/SPOT/SPOT4TAKE5_2015/NDVI_10m/"
  LF <- list.files(DIR.SPOT.NDVI, full=T)
  DIR.SPOT.NDSI <- "C:/Users/cholerp/Documents/SIG/SPOT/SPOT4TAKE5_2015/NDSI_10m/"
  LFS <- list.files(DIR.SPOT.NDSI, full=T)
    
  
  # OLD GRASS-derived topographic variables ----
  # SUIVRE strictement l'ordre !!
  
  execGRASS("r.in.gdal", flags="o", parameters=list(input="C:/Users/cholerp/Documents/PROJETS/SNOWBED/DATA_ANALYSIS/DEM10.tif", output="DEM10M",location="SNOB")) # create the location using parameters from .tif
  
  initGRASS("C:/GRASS GIS 7.4.1", gisDbase="C:/Users/cholerp/Documents/SIG/GRASS_DB",location="SNOB10Mn",mapset="DEM",override=TRUE) # no parameters
  gmeta() 
  
  # BEWARE : then need to change location ! and it works !
  execGRASS("r.slope.aspect",parameters=list(elevation="DEM10M",slope="SLO10M", aspect="ASP10M", pcurvature="PCUR10M",tcurvature="TCUR10M"))  # aspect, slope and curvature
  execGRASS("r.topidx",flags="overwrite",parameters=list(input="DEM10M",output="TPI10M"))
  execGRASS("r.flow",flags="overwrite",parameters=list(elevation="DEM10M",skip=as.integer(5),flowline="FLOWLINE10M",flowlength="FLOWLEN10M", flowaccumulation="FLOWACCU10M")) # flowlines
  
  # export to tiff
  TMP   <- raster(readRAST("SLO10M", ignore.stderr = T))
  writeRaster(TMP,"SLO10M.tif",overwrite=T,format="GTiff")
  
  TMP   <- raster(readRAST("ASP10M", ignore.stderr = T))
  writeRaster(TMP,"ASP10M.tif",overwrite=T,format="GTiff")
  
  TMP   <- raster(readRAST("PCUR10M", ignore.stderr = T))
  writeRaster(TMP,"PCUR10M.tif",overwrite=T,format="GTiff")
  
  TMP   <- raster(readRAST("TCUR10M", ignore.stderr = T))
  writeRaster(TMP,"TCUR10M.tif",overwrite=T,format="GTiff")
  
  TMP   <- raster(readRAST("TPI10M", ignore.stderr = T))
  writeRaster(TMP,"TPI10M.tif",overwrite=T,format="GTiff")
  
  
  # OLD DEMs ----
  # DEM at 1m fromBD_RGE aggregated -> DEM10.tif
  
  setwd("D:/BD_ALTI_SNOWBED/1m")
  
  TILE <- read.csv("./dalles73.csv")[,1]
  xt   <- as.numeric(substr(TILE,14,16))
  yt   <- as.numeric(substr(TILE,18,21))
  
  TAR <- sort(unique(as.vector(TILE[which(xt>=950 & xt<990 & yt>=6430 & yt<6465)])))
  
  for (i in 1:length(TAR)){
    print(i)
    NAME <- paste(TAR[i],".asc",sep="")
    GET(paste("ftp://osug@theresa.grenoble.cemagref.fr/BD_ALTI/1m/RGEALTI_MNT_1M_ASC_73/",NAME,sep=""), authenticate("osug", "#ign@osug38!"),write_disk(NAME,overwrite=T))
  }
  
  TILEf <- list.files("D:/BD_ALTI_SNOWBED/2M_TIF",full=T)
  TILE <- list.files("D:/BD_ALTI_SNOWBED/2M_TIF")
  xt   <- as.numeric(substr(TILE,14,16))
  yt   <- as.numeric(substr(TILE,18,21))
  
  TAR   <- TILEf[which(xt>=950 & xt<970 & yt>=6440 & yt<6465)] # NORTH-WEST
  TAR   <- TILEf[which(xt>=970 & xt<990 & yt>=6430 & yt<6465)] # NORTH-EAST
  
  setMethod('mosaic', signature(x='list', y='missing'), 
            function(x, y, fun, tolerance=0.05, filename=""){
              stopifnot(missing(y))
              args <- x
              if (!missing(fun)) args$fun <- fun
              if (!missing(tolerance)) args$tolerance<- tolerance
              if (!missing(filename)) args$filename<- filename
              do.call(mosaic, args)
            })
  FUN <- function(x) mean(x,na.rm=T)
  imgl  <- lapply(TAR, raster)
  # y     <- do.call(mosaic, args=list(x=imgl,fun=FUN,tolerance=0.5))  # issue
  y     <- do.call(mosaic, args=list(x=imgl,fun=mean,tolerance=0.5)) # bizarre values produced
  plot(y,colNA="black",col=DEM_pal2(100))
  writeRaster(round(y),"DEM_2m_EAST.tif",format="GTiff",overwrite=T)
  graphics.off();plot(raster("DEM_2m_EAST.tif"),col=DEM_pal2(100))
  graphics.off();plot(raster("DEM_2m_WEST.tif"),col=DEM_pal2(100))  
  
  # not to run
  LF  <- list.files("E:/BD_ALTI_SNOWBED/1m/",pattern=".asc",full=T)
  LFs <- list.files("E:/BD_ALTI_SNOWBED/1m/",pattern=".asc",full=F)
  
  setwd("E:/BD_ALTI_SNOWBED/2m_TIF/")
  for (i in 1:length(LF)){
    print(i)
    NAME <- gsub(".asc","_2m.tif",LFs[i])
    writeRaster(round(aggregate(raster(LF[i]),2,mean)),NAME,format="GTiff",overwrite=T)
  }
  
  WEST <- raster("DEM_2m_WEST.tif")
  EAST <- raster("DEM_2m_EAST.tif")
  
  DEM10     <- mosaic(WEST10 <- aggregate(WEST,fact=5,fun=mean),EAST10 <- aggregate(EAST,fact=5,fun=mean),fun=mean)
  DEM10@crs <- LAMB93
  
  writeRaster(round(DEM10),"DEM10.tif",format="GTiff",overwrite=T)
  
  # OLD DEM AUTRES ----
  # TPI from spatialEco
  DEM10 <- raster("DEM10.tif")
  
  TPIs   <- tpi(DEM10, scale = 49, win = "rectangle", normalize = FALSE)
  TPIs@crs<-LAMB93
  TPI25   <- projectRaster(TPIs,REF25)
  TPIs49t  <- TPI25*POI
  
  graphics.off();
  windows(10,8); plot(crop(TPIs5t,RNse),zlim=c(-6,6))
  windows(10,8); plot(crop(TPIs7t,RNse),zlim=c(-15,6))  
  windows(10,8); plot(crop(TPIs49t,RNse),zlim=c(-50,50)) 


# PART 2. OLD STUFF REMOTE SENSING 
  # OLD STUFF LANDSAT ----
  DIR.TEST        <- "K:/fyse/ALPAGES_SENTINELLES/ERL/TEST/"
  DIR.INPUT       <- "K:/fyse/ALPAGES_SENTINELLES/ERL/LANDSAT/L8/RAW/"
  DIR.OUTPUT      <- "K:/fyse/ALPAGES_SENTINELLES/ERL/LANDSAT/L8/BANDS/"
  # 1. Preparation des images LANDSAT 8 ? 30 m
  # voir le script landsat.r syr K:/fyse/alpages_sentinelles
  
  # Extraction , reprojection en lambert93 et correction de toutes les bandes des sc?nes Landsat intersectant l'arc alpin fran?ais
  
  DEM_30       <- raster("K:/fyse/ALPAGES_SENTINELLES/ERL/LANDSAT/DEMforLANDSAT.tif")
  dem_slopeasp <- slopeasp(as(DEM_30,  'SpatialGridDataFrame'))
  
  # rappel des bandes landsat
  # B1 = Blue
  # B2 = Green
  # B3 = Red
  # B4 = NIR
  # B5 = SWIR 1 (MIR)
  # B6 = Thermal IR with less saturation & lower sensitivity, 
  # http://landsat.usgs.gov/thermal_band_on_Landsat_7.php
  
  # select image
  
  # mi-juillet 2015 SCENE 195029
  DIR.INPUT <- "K:/fyse/ALPAGES_SENTINELLES/ERL/LANDSAT/L8/BANDS/"
  IMG.SEL   <- "LC81950292015194LGN00" # probablement trop tardive
  
  IMG.SEL   <- "LC81950292013172LGN00" # fin juin - pas mal
  B2        <- raster(paste(DIR.INPUT,IMG.SEL,"/",IMG.SEL,"_b2.tif",sep=""))
  B3        <- raster(paste(DIR.INPUT,IMG.SEL,"/",IMG.SEL,"_b3.tif",sep=""))
  B4        <- raster(paste(DIR.INPUT,IMG.SEL,"/",IMG.SEL,"_b4.tif",sep=""))
  B5        <- raster(paste(DIR.INPUT,IMG.SEL,"/",IMG.SEL,"_b5.tif",sep=""))
  B6        <- raster(paste(DIR.INPUT,IMG.SEL,"/",IMG.SEL,"_b6.tif",sep=""))
  
  NDVI      <- (B5 - B4) / (B4 + B5)
  NDSI      <- (B3 - B6) / (B3 + B6)
  
  
  # Selecting d'une zone test
  NDSIlaut <- crop(NDSI,LAUT)
  NDVIlaut <- crop(NDVI,LAUT)
  
  graphics.off()
  windows(10,20)
  par(mfrow=c(2,1))
  plot(NDSIlaut)
  plot(NDVIlaut)
  
  hist(NDVI05[])
  
  # seuillage NDSI
  THR        <-  0.6
  NDSIlaut.b <- reclassify(NDSIlaut, matrix(c(-Inf,THR,0,THR,Inf,1), ncol=3, byrow=TRUE), progress='text')
  
  # OLD STUFF SPOT4TAKE5 ----
  TMP     <- extent(869882.660460, 1016779.620920, 6385111.891312, 6508344.471111)
  TMP     <- c(TMP@xmin,TMP@xmax,TMP@ymin,TMP@ymax)/1000
  EXTspot <- extent(1000*c(floor(TMP[1]),ceiling(TMP[2]),floor(TMP[3]),ceiling(TMP[4])))
  # Use available raster from SPOT4TAKE5
  tmp <- raster("C:/Users/cholerp/Documents/SIG/SPOT4TAKE5_2015/NDVI_10m/SPOT_NDVI_2015-04-11.tif")
  EXTspotrast <- extent(tmp)
  
  # Images 2015 produced by Brad in the framework of SPOT4take5 project
  # Raster files are in LAMB93
  # Biblio : Dedieu & al. 2017 Remote Sensing
  
  ARA.EXT.300  <- extent(966900,967200,6446500,6446800)  # 300 m par 300 m
  ARA.EXT.400  <- extent(966850,967250,6446400,6446800)  # 400 m par 400 m
  
  RNse         <- raster::extent(292200, 298700, 4990000, 4994200)  # NEW in UTM
  RATIO        <- (RNse@ymax-RNse@ymin)/(RNse@xmax-RNse@xmin)
  RN10m<-raster(RNse, nrows=(RNse@ymax-RNse@ymin)/10, ncols=(RNse@xmax-RNse@xmin)/10, crs=UTM32)
  RN10_LAMB93  <- projectExtent(RN10m,LAMB93)
  RON.EXT      <- extent(RN10_LAMB93)
  
  FSFD   <- extract(FSFDrast,SIT_TOT[,3:4])
  NDVI04 <- extract(raster(LF[1]),SIT_TOT[,3:4])  
  NDVI05 <- extract(raster(LF[2]),SIT_TOT[,3:4])
  NDVI06 <- extract(raster(LF[3]),SIT_TOT[,3:4])
  
  TESTB  <- max(raster(LF[7]),raster(LF[8]),na.rm=T,progress='text')
  NDVI07 <- extract(TESTB,SIT_TOT[,3:4])
  
  TESTC  <- max(raster(LF[9]),raster(LF[10]),na.rm=T,progress='text')
  NDVI08 <- extract(TESTC,SIT_TOT[,3:4])
  
  NDSI04 <- extract(raster(LFS[1]),SIT_TOT[,3:4])  
  NDSI05 <- extract(raster(LFS[2]),SIT_TOT[,3:4])
  NDSI06 <- extract(raster(LFS[3]),SIT_TOT[,3:4])
  
  TESTSB  <- max(raster(LFS[7]),raster(LFS[8]),na.rm=T,progress='text')
  NDSI07 <- extract(TESTSB,SIT_TOT[,3:4])
  
  # TESTD  <- max(raster(LFS[10]),raster(LFS[11]),na.rm=T,progress='text')
  # NDSI08 <- extract(TESTD,SIT_TOT[,3:4])
  
  plot(FSFD,NDVI07)
  
  NC <- 12 # number of clusters
  LEGEND <- c( # to be changed each run !
    "1  : m?gaphorbaie",
    "7  : subalpin thermique tendance x?rophile Patzkea",
    "3  : subalpin thermique plus frais",
    "12  : Nardion - subalpin sup?rieur tendance nivale",
    "6  : for?t subalpine - m?aphorbaie",
    "1  : ?boulis alpin nival mat. fin",
    "10  : ?boulis alpin nival gros blocs",
    "8  : montagnard-subalpin - prairie fertile dense",
    "9  : montagnard-subalpin - pelouse thermique", 
    "10 : lande",
    "8 : alpin thermique",
    "2 : alpin nival dense"
  )
  
  COL <- rep("white",12);COL[4] <- "green"; COL[6] <- "blue"; COL[7] <- "orange"; COL[12]<-"lightblue"; COL[10] <- "brown"; COL[5] <- "gray" ; COL[2] <- "blue"; COL[c(1,10)]<-"blue"; COL[9] <- "red"; COL[11] <- "magenta"; 
  
  # Predictive maps of snowbeds
  # Condition #1. Under snow or right after snow on 5 juin (img 3)
  THR1  <- 0.05 # consistent with NDVI values in early June
  COND1 <- reclassify(raster(LF[3]), matrix(c(-1,THR1,1,THR1,1,0),ncol=3,byrow=T),progress='text')
  
  # Conditon #2 - peak standing NDVI>=0.5 to get more values 
  # Le faire sur une composite de fin de saison (end of August) : TESTC
  THR2  <- 0.5
  COND2 <- reclassify(TESTC, matrix(c(-1,THR2,0,THR2,1,1),ncol=3,byrow=T),progress='text')
  
  SNOWBET <- COND1*COND2
  plot(SNOWBET)
  
  writeRaster(SNOWBET,"SNOWBET_02_05.tif",format="GTiff")  # 363 282 cells / 
  writeRaster(SNOWBET,"SNOWBET_005_05.tif",format="GTiff") # 148 108 cells / 14 982 cel. de 100 m
  length(which(SNOWBET[]==1))                   
  length(which(aggregate(SNOWBET,10,max)[]==1))  # 29 805 cellules
  length(which(aggregate(SNOWBET,10,sum)[]>=60)) # 194 cellules
  
  # OLD STUFF HR SNOWCOVER EXTENT (sentinel + Landsat) ----
  # Products provided by Simon Gascoin on 07.2019
  # Number of days with snow per month (2016-2018)
  # Need to assemble two tiles (31TGK & 31TGL) for the AOI
  # produits en UTM31 (resolution 30m) to be reprojected at 25M
  
  tmp1    <- raster("C:/Users/cholerp/Documents/SIG/SNOW/SNOW_OCCURENCE/SNOW_OCCURENCE_T31TGK_20160701_20160731.tif")
  tmp2    <- raster("C:/Users/cholerp/Documents/SIG/SNOW/SNOW_OCCURENCE/SNOW_OCCURENCE_T31TGL_20160701_20160731.tif")
  SCE0716 <- merge(projectRaster(tmp1,REF25),projectRaster(tmp2,REF25))
  
  tmp1    <- raster("C:/Users/cholerp/Documents/SIG/SNOW/SNOW_OCCURENCE/SNOW_OCCURENCE_T31TGK_20170701_20170731.tif")
  tmp2    <- raster("C:/Users/cholerp/Documents/SIG/SNOW/SNOW_OCCURENCE/SNOW_OCCURENCE_T31TGL_20170701_20170731.tif")
  SCE0717 <- merge(projectRaster(tmp1,REF25),projectRaster(tmp2,REF25))
  
  tmp1    <- raster("C:/Users/cholerp/Documents/SIG/SNOW/SNOW_OCCURENCE/SNOW_OCCURENCE_T31TGK_20180701_20180731.tif")
  tmp2    <- raster("C:/Users/cholerp/Documents/SIG/SNOW/SNOW_OCCURENCE/SNOW_OCCURENCE_T31TGL_20180701_20180731.tif")
  SCE0718 <- merge(projectRaster(tmp1,REF25),projectRaster(tmp2,REF25))
  
  SCE07   <- mean(SCE0716,SCE0717,SCE0718,na.rm=T)
  writeRaster(SCE07,"SNOW_OCCU_07.tif",format="GTiff",overwrite=T)
  
  tmp1    <- raster("C:/Users/cholerp/Documents/SIG/SNOW/SNOW_OCCURENCE/SNOW_OCCURENCE_T31TGK_20160601_20160630.tif")
  tmp2    <- raster("C:/Users/cholerp/Documents/SIG/SNOW/SNOW_OCCURENCE/SNOW_OCCURENCE_T31TGL_20160601_20160630.tif")
  SCE0616 <- merge(projectRaster(tmp1,REF25),projectRaster(tmp2,REF25))
  
  tmp1    <- raster("C:/Users/cholerp/Documents/SIG/SNOW/SNOW_OCCURENCE/SNOW_OCCURENCE_T31TGK_20170601_20170630.tif")
  tmp2    <- raster("C:/Users/cholerp/Documents/SIG/SNOW/SNOW_OCCURENCE/SNOW_OCCURENCE_T31TGL_20170601_20170630.tif")
  SCE0617 <- merge(projectRaster(tmp1,REF25),projectRaster(tmp2,REF25))
  
  tmp1    <- raster("C:/Users/cholerp/Documents/SIG/SNOW/SNOW_OCCURENCE/SNOW_OCCURENCE_T31TGK_20180601_20180630.tif")
  tmp2    <- raster("C:/Users/cholerp/Documents/SIG/SNOW/SNOW_OCCURENCE/SNOW_OCCURENCE_T31TGL_20180601_20180630.tif")
  SCE0618 <- merge(projectRaster(tmp1,REF25),projectRaster(tmp2,REF25))
  
  SCE06   <- mean(SCE0616,SCE0617,SCE0618,na.rm=T)
  plot(SCE06)
  writeRaster(SCE06,"SNOW_OCCU_06.tif",format="GTiff",overwrite=T)
  
  
  

  
  # OLD STUFF BD ORTHO ----
  
  # UPDATE : la BD ortho est sur :
  DIR_ORTHO_IGN <- "K:/macroeco/GIS_DATA/Alpes/IGN/BD-ORTHO/D005/1_DONNEES_LIVRAISON_2010-10-00066/BDO_RVB_0M50_ECW_LAMB93_D05-ED03"   
  
  # Les archives de 1999 (en L2e), 2003 et 2009 sont sur :
  # serveur  : theresa.grenoble.cemagref.fr
  # login    : osug
  # pass     : #ign@osug38!
  
  # SOLVED ISSUE - direct download use httr -> see DEM at 1m
  
  # En local sur disque dur externe WD7
  # Images source de la base de donn?es ORTHOde l'IGN
  # le pavage est de 5x5km
  # Exemple de nom pour 2003 "05-2003-0980-6450-LA93.ecw", 
  
  # 2003 et 2009 2013 : convert in .tif & aggregate at 2m for a first look
  
  ORTHO        <- "D:/BD_ORTHO_SNOWBED/2013/"
  LF13         <- list.files(ORTHO,full=T,pattern=".ecw")
  for(i in 1:length(LF03)){
    INFILE     <- LF13[i]
    shell(cmd=paste('gdalinfo',INFILE,sep=' '))  
    OUTFILE    <- gsub(".ecw","_2m.tif", INFILE)
    shell(cmd=paste('gdalwarp -of GTiff -tr 2 2 ',INFILE,' ',OUTFILE,sep=''))
  }
  
  ORTHO        <- "D:/BD_ORTHO_SNOWBED/2003/"
  LF03         <- list.files(ORTHO,full=T,pattern=".ecw")
  for(i in 1:length(LF03)){
    INFILE     <- LF03[i]
    shell(cmd=paste('gdalinfo',INFILE,sep=' '))  
    OUTFILE     <- gsub(".ecw","_2m.tif", INFILE)
    shell(cmd=paste('gdalwarp -of GTiff -tr 2 2 ',INFILE,' ',OUTFILE,sep=''))
  }
  
  ORTHO       <- "D:/BD_ORTHO_SNOWBED/2009/"
  LF09        <- list.files(ORTHO,full=T,pattern=".ecw")
  for(i in 1:length(LF09)){
    INFILE     <- LF09[i]
    shell(cmd=paste('gdalinfo',INFILE,sep=' '))  
    OUTFILE     <- gsub(".ecw","_2m.tif", INFILE)
    shell(cmd=paste('gdalwarp -of GTiff -tr 2 2 ',INFILE,' ',OUTFILE,sep=''))
  }
  
  # 1999 : convert in .tif & lamb93 & aggregate at 2m for a first look
  ORTHO        <- "D:/BD_ORTHO_SNOWBED/1999/"
  LF99         <- list.files(ORTHO,full=T,pattern=".ecw")
  
  for(i in 1:length(LF99)){
    INFILE     <- LF99[i]
    shell(cmd=paste('gdalinfo',INFILE,sep=' '))  
    OUTFILE    <- gsub("LA2E-C10.ecw","LA93-C10_2m.tif", INFILE)
    shell(cmd=paste('gdalwarp -t_srs "+proj=lcc +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +lon_0=3 +x_0=700000.0 +y_0=6600000 +lat_0=46.5 +lat_1=49 +lat_2=44 +units=m +no_defs" -r bilinear -co INTERLEAVE=BAND -of GTiff -tr 2 2 ',INFILE,' ',OUTFILE,sep=''))
    shell(cmd=paste('gdalwarp -of GTiff -tr 2 2 ',INFILE,' ',OUTFILE,sep=''))
  } 
  
  # Mosaic & fenetrage selon les grilles lambert93 
  TEST  <- list.files(ORTHO,full=T,pattern="_2m.tif")
  setMethod('mosaic', signature(x='list', y='missing'), 
            function(x, y, fun, tolerance=0.05, filename=""){
              stopifnot(missing(y))
              args <- x
              if (!missing(fun)) args$fun <- fun
              if (!missing(tolerance)) args$tolerance<- tolerance
              if (!missing(filename)) args$filename<- filename
              do.call(mosaic, args)
            })
  
  imgl  <- lapply(TEST, brick)
  y     <- do.call(mosaic, args=list(x=imgl,fun=max,tolerance=0.5))
  plotRGB(y)
  
  LFtif09s <- list.files("D:/BD_ORTHO_SNOWBED/2009/",pattern="_2m.tif",full=F)
  ORTHO        <- "D:/BD_ORTHO_SNOWBED/1999/"
  for (i in 1:length(LF09)){
    print(i)
    NAME     <- gsub("-2009-","-1999-",LFtif09s[i])
    OUTFILE  <- paste(ORTHO,NAME,sep="")
    IMG      <- projectRaster(y,stack(LFtif09[i]))
    IMG      <- round(IMG)
    writeRaster(IMG,OUTFILE,format="GTiff",overwrite=T)
  }
  
  
  # synthese
  ORTHO           <- "D:/BD_ORTHO_SNOWBED/1999/"
  LFtif99         <- list.files(ORTHO,full=T,pattern="_2m.tif")
  ORTHO           <- "D:/BD_ORTHO_SNOWBED/2009/"
  LFtif09         <- list.files(ORTHO,full=T,pattern="_2m.tif")
  ORTHO           <- "D:/BD_ORTHO_SNOWBED/2003/"
  LFtif03         <- list.files(ORTHO,full=T,pattern="_2m.tif")
  
  
  # OLD STUFF THEIA-OSO ----
  # Nomenclature :culture ete:11; culture hiver:12; foret feuillus:31; foret coniferes:32; pelouses (meadows):34; landes ligneuses:36; urbain dense:41; urbain diffus:42; zones ind et com:43; surfaces routes:44; surfaces minerales:45; plages et dunes:46; eau:51; glaciers ou neige: 53; prairies (grasslands) : 211; vergers:221; vignes:222
  
  TMP    <- raster("C:/Users/cholerp/Documents/SIG/LAND_COVER/THEIA_LC2016_FR/LC_theia.tif") 
  SNOBoso <- projectRaster(TMP,REF25,method="ngb")
  TABoso <- table(SNOBoso[])
  
  OSOcol1   <-  c(rep("yellow",2),rep("darkgreen",2),"lightgreen","brown",rep("black",4),rep("gray",2),"blue","white","green")
  df1  <- data.frame(id=as.numeric(names(TABoso)),v=1:15)
  SNOBoso1  <- subs(SNOBoso, df1)
  plot(SNOBoso1,col=OSOcol1) 
  
  OSOcol2   <-  c("yellow","darkgreen","green","brown","black","gray","blue","white")
  COL <- adjustcolor(OSOcol2, alpha.f = 0.3) 
  df2 <- data.frame(id=as.numeric(names(TABoso)),v=c(rep(1,2),rep(2,2),3,4,rep(5,4),rep(6,2),7,8,3)) # with grouping of land cover classes
  SNOBoso2  <- subs(SNOBoso, df2)
  plot(SNOBoso2,col=COL) # works fine
  
  windows(10,8); plot(crop(SNOBoso2*POI,RNse),col=rainbow(12))
  windows(10,8); plot(crop(POI,RNse),col=rainbow(12))
  windows(10,8); plot(crop(SNOBoso2,RNse),col=rainbow(8))
  
  # Reclassify using Points Of Interest -> POI
  POItmp  <- SNOBoso2*MASKfsfd*MASKaoi
  graphics.off(); windows(10,8) ; plot(crop(SNOBoso2,RNse),col=rainbow(8))
  graphics.off(); windows(10,8) ; plot(crop(MASKaoi,RNse),col=rainbow(8))
  graphics.off(); windows(10,8) ; plot(crop(MASKfsfd,RNse),col=rainbow(8))
  df      <- data.frame(id=1:8, v=c(NA,NA,1,NA,NA,1,1,NA))
  POI     <- subs(POItmp,df)
  
  
  