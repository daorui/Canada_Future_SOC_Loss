
# 
#climate variables
unixtools::set.tempdir('/mnt/Fastrun/temp4r')
library(terra)

library(lubridate)

library(sf)

library(dplyr)

library(spData)

library(zoo) # rollmean

library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)
library(readxl)
library(hrbrthemes)
library(viridis)
library(PupillometryR) 
library(scales)
library(ggpmisc)

t1 <- list.files(pattern = 'nc','/mnt/DataSpace/Data_pool/CMIP6_rebuild',full.names = T)
var0 = 'tas'
model0='CanESM'
time0= 'historical'
#/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/001code/00CMIP2300_Preprocessing.R
model_2100 <- c("ACCESS-CM2", "CanESM5", "CMCC-ESM2", "EC-Earth3-Veg", "FIO-ESM-2-0", "GISS-E2-1-G", "GISS-E2-1-H", "HadGEM3-GC31-LL", "INM-CM4-8", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "UKESM1-0-LL", "ACCESS-ESM1-5", "CanESM5-CanOE", "CNRM-CM6-1", "CNRM-CM6-1-HR", "CNRM-ESM2-1", "EC-Earth3-Veg-LR", "INM-CM5-0")
model_gpp <- c("EC-Earth3-Veg-LR", "IPSL-CM6A-LR", "ACCESS-ESM1-5", "BCC-CSM2-MR",
               "CanESM5", "CanESM5-1", "CAS-ESM2-0", "CESM2-WACCM", 
               "CMCC-CM2-SR5", "EC-Earth3-CC", "INM-CM4-8", "INM-CM5-0",
               "TaiESM1", "EC-Earth3-Veg", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR")
model_2100_gpp <- setdiff(model_gpp,model_2100)[-c(4,5,7)]

model_2100_new <- c(model_2100[-c(6,7,8,15,17,18,19,20)],model_2100_gpp)

future_period <- data.frame(stardyear=c(1970,seq(2001,2100,10)),endyear=c(2000,seq(2010,2100,10)))

# calculate annual mean value of outputs 
models_output <- function (var0 , model0 ,time0 ) {
  ind00 <- grep(pattern = paste0('.*',var0,'_.*',model0,'_.*','historical'),x=t1)
  ind01 <- grep(pattern = paste0('.*',var0,'_.*',model0,'_.*',time0),x=t1) 
  ind00 <- unique(c(ind00,ind01)) 
  f2 <- rast(t1[ind00])
  
  times_span0 <- time(f2)
  #
  annual_month_average <- list()
  if(time0 == 'historical'){m=1} else {m=2:nrow(future_period)}
  for (k in m) {
    startdate <- future_period[k,1]
    enddate <- future_period[k,2]
    if(year(times_span0)[1] <= startdate) {  
      cal_year_ind <- which(year(times_span0) %in% startdate:enddate) 
      f3 <- f2[[cal_year_ind]] 
      mon_days <- list()
      for (i in startdate:enddate){
        all_doy <- seq(as.Date(paste0(i,'-1-1')),as.Date(paste0(i,'-12-31')),1)
        ind000 <- which(startdate:enddate==i)
        mon_days[[ind000]] <- data.frame(table(month(all_doy)))
      }
      mon_days_ind <- do.call(rbind,mon_days)[,2]
      #
      if(var0 == 'pr') {f40 <- f3*86400*mon_days_ind + 0.1} else {f40 <- f3 - 273.15}
      all_ind <- which(1:nrow(future_period) ==k)
      annual_month_average[[all_ind]] <- tapp(x = f40,index = 1:12,fun = 'mean')  ##median 20230726 change to mean:monthly average value for 1970-2000
    }
    models_output0 <- do.call(c,annual_month_average)
  }
  return(models_output0)
}

#
time <- c('historical','ssp245','ssp585')
var00 = c('pr','tas','tasmin','tasmax')
all_combin <- expand.grid(model_2100_new,var00,time)
#paste(all_combin[1,],sep="_")
#unlink("/mnt/DataSpace/Data_pool/CMIP6_rebuild/tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_196601-196612.nc")
#unlink("/mnt/DataSpace/Data_pool/CMIP6_rebuild/tas_ImonAnt_IPSL-CM6A-LR_historical_r1i1p1f1_gra_185001-201412.nc")
#unlink("/mnt/DataSpace/Data_pool/CMIP6_rebuild/tas_ImonGre_IPSL-CM6A-LR_historical_r1i1p1f1_grg_185001-201412.nc")
ii=146

for (ii in 1:nrow(all_combin)) {
  ind0000 <- grep(pattern = paste0('.*',all_combin[ii,2],'_.*',all_combin[ii,1],'_.*',all_combin[ii,3]), x = t1)
  f20 <- rast(t1[ind0000])
  times_span00 <- year(time(f20))
  assign(do.call(paste,c(all_combin[ii,],min(times_span00),max(times_span00),sep='_')) , 
         models_output(var0 = all_combin[ii,2], model0 = all_combin[ii,1],time0 = all_combin[ii,3]))
}

#pr ssp245 ACCESS-CM2
#
model <- model_2100_new
time00 <- c('historical','ssp245','ssp585')
var00=c('pr','tas','tasmin','tasmax')
# 
for (i in  1:length(model)){ 
  for(j in 1:4) {
    for (k in 2:3){
      hist00 <- grep(pattern = paste0(model[i],"_",var00[j],"_",'historical'),x=ls(),value = TRUE)
      future00 <- grep(pattern = paste0(model[i],"_",var00[j],"_",time00[k]),x=ls(),value = TRUE)
      if(var00[j] == 'pr') {
        assign(paste0(paste0(model[i],"_",var00[j],"_",time00[k],"_","bias")),lapply(X= get(future00),
                                                                                     FUN = function(x) {
                                                                                       (x - get(hist00))/(get(hist00) + 0.1) })) # changed based on Carlos Navarro-Racines #>sapply(lapply(X= get(future00),FUN = function(x) {x/get(hist00)}),tail,10)) #sapply 选取最后10个 即2100-2300
      } else {
        assign(paste0(paste0(model[i],"_",var00[j],"_",time00[k],"_","bias")),lapply(X= get(future00),FUN = function(x) { x - get(hist00)}))
      }
    }
  }
}

#
all_bias <- grep(pattern = 'bias',x=ls(),value=T)
#
# 
for (i in all_bias) {
  assign(i,do.call(c,tail(get(i),10))) 
}
# write them out
#dir.create('/mnt/DataSpace/Data_pool/02climate_bias_10/')
for (i in all_bias) {
  terra::writeRaster(get(i),paste0('/mnt/DataSpace/Data_pool/02climate_bias_10/',i,'.tif'),overwrite=TRUE)
}

# resample 
#
a0 <- list.files(pattern = 'bias.tif$','/mnt/DataSpace/Data_pool/02climate_bias_10',full.names = TRUE)
a_0 <- list.files(pattern = 'tas_.*.tif$','/mnt/DataSpace/Data_pool/02climate_bias_10',full.names = TRUE)
a_file <- setdiff(a0,a_0)
a_file0 <- grep('ssp585',a_file,value=T)

temp1 <- '/mnt/File0/DAAATAAA/world_climate/5mins/wc2.1_5m_prec/wc2.1_5m_prec_01.tif'
tmplt <- rast(temp1)

dir.create('/mnt/DataSpace/Data_pool/02_1climate_bias5min_10')
for( i in 1:length(a_file0)) {
  rf1 <- terra::rast(a_file0[i])
  wd5bias <- resample(x=rf1,y=tmplt,method='bilinear',threads=T)
  wd5bias1 <- round(wd5bias,2)
  writeRaster(wd5bias1,filename=paste0('/mnt/DataSpace/Data_pool/02_1climate_bias5min_10/0',basename(a_file0[i])),overwrite=T)
  rm(wd5bias,wd5bias1);gc();gc();gc();gc();gc();gc()
}


# 
library(terra)
library(sf)
unixtools::set.tempdir('/mnt/DataSpace/temp4r')
t3 <- st_read('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/2021digitallpr_000b21a_e/lpr_000b21a_e.shp')
t5 <- st_transform(t3,4326)
t6 <- vect(t5)

#20230719, 

biovariables00 <- function (x) {
  if (is.na(x[1])) 
    return(rep(NA, 19))
  x <- matrix(c(x, rep(NA, length(x))), ncol = 6)
  if (nrow(x) == 12) 
    x <- rbind(x, x[1:2, ])
  x[, 4] <- (x[, 2] + x[, 3])/2
  x[1:12, 5] <- caTools::runmean(x[, 1], 3, alg = "fast", endrule = "trim", 
                                 align = "left") * 3
  x[1:12, 6] <- caTools::runmean(x[, 4], 3, alg = "fast", endrule = "trim", 
                                 align = "left")
  x <- x[1:12, ]
  means <- colMeans(x)
  b1 <- means[4]
  b2 <- means[2] - means[3]
  b4 <- sd(x[, 4]) * 100  # based on dismo::biovars , we follow their instruction that multiple 100 
  b5 <- max(x[, 2])
  b6 <- min(x[, 3])
  b7 <- b5 - b6
  b3 <- b2/b7 * 100
  b10 <- max(x[, 6])
  b11 <- min(x[, 6])
  b12 <- means[1] * 12
  b13 <- max(x[, 1])
  b14 <- min(x[, 1])
  b16 <- max(x[, 5])
  b17 <- min(x[, 5])
  b15 <- sd(x[, 1])/(1 + means[1]) * 100
  b8 <- x[, 6][x[, 5] == b16][1]
  b9 <- x[, 6][x[, 5] == b17][1]
  b18 <- x[, 5][x[, 6] == b10][1]
  b19 <- x[, 5][x[, 6] == b11][1]
  return(c(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, 
           b13, b14, b15, b16, b17, b18, b19))
}
##################################################################

var2 <- c('prec','tmax','tmin')
var_cmip <- c('pr','tasmax','tasmin')

pr_hist0 <- list.files(pattern = 'tif', path = paste0('/mnt/File0/DAAATAAA/world_climate/5mins/wc2.1_5m_',var2[1]),full.names = TRUE)
tmax_hist0 <- list.files(pattern = 'tif', path = paste0('/mnt/File0/DAAATAAA/world_climate/5mins/wc2.1_5m_',var2[2]),full.names = TRUE)
tmin_hist0 <- list.files(pattern = 'tif', path = paste0('/mnt/File0/DAAATAAA/world_climate/5mins/wc2.1_5m_',var2[3]),full.names = TRUE)
#}
library(terra)
pr_hist1 <- rast(pr_hist0)
tmax_hist1 <- rast(tmax_hist0)
tmin_hist1 <- rast(tmin_hist0)
#
library(climatica)
library(caTools)
model_gpp <- c("EC-Earth3-Veg-LR", "IPSL-CM6A-LR", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CanESM5-1", "CAS-ESM2-0", "CESM2-WACCM", "CMCC-CM2-SR5", "EC-Earth3-CC", "INM-CM4-8", "INM-CM5-0", "TaiESM1", "EC-Earth3-Veg", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR")
model_2100 <- c("ACCESS-CM2", "CanESM5", "CMCC-ESM2", "EC-Earth3-Veg", "FIO-ESM-2-0", "GISS-E2-1-G", "GISS-E2-1-H", "HadGEM3-GC31-LL", "INM-CM4-8", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "UKESM1-0-LL", "ACCESS-ESM1-5", "CanESM5-CanOE", "CNRM-CM6-1", "CNRM-CM6-1-HR", "CNRM-ESM2-1", "EC-Earth3-Veg-LR", "INM-CM5-0")
model_2100_gpp <- setdiff(model_gpp,model_2100)[-c(4,5,7)]
model_2100_new <- c(model_2100[-c(6,7,8,15,17,18,19,20)],model_2100_gpp)

model <- model_2100_new[-c(2,4,7,11,12)] 

scaro <- c('ssp245','ssp585') 
model <- model_2100_new[c(2,4,7,11,12)] 
scaro <- c('ssp245')

for (j in 2:length(model)) {
  for (k in 1) { 
    f_future_bias <- list.files(pattern = paste0(model[j],'_.*',scaro[k],'.*tif$'),path = '/mnt/DataSpace/Data_pool/02_1climate_bias5min_10',full.names = TRUE)
    pr_future <- rast(f_future_bias[1])
    tmax_future <- rast(f_future_bias[2])
    tmin_future <- rast(f_future_bias[3])
    #
    ind03 <- data.frame(b0 = seq(1,12*10,12),e0 = seq(12,12*10,12)) 
    for (i in 1:10){ 
      pr_future1 <- (abs(pr_future[[ind03[i,1]:ind03[i,2]]]) + 1 ) * pr_hist1 +0.1 
      tmax_future1 <- tmax_future[[ind03[i,1]:ind03[i,2]]] + tmax_hist1
      tmin_future1 <- tmin_future[[ind03[i,1]:ind03[i,2]]] + tmin_hist1
      futur_bc_list <- terra::app(x=c(pr_future1,tmax_future1,tmin_future1),fun=biovariables00,cores=32)#
      nam1 <- rowMeans(data.frame(s1=seq(2000,2090,10),e1=seq(2010,2100,10)))[i]
      names(futur_bc_list) <- paste0('BIO',1:19)
      
      t9 <- crop(futur_bc_list,ext(t6) + 1) 
      t1 <- disagg(t9,fact=10,method='bilinear',
                   filename = paste0('/mnt/DataSpace/Data_pool/03cmip6_processed_bio19_5min_10/',model[j],"_",scaro[k],"_",nam1,"_CA.tif"),
                   overwrite=TRUE,gdal=c('COMPRESS=LZW','ZLEVEL=9','NUM_THREADS=64'),tempdir ='/mnt/temp4r/temp')
      rm(t1,t9);gc();gc();gc();gc();gc();Sys.sleep(20)
    }
  }
}
############# ########### ########### ########### ###########  2300 
########### ########### ########### ########### ########### ########### 
########### ########### ########### ########### ########### ########### 
t1 <- list.files(pattern = 'nc','/mnt/DataSpace/Data_pool/CMIP6_rebuild',full.names = T)

future_period <- data.frame(stardyear=c(1970,seq(2001,2300,10)),endyear=c(2000,seq(2010,2300,10)))
#NOTE: start from 2100 
model_2300 <- c("CanESM5", "EC-Earth3-Veg", "IPSL-CM6A-LR", "MRI-ESM2-0", "ACCESS-ESM1-5") 
time <- c('historical','ssp585')
var00=c('pr','tas','tasmin','tasmax')
all_combin <- expand.grid(model_2300,var00,time)
#
for (ii in 1:nrow(all_combin)) {
  ind0000 <- grep(pattern = paste0('.*',all_combin[ii,2],'_.*',all_combin[ii,1],'_.*',all_combin[ii,3]), x = t1)
  f20 <- rast(t1[ind0000])
  times_span00 <- year(time(f20))
  assign(do.call(paste,c(all_combin[ii,],min(times_span00),max(times_span00),sep='_')) , 
         models_output(var0 = all_combin[ii,2], model0 = all_combin[ii,1],time0 = all_combin[ii,3]))
}


#
model <- model_2300
time00 <- c('historical','ssp585')
var00=c('pr','tas','tasmin','tasmax')

for (i in  1:length(model)){ 
  for(j in 1:4) {
    for (k in 2){
      hist00 <- grep(pattern = paste0(model[i],"_",var00[j],"_",'historical'),x=ls(),value = TRUE)
      future00 <- grep(pattern = paste0(model[i],"_",var00[j],"_",time00[k]),x=ls(),value = TRUE)
      if(var00[j] == 'pr') {
        assign(paste0(paste0(model[i],"_",var00[j],"_",time00[k],"_","bias")),lapply(X= get(future00),
                                                                                     FUN = function(x) {
                                                                                       (x - get(hist00))/(get(hist00) + 0.1) })) # changed based on Carlos Navarro-Racines #>sapply(lapply(X= get(future00),FUN = function(x) {x/get(hist00)}),tail,10)) #sapply 选取最后10个 即2100-2300
      } else {
        assign(paste0(paste0(model[i],"_",var00[j],"_",time00[k],"_","bias")),lapply(X= get(future00),FUN = function(x) { x - get(hist00)}))
      }
    }
  }
}

#
all_bias <- grep(pattern = 'bias',x=ls(),value=T)
#
# 
for (i in all_bias) {
  assign(i,do.call(c,tail(get(i),30))) 
}
# 
for (i in all_bias) {
  terra::writeRaster(get(i),paste0('/mnt/DataSpace/Data_pool/02climate_bias_10/',i,'.tif'),overwrite=TRUE) 
}
###

a0 <- list.files(pattern = 'bias.tif$','/mnt/DataSpace/Data_pool/02climate_bias_10',full.names = TRUE)

a_file0 <- grep('ssp585',a0,value=T)
a_file1 <- do.call(c,lapply(model_2300,function(x) grep(paste0(x,'_'),a_file0,value = T)))
temp1 <- '/mnt/File0/DAAATAAA/world_climate/5mins/wc2.1_5m_prec/wc2.1_5m_prec_01.tif'
tmplt <- rast(temp1)

for( i in 1:length(a_file1)) {
  rf1 <- terra::rast(a_file1[i])
  wd5bias <- resample(x=rf1,y=tmplt,method='bilinear',threads=T)
  wd5bias1 <- round(wd5bias,2)
  writeRaster(wd5bias1,filename=paste0('/mnt/DataSpace/Data_pool/02_1climate_bias5min_10/0',basename(a_file1[i])),overwrite=T)
  rm(wd5bias,wd5bias1);gc();gc();gc();gc();gc();gc()
}
#
#
unixtools::set.tempdir('/mnt/Fastrun/temp4r')
library(terra)
library(sf)
t3 <- st_read('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/2021digitallpr_000b21a_e/lpr_000b21a_e.shp')
t5 <- st_transform(t3,4326)
t6 <- vect(t5)

#climatica::biovariables
biovariables00 <- function (x) {
  if (is.na(x[1])) 
    return(rep(NA, 19))
  x <- matrix(c(x, rep(NA, length(x))), ncol = 6)
  if (nrow(x) == 12) 
    x <- rbind(x, x[1:2, ])
  x[, 4] <- (x[, 2] + x[, 3])/2
  x[1:12, 5] <- caTools::runmean(x[, 1], 3, alg = "fast", endrule = "trim", 
                                 align = "left") * 3
  x[1:12, 6] <- caTools::runmean(x[, 4], 3, alg = "fast", endrule = "trim", 
                                 align = "left")
  x <- x[1:12, ]
  means <- colMeans(x)
  b1 <- means[4]
  b2 <- means[2] - means[3]
  b4 <- sd(x[, 4]) * 100
  b5 <- max(x[, 2])
  b6 <- min(x[, 3])
  b7 <- b5 - b6
  b3 <- b2/b7 * 100
  b10 <- max(x[, 6])
  b11 <- min(x[, 6])
  b12 <- means[1] * 12
  b13 <- max(x[, 1])
  b14 <- min(x[, 1])
  b16 <- max(x[, 5])
  b17 <- min(x[, 5])
  b15 <- sd(x[, 1])/(1 + means[1]) * 100
  b8 <- x[, 6][x[, 5] == b16][1]
  b9 <- x[, 6][x[, 5] == b17][1]
  b18 <- x[, 5][x[, 6] == b10][1]
  b19 <- x[, 5][x[, 6] == b11][1]
  return(c(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, 
           b13, b14, b15, b16, b17, b18, b19))
}
##################################################################

var2 <- c('prec','tmax','tmin')
var_cmip <- c('pr','tasmax','tasmin')
model <- c("CanESM5", "EC-Earth3-Veg", "IPSL-CM6A-LR", "MRI-ESM2-0", "ACCESS-ESM1-5")

pr_hist0 <- list.files(pattern = 'tif', path = paste0('/mnt/File0/DAAATAAA/world_climate/5mins/wc2.1_5m_',var2[1]),full.names = TRUE)
tmax_hist0 <- list.files(pattern = 'tif', path = paste0('/mnt/File0/DAAATAAA/world_climate/5mins/wc2.1_5m_',var2[2]),full.names = TRUE)
tmin_hist0 <- list.files(pattern = 'tif', path = paste0('/mnt/File0/DAAATAAA/world_climate/5mins/wc2.1_5m_',var2[3]),full.names = TRUE)

library(terra)
pr_hist1 <- rast(pr_hist0)
tmax_hist1 <- rast(tmax_hist0)
tmin_hist1 <- rast(tmin_hist0)
#
library(climatica)
library(caTools)

model <- c("CanESM5", "EC-Earth3-Veg", "IPSL-CM6A-LR", "MRI-ESM2-0", "ACCESS-ESM1-5")

scaro <- c('ssp585')
#scaro <- c('ssp585')
dir.create('/mnt/DataSpace/Data_pool/04_bios_variable_10')
for (j in 1:length(model)) {
  for (k in 1) {
    f_future_bias <- list.files(pattern = paste0(model[j],'_.*',scaro[k],'.*tif$'),path = '/mnt/DataSpace/Data_pool/02_1climate_bias5min_10',full.names = TRUE)
    pr_future <- rast(f_future_bias[1])
    tmax_future <- rast(f_future_bias[3])
    tmin_future <- rast(f_future_bias[4])
    #
    ind03 <- data.frame(b0 = seq(1,12*30,12),e0 = seq(12,12*30,12)) 
    for(i in 1:30){ 
      pr_future1 <- (abs(pr_future[[ind03[i,1]:ind03[i,2]]]) + 1 ) * pr_hist1 +0.1 
      tmax_future1 <- tmax_future[[ind03[i,1]:ind03[i,2]]] + tmax_hist1
      tmin_future1 <- tmin_future[[ind03[i,1]:ind03[i,2]]] + tmin_hist1
      futur_bc_list <- terra::app(x=c(pr_future1,tmax_future1,tmin_future1),fun=biovariables00,cores=96)
      
      nam1 <- rowMeans(data.frame(s1=seq(2000,2290,10),e1=seq(2010,2300,10)))[i]
      names(futur_bc_list) <- paste0('BIO',1:19)
      
      t9 <- crop(futur_bc_list,ext(t6) + 1) 
      t1 <- disagg(t9,fact=10,method='bilinear',
                   filename = paste0('/mnt/DataSpace/Data_pool/03cmip6_processed_bio19_5min_10/',model[j],"_",scaro[k],"_",nam1,"_CA.tif"),
                   overwrite=TRUE,gdal=c('COMPRESS=LZW','ZLEVEL=9','NUM_THREADS=16'),tempdir ='/mnt/temp4r/temp')
      rm(t1,t9);gc();gc();gc();gc();gc();Sys.sleep(20)
    }
    
  
  }
}
##############################################
##############################################
############################################
#2300 new 4models like access and giss
library(terra)
cmip6_ca <- list.files(pattern = 'ssp585.*2100_CA.tif$',path = '/mnt/DataSpace/Data_pool/03cmip6_processed_bio19_5min/',recursive = F,full.names = TRUE)

#
unixtools::set.tempdir('/mnt/DataSpace/temp4r')


periods0 <- data.frame(a=seq(2001,2300,10),b=seq(2010,2300,10))
ni <- function(x) {
  begin <- periods0[floor(x/19.1)+1,1] 
  end <- periods0[floor(x/19.1)+1,2]
  writeRaster(round(t1[[x]],1),filename = paste0('/mnt/DataSpace/Data_pool/04_bios_variable_10/BIO',
                                                 ceiling(x %% 19.1),'_',gsub('2020_2100_CA.tif','',basename(sources(t1))),begin,'-',end,'.tif'),overwrite=TRUE)
}

library(parallel)
for(i in 1:length(cmip6_ca)){
  t1 <- rast(cmip6_ca[i])
  mclapply(as.list(1:nlyr(t1)),FUN=ni,mc.cores = getOption("mc.cores", 4L))
}


############################.
#
unixtools::set.tempdir('/mnt/Fastrun/temp4r')
library(terra)
library(FAO56)
f1 <- list.files(pattern = 'nc','/mnt/DataSpace/Data_pool/CMIP6_rebuild',full.names = T)

model_pet <- c("ACCESS-CM2", "ACCESS-ESM1-5", "CanESM5", "CanESM5-1", "CAS-ESM2-0", "FGOALS-g3", "GFDL-ESM4", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "EC-Earth3-CC", "MPI-ESM1-2-LR")
extrass<- c("EC-Earth3-Veg-LR","INM-CM4-8", "INM-CM5-0"  )


#    function 
et_pm <- function(Delta = SlpSVPC(T_mean), T_mean = (T_min + T_max)/2, 
                  R_n = NULL, G = 0, gamma = PsyCon(AtmPres(elev)), u_2 = NULL, 
                  u_z = NULL, z = NULL, e_s = MSVP(T_max, T_min), T_dew = NULL, 
                  e_a = NULL, T_min = NULL, T_max = NULL, phi_deg = NULL, elev = NULL, 
                  date = NULL, n = NULL, N = NULL, a_s = 0.25, b_s = 0.5){
  (0.408 * Delta * (R_n - G) + gamma * (900/(T_mean + 273)) * 
     u_2 * (e_s - e_a))/(Delta + gamma * (1 + 0.34 * u_2))
}
#
PET0 <-
  et_pm(
    Delta = Delta0,
    T_mean = Tmean,
    R_n = Rn,
    G = 0,gamma = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = Delta = 
      gamma = gama,
    u_2 = u2,
    e_s = e_s,
    e_a = e_a
  )# two minutes

get(paste0('tas',0))
world_eve1 <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/Global_eve10km.tif')

####  
mod_use <- extrass
scenarios <- c('historical', 'ssp245', 'ssp585')
for (i in 1:length(mod_use)) {
  for (j in scenarios) {
    #
    hurs0 <- grep(paste0('hurs', '.*', mod_use[i], '_', j), f1, value = T)
    
    rsds0 <- grep(paste0('rsds', '.*', mod_use[i], '_', j), f1, value = T)
    
    sfcWind0 <- grep(paste0('sfcWind', '.*', mod_use[i], '_', j), f1, value = T)
    
    tas0 <- grep(paste0('tas_', '.*', mod_use[i], '_', j), f1, value = T)
    
    tasmax0 <- grep(paste0('tasmax', '.*', mod_use[i], '_', j), f1, value = T)
    
    tasmin0 <- grep(paste0('tasmin', '.*', mod_use[i], '_', j), f1, value = T)
    #
    Tmean <- rast(tas0)
    RelH <- rast(hurs0)
    Tmax <- rast(tasmax0)
    Tmin <- rast(tasmin0)
    u2 <- rast(sfcWind0)
    Rn <- rast(rsds0) * 0.0864 * 30.44 
    if (identical(
      nlyr(Tmean),
      nlyr(RelH),
      nlyr(Tmax),
      nlyr(Tmin),
      nlyr(u2),
      nlyr(Rn)
    )) {
      
      #
      Delta0 <- SlpSVPC(Tmean)
      eoT <- SatVP(Tmean)
      e_a <- eoT * RelH / 100
      e_s <- MSVP(T_max=Tmax, T_min=Tmin)
      #
      elev0 <- terra::extract(x=world_eve1, y= terra::as.points(rotate(Tmean[[1]])))[,2]
      elev <- Tmean[[1]] # as template
      values(elev) <- elev0
      gama <- PsyCon(AtmPres(elev))
      #
      PET0 <-
        et_pm(
          Delta = Delta0,
          T_mean = Tmean,
          R_n = Rn,
          G = 0,
          gamma = gama,
          u_2 = u2,
          e_s = e_s,
          e_a = e_a
        )# two minutes
      #
      assign(paste0(mod_use[i],'_',j),PET0)
    } else {
      message(mod_use[i],'Variables are not equal in layers')
    }
  }
}

ls()
############>>>>>>>>>>
fls <- list.files('tif',path = '/mnt/DataSpace/Data_pool/04_bios_variable/',full.names = T)
dir.create('/mnt/DataSpace/Data_pool/01PET_cmip/Processed_cmip6_pet_10')
dir.create('/mnt/DataSpace/Data_pool/02pet_use_10/')
tdaye <- rast(fls[1])
ca_covas <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global_re_align_2300_d1.tif') 
library(lubridate)
for (i in 1:length(mod_use)){
  time00 <- paste0(seq(2001,2291,10),'_',seq(2010,2300,10))
  historical_pet <- get(paste0(mod_use[i],'_historical'))
  ssp245_pet <- get(paste0(mod_use[i],'_ssp245'))
  ssp585_pet <- get(paste0(mod_use[i],'_ssp585'))
  if(i==6) {
    set.ext(historical_pet,ext(0,360,-90,90))
    set.ext(ssp245_pet,ext(0,360,-90,90))
    set.ext(ssp585_pet,ext(0,360,-90,90))
  }
  hist_ind <- which(time(historical_pet) %in% seq(as.Date('1970-1-1'),as.Date('2000-12-31'),by = 'day')) 
  historical_pet0 <- app(tapp(historical_pet[[hist_ind]],index = rep(1:(length(hist_ind)/12),rep(12,length(hist_ind)/12)),fun= sum),fun=mean)  
  ssp245_pet <- c(historical_pet,ssp245_pet)
  ssp585_pet <- c(historical_pet,ssp585_pet)
  ssp_ind <- which(time(ssp245_pet) %in% seq(as.Date('2001-1-1'),as.Date('2300-12-31'),by = 'day'))
  ssp245_pet0 <- tapp(tapp(ssp245_pet[[ssp_ind]],index = rep(1:(nlyr(ssp245_pet[[ssp_ind]])/12),rep(12,nlyr(ssp245_pet[[ssp_ind]])/12)),fun= sum),index = rep(1:(length(ssp_ind)/(10*12)),rep(10 ,length(ssp_ind)/(10*12))),fun=mean)  
  ssp585_pet0 <- tapp(tapp(ssp585_pet[[ssp_ind]],index = rep(1:(nlyr(ssp585_pet[[ssp_ind]])/12),rep(12,nlyr(ssp585_pet[[ssp_ind]])/12)),fun= sum),index = rep(1:(length(ssp_ind)/(10*12)),rep(10 ,length(ssp_ind)/(10*12))),fun=mean)  
  #
  if(ext(historical_pet0)[2] > 200) {historical_pet0 <- rotate(historical_pet0)}
  if(ext(ssp245_pet0)[2] > 200) {ssp245_pet0 <- rotate(ssp245_pet0)}
  if(ext(ssp585_pet0)[2] > 200) {ssp585_pet0 <- rotate(ssp585_pet0)}
  
  ssp245_bc <- crop(ssp245_pet0/historical_pet0,tdaye)
  ssp585_bc <- crop(ssp585_pet0/historical_pet0,tdaye)
  #
  diff_ssp245_1km <- resample(ssp245_bc,tdaye,method='bilinear',threads=TRUE)
  diff_ssp585_1km <- resample(ssp585_bc,tdaye,method='bilinear',threads=TRUE)
  #
  diff_ssp245_1km <- clamp(diff_ssp245_1km,upper=10,values=T)
  diff_ssp585_1km <- clamp(diff_ssp585_1km,upper=10,values=T)
  
  #
  future_ssp245 <- diff_ssp245_1km * ca_covas[[5]]
  future_ssp585 <- diff_ssp585_1km * ca_covas[[5]]
  writeRaster(historical_pet,filename = paste0('/mnt/DataSpace/Data_pool/01PET_cmip/Processed_cmip6_pet_10/',mod_use[i],'_','historical','_',year(min(time(ssp245_pet))),'-',year(max(time(ssp245_pet))),'.tif'))
  writeRaster(ssp245_pet,filename = paste0('/mnt/DataSpace/Data_pool/01PET_cmip/Processed_cmip6_pet_10/',mod_use[i],'_','ssp245','_',year(min(time(ssp585_pet))),'-',year(max(time(ssp585_pet))),'.tif')) 
  writeRaster(ssp585_pet,filename = paste0('/mnt/DataSpace/Data_pool/01PET_cmip/Processed_cmip6_pet_10/',mod_use[i],'_','ssp585','_',year(min(time(ssp245_pet))),'-',year(max(time(ssp245_pet))),'.tif'))
  #
  writeRaster(future_ssp245,filename = paste0('/mnt/DataSpace/Data_pool/02pet_use_10/',mod_use[i],'_','ssp245','_',year(min(time(ssp245_pet))),'-',year(max(time(ssp245_pet))),'.tif'))
  writeRaster(future_ssp585,filename = paste0('/mnt/DataSpace/Data_pool/02pet_use_10/',mod_use[i],'_','ssp585','_',year(min(time(ssp585_pet))),'-',year(max(time(ssp585_pet))),'.tif')) 
}

##. to each 
##. to each 
library(terra)
ff0 <- list.files('tif',path = '/mnt/DataSpace/Data_pool/02pet_use_10/',full.names = T)
# 
dir.create('/mnt/DataSpace/Data_pool/03pet_use_individual_10')
ff0_b =basename(ff0)

ff0_b1 = gsub('_1850-2100.tif|_1850-2300.tif','',ff0_b)
for(i in 1:length(ff0)){
  (time00 <- paste0(seq(2001,2291,10),'-',seq(2010,2300,10)))
  t1 <- rast(ff0[i])
  for(j in 1:nlyr(t1)){
    writeRaster(round(t1[[j]],0),filename = paste0('/mnt/DataSpace/Data_pool/03pet_use_individual_10/','pet_',ff0_b1[i],'_',time00[j],'.tif'),overwrite=T)
  }
}
###for extss
extbb <- list()
for(i in extrass){
  extbb[[i]] <- grep(i,ff0_b,value=T)
}
extff <- list()
for(i in extrass){
  extff[[i]] <- grep(i,ff0,value=T)
}
extff1 <- unlist(extff)
extbb1 <- unlist(extbb)
ff0_b2 = gsub('_1850-2100.tif|_1850-2300.tif','',extbb1)
for(i in 1:length(extff1)){
  (time00 <- paste0(seq(2001,2291,10),'-',seq(2010,2300,10)))
  t1 <- rast(extff1[i])
  for(j in 1:nlyr(t1)){
    writeRaster(round(t1[[j]],0),filename = paste0('/mnt/DataSpace/Data_pool/03pet_use_individual_10/','pet_',ff0_b2[i],'_',time00[j],'.tif'),overwrite=T)
  }
}


bc_ssp <- crop(diff_ssp,tdaye)
diff_ssp_1km <- resample(bc_ssp,tdaye,method='bilinear',threads=TRUE)

##########################################################
########################################################################################################################################################
############# water

f1 <- list.files(pattern = 'nc','/mnt/DataSpace/Data_pool/CMIP6_rebuild',full.names = T)
mrsos0 <- grep('mrsos_',f1,value=T)
mrso0 <- grep('mrso_',f1,value=T)
mod1 <- unique(unlist(lapply(strsplit(mrsos0,'/|_'),function(x) x[10])))
mod2 <- unique(unlist(lapply(strsplit(mrso0,'/|_'),function(x) x[10])))

mod3 <- intersect(mod1,mod2)

tt1 <- data.frame(model = unlist(lapply(strsplit(mrsos0,'/|_'),function(x) x[10])),
                  period = unlist(lapply(strsplit(mrsos0,'/|_'),function(x) x[11]))
) %>% mutate(mod_period = paste0(model,'-',period))

tt4 <- table(tt1$model[unlist(lapply(unique(tt1$mod_period),function(x) which(tt1$mod_period == x)[1]))])
w_top <- tt4[which(tt4==3)]

tt2 <- data.frame(model = unlist(lapply(strsplit(mrso0,'/|_'),function(x) x[10])),
                  period = unlist(lapply(strsplit(mrso0,'/|_'),function(x) x[11]))
) %>% mutate(mod_period = paste0(model,'-',period))

tt3 <- table(tt2$model[unlist(lapply(unique(tt2$mod_period),function(x) which(tt2$mod_period == x)[1]))])
w_all<- tt3[which(tt3==3)]
mod3 <- intersect(names(w_top),names(w_all))
# refer to /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/001code/00_1complementary_experiments.R
library(terra)

library(lubridate)

library(sf)

library(dplyr)


f1 <- list.files(pattern = 'nc','/mnt/DataSpace/Data_pool/CMIP6_rebuild',full.names = T)
ca_shap1 <- vect('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/2021digitallpr_000b21a_e/lpr_000b21a_e.shp')
ca_shp2 <- terra::project(ca_shap1,'epsg:4326')

ca_covas <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global_re_align_2300_d1.tif') 
tdaye <- ca_covas[[1]]
model = mod3[7]
var1 = 'mrsos'
scenario = 'ssp245'

dir.create('/mnt/DataSpace/Data_pool/05_smap_10')
get_bc <- function(model,scenario) {
  time00 <- paste0(seq(2001,2291,10),'-',seq(2010,2300,10))
 
  pure_hist0 <- rast(grep(paste0(model,'_historical'),grep('mrsos_',f1,value = T),value = T)) #
  pure_hist_all <- rast(grep(paste0(model,'_historical'),grep('mrso_',f1,value = T),value = T)) # 
  if(model == 'EC-Earth3-Veg') {
    pure_hist1 <- pure_hist_all[[1:1980]] - pure_hist0 
  } else {
    pure_hist1 <- pure_hist_all - pure_hist0 
  }
  
  crs(pure_hist1) <- crs(tdaye)
  crs(pure_hist0) <- crs(tdaye)
  #
  if(xmax(pure_hist0) > 200) {pure_hist_t2 <- rotate(pure_hist0)} else {pure_hist_t2 <- pure_hist0}
  if(xmax(pure_hist1) > 200) {pure_hist_t2_2 <- rotate(pure_hist1)} else {pure_hist_t2_2 <- pure_hist1}
  
  # 
  hist0 <- rast(grep(paste0(model,'_ssp245'),grep('mrsos_',f1,value = T),value = T)) 
  hist_all <- rast(grep(paste0(model,'_ssp245'),grep('mrso_',f1,value = T),value = T)) # 
  if(model == 'EC-Earth3-Veg') {
    hist1 <- hist_all[[1:1032]] - hist0 
  } else {
    hist1 <- hist_all - hist0 
  }
  crs(hist1) <- crs(tdaye)
  crs(hist0) <- crs(tdaye)
  if(model=='BCC-CSM2-MR') {time(hist0) <- seq(as.Date('2015-01-16'),as.Date('2100-12-16'),by='month')}
  if(model=='BCC-CSM2-MR') {time(hist1) <- seq(as.Date('1850-01-16'),as.Date('2014-12-16'),by='month')}
  #
  if(xmax(hist0) > 200) {hist_t2 <- rotate(hist0)} else {hist_t2 <- hist0}
  if(xmax(hist1) > 200) {hist_t2_2 <- rotate(hist1)} else {hist_t2_2 <- hist1}
  
  hist_ind <- which(time(hist_t2) %in% seq(as.Date('2016-01-01'),as.Date('2021-12-31'),by='day'))
  hist_t3 <- hist_t2[[hist_ind]] # for mrsos
  hist_t4 <- app(tapp(hist_t3,index= rep(1:(length(hist_ind)/12),rep(12,length(hist_ind)/12)), fun= sum),fun=mean)
  
  hist_ind_2 <- which(time(hist_t2_2) %in% seq(as.Date('2016-01-01'),as.Date('2021-12-31'),by='day'))
  hist_t3_2 <- hist_t2_2[[hist_ind_2]] # for mrsumap
  hist_t4_2 <- app(tapp(hist_t3_2,index= rep(1:(length(hist_ind_2)/12),rep(12,length(hist_ind_2)/12)), fun= sum),fun=mean)
  
  # 
  #
  f_ssp <- grep(paste0(model,'_',scenario),grep('mrsos_',f1,value = T),value = T)
  
  if(length(f_ssp) < 1 ) {diff_ssp <- NULL} else {
    ssp0 <- rast(f_ssp)
    crs(ssp0) <- crs(tdaye)
    if(model=='BCC-CSM2-MR') {time(ssp0) <- seq(as.Date('2015-01-16'),as.Date('2100-12-16'),by='month')}
    if(xmax(ssp0) > 200) {ssp_t2 <- rotate(ssp0)} else {ssp_t2 <- ssp0}
    
    ssp_t2 = c(pure_hist_t2,ssp_t2) # 
    ind <- which(time(ssp_t2) %in% seq(as.Date('2001-01-01'),as.Date('2300-12-31'),by='day')) # 
    ssp_t3 <- ssp_t2[[ind]]
    ssp_t4 <- tapp(tapp(ssp_t3,index = rep(1:(length(ind)/(12)),rep(12,length(ind)/(12))),fun= sum),index= rep(1:(length(ind)/(12*10)),rep(10,length(ind)/(12*10))),fun=mean) # 
    diff_ssp <- ssp_t4 / (hist_t4 + 1*10^-20)  #  
    bc_ssp <- crop(x= diff_ssp,y = ext(ca_shp2) + 1)
    bc_ssp1 <- clamp(bc_ssp,upper=10,values=T) # 
    diff_ssp_1km <- resample(bc_ssp1,tdaye,method='bilinear',threads=TRUE)
    f_ssp_1km  <- diff_ssp_1km * ca_covas[[50]]
  }
  if(!is.null(diff_ssp)) {
    lapply(1:nlyr(f_ssp_1km), function(x) writeRaster(f_ssp_1km[[x]],filename = paste0('/mnt/DataSpace/Data_pool/05_smap_10/','SMAP','_',model,'_',scenario,'_',time00[x],'.tif'),overwrite=T))
  }
  
  f_ssp_1 <- grep(paste0(model,'_',scenario),grep('mrso_',f1,value = T),value = T)
  
  if(length(f_ssp_1) < 1 ) {diff_ssp <- NULL} else {
    ssp1 <- rast(f_ssp_1)
    ssp10 <- ssp1 - ssp0 # 
    crs(ssp10) <- crs(tdaye)
    if(model=='BCC-CSM2-MR') {time(ssp10) <- seq(as.Date('2015-01-16'),as.Date('2100-12-16'),by='month')}
    if(xmax(ssp10) > 200) {ssp_t20 <- rotate(ssp10)} else {ssp_t20 <- ssp10}
    
    ssp_t20 <- c(pure_hist_t2_2,ssp_t20) # 
    ind <- which(time(ssp_t20) %in% seq(as.Date('2001-01-01'),as.Date('2300-12-31'),by='day')) 
    ssp_t30 <- ssp_t20[[ind]]
    ssp_t40 <- tapp(tapp(ssp_t30,index = rep(1:(length(ind)/(12)),rep(12,length(ind)/(12))),fun= sum),index= rep(1:(length(ind)/(12*10)),rep(10,length(ind)/(12*10))),fun=mean) # 
    diff_ssp0 <- ssp_t40 / (hist_t4_2 + 1*10^-20)  # 
    bc_ssp0 <- crop(x= diff_ssp0,y = ext(ca_shp2) + 1)
    bc_ssp1_2 <- clamp(bc_ssp0,upper=10,values=T) # 
    diff_ssp_1km_2 <- resample(bc_ssp1_2,tdaye,method='bilinear',threads=TRUE)
    f_ssp_1km_2  <- diff_ssp_1km_2 * ca_covas[[51]]
  }
  if(!is.null(diff_ssp0)) {
    lapply(1:nlyr(f_ssp_1km_2), function(x) writeRaster(f_ssp_1km_2[[x]],filename = paste0('/mnt/DataSpace/Data_pool/05_smap_10/','SUMAP','_',model,'_',scenario,'_',time00[x],'.tif'),overwrite=T))
  }
}

#16 models and output done done done 

# mod3[5] ssp245 ssp585 's mrso all to 2300 but mrsos to 2100, so SUMAP 2100- 2300 were not comparable to SMAP, then we delete 2100-2300 SUMAP in folder 
for(m in mod3[-c(7:9)]) {
  for (j in c('ssp245', 'ssp585')) {
    get_bc(model=m,scenario = j)
  }
}

#
grep(glob2rx(paste0('*mrso_*',mod3[11],'*')),f1,value=T)


