## generate future GPP file based on EC 



library(terra)
library(sf)
library(spData)
library(dplyr)
library(tidyr)
# COX et al. 2013
linear_reg <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  
  xm <- mean(x)
  ym <- mean(y)
  
  x2 <- x * x
  y2 <- y * y
  xy <- x * y
  
  ssxx <- sum(x2) - nx * xm^2
  ssyy <- sum(y2) - ny * ym^2
  ssxy <- sum(xy) - ny * xm * ym
  
  b <- ssxy / ssxx
  a <- ym - b * xm
  
  yf <- a + b * x
  
  r2 <- ssxy^2 / (ssxx * ssyy)
  
  e2 <- (y - yf)^2
  s2 <- sum(e2) / (nx - 2)
  
  s <- sqrt(s2)
  
  da <- s * sqrt(1.0 / nx + xm^2 / ssxx)
  db <- s / sqrt(ssxx)
  
  #
  minx <- min(x) - 0.1 * (max(x) - min(x))
  maxx <- max(x) + 0.1 * (max(x) - min(x))
  nfit <- 200
  dx <- (maxx - minx) / nfit
  xfit <- minx + seq(0, nfit - 1) * dx
  yfit <- a + b * xfit
  yband <- numeric(nfit)
  
  #
  for (n in seq_along(xfit)) {
    yband[n] <- sqrt(s2 * (1.0 + 1.0 / nx + (xfit[n] - xm)^2 / ssxx))
  }
  
  return(list(
    yf = yf,
    a = a,
    b = b,
    da = da,
    db = db,
    xfit = xfit,
    yfit = yfit,
    yband = yband
  ))
}

EC_method1 <- readRDS('/mnt/DataSpace/Projects/Canada_C/processed_output/ec_method_ssp245_test.rds')
EC_method1 <- readRDS('/mnt/DataSpace/Projects/Canada_C/processed_output/ec_method_ssp585_test.rds')

time0 <- paste0(seq(2000, 2290, 10), '-', seq(2010, 2300, 10))
par(mfrow = c(3, 4))
ratio = list()
for (i in 1:30) {
  #>
  tt1 <- unlist(lapply(EC_method1, function(x)
    x[[1]][i, 4]))
  tt2 <- unlist(lapply(EC_method1, function(x)
    x[[5]]))
  #
  if (i > 10) {
    tt1 <- c(tt1[c(1, 2)],
             mean(tt1[3:4], na.rm = T),
             tt1[5],
             mean(tt1[6:8], na.rm = T),
             mean(tt1[9:10], na.rm = T),
             tt1[c(11, 12)])[c(3, 4, 5, 7)]
    tt2 <- c(tt2[c(1, 2)],
             mean(tt2[3:4], na.rm = T),
             tt2[5],
             mean(tt2[6:8], na.rm = T),
             mean(tt2[9:10], na.rm = T),
             tt2[c(11, 12)])[c(3, 4, 5, 7)]
    tt2 <- log(tt2 * 10000)
  } else{
    tt1 <- c(tt1[c(1, 2)],
             mean(tt1[3:4], na.rm = T),
             tt1[5],
             mean(tt1[6:8], na.rm = T),
             mean(tt1[9:10], na.rm = T),
             tt1[c(11, 12)])
    tt2 <- c(tt2[c(1, 2)],
             mean(tt2[3:4], na.rm = T),
             tt2[5],
             mean(tt2[6:8], na.rm = T),
             mean(tt2[9:10], na.rm = T),
             tt2[c(11, 12)])
    tt2 <- log(tt2 * 10000)
  }
  #>
  (y <- tt1)
  (x <- tt2)
  lm0 <- lm(y ~ x)
  t2 <- predict(lm0)
  (mn_pr <- mean(y))
  (std_pr <- sd(y))
  #
  #(x_obs <- 8.976 )# * 2.9171) # sensitivity * gdd slope
  #(dx_obs <- 1.649 )# * 0.4748) #
  
  #if log tt2 based on X ~ N(u,q^2),then log(X) ~ N(log(u)-1/2q^2, q^2/u^2)
  (x_0 <- log(8.976 / 10000))#
  (dx_obs <- 1.649 / 8.976)#
  (x_obs <- x_0 - 0.5 * dx_obs^2 + log(10000))
  #x_obs <- log(8.976)
  #
  mn <- x_obs
  std <- dx_obs
  xbest <- x_obs
  xlo <- x_obs - dx_obs
  xhi <- x_obs + dx_obs
  
  #
  (yr1 <- min(y) - 0.1 * (max(y) - min(y)))
  (yr2 <- max(y) + 0.1 * (max(y) - min(y)))
  (xr1 <- min(x) - 0.1 * (max(x) - min(x)))
  (xr2 <- max(x) + 0.1 * (max(x) - min(x)))
  
  #
  fit <- linear_reg(x, y)
  (yf <- fit$yf)
  (a <- fit$a)
  (b <- fit$b)
  (da <- fit$da)
  (db <- fit$db)
  (xfit <- fit$xfit)
  (yfit <- fit$yfit)
  (yband <- fit$yband)
  
  #
  x2 <- xfit
  nfitx <- length(xfit)
  dx <- diff(x2)[1] #
  Px <- x2
  Pi <- pi
  Px <- 1 / sqrt(2 * Pi * std^2) * exp(-((x2 - mn) / (sqrt(2) * std))^2)
  #
  Px0 <- dnorm(x2, mean = mn, sd = std)
  
  #
  miny <- mn_pr - 5 * std_pr
  maxy <- mn_pr + 5 * std_pr
  mfity <- 2000
  dy <- (maxy - miny) / mfity #
  y2 <- seq(miny, maxy, length.out = mfity) #
  
  #
  Py_pr <- dnorm(y2, mean = mn_pr, sd = std_pr)
  # plot(y2,Py_pr,col='red')
  # #Py_r <- dnorm(y2, mean = 0.37, sd = 0.14)
  
  #
  Pxy <- matrix(0, nrow = nfitx, ncol = mfity)
  Pyx <- matrix(0, nrow = mfity, ncol = nfitx)
  Py <- rep(0, mfity)
  Py_norm <- 0
  pxy_test <- vector()
  py_test <- vector()
  py_norm_test <- vector()
  for (m in seq_len(mfity)) {
    for (n in seq_len(nfitx)) {
      Py_given_x <- dnorm(y2[m], mean = yfit[n], sd = yband[n]) #
      pxy_test[n] <- Py_given_x #plot(x2,pxy_test,col='blue')#plot(x2,Px*pxy_test,type='l')
      Pxy[n, m] <- Px[n] * Py_given_x #
      Pyx[m, n] <- Pxy[n, m] #
      #
      Py[m] <- Py[m] + Pxy[n, m] * dx
      py_test[n] <- Py[m] #plot(x2,py_test)
    }
    Py_norm <- Py_norm + Py[m] * dy #
    py_norm_test[m] <- Py_norm
  }
  
  #
  # ttpy<- dnorm(y2,mean= a+b*mn,sd = sqrt(b^2*dx_obs + std_pr^2))
  # lines(y2,ttpy,type='l',col='black')
  
  Py <- Py / Py_norm #
  
  plot(
    y2,
    Py,
    type = "l",
    col = "blue",
    lwd = 5,
    xlab = "Log based sensitivity",
    ylab = "Probability Density Per PgC"
  )
  lines(y2, Py_pr, col = "red", lwd = 5)
  
  
  #>
  (mean_post <- sum(y2 * Py) * dy)
  (mean_pr <- sum(y2 * Py_pr) * dy)
  
  (sd_pr <- sqrt(sum((y2 - mean_pr)^2 * Py_pr) * dy))
  (sd_post <- sqrt(sum((y2 - mean_post)^2 * Py) * dy))
  
  ###############################
  ###############################
  #
  quantile <- pnorm(y, mean = mean_pr, sd = sd_pr)
  quantile
  ##############################
  #
  p <- quantile
  #
  value <- qnorm(p, mean = mean_post, sd = sd_post)
  value
  #
  if (i > 10) {
    ratio[[i]] = data.frame(
      model_0 = c(
        "ACCESS-ESM1-5",
        "BCC-CSM2-MR",
        "CanESM5",
        "CESM2-WACCM",
        "EC-Earth3",
        "INM-CM",
        "IPSL-CM6A-LR",
        "TaiESM1"
      )[c(3:5, 7)],
      ratio = value / y
    )
  } else {
    ratio[[i]] = data.frame(
      model_0 = c(
        "ACCESS-ESM1-5",
        "BCC-CSM2-MR",
        "CanESM5",
        "CESM2-WACCM",
        "EC-Earth3",
        "INM-CM",
        "IPSL-CM6A-LR",
        "TaiESM1"
      ),
      ratio = value / y
    )
  }
}
names(ratio) = time0[1:length(ratio)]
ratio_df <- do.call(rbind, ratio) %>% mutate(time = data.frame(matrix(
  unlist(strsplit(rownames(.), split = '\\.')), ncol = 2, byrow = T
))[, 1])
# head(ratio_df)
# ratio_df %>% filter(model_0 %in% model)%>% pull(ratio)
##########>################>################>###ratio##########>################>################>################>################>################>################>######
##########>##########>################>################>################>################>################>################>######
##########>##########>################>################>################>################>################>################>######
##########
##########>################>################>################>################>################>######
##########>##########>################>################>################>################>######
##########>##########>################>################>################>################>######
##########>
# > from /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/003_code_lehigh/09_002_20240414updated_version_of_GPP.R
library(terra)

library(lubridate)

library(sf)

library(dplyr)


f1 <- list.files(
  pattern = 'nc',
  '/mnt/DataSpace/Data_pool/CMIP6_rebuild',
  full.names = T
)

ca_shap1 <- vect(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/2021digitallpr_000b21a_e/lpr_000b21a_e.shp'
)
ca_shp2 <- terra::project(ca_shap1, 'epsg:4326')

#ca_covas <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global_re_align_2300_d1.tif') # digital = 1 ,this version meet the extent to 2300
ca_covas <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global1214.tif')
gpp_fluxcom_hist <- rast(
  '/mnt/DataSpace/Projects/Canada_C/processed_output/GPP_historical_from_fluxcom_1981_2000_gCm2year.tif'
)
tdaye <- ca_covas[[1]]
model = 'EC-Earth3-Veg-LR'
var1 = 'gpp'
scenario = 'ssp245'
##### this function combine function from /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/005_python_EC/0006_GPP_EC_concise_and_used.R and ##### this function combine function from /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/005_python_EC/0006_GPP_EC_concise_and_used.R and /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/003_code_lehigh/09_002_20240414updated_version_of_GPP.R
##### this function combine function from /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/005_python_EC/0006_GPP_EC_concise_and_used.R and ##### this function combine function from /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/005_python_EC/0006_GPP_EC_concise_and_used.R and /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/003_code_lehigh/09_002_20240414updated_version_of_GPP.R
###### and /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/004_using_GPP_fVegLitter/001_BaseGPP_future_file_creation_fromCMIP.R
get_bc <- function(var1, model, scenario) {
  cat('\n', model, '+', scenario, '\n')
  mod_inds <- unlist(lapply(ratio_df$model_0, function(x)
    grepl(x, model)))
  EC_ratio <- ratio_df %>% filter(mod_inds) %>% pull(ratio)
  cat(EC_ratio, '\n')
  time00 <- paste0(seq(2001, 2291, 10), '-', seq(2010, 2300, 10))
  hist0 <- rast(grep(
    paste0(model, '_historical'),
    grep('gpp', f1, value = T),
    value = T
  ))
  crs(hist0) <- "EPSG:4326" #crs(tdaye)
  if (model == 'BCC-CSM2-MR') {
    terra::time(hist0) <- seq(as.Date('1850-01-16'), as.Date('2014-12-16'), by =
                                'month')
  }
  if (xmax(hist0) > 200) {
    hist_t2 <- rotate(hist0)
  } else {
    hist_t2 <- hist0
  }
  ###
  ###
  
  library(spData)
  data("world")
  ca_shp <- world %>% filter(name_long == 'Canada')
  ###
  if (var1 == 'rsds') {
    hist_ind <- which(time(hist_t2) %in% seq(as.Date('1970-01-01'), as.Date('2000-12-31'), by =
                                               'day'))
  } else {
    #hist_ind <- which(time(hist_t2) %in% seq(as.Date('2004-01-01'),as.Date('2014-12-31'),by='day'))
    hist_ind <- which(time(hist_t2) %in% seq(as.Date('1981-01-01'), as.Date('2000-12-31'), by =
                                               'day'))  #20240323 changed from 2004 - 2014 to 1970-2000
  }
  hist_t3 <- hist_t2[[hist_ind]]
  hist_t4 <- app(
    tapp(
      hist_t3,
      index = rep(1:(length(hist_ind) / 12), rep(12, length(hist_ind) / 12)),
      fun = function(x)
        sum(x, na.rm = T)
    ),
    fun = function(x)
      mean(x, na.rm = T)
  ) # mean value of historical
  #
  f_ssp <- grep(paste0(model, '_', scenario), grep('gpp', f1, value = T), value = T)
  if (length(f_ssp) < 1) {
    diff_ssp <- NULL
  } else {
    ssp0 <- rast(f_ssp)
    crs(ssp0) <- "EPSG:4326" #crs(tdaye)
    if (model == 'BCC-CSM2-MR') {
      terra::time(ssp0) <- seq(as.Date('2015-01-16'), as.Date('2100-12-16'), by =
                                 'month')
    }
    if (xmax(ssp0) > 200) {
      ssp_t2_1 <- rotate(ssp0)
    } else {
      ssp_t2_1 <- ssp0
    }
    # ind <- which(time(ssp_t2) %in% seq(as.Date('2021-01-01'),as.Date('2300-12-31'),by='day')) # some febrary is 2-15  and others are 01-16
    if (model == 'CESM2-WACCM' & scenario == 'ssp585') {
      ssp_t2 <- c(hist_t2, ssp_t2_1) # 20240325 newly added
      ind <- which(time(ssp_t2) %in% seq(as.Date('2001-01-01'), as.Date('2300-12-31'), by =
                                           'day'))## 20240325 newly updated 2021 to 2001
      ind <- c(ind, tail(ind, 12)) # 2299 + 2299 12 months
      ssp_t3 <- ssp_t2[[ind]]
      terra::time(ssp_t3) <- seq(as.Date('2001-01-15'), as.Date('2300-12-15'), by =
                                   'month')
      
      # ind_99 <- which(time(ssp_t2) %in% seq(as.Date('1970-01-01'),as.Date('2300-12-31'),by='day'))
      # ssp_t3_99 <- ssp_t2[[c(ind_99,tail(ind_99,12))]]
    } else if (model == 'TaiESM1' &
               scenario == 'ssp585') {
      #becaused from ssp585 taiesm1 as starts 2015-02-15
      ssp_t2 <- c(hist_t2, ssp_t2_1) # 20240325 newly added
      ind <- which(time(ssp_t2) %in% seq(as.Date('2001-01-01'), as.Date('2300-12-31'), by =
                                           'day'))## 20240325 newly updated 2021 to 2001
      insert_point <- which(ind == which(time(ssp_t2) == as.Date('2015-02-15'))) # for feb is 02-15
      insert_point_preyear <- which(ind == which(time(ssp_t2) == as.Date('2014-01-16'))) #the rest is 16th day of the month
      ind <- c(ind[1:c(insert_point - 1)], ind[insert_point_preyear], ind[c(insert_point):length(ind)]) # 2015-2-16 2015-1-16
      ssp_t3 <- ssp_t2[[ind]]
      terra::time(ssp_t3) <- seq(as.Date('2001-01-15'), as.Date('2100-12-15'), by =
                                   'month')
      
      # ind_99 <- which(time(ssp_t2) %in% seq(as.Date('1970-01-01'),as.Date('2300-12-31'),by='day'))
      # insert_point0 <- which(ind_99 == which(time(ssp_t2) == as.Date('2015-02-15')))
      # insert_point0_preyear <- which(ind_99 == which(time(ssp_t2) == as.Date('2014-01-16')))
      # ind_99 <- c(ind_99[1:c(insert_point0-1)],ind_99[insert_point0_preyear],ind_99[insert_point0:length(ind_99)]) # 2015-2-16 2015-1-16
      # ssp_t3_99 <- ssp_t2[[ind_99]]
    } else {
      ssp_t2 <- c(hist_t2, ssp_t2_1) # 20240325 newly added
      ind <- which(time(ssp_t2) %in% seq(as.Date('2001-01-01'), as.Date('2300-12-31'), by =
                                           'day'))
      ssp_t3 <- ssp_t2[[ind]]
      # ind_99 <- which(time(ssp_t2) %in% seq(as.Date('1970-01-01'),as.Date('2300-12-31'),by='day'))
      # ssp_t3_99 <- ssp_t2[[ind_99]]
    }
    ssp_t4 <- tapp(tapp(ssp_t3, index = rep(1:(
      length(ind) / (12)
    ), rep(
      12, length(ind) / (12)
    )), fun = sum),
    index = rep(1:(length(ind) / (12 * 10)), rep(10, length(ind) / (12 * 10))),
    fun = mean) # "ACCESS-ESM1-5" 11 months per year, wierd
    diff_ssp <- ssp_t4 / (hist_t4 + 1 * 10^-20)  # quantile(values(crop(hist_t4,ext(ca_shp2) +1 )),0.5,na.rm=T) as evaporation related to precipitation , so bias correction should use ratio
    bc_ssp <- crop(x = diff_ssp, y = ext(ca_shp2) + 1)
    bc_ssp1 <- clamp(bc_ssp, upper = 10, values = T) # set the largest increase as 10 , as it will divided by 10000 evently
    diff_ssp_1km <- resample(bc_ssp1, tdaye, method = 'bilinear', threads =
                               TRUE)
    if (var1 == 'gpp') {
      f_ssp_1km  <- diff_ssp_1km * gpp_fluxcom_hist * EC_ratio #* 365 * 24 * 3600 * 10000 as we just need multiplicator, so dont need to get the original value.   scale 0.0001 see gee from second to annual #ssp_t4 * 365 * 24 * 3600 # unit is s, our gpp from GEE is annual total amount
      #writeRaster(f_ssp_1km,paste0('/mnt/DataSpace/Data_pool/cmip6_complementary/processed_output/',model,'_',scenario,'_',var1,'_',year(min(time(ssp0))),'-',year(max(time(ssp0))),'.tif'))
    } else if (var1 == 'npp') {
      f_ssp_1km  <- diff_ssp_1km * ca_covas[[3]] #365 * 24 * 3600 * 10000
      #writeRaster(f_ssp_1km,paste0('/mnt/DataSpace/Data_pool/cmip6_complementary/processed_output/',model,'_',scenario,'_',var1,'_',year(min(time(ssp0))),'-',year(max(time(ssp0))),'.tif'))
    } else if (var1 == 'rsds') {
      f_ssp_1km  <- diff_ssp_1km  * ca_covas[[26]] #* 86.4242 ssp w m-2 , original kj m-2 d-1, 1w/m2 = 24/0.2777=86.4242 kJ m-2 day-1https://electronics.stackexchange.com/questions/538921/how-to-convert-calculate-textkj-textm2-to-textw-textm2
      #writeRaster(f_ssp_1km,paste0('/mnt/DataSpace/Data_pool/cmip6_complementary/processed_output/',model,'_',scenario,'_',var1,'_',year(min(time(ssp0))),'-',year(max(time(ssp0))),'.tif'))
    }
  }
  if (!is.null(diff_ssp)) {
    #lapply(1:nlyr(f_ssp_1km), function(x) writeRaster(f_ssp_1km[[x]],filename = paste0('/mnt/DataSpace/Data_pool/01cmip6_processed_output/',var1,'_',model,'_',scenario,'_',time00[x],'.tif'),overwrite=T))
    lapply(1:nlyr(f_ssp_1km), function(x)
      writeRaster(
        f_ssp_1km[[x]],
        filename = paste0(
          '/mnt/DataSpace/Data_pool/01cmip6_processed_output_EC_GPP/',
          var1,
          '_',
          model,
          '_',
          scenario,
          '_',
          time00[x],
          '.tif'
        ),
        overwrite = T
      )) # change the folder to 01cmip6_processed_output_updata_npp
  }
}


# ff <- list.files('/mnt/DataSpace/Data_pool/01cmip6_processed_output_EC_GPP/',full.names = T)
# fft <- rast(ff[2])
# plot(fft)
# #> Area
# ca_pure_area84 <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/ca_cmip6_area_WGS84.tif') #new created
# fft1 <- fft * ca_pure_area84
# global(fft1,function(x) sum(x,na.rm = T))/10^15
#
# ff0 <- list.files('/mnt/DataSpace/Data_pool/01cmip6_processed_output_update_GPP_NPP_from_fluxcom_20240423',full.names = T)
# fft0 <- rast(grep('IPSL-CM6A-LR_ssp245',ff0,value=T)[1])
# plot(fft0)
# fft2 <- fft0 * ca_pure_area84
# global(fft2,function(x) sum(x,na.rm = T))/10^15

#16 models and output done done done
model_gpp <- c(
  "EC-Earth3-Veg-LR",
  "IPSL-CM6A-LR",
  "ACCESS-ESM1-5",
  "BCC-CSM2-MR",
  "CanESM5",
  "CanESM5-1",
  "CAS-ESM2-0",
  "CESM2-WACCM",
  "CMCC-CM2-SR5",
  "EC-Earth3-CC",
  "INM-CM4-8",
  "INM-CM5-0",
  "TaiESM1",
  "EC-Earth3-Veg",
  "MPI-ESM1-2-HR",
  "MPI-ESM1-2-LR"
)
mod_u <- c(
  'ACCESS-ESM1-5',
  'BCC-CSM2-MR',
  'CanESM5',
  'CanESM5-1',
  'CAS-ESM2-0',
  'CESM2-WACCM',
  'EC-Earth3-CC',
  'EC-Earth3-Veg',
  'EC-Earth3-Veg-LR',
  'INM-CM4-8',
  'INM-CM5-0',
  'IPSL-CM6A-LR',
  'TaiESM1'
)
mod_u <- mod_u[-5]
mod_u <- setdiff(mod_u, model_gpp[1:6])
for (m in mod_u) {
  for (i in c('gpp')) {
    for (j in c('ssp585')) {
      get_bc(var1 = i,
             model = m,
             scenario = j)
    }
  }
}



############fVegLitter



library(terra)
library(spData)
library(sf)
library(dplyr)
data("world")
ca_shap2 <- world %>% filter(name_long == 'Canada')
##########



#
fnpp <- list.files(
  pattern = glob2rx('*gpp_*.nc'),
  '/mnt/DataSpace/Data_pool/CMIP6_rebuild',
  full.names = T
)
freq0 <- unique(unlist(lapply(strsplit(fnpp, split = "_"), function(x)
  x[4])))
model_list <- lapply(strsplit(fnpp, split = "_"), function(x)
  x[5])
scenario <- lapply(strsplit(fnpp, split = "_"), function(x)
  x[6])
time <- lapply(strsplit(fnpp, split = "_"), function(x)
  x[9])
begin <- lapply(lapply(time, function(x)
  strsplit(x, split = '-|\\.')), function(x)
    x[[1]][[1]][1])
end <- lapply(lapply(time, function(x)
  strsplit(x, split = '-|\\.')), function(x)
    x[[1]][2])

fnpp_info <- data.frame(
  model = unlist(model_list),
  scenario = unlist(scenario),
  time_begin = unlist(begin),
  time_end = unlist(end)
)

library(dplyr)
library(tidyr)
fnpp_info2 <- fnpp_info %>% group_by(model, scenario) %>% summarise(
  model = model[1],
  scenario = scenario[1],
  begin = min(time_begin),
  end = max(time_end)
)
length(unique(fnpp_info$model))
model_use <- names(which(table(fnpp_info2$model) == 3)) # 14models

#
use_ind <- which(fnpp_info$model %in% model_use)
#
#get files that we need from all downloaded files
fnpp_file_use <- fnpp[use_ind]
npp_model_use <- model_use

#
#
#
#
#

flittersoil <- list.files(
  pattern = glob2rx('*fVegLitter_*.nc'),
  '/mnt/DataSpace/Data_pool/CMIP6_rebuild_flittersoil',
  full.names = T
)

#
rast(flittersoil[1])
freq0 <- unique(unlist(lapply(strsplit(flittersoil, split = "_"), function(x)
  x[5])))
model_list <- lapply(strsplit(flittersoil, split = "_"), function(x)
  x[6])
scenario <- lapply(strsplit(flittersoil, split = "_"), function(x)
  x[7])
time <- lapply(strsplit(flittersoil, split = "_"), function(x)
  x[10])
begin <- lapply(lapply(time, function(x)
  strsplit(x, split = '-|\\.')), function(x)
    x[[1]][[1]][1])
end <- lapply(lapply(time, function(x)
  strsplit(x, split = '-|\\.')), function(x)
    x[[1]][2])

flittersoil_info <- data.frame(
  model = unlist(model_list),
  scenario = unlist(scenario),
  time_begin = unlist(begin),
  time_end = unlist(end)
)

library(dplyr)
library(tidyr)
flittersoil_info2 <- flittersoil_info %>% group_by(model, scenario) %>% summarise(
  model = model[1],
  scenario = scenario[1],
  begin = min(time_begin),
  end = max(time_end)
)
length(unique(flittersoil_info$model))
model_use <- names(which(table(flittersoil_info2$model) == 3)) # 14models

#
use_ind <- which(flittersoil_info$model %in% model_use)
#
#get files that we need from all downloaded files
flittersoil_file_use <- flittersoil[use_ind]
flittersoil_model_use <- model_use
#
item = 'fVegLitter'
x = model_use[2]
x = 'EC-Earth3-Veg-LR'
scenario = 'ssp585'


#
#dir.create('/mnt/DataSpace/Data_pool/NPP_litter_ratio')
ca_soc_stock <- rast(
  '/mnt/DataSpace/Projects/Canada_C/processed_output/base_soc_kgm2_WGS84_1214.tif'
)
ref1 <- ca_soc_stock[[1]]
get_ratio <- function(x,
                      item = 'fVegLitter',
                      freq = 'Lmon',
                      begin = 2000,
                      end = 2300,
                      scenario) {
  litter_mods_hist <- grep(paste0(item, '_', 'Lmon', '_', x, '_', 'historical'),
                           x = flittersoil_file_use,
                           value = TRUE)
  litter_mods_ssp <- grep(paste0(item, '_', 'Lmon', '_', x, '_', scenario),
                          x = flittersoil_file_use,
                          value = TRUE)
  litter_mods0 <- c(litter_mods_hist, litter_mods_ssp)
  #litter_mods0 <- grep(paste0(item,'_',freq,'_',x,'_',scenario),x= if(scenario =='ssp245') {cmip_soil_245} else {cmip_soil_585},value=TRUE) #paste0 like cesm2 cesm2-waccm etc
  if (length(litter_mods0) > 0) {
    #if(scenario=='ssp585' & x=='CESM2-WACCM' & item=='cSoil' & freq=='Eyr') {ra0 <- rast(litter_mods0)[[-87]]} else {ra0 <- rast(litter_mods0)} # 2101 duplicate
    litter_ra0 <- rast(litter_mods0)
    if (x == c('CESM2-WACCM') &
        scenario == 'ssp585') {
      litter_ra0 <- c(litter_ra0, litter_ra0[[(nlyr(litter_ra0) - 11):nlyr(litter_ra0)]])
    } #extend one year
    if (x == c('TaiESM1') &
        scenario == 'ssp585') {
      litter_ra0 <- c(litter_ra0[[1:1980]], litter_ra0[[1981]], litter_ra0[[1981:nlyr(litter_ra0)]])
    } #extend one year 20150101
    
    if (ext(litter_ra0)[2] > 200) {
      litter_ra01 <- rotate(litter_ra0)
    } else {
      litter_ra01 <- litter_ra0
    }
    #litter_ra1 <- crop(litter_ra0,ext(-143,-53,43,86))
    crs(litter_ra01) <- 'epsg:4326'
    ####
    ####
    ####
    #
    litter_ra1 <- terra::crop(litter_ra01, ext(ca_shap2) + 3) #
    litter_ra2 <- mask(litter_ra1, ca_shap2)
    # litter_ra2 <- clamp(litter_ra2,lower = 10^-26) # newly added
    ind_use <- which(
      time(litter_ra2) >= as.Date('2001-01-01') &
        time(litter_ra2) <= as.Date('2300-12-31')
    )
    litter_ra3 <- litter_ra2[[ind_use]]
    litter_use <- tapp(litter_ra3, index = rep(1:(nlyr(litter_ra3) / (12 *
                                                                        10)), times = rep(12 * 10, nlyr(litter_ra3) / (12 * 10))), function(x)
                                                                          sum(x, na.rm = T))# * 30 * 3600 *24/20) # 20 year and 10000 scale,but here it doesn't matter
  } else {
    litter_use = 'NA'
  }
  
  
  npp_mods_hist <- grep(paste0('gpp_', 'Lmon', '_', x, '_', 'historical'),
                        x = fnpp_file_use,
                        value = TRUE)
  npp_mods_ssp <- grep(paste0('gpp_', 'Lmon', '_', x, '_', scenario),
                       x = fnpp_file_use,
                       value = TRUE)
  npp_mods0 <- c(npp_mods_hist, npp_mods_ssp)
  #npp_mods0 <- grep(paste0(item,'_',freq,'_',x,'_',scenario),x= if(scenario =='ssp245') {cmip_soil_245} else {cmip_soil_585},value=TRUE) #paste0 like cesm2 cesm2-waccm etc
  if (length(npp_mods0) > 0) {
    #if(scenario=='ssp585' & x=='CESM2-WACCM' & item=='cSoil' & freq=='Eyr') {ra0 <- rast(npp_mods0)[[-87]]} else {ra0 <- rast(npp_mods0)} # 2101 duplicate
    npp_ra0 <- rast(npp_mods0)
    if (x == c('CESM2-WACCM') &
        scenario == 'ssp585') {
      npp_ra0 <- c(npp_ra0, npp_ra0[[(nlyr(npp_ra0) - 11):nlyr(npp_ra0)]])
    } #lack one year
    if (x == c('TaiESM1') &
        scenario == 'ssp585') {
      npp_ra0 <- c(npp_ra0[[1:1980]], npp_ra0[[1981]], npp_ra0[[1981:nlyr(npp_ra0)]])
    } #lack 2015 0101
    if (ext(npp_ra0)[2] > 200) {
      npp_ra01 <- rotate(npp_ra0)
    } else {
      npp_ra01 <- npp_ra0
    }
    #npp_ra1 <- crop(npp_ra0,ext(-143,-53,43,86))
    crs(npp_ra01) <- 'epsg:4326'
    ####
    ####
    ####
    
    npp_ra1 <- terra::crop(npp_ra01, ext(ca_shap2) + 3) #
    npp_ra2 <- mask(npp_ra1, ca_shap2)
    npp_ra2 <- clamp(npp_ra2, lower = 10^-15) #
    npp_ind_use <- which(time(npp_ra2) >= as.Date('2001-01-01') &
                           time(npp_ra2) <= as.Date('2300-12-31'))
    npp_ra3 <- npp_ra2[[npp_ind_use]]
    npp_use <- tapp(npp_ra3, index = rep(1:(nlyr(npp_ra3) / (12 * 10)), times = rep(12 *
                                                                                      10, nlyr(npp_ra3) / (12 * 10))), function(x)
                                                                                        sum(x, na.rm = T))# * 30 * 3600 *24/20) # 20 year and 10000 scale,but here it doesn't matter
  } else {
    npp_use = 'NA'
  }
  
  time00 <- seq(2000, 2300, 10)
  
  litter_use1 = ifel(litter_use > npp_use, npp_use, litter_use)
  
  npp_soil_ratio <- litter_use1 / npp_use
  npp_soil_ratio1 <- crop(npp_soil_ratio, ref1)
  npp_soil_ratio1 <- focal(npp_soil_ratio1, w = 3, function(x)
    min(x, na.rm = T), na.policy = 'only')
  npp_soil_ratio2 <- resample(npp_soil_ratio1, ref1, threads = TRUE)
  npp_soil_ratio3 <- mask(npp_soil_ratio2, ref1)
  writeRaster(
    npp_soil_ratio3,
    paste0(
      '/mnt/DataSpace/Data_pool/fveglitter_gpp_ratio_from_fluxcom/',
      x,
      '_',
      scenario,
      '_2000-',
      time00[nlyr(npp_soil_ratio) + 1],
      '.tif'
    ),
    overwrite = T
  )
  
}


model_gpp <- c(
  "EC-Earth3-Veg-LR",
  "IPSL-CM6A-LR",
  "ACCESS-ESM1-5",
  "BCC-CSM2-MR",
  "CanESM5",
  "CanESM5-1",
  "CAS-ESM2-0",
  "CESM2-WACCM",
  "CMCC-CM2-SR5",
  "EC-Earth3-CC",
  "INM-CM4-8",
  "INM-CM5-0",
  "TaiESM1",
  "EC-Earth3-Veg",
  "MPI-ESM1-2-HR",
  "MPI-ESM1-2-LR"
)

model_pet <- c(
  "ACCESS-CM2",
  "EC-Earth3-Veg-LR",
  "ACCESS-ESM1-5",
  "CanESM5",
  "CanESM5-1",
  "CAS-ESM2-0",
  "FGOALS-g3",
  "GFDL-ESM4",
  "IPSL-CM6A-LR",
  "MIROC6",
  "MPI-ESM1-2-HR",
  "MRI-ESM2-0",
  "EC-Earth3-CC",
  "MPI-ESM1-2-LR"
)
extrass <- c("EC-Earth3-Veg-LR", "INM-CM4-8", "INM-CM5-0")

model_2100_new <- c(
  "ACCESS-CM2",
  "CanESM5",
  "CMCC-ESM2",
  "EC-Earth3-Veg",
  "FIO-ESM-2-0",
  "INM-CM4-8",
  "IPSL-CM6A-LR",
  "MIROC6",
  "MPI-ESM1-2-HR",
  "MPI-ESM1-2-LR",
  "MRI-ESM2-0",
  "ACCESS-ESM1-5",
  "EC-Earth3-Veg-LR",
  "INM-CM5-0" ,
  "BCC-CSM2-MR",
  "CanESM5-1",
  "CAS-ESM2-0",
  "EC-Earth3-CC"
)
model_2300_clim <- c("CanESM5",
                     "EC-Earth3-Veg",
                     "IPSL-CM6A-LR",
                     "MRI-ESM2-0",
                     "ACCESS-ESM1-5")
model_soil_water <- c(
  "ACCESS-CM2",
  "ACCESS-ESM1-5",
  "CanESM5"  ,
  "CMCC-ESM2"  ,
  "EC-Earth3-Veg" ,
  "EC-Earth3-Veg-LR" ,
  "IPSL-CM6A-LR",
  "MIROC6",
  "MPI-ESM1-2-HR" ,
  "MPI-ESM1-2-LR",
  "MRI-ESM2-0"
)
model_solar_2300 <- c('IPSL-CM6A-LR',
                      'ACCESS-ESM1-5',
                      'CanESM5',
                      'CESM2-WACCM',
                      'EC-Earth3-Veg') #
fVegLitter_m <- c(
  'BCC-CSM2-MR',
  'CanESM5',
  'CanESM5-1',
  'CESM2-WACCM',
  'CMCC-CM2-SR5',
  'CMCC-ESM2',
  'EC-Earth3-CC',
  'EC-Earth3-Veg',
  'EC-Earth3-Veg-LR',
  'IPSL-CM6A-LR',
  'KIOST-ESM',
  'NorESM2-LM',
  'NorESM2-MM',
  'TaiESM1'
)

model_fVegLitter <- intersect(model_gpp, fVegLitter_m)


final_model <- model_fVegLitter
#
for (m in final_model) {
  for (w in c('ssp245', 'ssp585')) {
    get_ratio(x = m, scenario = w)
  }
}

#
library(terra)
f1 <- list.files('/mnt/DataSpace/Data_pool/fveglitter_gpp_ratio_from_fluxcom/',
                 full.names = T)
f2 <- rast(f1[3])
plot(f2[[1]])
library(terra)
ff <- list.files('/mnt/DataSpace/Data_pool/01cmip6_processed_output_EC_GPP/',
                 full.names = T)
landcover <- rast(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/LAND_cover/land_cover_mask_lt15_keep_origin_1km_w84_1214.tif'
) #> 15  #1214 with new extend as same as bd
ca_pure_area84 <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/ca_cmip6_area_WGS84.tif') #20231214 new created

ff <- grep('CanESM5', ff, value = T)
fft <- rast(ff[1])
plot(fft)

twsd <- fft * f2[[1]]
par(mfrow = c(1, 3))
plot(twsd, main = 'fgpp')
#title(main='fgpp')
t1 <- mask(basedata0[[3]], mask = ca_soc_use)
t4 <- mask(basedata0[[2]], mask = ca_soc_use)
plot(t1, main = 'card')
plot(t4, main = 'lu')


sdf <- mask(twsd, mask = landcover) * ca_pure_area84
sdf <- twsd * ca_pure_area84
global(sdf, function(x)
  sum(x, na.rm = T))
























