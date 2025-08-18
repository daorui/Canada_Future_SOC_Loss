

# rsds
# refer to /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/001code/00_1complementary_experiments.R
library(terra)

library(lubridate)

library(sf)

library(dplyr)
unlink(
  "/mnt/DataSpace/Data_pool/CMIP6_rebuild/rsds_ImonAnt_CESM2-WACCM_historical_r1i1p1f1_gn_185001-201412.nc"
)
unlink(
  "/mnt/DataSpace/Data_pool/CMIP6_rebuild/rsds_ImonGre_CESM2-WACCM_historical_r1i1p1f1_gn_185001-201412.nc"
)
f1 <- list.files(
  pattern = 'nc',
  '/mnt/DataSpace/Data_pool/CMIP6_rebuild',
  full.names = T
)
ca_shap1 <- vect(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/2021digitallpr_000b21a_e/lpr_000b21a_e.shp'
)
ca_shp2 <- terra::project(ca_shap1, 'epsg:4326')

ca_covas <- rast(
  '/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global_re_align_2300_d1.tif'
) # digital = 1 ,this version meet the extent to 2300
tdaye <- ca_covas[[1]]
model = 'CESM2-WACCM'
var1 = 'rsds'
scenario = 'ssp585'
dir.create('/mnt/DataSpace/Data_pool/01_rsds_output_10')
get_bc <- function(var1, model, scenario) {
  time00 <- paste0(seq(2001, 2291, 10), '-', seq(2010, 2300, 10))#
  hist0 <- rast(grep(paste0(model, '_historical'), grep(var1, f1, value = T), value = T))
  crs(hist0) <- crs(tdaye)
  if (model == 'BCC-CSM2-MR') {
    time(hist0) <- seq(as.Date('1850-01-16'), as.Date('2014-12-16'), by = 'month')
  }
  if (xmax(hist0) > 200) {
    hist_t2 <- rotate(hist0)
  } else {
    hist_t2 <- hist0
  }
  if (var1 == 'rsds') {
    hist_ind <- which(time(hist_t2) %in% seq(as.Date('1970-01-01'), as.Date('2000-12-31'), by =
                                               'day'))
  } else {
    hist_ind <- which(time(hist_t2) %in% seq(as.Date('2004-01-01'), as.Date('2014-12-31'), by =
                                               'day'))
  }
  hist_t3 <- hist_t2[[hist_ind]]
  hist_t4 <- app(tapp(hist_t3, index = rep(1:(
    length(hist_ind) / 12
  ), rep(
    12, length(hist_ind) / 12
  )), fun = sum), fun = mean)
  #
  f_ssp <- grep(paste0(model, '_', scenario), grep(var1, f1, value = T), value = T)
  if (length(f_ssp) < 1) {
    diff_ssp <- NULL
  } else {
    ssp0 <- rast(f_ssp)
    crs(ssp0) <- crs(tdaye)
    if (model == 'BCC-CSM2-MR') {
      time(ssp0) <- seq(as.Date('2015-01-16'), as.Date('2100-12-16'), by = 'month')
    }
    if (xmax(ssp0) > 200) {
      ssp_t2 <- rotate(ssp0)
    } else {
      ssp_t2 <- ssp0
    }
    ssp_t2 <- c(hist_t2, ssp_t2)
    # ind <- which(time(ssp_t2) %in% seq(as.Date('2021-01-01'),as.Date('2300-12-31'),by='day')) # some febrary is 2-15  and others are 01-16
    if (model == 'CESM2-WACCM' & scenario == 'ssp585') {
      ind <- which(time(ssp_t2) %in% seq(as.Date('2001-01-01'), as.Date('2300-12-31'), by =
                                           'day'))
      ind <- c(ind, tail(ind, 12)) #
      ssp_t3 <- ssp_t2[[ind]]
      time(ssp_t3) <- seq(as.Date('2001-01-15'), as.Date('2300-12-15'), by =
                            'month')
    } else {
      ind <- which(time(ssp_t2) %in% seq(as.Date('2001-01-01'), as.Date('2300-12-31'), by =
                                           'day'))
      ssp_t3 <- ssp_t2[[ind]]
    } #
    ssp_t4 <- tapp(tapp(ssp_t3, index = rep(1:(
      length(ind) / (12)
    ), rep(
      12, length(ind) / (12)
    )), fun = sum),
    index = rep(1:(length(ind) / (12 * 10)), rep(10, length(ind) / (12 * 10))),
    fun = mean) #
    diff_ssp <- ssp_t4 / (hist_t4 + 1 * 10^-20)  # quantile(values(crop(hist_t4,ext(ca_shp2) +1 )),0.5,na.rm=T)
    bc_ssp <- crop(x = diff_ssp, y = ext(ca_shp2) + 1)
    bc_ssp1 <- clamp(bc_ssp, upper = 10, values = T) #
    diff_ssp_1km <- resample(bc_ssp1, tdaye, method = 'bilinear', threads =
                               TRUE)
    if (var1 == 'gpp') {
      f_ssp_1km  <- diff_ssp_1km * ca_covas[[4]]
      #writeRaster(f_ssp_1km,paste0('/mnt/DataSpace/Data_pool/cmip6_complementary/processed_output/',model,'_',scenario,'_',var1,'_',year(min(time(ssp0))),'-',year(max(time(ssp0))),'.tif'))
    } else if (var1 == 'npp') {
      f_ssp_1km  <- diff_ssp_1km * ca_covas[[3]]
      #writeRaster(f_ssp_1km,paste0('/mnt/DataSpace/Data_pool/cmip6_complementary/processed_output/',model,'_',scenario,'_',var1,'_',year(min(time(ssp0))),'-',year(max(time(ssp0))),'.tif'))
    } else if (var1 == 'rsds') {
      f_ssp_1km  <- diff_ssp_1km  * ca_covas[[26]]
      #writeRaster(f_ssp_1km,paste0('/mnt/DataSpace/Data_pool/cmip6_complementary/processed_output/',model,'_',scenario,'_',var1,'_',year(min(time(ssp0))),'-',year(max(time(ssp0))),'.tif'))
    }
  }
  if (!is.null(diff_ssp)) {
    lapply(1:nlyr(f_ssp_1km), function(x)
      writeRaster(
        f_ssp_1km[[x]],
        filename = paste0(
          '/mnt/DataSpace/Data_pool/01_rsds_output_10/',
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
      ))
  }
}

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

for (m in model_gpp) {
  for (i in c('rsds')) {
    for (j in c('ssp245', 'ssp585')) {
      get_bc(var1 = i,
             model = m,
             scenario = j)
    }
  }
}
