#
#

## step 1 spline
## step 1 spline
## step 1 spline
## step 1 spline
## step 1 spline
## step 1 spline


library(readr)
library(dplyr)
library(sp)
library(sf)
library(data.table)
library(aqp)
#
if (R.Version()$os == "linux-gnu") {
  path00 <- "/mnt/File0"
} else {
  path00 <- "R:"
}
# spline to specific layers
datadirect_name <- '/DAAATAAA/Data_collection/MY_global_Temp_resilence_C/00Rawdata'
use_layer3 <- readr::read_csv(paste0(
  path00,
  datadirect_name,
  '/global_organic_carbon_with_extension1.csv'
))#
length(unique(use_layer3$profile_id))
#delete 0 value
use_layer3 <- use_layer3 %>% filter(orgc_value_avg != 0)
length(unique(use_layer3$profile_id))
#
t0 <- data.frame(table(use_layer3$profile_id))
use_layer3_0 <- use_layer3 %>% filter(profile_id %in% t0[which(t0$Freq >= 3), ]$Var1)
aa0 <- use_layer3_0 %>% group_by(profile_id) %>% summarise(x = sd(orgc_value_avg))
aaa0 <- aa0[which(aa0$x == 0), 1]
use_layer3_1 <- use_layer3_0 %>% filter(!profile_id %in% aaa0$profile_id)

use_layer3 <- use_layer3_1




library(aqp)
depths(use_layer3) <- profile_id ~ upper_depth + lower_depth


library(ithir)
use_data0 <- ithir::ea_spline(
  use_layer3,
  var.name = 'orgc_value_avg' ,
  lam = 0.1,
  d = c(0, 5, 10, 15, seq(20, 100, 10), seq(120, 300, 30))
) # ** 20241212 update 10 cm interval

#
t1 <- dplyr::distinct(use_layer3_1[, c(1:6)])
############---------------------------------------------
num1 <- list()
dup_datase_id <- names(which(table(t1$profile_id) > 1)) # more than one choose out
for (i in dup_datase_id) {
  num1[[which(dup_datase_id %in% i)]] <- which(t1$profile_id %in% i)[-1]
}
unlist(num1)
distinct_data0 <- t1[-unlist(num1), ]
#
use_data1 <- dplyr::left_join(use_data0[[1]], distinct_data0, by = c('id' =
                                                                       'profile_id')) #
saveRDS(
  use_data1,
  paste0(
    path00,
    datadirect_name,
    '/third_splined_data_to3m_use_20241212_10cmintervals.rds'
  )
) #


# step 2
# step 2
# step 2
# step 2
# step 2


library(sp)
library(rgdal)
library(gdalUtils)
library(sf)
library(raster)
library(terra)
library(dplyr)
library(spData)
library(readr)
library(data.table)
#
if (R.Version()$os == "linux-gnu") {
  #.Platform$OS.type =='unix'
  path00 <- "/mnt/File0"
} else {
  path00 <- "R:"
}
#
datadirect_name <- '/DAAATAAA/Data_collection/MY_global_Temp_resilence_C/00Rawdata'

use_data1 <- readr::read_rds(
  paste0(
    path00,
    datadirect_name,
    '/third_splined_data_to3m_use_20241212_10cmintervals.rds'
  )
) #



data_use00 <- fread(paste0(
  path00,
  datadirect_name,
  '/global_organic_carbon_with_extension1.csv'
))
who0 <- data_use00 %>% filter(latitude %in% no2$latitude)

# make a triangle simple feature
s2 <- list(
  c(0, 0),
  c(0, 10),
  c(10, 0),
  c(10, 10),
  c(0, 20),
  c(10, 20),
  c(0, 30),
  c(10, 30),
  c(0, 60),
  c(10, 60),
  c(0, 70),
  c(10, 70),
  c(0, 88),
  c(10, 88)
)
t1 <- lapply(s2, function(x)
  st_point(x))
s3.sf <- st_sfc(t1, crs = 4326)
s2.pts = st_sf(ID = "tr", s3.sf)
s3.pts <- st_transform(s2.pts, crs = '+proj=aeqd + lat_0=0 + lon_0=0 + x_0=0 + y_0=0')
s4.pts <- st_transform(s2.pts, crs = 4087) #
data("world")

w84 <- st_distance(s2.pts)
aeqd <- st_distance(s3.pts)
eqcl <- st_distance(s4.pts)


#
#aggregate
glob_raster <- rast(gsub(
  "\\\\",
  "/",
  paste0(
    path00,
    "\\DAAATAAA\\Global aridity and PET\\7504448 (1)\\Global-AI_ET0_annual_v3\\Global-AI_ET0_v3_annual\\ai_v3_yr.tif"
  )
))
head(use_data1)

use_data2 <- st_as_sf(use_data1,
                      coords = c('longitude', 'latitude'),
                      crs = st_crs(4326))

bg_extract <- extract(glob_raster, use_data2, cells = T, xy = T)
length(which(table(bg_extract$cell) > 1))  #>2 #8110
use_data3 <- cbind(use_data1, bg_extract[, c(1, 3:5)])
library(dplyr)

use_data4 <- use_data3 %>% group_by(cell) %>% summarise_at(vars('0-5 cm':'270-300 cm', 'latitude', 'longitude', 'x', 'y'),
                                                           mean,
                                                           na.rm = T)

use_data5 <- use_data3 %>% group_by(cell) %>% summarise_at(vars('country_id', 'country_name'), function(x)
  x[1])
use_data6 <- left_join(use_data4, use_data5, by = 'cell')
use_data6$ID_no <- 1:nrow(use_data6)

write.csv(
  use_data6,
  paste0(
    path00,
    datadirect_name,
    '/aggregated_20241212_3m_10cminterval.csv'
  )
)



#########step3
#########step3
#########step3
#########step3
#########step3
#########step3.


#load data
use_data6 <- read.csv(paste0(
  path00,
  datadirect_name,
  '/aggregated_20241212_3m_10cminterval.csv'
))
names(use_data6)[3:21] <- c('X5cm', 'X10cm', 'X15cm', 'X20cm', 'X30cm', paste0('X', c(seq(40, 100, 10), seq(120, 300, 30)), 'cm'))
use_data7 <- st_as_sf(use_data6,
                      coords = c('longitude', 'latitude'),
                      crs = st_crs(4326))
library(spData)
data('world')

system.time(extractdata0 <- terra::extract(covs_list, use_data7))
# for reference
use_data18 <- readr::read_csv(
  '/mnt/File0/DAAATAAA/Data_collection/MY_global_Temp_resilence_C/01Output/04base4model_all.csv'
)
#

data("world")
world_sp <- as(world, 'Spatial')
pts_sp <- as(use_data7, 'Spatial')
continent0 <- sp::over(pts_sp, world_sp['continent'])

use_data10 <- cbind(use_data6, continent0)
#
which(is.na(use_data10$continent))
use_data10$continent <- world$continent[st_nearest_feature(use_data7, world['continent'])]
use_data11 <- cbind(use_data10, extractdata0)
names(use_data11)[31:88] <- names(use_data18)[c(15:65, 67:73)]
use_data12 <- use_data11[, -30]
write.csv(
  use_data12,
  '/mnt/File0/DAAATAAA/Data_collection/MY_global_Temp_resilence_C/01Output/07base_for_model_all_20241212_UPDATE10CMinterval.csv'
)

#########step4
#########step4
#########step4
#########step4
#########step4
library(sp)
library(rgdal)
library(gdalUtils) #devtools:::install_github("gearslaboratory/gdalUtils")
library(raster)
library(terra)
library(spDataLarge)
library(spData)
library(sf)
library(s2)
library(dplyr)
library(nngeo)
library(recipes)
library(data.table)
library(tidyr)
library(xgboost)
library(caret)
library(Matrix)
library(ranger)

setwd('/mnt/File0/DAAATAAA/Data_collection/MY_global_Temp_resilence_C')
data_use <- read.csv('./01Output/07base_for_model_all_20241212_UPDATE10CMinterval.csv')
ca_data_use <- data_use %>% filter(country_name == 'Canada')


ca_data_sf0 <- ca_data_use %>%
  st_as_sf(., coords = c('longitude', 'latitude'), crs = 4326)


library(terra)
library(sf)

########  SP CV

library(spatialsample)
ca_data_u <- ca_data_sf0
#
set.seed(5234)
folds <- spatial_clustering_cv(ca_data_u, v = 100)


fold_n <- function(x) {
  assessment(x)$X.1
}

fold_n1 <- purrr::map(folds$splits, fold_n)

#
new_group <- 1:nrow(ca_data_u)
for (i in 1:nrow(folds)) {
  new_group[which(ca_data_u$X.1 %in% fold_n1[[i]])] <- i
}
#

ca_data_u$fold_col <- new_group


library(dplyr)
library(sf)
library(spData)
data("world")
plot(world %>% filter(name_long == 'Canada') %>% st_geometry())


set.seed(7180)
test4 <- sample(1:100, 5) #
plot(st_geometry(ca_data_u[which(new_group %in% sample(1:100, 5, replace =
                                                         F)), ]), col = 'red', add = T)
plot(st_geometry(ca_data_u), add = T)
##---------------------------------------------------
d1 <- ca_data_u[!(new_group %in% test4), ]
testdata_u <- setdiff(ca_data_u, d1)
saveRDS(
  testdata_u,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/ca_soc_5per_4test_20241212.rds'
)
saveRDS(
  d1,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/ca_soc_95per_4train_20241212.rds'
)
# L143


##########
##step 5
##step 5
##step 5
##step 5
##step 5
##step 5
ca_covas <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global1214.tif')
gpp_hist <- rast(
  '/mnt/DataSpace/Projects/Canada_C/processed_output/GPP_historical_from_fluxcom_1981_2000_gCm2year.tif'
)
ca_covas[[4]] <- gpp_hist
names(ca_covas)[4] <- 'GPP4'
#

ca_peat_sf <- readRDS(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/Sasha_peat_20240325_canada_sf_peatland_SOC_20241212update_10cminterval.rds'
)
ca_peat_cov <- terra::extract(ca_covas, ca_peat_sf)
names(ca_peat_cov)[29] <- 'Pop_den28'
ca_peat_use <- cbind(ca_peat_sf, ca_peat_cov)
ca_peat_use$Pop_den28 <- as.integer(ca_peat_use$Pop_den28)
ca_peat_use1 <- ca_peat_use[, c(2:20, 23:81)]
names(ca_peat_use1)[1:19] <- c('X5cm', 'X10cm', 'X15cm', 'X20cm', 'X30cm', paste0('X', c(seq(40, 100, 10), seq(120, 300, 30)), 'cm'))

######
d5 <- readRDS(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/ca_soc_train_wosis_and_only_cheryl20241212_updateinterval.rds'
)
d5$Peat <- c(rep(0, 5088), rep(1, 20))
ca_peat_use1$Peat <- 1
names(ca_peat_use1) <- names(d5)
d6 <- rbind(d5, ca_peat_use1)
saveRDS(
  d6,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/ca_soc_train_wosis_and_cheryl_and_sasha_20241212_updateinterval.rds'
)

############################################

###step 6
###step 6
###step 6
###step 6
###step 6
#
#########
library(spatialsample)
library(dplyr)
library(tidyr)
library(h2o)
library(sf)
library(terra)

#
d6 <- readRDS(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/ca_soc_train_wosis_and_cheryl_and_sasha_20241212_updateinterval.rds'
) # wosis - 5% + cheryl + peatcan
#
ca_test <- readRDS(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/ca_soc_5per_4test_20241212.rds'
)
ca_test$Peat <- 0
t_test <- ca_test[, c(5:23, 30:87, 90, 89)]
identical(names(d6), names(t_test))
d7 <- rbind(d6, t_test)
names(d7)[1:19] <- c('X5cm', 'X10cm', 'X15cm', 'X20cm', 'X30cm', paste0('X', c(seq(40, 100, 10), seq(120, 300, 30)), 'cm')) #
d6 <- d7
plot(world %>% filter(name_long == 'Canada') %>% st_geometry())
plot(st_geometry(d6),
     add = T,
     pch = 16,
     cex = 0.6)
gpp_hist <- rast(
  '/mnt/DataSpace/Projects/Canada_C/processed_output/GPP_historical_from_fluxcom_1981_2000_gCm2year.tif'
)
gpp_new <- terra::extract(x = gpp_hist, y = d6)
d6$GPP4 <- gpp_new$lyr.1 #
d1 <- d6
d1$X.1 = 1:nrow(d1)

test_ind <- seq.int(1, nrow(d1), by = ceiling(nrow(d1) / 150))
data_test <- d1[test_ind, ]
data_train <- d1[-test_ind, ]

#


saveRDS(
  data_train,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/ca_soc_traindataset_20241212.rds'
)
saveRDS(
  data_test,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/ca_soc_testdataset_20241212.rds'
)


library(spData)
data('world')
plot(world %>% filter(name_long == 'Canada') %>% st_geometry())
plot(data_train[0], add = T)
plot(data_test[0],
     col = 'red',
     pch = 15,
     add = T)


###train data prep
###train data prep
###train data prep
path = '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/'
d100 <- readRDS(paste0(path, 'ca_soc_traindataset_20241212.rds'))
d6 <- d100[, 1:78]
d6[1:19] <- lapply(d6[1:19], function(x) {
  x[x == 0 | x == -9999] <- NA
  return(x)
})
#
gpp_hist <- rast(
  '/mnt/DataSpace/Projects/Canada_C/processed_output/GPP_historical_from_fluxcom_1981_2000_gCm2year.tif'
)
gpp_new <- terra::extract(x = gpp_hist, y = d100)
d6$GPP4 <- gpp_new$lyr.1 #


peat_ca <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/peatland_area1214.tif') #change name from peatland_area to peatland1214
peat_ca <- ifel(is.na(peat_ca), 0, peat_ca)
peat_new <- terra::extract(x = peat_ca, y = d100)
ind_wosis <- which(d6$Peat == 0)
ind_peat <- which(d6$Peat == 1)
d6$Peat[ind_wosis] <- peat_new$PEAT_PER[ind_wosis]
d6$Peat[ind_peat] <- peat_new$PEAT_PER[ind_peat] #

#
for (i in 1:nrow(d6)) {
  #
  row <- d6[i, 1:19]
  if (is.data.frame(row)) {
    row <- as.numeric(unlist(row))[-c(20:21)]
  }
  
  non_na_indices <- which(!is.na(row))
  
  if (length(non_na_indices) > 0) {
    last_index <- tail(non_na_indices, n = 1)
    d6[i, last_index] <- NA
  }
}





d6$Peat[ind_peat[which(d6$Peat[ind_peat] < 5 &
                         d6$X20cm[ind_peat] > 300 &
                         d6$X10cm[ind_peat] > 300 &
                         d6$X50cm[ind_peat] > 300 &
                         st_coordinates(d6)[ind_peat, 2] < 65)]] <- 99

d6$Peat[ind_peat[which(d6$Peat[ind_peat] > 60 &
                         (is.na(d6$X40cm[ind_peat]) |
                            d6$X40cm[ind_peat] < 300))]] <- 50
indwosis_no1peat <- which(d6$Peat[ind_wosis] > 60 &
                            (is.na(d6$X40cm[ind_wosis]) |
                               d6$X15cm[ind_wosis] < 300 |
                               d6$X40cm[ind_wosis] < 300)) #
d6$Peat[ind_wosis[indwosis_no1peat]] <- 45
indwosis_no2peat <- which(d6$Peat[ind_wosis] > 30 &
                            d6$Peat[ind_wosis] <= 60 &
                            d6$X40cm[ind_wosis] < 150) #
d6$Peat[ind_wosis[indwosis_no2peat]] <- 25


####

d6[ind_peat, 1:19] <- lapply(d6[ind_peat, 1:19], function(x) {
  if (is.numeric(x)) {
    ifelse(x < 150, NA, x)
  } else {
    x
  }  #
})

peat_peat <- which(d6$Peat > 60)
d6[peat_peat, 1:19] <- lapply(d6[peat_peat, 1:19], function(x) {
  if (is.numeric(x)) {
    ifelse(x < 300, NA, x)
  } else {
    x
  }
})

d6$Peat <- cut(
  d6$Peat,
  breaks = c(0, 5, 30, 60, 100),
  labels = c(5, 30, 60, 100),
  include.lowest = TRUE,
  right = TRUE
)
table(d6$Peat)
# 5   30   60  100
#2892 2614  503  425
reclass_table <- matrix(c(0, 5, 5, 5, 30, 30, 30, 60, 60, 60, 100, 100),
                        ncol = 3,
                        byrow = TRUE)


peat_reclassified <- classify(peat_ca, reclass_table)



d1 <- d6
d1$X.1 = 1:nrow(d1)
#end from all #

#
library(spatialsample)
set.seed(1239)
folds <- spatial_clustering_cv(d1, v = 8)
fold_n <- function(x) {
  assessment(x)$X.1
}
fold_n1 <- purrr::map(folds$splits, fold_n)

new_group <- 1:nrow(d1)
for (i in 1:nrow(folds)) {
  new_group[which(d1$X.1 %in% fold_n1[[i]])] <- i
}
d1$fold_col <- new_group
unique(d1$fold_col)
autoplot(folds, show.legend = T)

#

#
item0 <- paste0("X", c(5, 10, 15, 20, 30, seq(40, 100, 10), seq(120, 300, 30)), "cm")

covs0 <- c(
  'NDVI1',
  'EVI2',
  'NPP3',
  'GPP4',
  'PET5',
  'AI6',
  'BIO1',
  'BIO10',
  'BIO11',
  'BIO12',
  'BIO13',
  'BIO14',
  'BIO15',
  'BIO16',
  'BIO17',
  'BIO18',
  'BIO19',
  'BIO2',
  'BIO3',
  'BIO4',
  'BIO5',
  'BIO6',
  'BIO7',
  'BIO8',
  'BIO9',
  'Solar26',
  'So_seasonal27',
  'Pop_den28',
  'HFP29',
  'Vrm30',
  'Tcurv31',
  'Pcurv32',
  'Ele33',
  'Slope34',
  'Asp_c35',
  'Asp_s36',
  'Est37',
  'Nor38',
  'Roug39',
  'TPI40',
  'TRI41',
  'Dx42',
  'Dxx43',
  'Dy44',
  'Dyy45',
  'Maj48',
  'lith.8',
  'lith.24',
  'WTD85',
  'smap',
  'SUMAP',
  'Ref_band1',
  'ref_band2',
  'ref_band3',
  'ref_band4',
  'ref_band5',
  'ref_band6',
  'ref_band7',
  'Peat'
)
added_cov <- 'fold_col'


d10 <-
  d1 %>% dplyr::select (all_of(item0), all_of(covs0), all_of(added_cov)) %>% #
  tidyr::pivot_longer(
    cols = paste0("X", c(
      5, 10, 15, 20, 30, seq(40, 100, 10), seq(120, 300, 30)
    ), "cm"),
    names_to = 'Depth',
    names_transform = list(
      Depth = function(x)
        gsub('X|cm', '', x) %>%
        as.integer()
    ),
    values_drop_na = TRUE,
    values_to = 'SOC'
  ) %>% mutate(
    Depth = case_when(
      Depth == '5' ~ 2.5,
      Depth == '10' ~ 7.5,
      Depth == '15' ~ 12.5,
      Depth == '20' ~ 17.5,
      Depth == '30' ~ 25,
      Depth == '40' ~ 35,
      Depth == '50' ~ 45,
      Depth == '60' ~ 55,
      Depth == '70' ~ 65,
      Depth == '80' ~ 75,
      Depth == '90' ~ 85,
      Depth == '100' ~ 95,
      Depth == '120' ~ 110,
      Depth == '150' ~ 135,
      Depth == '180' ~ 165,
      Depth == '210' ~ 195,
      Depth == '240' ~ 225,
      Depth == '270' ~ 255,
      Depth == '300' ~ 285,
    )
  ) %>%
  mutate_all(list( ~ replace(., . == c(-9999), NA))) %>%
  filter(
    !!as.name('SOC') != 0 &
      !!as.name('SOC') < 580 & !!as.name('SOC') > 0.01 &
      !is.na(!!as.name('SOC'))
  ) %>%
  mutate(
    SOC = log(SOC * 100)^2 ,
    NDVI1 = NDVI1 / 10000,
    EVI2 = EVI2 / 10000,
    NPP3 = NPP3 / 10000,
    GPP4 = GPP4 / 1000,
    PET5 = log(PET5),
    #
    AI6 = AI6 / 10000,
    BIO4 = BIO4 / 100,
    Solar26 = Solar26 / 10000,
    So_seasonal27 = So_seasonal27 / 10000,
    Vrm30 = Vrm30 * 10^4,
    #
    Tcurv31 = Tcurv31 * 10^4,
    #
    Pcurv32 = Pcurv32 * 10^4,
    #
    Ref_band1 = Ref_band1 / 10000,
    ref_band2 = ref_band2 / 10000,
    ref_band3 = ref_band3 / 10000,
    ref_band4 = ref_band4 / 10000,
    ref_band5 = ref_band5 / 10000,
    ref_band6 = ref_band6 / 10000,
    ref_band7 = ref_band7 / 10000,
    Maj48 = factor(Maj48, levels = 1:10),
    
    lith.8 = factor(lith.8, levels = 1:7),
    lith.24 = factor(lith.24, levels = 1:15)
  ) %>% na.omit()
#
#
d10[, 10:17] <- log(st_drop_geometry(d10[, 10:17]))
d11 <- d10
d11$Pop_den28 <- as.integer(d11$Pop_den28)
str(d11)

coord00 <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_coordinates.tif')
sedimentary_thick <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/Canada_sedimentary_thickness.tif')
xy_thick <- terra::extract(x = c(coord00, sedimentary_thick), y = d10)
d11$lon <- xy_thick$lon / 10^6
d11$lat <- xy_thick$lat / 10^6
d11$sed_thick <- xy_thick$`average_soil_and_sedimentary-deposit_thickness`
## 20240726 add xy ## 20240726 add xy ## 20240726 add xy

d12 <- sf::st_drop_geometry(d11)
saveRDS(
  d11,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/SOC_ca_distribute1214_updated_20240618_add_sasha20241212.rds'
)
saveRDS(
  d12,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/SOC_ca_distribute1214_updated_20240730_add_sasha20241212.rds'
)

#before

library(sf)
library(terra)
library(dplyr)
library(tidyr)
d110 <- readRDS(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/SOC_ca_distribute1214_updated_20240618_add_sasha20241212.rds'
)
# select data
d11 <- d110
library(h2o)
h2o.init(nthreads = 36, max_mem_size = '300G')
d1_hdf <- as.h2o(x = sf::st_drop_geometry(d11))
set.seed(1199)
system.time(
  h_auto <- h2o.automl(
    x = setdiff(names(d1_hdf), c('SOC', 'fold_col', 'NDVI1', 'EVI2', 'NPP3')) ,
    y = 'SOC',
    # 20240618 remove NDVI EVI NPP
    training_frame = d1_hdf,
    exclude_algos = c("DeepLearning", "GLM"),
    max_models = NULL,
    fold_column = 'fold_col',
    stopping_metric =  'RMSE',
    keep_cross_validation_predictions = TRUE,
    #
    seed = 1807
  )
)   #done 3580.671s



#  next
lb <- h2o.get_leaderboard(object = h_auto, extra_columns = "ALL")
print(lb, n = nrow(lb))
nidaye <- as.data.frame(lb)
saveRDS(
  nidaye,
  '/mnt/DataSpace/Projects/Canada_C/models/soc_h2o_leaderboad0619_add_wosis_cheryl_sasha_remove_npp_ndvi_evi_add_xy_sed_thick_20241215_EXCLUDE_Deepl_GLM_alldata_minus_nospcv_20fcv_with_5factorPeatall_reclean_fv300.rds'
)#2024.1219 add some peatcan to 100% peat

nidaye0 <- readRDS(
  '/mnt/DataSpace/Projects/Canada_C/models/soc_h2o_leaderboad0619_add_wosis_cheryl_sasha_remove_npp_ndvi_evi_add_xy_sed_thick.rds'
)
unique(nidaye$algo)

print(lb, n = nrow(lb))
m1 <- h_auto@leader
# this is equivalent to
# m <- h2o.get_best_model(h_auto)

# Get the best model using a non-default metric
#m <- h2o.get_best_model(h_auto, criterion = "RMSE")

# Get the best XGBoost model using default sort metric
xgb2 <- h2o.get_best_model(h_auto, algorithm = "xgboost")

# Get the best GBM model, ranked by logloss
gbm1 <- h2o.get_best_model(h_auto, algorithm = "GBM", criterion = "RMSE")

# Get the best randomforest model, ranked by logloss
rf3 <- h2o.get_best_model(h_auto, algorithm = "DRF", criterion = "RMSE")

# Get the best deepleardning model, ranked by logloss
deeplearn4 <- h2o.get_best_model(h_auto, algorithm = "DeepLearning", criterion = "RMSE")
glm <- h2o.get_best_model(h_auto, algorithm = 'GLM', criterion = 'RMSE')
# # Get the best stackensemble model, ranked by ,also by what we want , best of Family
stacken_best_family <- h2o.get_best_model(h_auto, algorithm = "StackedEnsemble", criterion = "RMSE")
dir.create('/mnt/DataSpace/Projects/Canada_C/models_after_EC/')
h2o.saveModel(
  object = m1,
  '/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble',
  export_cross_validation_predictions = TRUE
) #
h2o.saveModel(
  object = xgb2,
  '/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb',
  export_cross_validation_predictions = TRUE
)
h2o.saveModel(
  object = gbm1,
  '/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm',
  export_cross_validation_predictions = TRUE
)
h2o.saveModel(
  object = rf3,
  '/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf',
  export_cross_validation_predictions = TRUE
)
h2o.saveModel(
  object = deeplearn4,
  '/mnt/DataSpace/Projects/Canada_C/models_after_EC/deeplearn',
  export_cross_validation_predictions = TRUE
)
h2o.saveModel(
  object = stacken_best_family,
  '/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble',
  export_cross_validation_predictions = TRUE
) #
h2o.saveModel(
  object = glm,
  '/mnt/DataSpace/Projects/Canada_C/models_after_EC/glm',
  export_cross_validation_predictions = TRUE
) #



###
##20241219 add SOC + real peat
stack_m <- h2o.loadModel(
  "/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_5_AutoML_2_20241219_120448"
)
xgb2 <- h2o.loadModel(
  "/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_2_20241219_120448_model_1"
)
gbm1 <- h2o.loadModel(
  "/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_2_AutoML_2_20241219_120448"
)
rf3 <- h2o.loadModel(
  "/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/XRT_1_AutoML_2_20241219_120448"
)


# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_3_AutoML_1_20241219_173724"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_1_20241219_173724_model_4"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_2_AutoML_1_20241219_173724"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/XRT_1_AutoML_1_20241219_173724"
#
# # nospcv 20 folder
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_2_20241219_200021"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_2_20241219_200021_model_3"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_2_20241219_200021"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_2_20241219_200021"
#
# #### 20 nfolder  all data peat in 5 factor
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241219_233834"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_1_20241219_233834_model_3"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241219_233834"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/XRT_1_AutoML_1_20241219_233834"
#
#
# # reclearned wosis data to peat and 20 nfolder no spcv
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241220_162547"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_1_20241220_162547_model_3"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241220_162547"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/XRT_1_AutoML_1_20241220_162547"
#
# ### recleaned wosis with 4 factor peat and 20 nfolder no spcv instead of 5 above
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241220_194102"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_1_20241220_194102_model_3"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241220_194102"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20241220_194102"
#
# #####1221recleaned wosis with 5 factor peat and 20 nfolder no spcv ,keep pure of 1st class and let other extract from peat ca
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241221_154712"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_1_AutoML_1_20241221_154712"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241221_154712"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20241221_154712"
#
#
# # give some <5 to 1st class
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241221_224322"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_1_AutoML_1_20241221_224322"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241221_224322"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20241221_224322"
#
#
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241222_132638"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_1_AutoML_1_20241222_132638"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241222_132638"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20241222_132638"
#
#
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241222_203522"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_1_AutoML_1_20241222_203522"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241222_203522"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20241222_203522"
#
#
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241228_120929"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_1_AutoML_1_20241228_120929"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241228_120929"
# [1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20241228_120929"
##step 6
##step 6
##step 6 base map
#1 load covars
unixtools::set.tempdir('/mnt/Fastrun/temp4r')
library(terra)
library(dplyr)
library(tidyr)
library(data.table)

ca_covas <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global1214.tif')
gpp_hist <- rast(
  '/mnt/DataSpace/Projects/Canada_C/processed_output/GPP_historical_from_fluxcom_1981_2000_gCm2year.tif'
)
ca_covas[[4]] <- gpp_hist

# add xy sedimentary
coord00 <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_coordinates.tif')
sedimentary_thick <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/Canada_sedimentary_thickness.tif')
ca_covas <- c(ca_covas, coord00, sedimentary_thick)

#add peat 20240404
peat_ca <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/peatland_area1214.tif')

peat_ca <- ifel(is.na(peat_ca), 0, peat_ca)
ca_covas[['Peat']] <- peat_ca


reclass_table <- matrix(c(0, 5, 5, 5, 25, 25, 25, 50, 50, 50, 60, 60, 60, 100, 100),
                        ncol = 3,
                        byrow = TRUE)


peat_reclassified <- classify(peat_ca, reclass_table)

# add sedimentary depth
# f2 <- list.files(pattern = 'tif','/mnt/File0/DAAATAAA/soil_thickness_and_sediments/Global_Soil_Regolith_Sediment_1304/data',full.names = T)
# f_sed <- rast(f2[1])# https://daac.ornl.gov/SOILS/guides/Global_Soil_Regolith_Sediment.html
# f_sed2 <- crop(f_sed,ext(ca_covas))
# f_sed3 <- ifel(f_sed2 >50,0,f_sed2)
# writeRaster(f_sed3,'/mnt/File0/DAAATAAA/forfinalmap_1km/Canada_sedimentary_thickness.tif')
# ca_covas <- c(ca_covas,f_sed3)
#writeRaster(ca_covas,'/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global1214_fluxcomgpp.tif',overwrite=T)#20240619 add fluxcom 1981-2000 GPP
#ca_covas <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global_re_align_with_wrb.tif')

names(ca_covas) <- c(
  'NDVI1',
  'EVI2',
  'NPP3',
  'GPP4',
  'PET5',
  'AI6',
  'BIO1',
  'BIO10',
  'BIO11',
  'BIO12',
  'BIO13',
  'BIO14',
  'BIO15',
  'BIO16',
  'BIO17',
  'BIO18',
  'BIO19',
  'BIO2',
  'BIO3',
  'BIO4',
  'BIO5',
  'BIO6',
  'BIO7',
  'BIO8',
  'BIO9',
  'Solar26',
  'So_seasonal27',
  'Pop_den28',
  'HFP29',
  'Vrm30',
  'Tcurv31',
  'Pcurv32',
  'Ele33',
  'Slope34',
  'Asp_c35',
  'Asp_s36',
  'Est37',
  'Nor38',
  'Roug39',
  'TPI40',
  'TRI41',
  'Dx42',
  'Dxx43',
  'Dy44',
  'Dyy45',
  'Maj48',
  'lith.8',
  'lith.24',
  'WTD85',
  'smap',
  'SUMAP',
  'Ref_band1',
  'ref_band2',
  'ref_band3',
  'ref_band4',
  'ref_band5',
  'ref_band6',
  'ref_band7',
  'lon',
  'lat',
  'Sed_thick',
  'Peat'
) #20240730 newly added lon lat sed_thick
##'Histel','Histosol','WRB'
#
#writeRaster(ca_covas,'/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global1214_fluxcomgpp_0730_lonlat_sedi.tif')

writeRaster(
  ca_covas,
  '/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global1214_fluxcomgpp_0730_lonlat_sedi_peat.tif'
)


#20231214
preproce <- function(cov_raster, dep0) {
  cov_raster$Depth <- dep0
  data0 = as.data.frame(as.matrix(cov_raster))
  new_d <- data0 %>%
    mutate_all(list( ~ replace(., . == c(-9999), NA))) %>%
    mutate(
      #Depth = factor(Depth),
      NDVI1 = NDVI1 / 10000,
      EVI2 = EVI2 / 10000,
      NPP3 = NPP3 / 10000,
      GPP4 = GPP4 / 1000,
      ##
      PET5 = log(PET5),
      #
      AI6 = AI6 / 10000,
      BIO4 = BIO4 / 100,
      #
      BIO12 = ifelse(BIO12 == 0, -1, log(BIO12)),
      #2024 0111 thursday
      BIO13 = ifelse(BIO13 == 0, -1, log(BIO13)),
      #
      BIO14 = ifelse(BIO14 == 0, -1, log(BIO14)),
      #
      BIO15 = ifelse(BIO15 == 0, -1, log(BIO15)),
      #
      BIO16 = ifelse(BIO16 == 0, -1, log(BIO16)),
      #
      BIO17 = ifelse(BIO17 == 0, -1, log(BIO17)),
      #
      BIO18 = ifelse(BIO18 == 0, -1, log(BIO18)),
      #
      BIO19 = ifelse(BIO19 == 0, -1, log(BIO19)),
      #
      Solar26 = Solar26 / 10000,
      So_seasonal27 = So_seasonal27 / 10000,
      Vrm30 = Vrm30 * 10^4,
      #
      Tcurv31 = Tcurv31 * 10^4,
      #
      Pcurv32 = Pcurv32 * 10^4,
      #
      Ref_band1 = Ref_band1 / 10000,
      ref_band2 = ref_band2 / 10000,
      ref_band3 = ref_band3 / 10000,
      ref_band4 = ref_band4 / 10000,
      ref_band5 = ref_band5 / 10000,
      ref_band6 = ref_band6 / 10000,
      ref_band7 = ref_band7 / 10000,
      Maj48 = factor(Maj48, levels = 1:10),
      lith.8 = factor(lith.8, levels = 1:7),
      lith.24 = factor(lith.24, levels = 1:15),
      
      lon = lon / 10^6,
      lat = lat / 10^6,
      Peat = cut(
        Peat,
        breaks = c(0, 5, 30, 60, 100),
        labels = c(5, 30, 60, 100),
        include.lowest = TRUE,
        right = TRUE
      )
    )
  return(new_d)
}

ca_base_new <- preproce(ca_covas, dep0 = 0)
ind <- rowSums(is.na(ca_base_new)) == 0
length(which(ind))
niubi <- ca_base_new[ind, ] %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

fwrite(
  ca_base_new,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ca_base_new_use_for_base_only.csv'
)
fwrite(
  niubi,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ca_base_new_df20240730.csv'
) #lon and lat and sedimeantry
saveRDS(
  ind,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ind20240730.rds'
)


#############################note

ca_base_new <- fread(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ca_base_new_use_for_base_only.csv'
)
niubi <- fread(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ca_base_new_df20240730.csv'
) #lon and lat and sedimeantry
ind <- readRDS(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ind20240730.rds'
)

niubi$Maj48 <- factor(niubi$Maj48, levels = c(1:10))
niubi$lith.8 <- factor(niubi$lith.8, levels = c(1:7))
niubi$lith.24 <- factor(niubi$lith.24, levels = c(1:15))
niubi$Peat <- factor(niubi$Peat, levels = c(5, 30, 60, 100))
# # ca_base_new <- fread('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ca_base_new_df2024.csv')
# # ind <- readRDS('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ind2024.rds')
# file.copy(from = '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ca_base_new_df2024.csv',
#           to = '/mnt/Fastrun/Canada_C/rawdata/004_sensitive_analysis/ca_base_new_df2024.csv',overwrite = T)
# file.copy(from = '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ca_base_new_df20240110.csv',
#           to = '/mnt/Fastrun/Canada_C/rawdata/004_sensitive_analysis/ca_base_new_df20240110.csv',overwrite = T)
#
# file.copy(from = '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ind2024.rds',
#           to = '/mnt/Fastrun/Canada_C/rawdata/004_sensitive_analysis/ind2024.rds')


#2 load models
library(h2o)

h2o.init(nthreads = 64,
         max_mem_size = '300G',
         port = 15929)


#20231214
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models/bestfamiliy_ensemble/StackedEnsemble_BestOfFamily_6_AutoML_1_20231214_210359")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models/gbm/GBM_lr_annealing_selection_AutoML_1_20231214_210359_select_model")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models/xgb/XGBoost_lr_search_selection_AutoML_1_20231214_210359_select_grid_model_4")
# deepl4 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models/deeplearn/DeepLearning_grid_1_AutoML_1_20231214_210359_model_3")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models/rf/DRF_1_AutoML_1_20231214_210359")
# #
# 20240320 add lehigh
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bestfamiliy_ensemble/StackedEnsemble_BestOfFamily_3_AutoML_2_20240320_221513")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/gbm/GBM_grid_1_AutoML_2_20240320_221513_model_16")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/xgb/XGBoost_grid_1_AutoML_2_20240320_221513_model_18")
# deepl4 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/deeplearn/DeepLearning_grid_3_AutoML_2_20240320_221513_model_1")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/rf/XRT_1_AutoML_2_20240320_221513")
# 20240321 add sheryl
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bestfamiliy_ensemble/StackedEnsemble_BestOfFamily_4_AutoML_4_20240321_174522")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/gbm/GBM_grid_1_AutoML_4_20240321_174522_model_16")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/xgb/XGBoost_grid_1_AutoML_4_20240321_174522_model_52")
# deepl4 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/deeplearn/DeepLearning_grid_1_AutoML_4_20240321_174522_model_1")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/rf/XRT_1_AutoML_4_20240321_174522")
# 20240321 only sheryl
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bestfamiliy_ensemble/StackedEnsemble_AllModels_6_AutoML_1_20240321_214650")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/gbm/GBM_5_AutoML_1_20240321_214650")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/xgb/XGBoost_lr_search_selection_AutoML_1_20240321_214650_select_grid_model_3")
# deepl4 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/deeplearn/DeepLearning_grid_1_AutoML_1_20240321_214650_model_3")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/rf/XRT_1_AutoML_1_20240321_214650")
#
# 20240325 wosis cheryl sasha
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bestfamiliy_ensemble/StackedEnsemble_BestOfFamily_6_AutoML_1_20240325_161549")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/gbm/GBM_grid_1_AutoML_1_20240325_161549_model_3")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/xgb/XGBoost_grid_1_AutoML_1_20240325_161549_model_13")
# deepl4 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/deeplearn/DeepLearning_1_AutoML_1_20240325_161549")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/rf/DRF_1_AutoML_1_20240325_161549")


# #20240618 exclude NPP NDVI EVi:  wosis cheryl sasha
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_BestOfFamily_1_AutoML_1_20240618_110330")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/GBM_grid_1_AutoML_1_20240618_110330_model_34")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_1_20240618_110330_model_11")
# deepl4 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/deeplearn/DeepLearning_1_AutoML_1_20240618_110330")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20240618_110330")

#20240619 change GPP4/100 to GPP4/1000
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_BestOfFamily_5_AutoML_1_20240619_130612")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_grid_1_AutoML_1_20240619_130612_model_23")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_1_20240619_130612_model_1")
# deepl4 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/deeplearn/DeepLearning_grid_1_AutoML_1_20240619_130612_model_1")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20240619_130612")

#20240730
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_BestOfFamily_4_AutoML_1_20240730_172301")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_grid_1_AutoML_1_20240730_172301_model_1")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_1_20240730_172301_model_22")
# deepl4 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/deeplearn/DeepLearning_grid_1_AutoML_1_20240730_172301_model_1")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20240730_172301")

# #20241212 best = GBM > stack > xgb > rf > glm
# gbm1 <- h2o.loadModel( "/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/GBM_grid_1_AutoML_2_20241211_232609_model_4")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_2_20241211_232609_model_10")
# #[1] "/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_grid_1_AutoML_2_20241211_232609_model_4"
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/XRT_1_AutoML_2_20241211_232609")
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_BestOfFamily_6_AutoML_2_20241211_232609")
# deepl4 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/glm/GLM_1_AutoML_2_20241211_232609")
#
# # with deepl  20241212
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_1_20241212_75350_model_4")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_grid_1_AutoML_1_20241212_75350_model_4")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/XRT_1_AutoML_1_20241212_75350")
# deepl4 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/deeplearn/DeepLearning_1_AutoML_1_20241212_75350")
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_BestOfFamily_5_AutoML_1_20241212_75350")
# glm <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/glm/GLM_1_AutoML_1_20241212_75350")
#
# # no spcv
#
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_1_AutoML_1_20241212_120136")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241212_120136")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20241212_120136")
# deepl4 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/deeplearn/DeepLearning_grid_1_AutoML_1_20241212_120136_model_1")
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_2_AutoML_1_20241212_120136")
# glm <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/glm/GLM_1_AutoML_1_20241212_120136")

#
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241212_181356")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_1_20241212_181356_model_2")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_grid_1_AutoML_1_20241212_181356_model_2")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/XRT_1_AutoML_1_20241212_181356")


# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/GBM_grid_1_AutoML_1_20241213_125618_model_5")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_3_AutoML_1_20241213_125618")
# #h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_grid_1_AutoML_1_20241213_125618_model_5")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/XRT_1_AutoML_1_20241213_125618")
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_BestOfFamily_4_AutoML_1_20241213_125618")

### add peat now

# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_BestOfFamily_6_AutoML_1_20241216_155555")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_1_20241216_155555_model_2")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_5_AutoML_1_20241216_155555")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/XRT_1_AutoML_1_20241216_155555")

### real peat + depth < 150
# "StackedEnsemble" "GBM"             "XGBoost"         "DRF"
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_5_AutoML_2_20241219_120448")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_2_20241219_120448_model_1")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_2_AutoML_2_20241219_120448")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/XRT_1_AutoML_2_20241219_120448")
#
# stack_m <-h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_3_AutoML_1_20241219_173724")
#  xgb2 <-h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_1_20241219_173724_model_4")
#  gbm1 <-h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_2_AutoML_1_20241219_173724")
#  rf3 <-h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/XRT_1_AutoML_1_20241219_173724")
#
# # no spcv 20 cv folder
# stack_m <-h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_2_20241219_200021")
# xgb2 <-h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_2_20241219_200021_model_3")
# gbm1 <-h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_2_20241219_200021")
# rf3 <-h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_2_20241219_200021")


# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241219_233834")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_1_20241219_233834_model_3")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241219_233834")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/XRT_1_AutoML_1_20241219_233834")
#
#
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241220_162547")
#  xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_1_20241220_162547_model_3")
#  gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241220_162547")
#  rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/XRT_1_AutoML_1_20241220_162547")

# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241220_194102")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_grid_1_AutoML_1_20241220_194102_model_3")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241220_194102")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20241220_194102")
#

# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241221_154712")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_1_AutoML_1_20241221_154712")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241221_154712")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20241221_154712")
#
#
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241221_224322")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_1_AutoML_1_20241221_224322")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241221_224322")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20241221_224322")
#

# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241222_132638")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_1_AutoML_1_20241222_132638")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241222_132638")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20241222_132638")

# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241222_203522")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_1_AutoML_1_20241222_203522")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241222_203522")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20241222_203522")

stack_m <- h2o.loadModel(
  "/mnt/DataSpace/Projects/Canada_C/models_after_EC/bestfamiliy_ensemble/StackedEnsemble_AllModels_4_AutoML_1_20241228_120929"
)
xgb2 <- h2o.loadModel(
  "/mnt/DataSpace/Projects/Canada_C/models_after_EC/xgb/XGBoost_1_AutoML_1_20241228_120929"
)
gbm1 <- h2o.loadModel(
  "/mnt/DataSpace/Projects/Canada_C/models_after_EC/gbm/GBM_4_AutoML_1_20241228_120929"
)
rf3 <- h2o.loadModel(
  "/mnt/DataSpace/Projects/Canada_C/models_after_EC/rf/DRF_1_AutoML_1_20241228_120929"
)
models <- list(stack_m, gbm1, xgb2, rf3)

#
dep00 <- c(0, 5, 15, 30, 60, 100, 200, 300)
#dep00 <- c(2.5,10,22.5,45,80,150,250)
#ca_base_new <- preproce(ca_covas,dep0 = 0)
#ind = rowSums(is.na(ca_base_new)) == 0;length(which(ind))
#
#test <- as.h2o(x = ca_base_new[ind, ])

pred_fun0 <- function(model, t_data) {
  maps0 <- h2o.predict(model, t_data)
  ca_base_new[ind, 'pred'] <- as.data.frame(maps0)['predict']
  pred_r <- ca_covas$BIO1
  pred_r$pred = ca_base_new$pred
  return(pred_r$pred)
}
#

lt1000 <- function(x) {
  x[which(x > 120.3014)] <- 120.3014  # log(SOC *100) ^2 = log(580 *100) ^2 = 120.3014  som to soc log10(580) ; 1000/1.724
  return(x)
}
#
for (i in dep00[c(5, 4, 1, 2, 3, 6, 7, 8)]) {
  niubi$Depth = i
  test <- as.h2o(x = niubi)
  maps <- lapply(models, pred_fun0, t_data = test)
  maps0 <- do.call(c, maps)
  maps1 <- terra::clamp(maps0, upper = 120.3014, values = T)
  maps1 <- exp(sqrt(maps1)) / 100
  maps1$sd <- app(maps1[[2:4]], cores = 80, fun = 'sd')
  maps1$mean <- app(maps1[[2:4]], cores = 80, fun = 'mean')
  maps1$sd <- subst(maps1$sd, 0, minmax(maps1$sd)[2])
  maps1$cv <- maps1$sd / maps1$mean
  
  names(maps1) <- c('Stack', 'GBM', 'XGBoost', 'RF', 'sd', 'mean', 'cv')
  terra::writeRaster(
    maps1,
    paste0(
      '/mnt/DataSpace/Projects/Canada_C/SOC_maps/Base20240730_without_ndvievinpp_wosis_sasa_cheryl_',
      i,
      'cm_automl1213_nodeepglm_4categoricalpeat_nospcv_alldepth_recleaned_wosis_v3.tif'
    ),
    overwrite = TRUE
  )
  
}

socs <- list.files(
  pattern = 'tif',
  '/mnt/DataSpace/Projects/Canada_C/SOC_maps/',
  full.names = T
)
t1 <- grep('cm_automl1213_nodeepglm.tif', socs, value =  T)[c(1, 7, 3, 6, 8, 2)]
t1 <- grep('cm_automl1213_nodeepglm_peat.tif', socs, value = T)[c(1, 7, 3, 6, 8, 2)]
t2 <- rast(t1)
plot(t2[[seq(1, 42, 7)]], breaks = c(0, 5, 10, 20, 30, 45, 60, 80, 100, 200, 300, 400, 600))
plot(t2[[seq(1, 42, 7)]])
library(RColorBrewer)
#display.brewer.all()
col0 <- colorRampPalette(rev(brewer.pal(n = 8, name = 'Spectral')))
breaks0 = c(0, 5, 20, 45, 60, 80, 100, 200, 300, 400, 600)
par(mfrow = c(2, 3))
#cherylonly
plot(rast(socs[41])[[1]], breaks = breaks0, col = col0(length(breaks0)))
#cheryl+lehigh
plot(rast(socs[17])[[1]], breaks = breaks0, col = col0(length(breaks0)))
breaks0 = c(-300, -60, -30, -5, 0, 5, 30, 60, 600)
plot(rast(socs[41])[[1]] - rast(socs[33])[[1]],
     breaks = breaks0,
     col = col0(length(breaks0)))
#add lon lat sedimentary
for (i in 1:6) {
  plot(
    rast(socs[65])[[i]],
    breaks = breaks0,
    col = col0(length(breaks0)),
    main = names(rast(socs[65]))[i]
  )
}

###>>>>>>>>>>>>>>>>>>
##################
for (i in 1:4) {
  bm <- t2[[seq(i, 42, 7)]]
  gb <- (bm[[1:5]] + bm[[2:6]]) / 2
  ca_pure_area84 <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/ca_cmip6_area_WGS84.tif') #new created
  
  BD0 <- 0.08403 + 1.34364 * exp(-0.00874 * gb)
  socs_ml <- gb * BD0 * c(5, 10, 15, 30, 40)
  socs_ml_1m <- sum(socs_ml[[1:5]])
  socs_global_1m <- socs_ml_1m * ca_pure_area84
  ttc <- global(socs_global_1m, function(x)
    sum(x, na.rm = T)) / 10^14
  cat(ttc[[1]])
}






############################################################################
#####for bd bulk density ##########
#####for bd bulk density ##########
#####for bd bulk density ##########
#####for bd bulk density ##########
#####for bd bulk density ##########
#####for bd bulk density ##########
#####for bd bulk density ##########
#####for bd bulk density ##########
#####for bd bulk density ##########
#####for bd bulk density ##########
#####for bd bulk density ##########
#####for bd bulk density ##########

soc_ca <- use_layer3 %>% filter(country_name == 'Canada')

bd_ca <- read.csv(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/bd_ca_for_jitter_visual.csv'
)
bd_ca0 <- bd_ca[, -1]
bd_ca0[which(!complete.cases(bd_ca0)), ]

library(aqp)
depths(bd_ca0) <- profile_id ~ upper_depth + lower_depth
library(ithir)
use_bd0 <- ithir::ea_spline(
  bd_ca0,
  var.name = 'bd_use' ,
  lam = 0.1,
  d = c(0, 5, 10, 15, seq(20, 100, 10), seq(120, 300, 30))
) #20241213 try it again to predict bd again before # 20230615 down to 3m
#
saveRDS(use_bd0, file = '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/bd_use_splined_updated_intervalto10cm20241213.rds')
#
bd_ca1 <- bd_ca[, -1]
t1 <- dplyr::distinct(bd_ca1[, c(1:6)])

use_bd0[[1]]$id <- as.numeric(use_bd0[[1]]$id)
use_data1 <- dplyr::left_join(use_bd0[[1]], t1, by = c('id' = 'profile_id')) #
#
write_csv(
  use_data1,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/bd_use_ca_with_xy_updated_interval20241213.csv'
)
#
library(terra)
#aggregate
glob_raster <- rast(gsub(
  "\\\\",
  "/",
  paste0(
    '\\mnt\\File0',
    "\\DAAATAAA\\Global aridity and PET\\7504448 (1)\\Global-AI_ET0_annual_v3\\Global-AI_ET0_v3_annual\\ai_v3_yr.tif"
  )
))
head(use_data1)
#
use_data1[use_data1 == c(-9999)] <- NA

use_data2 <- st_as_sf(use_data1,
                      coords = c('longitude', 'latitude'),
                      crs = st_crs(4326))

bg_extract <- terra::extract(x = glob_raster,
                             y = use_data2,
                             cells = T,
                             xy = T)
length(which(table(bg_extract$cell) > 1))  #>2 #213
use_data3 <- cbind(use_data1, bg_extract[, c(1, 3:5)])
library(dplyr)

#20230504 for 3m depth
use_data4 <- use_data3 %>% group_by(cell) %>% summarise_at(vars('0-5 cm':'270-300 cm', 'latitude', 'longitude', 'x', 'y'),
                                                           mean,
                                                           na.rm = T)

use_data5 <- use_data3 %>% group_by(cell) %>% summarise_at(vars('country_id', 'country_name'), function(x)
  x[1])
use_data6 <- left_join(use_data4, use_data5, by = 'cell') # we now stilll use longi and lati for distance
use_data6$ID_no <- 1:nrow(use_data6)
write.csv(
  use_data6,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/bd_ca_aggregated1km_20230615_20241213updated.csv'
)

##joint peat bd

ca_covas <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global1214.tif')

names(ca_covas) <- c(
  'NDVI1',
  'EVI2',
  'NPP3',
  'GPP4',
  'PET5',
  'AI6',
  'BIO1',
  'BIO10',
  'BIO11',
  'BIO12',
  'BIO13',
  'BIO14',
  'BIO15',
  'BIO16',
  'BIO17',
  'BIO18',
  'BIO19',
  'BIO2',
  'BIO3',
  'BIO4',
  'BIO5',
  'BIO6',
  'BIO7',
  'BIO8',
  'BIO9',
  'Solar26',
  'So_seasonal27',
  'Pop_den28',
  'HFP29',
  'Vrm30',
  'Tcurv31',
  'Pcurv32',
  'Ele33',
  'Slope34',
  'Asp_c35',
  'Asp_s36',
  'Est37',
  'Nor38',
  'Roug39',
  'TPI40',
  'TRI41',
  'Dx42',
  'Dxx43',
  'Dy44',
  'Dyy45',
  'Maj48',
  'lith.8',
  'lith.24',
  'WTD85',
  'smap',
  'SUMAP',
  'Ref_band1',
  'ref_band2',
  'ref_band3',
  'ref_band4',
  'ref_band5',
  'ref_band6',
  'ref_band7'
)

bd_ca <- read.csv(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/bd_ca_aggregated1km_20230615_20241213updated.csv'
)

#
######################20240325
peat_use <- readRDS(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/Sasha_peat_20240325_fortrian.rds'
)
peat_use1 <- peat_use[, c(1:5, 17, 34, 35)]

peat_use1$SAMP_THICK[which(is.na(peat_use1$SAMP_THICK))] <- 3
# create lower depth
library(dplyr)
peat_use2 <- peat_use1 %>% mutate(lower_depth = UPPER_SAMP_DEPTH + SAMP_THICK) %>% filter(complete.cases(.))
names(peat_use2) <- c('ID',
                      'Sub_ID',
                      'upper',
                      'thick',
                      'BD',
                      'SOC',
                      'lat',
                      'lon',
                      'lower')
peat_use2$SOC <- peat_use2$SOC * 10  # % to g/kg
library(aqp)
depths(peat_use2) <- ID ~ upper + lower

library(ithir)
thd_peat0 <- ithir::ea_spline(
  peat_use2,
  var.name = 'BD' ,
  lam = 0.1,
  d = c(0, 5, 10, 15, seq(20, 100, 10), seq(120, 300, 30))
) #
peat_use3 <- thd_peat0[[1]]

#canada all peat with geo info
all_peat4 <- left_join(peat_use3,
                       peat_use1[, c(1, 7, 8)],
                       by = c('id' = 'CORE_ID'),
                       multiple = 'first')
#
all_peat4$geo <- paste0(all_peat4$LONGITUDE, '_', all_peat4$LATITUDE)
length(unique(all_peat4$geo))
t1 <- table(all_peat4$id)
t2 <- table(all_peat4$geo)
#
dup0 <- all_peat4[which(all_peat4$geo %in% names(which(t2 > 1))), ] %>% select(id, LATITUDE, LONGITUDE)
dup1 <- st_as_sf(dup0,
                 coords = c('LONGITUDE', 'LATITUDE'),
                 crs = st_crs(4326))
all_peat5 <- all_peat4 %>% group_by(geo) %>% summarise(id = min(id), across(where(is.numeric), function(x)
  mean(x, na.rm = T))) %>% arrange(id) %>% select(-geo)
thd_peat3 <- all_peat5
bd_peat7 <- thd_peat3[, c(1:20, 22, 23)]
names(bd_peat7) <- names(bd_ca)[c(1, 3:23)]
bd_use1 <- rbind(bd_ca[, c(1, 3:23)], bd_peat7) #
bd_peat8 <- st_as_sf(bd_peat7,
                     coords = c('longitude', 'latitude'),
                     crs = st_crs(4326))
write_sf(
  bd_peat8,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/bd_and_soc_fromsasha_updated10interval20241213.shp'
)
unlink(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/bd_peatonly_fromsasha_updated10interval20241213.shp'
)

nrow(bd_peat8)

####################add SOC to train the model , no soc using empiricall relationship
### add SOC version
#wosis
soc_spline <- read_rds(
  '/mnt/File0/DAAATAAA/Data_collection/MY_global_Temp_resilence_C/00Rawdata/third_splined_data_to3m_use_20241212_10cmintervals.rds'
)
bd_spline <- read_csv(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/bd_use_ca_with_xy_updated_interval20241213.csv'
)
names(soc_spline)
soc_ca_spline <- soc_spline %>% filter(country_name == 'Canada')
soc_ca_spline$id <- as.numeric(soc_ca_spline$id)
bd_soc <- left_join(x = bd_spline,
                    y = soc_ca_spline,
                    by = c('id', 'soil depth'))
use_data1 <- bd_soc # use code above
#
library(terra)
#aggregate
glob_raster <- rast(gsub(
  "\\\\",
  "/",
  paste0(
    '\\mnt\\File0',
    "\\DAAATAAA\\Global aridity and PET\\7504448 (1)\\Global-AI_ET0_annual_v3\\Global-AI_ET0_v3_annual\\ai_v3_yr.tif"
  )
))
head(use_data1)
#
use_data1[use_data1 == c(-9999)] <- NA

use_data2 <- st_as_sf(use_data1,
                      coords = c('longitude.x', 'latitude.x'),
                      crs = st_crs(4326))

bg_extract <- terra::extract(x = glob_raster,
                             y = use_data2,
                             cells = T,
                             xy = T)
length(which(table(bg_extract$cell) > 1))  #>2 #213
use_data3 <- cbind(use_data1, bg_extract[, c(1, 3:5)])
library(dplyr)

#20230504
use_data4 <- use_data3 %>% group_by(cell) %>% summarise_at(
  vars(
    '0-5 cm.x':'270-300 cm.x',
    '0-5 cm.y':'270-300 cm.y',
    'latitude.x',
    'longitude.y',
    'x',
    'y'
  ),
  mean,
  na.rm = T
)

use_data5 <- use_data3 %>% group_by(cell) %>% summarise_at(vars('country_id.x', 'country_name.y'), function(x)
  x[1])
use_data6 <- left_join(use_data4, use_data5, by = 'cell')
use_data6$ID_no <- 1:nrow(use_data6)
write.csv(
  use_data6,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/bd_ca_aggregated1km_20230615_20241213updated_bd_add_soc_spline.csv'
)
use_data6 <- read.csv(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/bd_ca_aggregated1km_20230615_20241213updated_bd_add_soc_spline.csv'
)[, -1]

###> can peat
soct <- readRDS(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/Sasha_peat_20240325_canada_sf_peatland_SOC_20241212update_10cminterval.rds'
)

bdt <- read_sf(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/bd_and_soc_fromsasha_updated10interval20241213.shp'
)

can_peat_soc_bd <- left_join(
  x = data.frame(st_drop_geometry(bdt)),
  y = data.frame(st_drop_geometry(soct)),
  by = c('X' = 'id')
)
can_peat_soc_bd1 <- cbind(can_peat_soc_bd, st_coordinates(bdt))
wosis_t1 <- use_data6[, c(2:39, 42:43)]
peat_t2 <- can_peat_soc_bd1[, c(2:39, 41, 42)]
names(peat_t2) <- names(wosis_t1)
combine_soc_bd <- rbind(wosis_t1, peat_t2)
saveRDS(
  combine_soc_bd,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/soc_and_bd_for_bdtrain_20241213.rds'
)
names(combine_soc_bd)[39:40] <- c('lon', 'lat')
names(peat_t2)[39:40] <- c('lon', 'lat')
names(wosis_t1)[39:40] <- c('lon', 'lat')

t10 <-
  combine_soc_bd %>% pivot_longer(
    cols = matches("cm\\.x|cm\\.y"),
    names_to = c("depth", ".value"),
    names_pattern = "(.*) cm\\.(x|y)"
  ) %>%
  rename(bd = x, SOC = y)

t10$Peat <- c(rep(0, 39026), rep(1, 16340))
use_ind_bd <- which(is.na(t10$bd))
use_ind_bd1 <- which(t10$bd <= 0)

ttsoc <- t10[which(t10$bd > c(1.34364 + 0.08403)) , ] # threshold of BD get negtive SOC
summary(ttsoc$SOC)
ttbd <- t10[which(t10$bd < c(0.08403)) , ] # threshold of BD get negtive SOC
summary(ttbd$SOC)



t11 <- t10[-c(use_ind_bd, use_ind_bd1), ] %>%
  mutate(SOC = ifelse(is.na(SOC) |
                        SOC == -9999, log((bd - 0.08403) / 1.34364) / (-0.00874), SOC)) %>%
  mutate(SOC = ifelse(is.na(SOC), 480.91, SOC)) %>%
  mutate(SOC = ifelse(SOC <= 0, 1.29, SOC)) %>%
  mutate(SOC = ifelse(SOC >= 580, 580, SOC)) %>%
  mutate(bd = ifelse(bd > 2, 2, bd))

summary(t10$bd)
summary(t10$SOC)
summary(t11$bd)
summary(t11$SOC)
t11[which(is.na(t11$SOC)), ]
t11[which(t11$SOC < 0) , ]


length(which(!is.na(t10$bd))) / 19
saveRDS(
  t11,
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/soc_and_bd_for_bdtrain_20241214_longformat.rds'
)

#########################model build #########################
#########################model build #########################
#########################model build #########################
#########################model build #########################
#########################model build #########################
#########################model build #########################
#########################model build #########################

# #new add SOC
# library(sf)
# library(terra)
# library(dplyr)
# library(tidyr)
# bd_ca <- t11
# bd_ca_sf0 <- st_as_sf(bd_ca,coords= c('lon','lat'),crs=st_crs(4326))
#
# ca_covas <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global1214.tif')
# gpp_hist <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/GPP_historical_from_fluxcom_1981_2000_gCm2year.tif' )
# ca_covas[[4]] <- gpp_hist
# coord00 <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_coordinates.tif')
# sedimentary_thick <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/Canada_sedimentary_thickness.tif')
# ca_covas <- c(ca_covas,coord00,sedimentary_thick)
#
# #ca_covas <- rast('/mnt/DataSpace/nidayegedan.tif')
# names(ca_covas) <- c('NDVI1','EVI2','NPP3','GPP4','PET5','AI6',
#                      'BIO1','BIO10','BIO11','BIO12','BIO13','BIO14','BIO15','BIO16','BIO17','BIO18','BIO19','BIO2','BIO3','BIO4','BIO5','BIO6','BIO7','BIO8','BIO9',
#                      'Solar26','So_seasonal27','Pop_den28','HFP29','Vrm30','Tcurv31','Pcurv32','Ele33','Slope34','Asp_c35','Asp_s36','Est37','Nor38','Roug39','TPI40','TRI41','Dx42','Dxx43','Dy44','Dyy45',
#                      'Maj48','lith.8','lith.24','WTD85','smap','SUMAP','Ref_band1','ref_band2','ref_band3','ref_band4','ref_band5','ref_band6','ref_band7',
#                      'lon','lat','Sed_thick')
# ####add peat
# #add peat 20240404
# peat_ca <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/peatland_area1214.tif') #change name from peatland_area to peatland1214
# #first try using peatfraction 20240404
# # peat_ca <- ifel(is.na(peat_ca),0,peat_ca)
# # ca_covas[['Peat']] <- peat_ca
# #second try using peat 0 1 20240405
# peat_ca <- ifel(is.na(peat_ca),0,peat_ca)
# peat_ca <- ifel(peat_ca >= 25,1,0)
# ca_covas[['Peat']] <- peat_ca
#
# plot(peat_ca)
#
# #
# # ttt <- terra::extract(ca_covas,vect(bd_ca_sf0))
# # nrow(ttt)
# # bd_ca_sf1 <- cbind(bd_ca_sf0,ttt)
# # bd_ca_sf1$Pop_den28 <- as.integer(bd_ca_sf1$Pop_den28)
# #
#
#
# #
# ttt <- terra::extract(ca_covas,vect(bd_ca_sf0))
# names(ttt)
# bd_ca_sf1 <- cbind(bd_ca_sf0,ttt)
# bd_ca_sf1$Pop_den28 <- as.integer(bd_ca_sf1$Pop_den28)
# bd_ca_sf1$Peat <- factor(bd_ca_sf1$Peat,levels = c(0,1))
# plot(SOC~bd,data= bd_ca_sf1,pch = 16,cex = 0.1)

t11 <- readRDS(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/000rawdata/soc_and_bd_for_bdtrain_20241214_longformat.rds'
)
#delete last row of each site
t111 <- t11 %>% mutate(site = paste0(lon, '_', lat)) %>% group_by(site) %>% filter(row_number() != n()) %>% filter(n() > 3) %>% ungroup()
#new add SOC
library(sf)
library(terra)
library(dplyr)
library(tidyr)
bd_ca <- t111
bd_ca_sf0 <- st_as_sf(bd_ca, coords = c('lon', 'lat'), crs = st_crs(4326))
ca_covas <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global1214.tif')
gpp_hist <- rast(
  '/mnt/DataSpace/Projects/Canada_C/processed_output/GPP_historical_from_fluxcom_1981_2000_gCm2year.tif'
)
ca_covas[[4]] <- gpp_hist
coord00 <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_coordinates.tif')
sedimentary_thick <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/Canada_sedimentary_thickness.tif')
ca_covas <- c(ca_covas, coord00, sedimentary_thick)

names(ca_covas) <- c(
  'NDVI1',
  'EVI2',
  'NPP3',
  'GPP4',
  'PET5',
  'AI6',
  'BIO1',
  'BIO10',
  'BIO11',
  'BIO12',
  'BIO13',
  'BIO14',
  'BIO15',
  'BIO16',
  'BIO17',
  'BIO18',
  'BIO19',
  'BIO2',
  'BIO3',
  'BIO4',
  'BIO5',
  'BIO6',
  'BIO7',
  'BIO8',
  'BIO9',
  'Solar26',
  'So_seasonal27',
  'Pop_den28',
  'HFP29',
  'Vrm30',
  'Tcurv31',
  'Pcurv32',
  'Ele33',
  'Slope34',
  'Asp_c35',
  'Asp_s36',
  'Est37',
  'Nor38',
  'Roug39',
  'TPI40',
  'TRI41',
  'Dx42',
  'Dxx43',
  'Dy44',
  'Dyy45',
  'Maj48',
  'lith.8',
  'lith.24',
  'WTD85',
  'smap',
  'SUMAP',
  'Ref_band1',
  'ref_band2',
  'ref_band3',
  'ref_band4',
  'ref_band5',
  'ref_band6',
  'ref_band7',
  'lon',
  'lat',
  'Sed_thick'
)

#add peat 20240404
peat_ca <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/peatland_area1214.tif') #
#first try using peatfraction 20240404
# peat_ca <- ifel(is.na(peat_ca),0,peat_ca)
# ca_covas[['Peat']] <- peat_ca
#second try using peat 0 1 20240405
peat_ca <- ifel(is.na(peat_ca), 0, peat_ca)
reclass_table <- matrix(c(0, 5, 5, 5, 30, 30, 30, 60, 60, 60, 100, 100),
                        ncol = 3,
                        byrow = TRUE)


peat_reclassified <- classify(peat_ca,
                              reclass_table,
                              include.lowest = T,
                              right = T)
peat_reclassified <- terra::focal(peat_reclassified,
                                  w = 7,
                                  fun = 'modal',
                                  na.policy = 'only')

ca_covas[['Peat']] <- peat_reclassified  # peat_ca

plot(peat_ca)

#



#
ttt <- terra::extract(ca_covas, vect(bd_ca_sf0))
names(ttt)
bd_ca_sf1 <- cbind(bd_ca_sf0[, 1:3], ttt)
bd_ca_sf1$Pop_den28 <- as.integer(bd_ca_sf1$Pop_den28)
bd_ca_sf1$Peat <- factor(bd_ca_sf1$Peat, levels = c(5, 30, 60, 100))
str(bd_ca_sf1)
plot(SOC ~ bd,
     data = bd_ca_sf1,
     pch = 16,
     cex = 0.1)

reclass_table <- matrix(c(0, 5, 5, 5, 30, 30, 30, 45, 45, 45, 100, 100),
                        ncol = 3,
                        byrow = TRUE)


peat_reclassified22 <- classify(peat_ca,
                                reclass_table,
                                include.lowest = T,
                                right = T)

peat000 <- terra::extract(peat_reclassified22, vect(bd_ca_sf0))
bd_ca_sf_test <- cbind(bd_ca_sf0[, 1:5], peat000[, 2])
bd_ca_sf_test$id <- 1:nrow(bd_ca_sf_test)
names(bd_ca_sf_test)[6] <- 'Peat1'
site100_0 <- bd_ca_sf_test %>% filter(Peat1 == 100 &
                                        depth == '30-40' &
                                        SOC < 300) %>% summarise(site0 = unique(site))
site100_1 <- bd_ca_sf_test %>% filter(Peat1 == 100 &
                                        depth == '20-30' &
                                        SOC < 300) %>% summarise(site0 = unique(site))
site100_3 <- bd_ca_sf_test %>% filter(Peat1 == 100 &
                                        depth == '5-10' &
                                        SOC < 300) %>% summarise(site0 = unique(site))
site100_2 <- bd_ca_sf_test %>% filter(Peat1 == 100 &
                                        depth == '15-20' &
                                        SOC < 300) %>% summarise(site0 = unique(site))
site100_4 <- bd_ca_sf_test %>% filter(Peat1 == 100 &
                                        depth == '10-15' &
                                        SOC < 300) %>% summarise(site0 = unique(site))

site100_all <- unique(
  c(
    site100_0$site0,
    site100_1$site0,
    site100_2$site0,
    site100_3$site0,
    site100_4$site0
  )
)
site100_ind <- which(bd_ca_sf0$site %in% site100_all)

#final bd use
bd_ca_sf1 <- bd_ca_sf1[-site100_ind, ]



item0 <- paste0("X", c(5, 10, 15, 20, 30, seq(40, 100, 10), seq(120, 300, 30)), "cm")

# covs0 <- c('NDVI1','EVI2','NPP3','GPP4','PET5','AI6','BIO1','BIO10','BIO11','BIO12','BIO13','BIO14','BIO15','BIO16','BIO17','BIO18','BIO19','BIO2',
#            'BIO3','BIO4','BIO5','BIO6','BIO7','BIO8','BIO9','Solar26','So_seasonal27','Pop_den28','HFP29','Vrm30','Tcurv31','Pcurv32','Ele33',
#            'Slope34','Asp_c35','Asp_s36','Est37','Nor38','Roug39','TPI40','TRI41','Dx42','Dxx43','Dy44','Dyy45','Maj48','lith.8','lith.24','WTD85',
#            'smap','SUMAP','Ref_band1','ref_band2','ref_band3','ref_band4','ref_band5',
#            'ref_band6','ref_band7')
# added_cov <- 'fold_col'
#before
#added_cov <- c('Histosol','Histel','WRB','fold_col')
#ca_data_sf0 from 01_data_and_model_train.r   # only load dplyr tidyr and sf, it works, but load other packages, it will not work.

# after 20231214
names(bd_ca_sf1)[1] <- 'Depth'
d10 <-
  bd_ca_sf1 %>% select(c(1:3, 5:66)) %>% mutate(
    Depth = case_when(
      Depth == '0-5' ~ 2.5,
      Depth == '5-10' ~ 7.5,
      Depth == '10-15' ~ 12.5,
      Depth == '15-20' ~ 17.5,
      Depth == '20-30' ~ 25,
      Depth == '30-40' ~ 35,
      Depth == '40-50' ~ 45,
      Depth == '50-60' ~ 55,
      Depth == '60-70' ~ 65,
      Depth == '70-80' ~ 75,
      Depth == '80-90' ~ 85,
      Depth == '90-100' ~ 95,
      Depth == '100-120' ~ 110,
      Depth == '120-150' ~ 135,
      Depth == '150-180' ~ 165,
      Depth == '180-210' ~ 195,
      Depth == '210-240' ~ 225,
      Depth == '240-270' ~ 255,
      Depth == '270-300' ~ 285,
    )
  ) %>%
  mutate(
    bd = log(bd * 100)^2 ,
    NDVI1 = NDVI1 / 10000,
    EVI2 = EVI2 / 10000,
    NPP3 = NPP3 / 10000,
    GPP4 = GPP4 / 1000,
    #
    PET5 = log(PET5),
    #
    AI6 = AI6 / 10000,
    BIO4 = BIO4 / 100,
    #++++++
    # BIO12 = log10(BIO12),#
    # BIO18 = log10(BIO18), #
    Solar26 = Solar26 / 10000,
    So_seasonal27 = So_seasonal27 / 10000,
    Vrm30 = Vrm30 * 10^4,
    #
    Tcurv31 = Tcurv31 * 10^4,
    #
    Pcurv32 = Pcurv32 * 10^4,
    #
    Ref_band1 = Ref_band1 / 10000,
    ref_band2 = ref_band2 / 10000,
    ref_band3 = ref_band3 / 10000,
    ref_band4 = ref_band4 / 10000,
    ref_band5 = ref_band5 / 10000,
    ref_band6 = ref_band6 / 10000,
    ref_band7 = ref_band7 / 10000,
    Maj48 = factor(Maj48, levels = 1:10),
    lith.8 = factor(lith.8, levels = 1:7),
    lith.24 = factor(lith.24, levels = 1:15),
    Peat = factor(Peat, levels = c(5, 30, 60, 100)),
    lon = lon / 10^6,
    lat = lat / 10^6
  ) %>% na.omit()
#
#
d10[, paste0('BIO', 12:19)] <- log(st_drop_geometry(d10[, paste0('BIO', 12:19)]))
d11 <- d10
d11$Pop_den28 <- as.integer(d11$Pop_den28)
str(d11)


library(h2o)
h2o.init(nthreads = 32, max_mem_size = '100G')
d1_hdf <- as.h2o(x = sf::st_drop_geometry(d11))
set.seed(23)
system.time(
  h_auto <- h2o.automl(
    x = setdiff(names(d1_hdf), c('NDVI1', 'EVI2', 'NPP3', 'bd')) ,
    y = 'bd',
    training_frame = d1_hdf,
    exclude_algos = c("GLM", "DeepLearning"),
    max_models = NULL,
    nfolds = 10,
    stopping_metric = 'RMSE',
    # or AUTO
    keep_cross_validation_predictions = TRUE,
    #
    seed = 927
  )
)   #done 3571.560



#
lb <- h2o.get_leaderboard(object = h_auto, extra_columns = "ALL")
print(lb, n = nrow(lb))  #
nidaye <- as.data.frame(lb)
saveRDS(
  nidaye,
  '/mnt/DataSpace/Projects/Canada_C/models/bd_h2o_leaderboad1214_log_precip_add_sasha_20241214_addSOC.rds'
)
saveRDS(
  nidaye,
  '/mnt/DataSpace/Projects/Canada_C/models/bd_h2o_leaderboad1214_log_precip_add_sasha_20241214_addSOC_peatland.rds'
)

saveRDS(
  nidaye,
  '/mnt/DataSpace/Projects/Canada_C/models/bd_h2o_leaderboad1214_log_precip_add_sasha_20241214_addSOC_peatland_realpeat.rds'
)#25% peat version
saveRDS(
  nidaye,
  '/mnt/DataSpace/Projects/Canada_C/models/bd_h2o_leaderboad1214_log_precip_add_sasha_20241214_addSOC_peatland_realpeat_v2_4category.rds'
) # 4category peat version

nidaye <- readRDS(
  '/mnt/DataSpace/Projects/Canada_C/models/bd_h2o_leaderboad1214_log_precip_add_sasha_20241214_addSOC_peatland_realpeat.rds'
)
unique(nidaye$algo)

m0 <- h_auto@leader



# Get the best XGBoost model using default sort metric
xgb3 <- h2o.get_best_model(h_auto, algorithm = "xgboost")

# Get the best GBM model, ranked by logloss
gbm1 <- h2o.get_best_model(h_auto, algorithm = "GBM", criterion = "RMSE")

rf4 <- h2o.get_best_model(h_auto, algorithm = "DRF", criterion = "RMSE")

# Get the best deepleardning model, ranked by logloss

# Get the best stackensemble model, ranked by ,also by what we want , best of Family
stacken_best_family <- h2o.get_best_model(h_auto, algorithm = "StackedEnsemble", criterion = "RMSE")
#


#
h2o.saveModel(
  object = m0,
  '/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm',
  export_cross_validation_predictions = TRUE
)
h2o.saveModel(
  object = gbm1,
  '/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm',
  export_cross_validation_predictions = TRUE
)
#h2o.saveModel(object= deeplearn2,'/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm',export_cross_validation_predictions = TRUE)
h2o.saveModel(
  object = xgb3,
  '/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm',
  export_cross_validation_predictions = TRUE
)
h2o.saveModel(
  object = rf4,
  '/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm',
  export_cross_validation_predictions = TRUE
)
h2o.saveModel(
  object = stacken_best_family,
  '/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm',
  export_cross_validation_predictions = TRUE
)

#from extract peat from peat map

[1] "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/StackedEnsemble_AllModels_4_AutoML_1_20241215_183014"
> h2o.saveModel(
  object = gbm1,
  '/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm',
  export_cross_validation_predictions = TRUE
)
[1] "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/GBM_grid_1_AutoML_1_20241215_183014_model_4"
> #h2o.saveModel(object= deeplearn2,'/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm',export_cross_validation_predictions = TRUE)
  > h2o.saveModel(
    object = xgb3,
    '/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm',
    export_cross_validation_predictions = TRUE
  )
[1] "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/XGBoost_grid_1_AutoML_1_20241215_183014_model_16"
> h2o.saveModel(
  object = rf4,
  '/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm',
  export_cross_validation_predictions = TRUE
)
[1] "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/DRF_1_AutoML_1_20241215_183014"

###add soc peat
[1] "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/StackedEnsemble_AllModels_4_AutoML_1_20241215_205434"
[1] "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/GBM_grid_1_AutoML_1_20241215_205434_model_13"
[1] "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/XGBoost_grid_1_AutoML_1_20241215_205434_model_20"
[1] "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/DRF_1_AutoML_1_20241215_205434"
###add soc peat 4 category
[1] "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/StackedEnsemble_AllModels_4_AutoML_5_20241231_120738"
[1] "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/GBM_4_AutoML_5_20241231_120738"
[1] "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/XGBoost_grid_1_AutoML_5_20241231_120738_model_20"
[1] "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/DRF_1_AutoML_5_20241231_120738"

################# mapping

unixtools::set.tempdir('/mnt/Fastrun/temp4r')
library(terra)
library(dplyr)
library(tidyr)
library(data.table)

ca_covas <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global1214.tif')
gpp_hist <- rast(
  '/mnt/DataSpace/Projects/Canada_C/processed_output/GPP_historical_from_fluxcom_1981_2000_gCm2year.tif'
)
ca_covas[[4]] <- gpp_hist

#
coord00 <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_coordinates.tif')
sedimentary_thick <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/Canada_sedimentary_thickness.tif')
#
peat_ca <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/peatland_area1214.tif')
#first try using peatfraction 20240404
# peat_ca <- ifel(is.na(peat_ca),0,peat_ca)
# ca_covas[['Peat']] <- peat_ca
#second try using peat 0 1 20240405

peat_ca <- ifel(is.na(peat_ca), 0, peat_ca)
reclass_table <- matrix(c(0, 5, 5, 5, 30, 30, 30, 60, 60, 60, 100, 100),
                        ncol = 3,
                        byrow = TRUE)


peat_reclassified <- classify(peat_ca,
                              reclass_table,
                              include.lowest = T,
                              right = T)

ca_covas[['Peat']] <- peat_reclassified


ca_covas <- c(ca_covas, coord00, sedimentary_thick)


names(ca_covas) <- c(
  'NDVI1',
  'EVI2',
  'NPP3',
  'GPP4',
  'PET5',
  'AI6',
  'BIO1',
  'BIO10',
  'BIO11',
  'BIO12',
  'BIO13',
  'BIO14',
  'BIO15',
  'BIO16',
  'BIO17',
  'BIO18',
  'BIO19',
  'BIO2',
  'BIO3',
  'BIO4',
  'BIO5',
  'BIO6',
  'BIO7',
  'BIO8',
  'BIO9',
  'Solar26',
  'So_seasonal27',
  'Pop_den28',
  'HFP29',
  'Vrm30',
  'Tcurv31',
  'Pcurv32',
  'Ele33',
  'Slope34',
  'Asp_c35',
  'Asp_s36',
  'Est37',
  'Nor38',
  'Roug39',
  'TPI40',
  'TRI41',
  'Dx42',
  'Dxx43',
  'Dy44',
  'Dyy45',
  'Maj48',
  'lith.8',
  'lith.24',
  'WTD85',
  'smap',
  'SUMAP',
  'Ref_band1',
  'ref_band2',
  'ref_band3',
  'ref_band4',
  'ref_band5',
  'ref_band6',
  'ref_band7',
  'Peat',
  'lon',
  'lat',
  'Sed_thick'
) #20240730 newly added lon lat sed_thick

##Add SOC
socs <- list.files(
  pattern = 'tif',
  '/mnt/DataSpace/Projects/Canada_C/SOC_maps/',
  full.names = T
)

t1 <- grep('v3.tif', socs, value =  T)[c(1, 7, 3, 6, 8, 2, 4, 5)]

t2 <- rast(t1)[[seq(1, 56, 7)]]
ind <- readRDS(
  '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ind20240730.rds'
)

SOC <- as.data.frame(as.matrix(t2))
SOC_pred <- SOC[ind, ]
names(SOC_pred) <- c(0, 5, 15, 30, 60, 100, 200, 300)
#20231214
preproce <- function(cov_raster, dep0) {
  cov_raster$Depth <- dep0
  data0 = as.data.frame(as.matrix(cov_raster))
  new_d <- data0 %>%
    mutate_all(list( ~ replace(., . == c(-9999), NA))) %>%
    mutate(
      #Depth = factor(Depth),
      NDVI1 = NDVI1 / 10000,
      EVI2 = EVI2 / 10000,
      NPP3 = NPP3 / 10000,
      GPP4 = GPP4 / 1000,
      ##
      PET5 = log(PET5),
      #
      AI6 = AI6 / 10000,
      BIO4 = BIO4 / 100,
      #
      BIO12 = ifelse(BIO12 == 0, -1, log(BIO12)),
      #2024 0111 thursday
      BIO13 = ifelse(BIO13 == 0, -1, log(BIO13)),
      #
      BIO14 = ifelse(BIO14 == 0, -1, log(BIO14)),
      #
      BIO15 = ifelse(BIO15 == 0, -1, log(BIO15)),
      #
      BIO16 = ifelse(BIO16 == 0, -1, log(BIO16)),
      #
      BIO17 = ifelse(BIO17 == 0, -1, log(BIO17)),
      #
      BIO18 = ifelse(BIO18 == 0, -1, log(BIO18)),
      #
      BIO19 = ifelse(BIO19 == 0, -1, log(BIO19)),
      #
      Solar26 = Solar26 / 10000,
      So_seasonal27 = So_seasonal27 / 10000,
      Vrm30 = Vrm30 * 10^4,
      #
      Tcurv31 = Tcurv31 * 10^4,
      #
      Pcurv32 = Pcurv32 * 10^4,
      #
      Ref_band1 = Ref_band1 / 10000,
      ref_band2 = ref_band2 / 10000,
      ref_band3 = ref_band3 / 10000,
      ref_band4 = ref_band4 / 10000,
      ref_band5 = ref_band5 / 10000,
      ref_band6 = ref_band6 / 10000,
      ref_band7 = ref_band7 / 10000,
      Maj48 = factor(Maj48, levels = 1:10),
      lith.8 = factor(lith.8, levels = 1:7),
      lith.24 = factor(lith.24, levels = 1:15),
      Peat = factor(Peat, levels = c(5, 30, 60, 100)),
      lon = lon / 10^6,
      lat = lat / 10^6 # 20240730 newly added
      
    )
  return(new_d)
}

ca_base_new <- preproce(ca_covas, dep0 = 0)
ind <- rowSums(is.na(ca_base_new)) == 0
length(which(ind))
niubi <- ca_base_new[ind, ] %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

#2 load models
library(h2o)

h2o.init(nthreads = 64,
         max_mem_size = '300G',
         port = 15929)
Sys.sleep(20)
# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/StackedEnsemble_AllModels_4_AutoML_1_20241214_211305")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/GBM_4_AutoML_1_20241214_211305")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/XGBoost_grid_1_AutoML_1_20241214_211305_model_16")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/DRF_1_AutoML_1_20241214_211305")

# stack_m <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/StackedEnsemble_AllModels_4_AutoML_1_20241215_205434")
# gbm1 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/GBM_grid_1_AutoML_1_20241215_205434_model_13")
# xgb2 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/XGBoost_grid_1_AutoML_1_20241215_205434_model_20")
# rf3 <- h2o.loadModel("/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/DRF_1_AutoML_1_20241215_205434")

stack_m <- h2o.loadModel(
  "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/StackedEnsemble_AllModels_4_AutoML_5_20241231_120738"
)
gbm1 <- h2o.loadModel(
  "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/GBM_4_AutoML_5_20241231_120738"
)
xgb2 <- h2o.loadModel(
  "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/XGBoost_grid_1_AutoML_5_20241231_120738_model_20"
)
rf3 <- h2o.loadModel(
  "/mnt/DataSpace/Projects/Canada_C/models_add_lehigh/bd_gbm/DRF_1_AutoML_5_20241231_120738"
)



models <- list(stack_m, gbm1, xgb2, rf3)

#
dep00 <- c(0, 5, 15, 30, 60, 100, 200, 300)

pred_fun0 <- function(model, t_data) {
  maps0 <- h2o.predict(model, t_data)
  ca_base_new[ind, 'pred'] <- as.data.frame(maps0)['predict']
  pred_r <- ca_covas$BIO1
  pred_r$pred = ca_base_new$pred
  return(pred_r$pred)
}
#
lt1000 <- function(x) {
  x[which(x > 28.07217)] <- 28.07217
  return(x)
}
#
Sys.sleep(200)
for (i in dep00) {
  niubi$Depth = i
  niubi$SOC = SOC_pred[, as.character(i)]
  test <- as.h2o(x = niubi)
  maps <- lapply(models, pred_fun0, t_data = test)
  maps0 <- do.call(c, maps)
  maps1 <- terra::clamp(maps0, upper = 28.07217, values = T)
  maps1 <- exp(sqrt(maps1)) / 100
  maps1$sd <- app(maps1[[2:4]], cores = 80, fun = 'sd')
  maps1$mean <- app(maps1[[2:4]], cores = 80, fun = 'mean')
  maps1$sd <- subst(maps1$sd, 0, minmax(maps1$sd)[2])
  maps1$cv <- maps1$sd / maps1$mean
  names(maps1) <- c('Stack', 'GBM', 'XGBoost', 'RF', 'sd', 'mean', 'cv')
  terra::writeRaster(
    maps1,
    paste0(
      '/mnt/DataSpace/Projects/Canada_C/SOC_maps/base_BD_addSOC_',
      i,
      'cm_automl1215_v3_soc_peatv2.tif'
    ),
    overwrite = TRUE
  )
}


f1 <- list.files(
  pattern = 'tif',
  '/mnt/DataSpace/Projects/Canada_C/SOC_maps/',
  full.names = T
)
f2 <- grep('automl1215_v3_soc_peat', f1, value = T)[c(1, 7, 3, 6, 8, 2, 4, 5)]
f3 <- lapply(f2, function(x)
  rast(x)[[1]]) # "StackedEnsemble" "GBM"             "XGBoost"         "DRF"   [5] "sd"    [6]  "mean"    "cv"
f4 <- do.call(c, f3)
f5 <- (f4[[1:7]] + f4[[2:8]]) / 2
plot(
  f5[[1:5]],
  breaks = c(0, 0.1, 0.3, 0.6, 1.5, 2),
  col = c('red', 'red', 'blue', 'black', 'black')
)
writeRaster(
  f5,
  '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/BD_basedon_ml20241214.tif',
  overwrite = T
) # just add SOC

####
socs <- list.files(
  pattern = 'tif',
  '/mnt/DataSpace/Projects/Canada_C/SOC_maps/',
  full.names = T
)

all_maps <- grep('v3.tif', socs, value = T)[c(1, 7, 3, 6, 8, 2, 4, 5)]

all_maps2 <- lapply(all_maps, function(x)
  rast(x)[[1]]) # Stack,        GBM,     XGBoost,          RF,
all_maps3 <- do.call(c, all_maps2)
all_maps4 <- (all_maps3[[1:7]] + all_maps3[[2:8]]) / 2

f4 <- all_maps4
BD0 <- 0.08403 + 1.34364 * exp(-0.00874 * f4)
writeRaster(
  BD0,
  '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/BD_based_on_SOC_20250107_newlyasused.tif'
)
socs_ml <- f4 * BD0 * c(5, 10, 15, 30, 40)
socs_ml_1m <- sum(socs_ml[[1:5]]) / 100
writeRaster(
  socs_ml_1m,
  '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/canada_SOCs_1m_kgm2_expoential_relationship_SOC_BD_with_originaldata_20250107.tif'
)

for (i in 2:4) {
  socs <- list.files(
    pattern = 'tif',
    '/mnt/DataSpace/Projects/Canada_C/SOC_maps/',
    full.names = T
  )
  all_maps <- grep('v3.tif', socs, value = T)[c(1, 7, 3, 6, 8, 2, 4, 5)]
  modname <- c('Stack', 'GBM', 'XGBoost', 'RF')[i]
  all_maps2 <- lapply(all_maps, function(x)
    rast(x)[[i]]) # Stack,        GBM,     XGBoost,          RF,
  all_maps3 <- do.call(c, all_maps2)
  all_maps4 <- (all_maps3[[1:7]] + all_maps3[[2:8]]) / 2
  f4 <- all_maps4
  BD0 <- 0.08403 + 1.34364 * exp(-0.00874 * f4)
  socs_ml <- f4 * BD0 * c(5, 10, 15, 30, 40)
  socs_ml_1m <- sum(socs_ml[[1:5]]) / 100
  writeRaster(
    socs_ml_1m,
    paste0(
      '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/canada_SOCs_1m_kgm2_expoential_relationship_SOC_BD_with_originaldata_20250107_',
      modname,
      '.tif'
    )
  )
}
writeRaster(
  socs_ml_1m,
  '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/canada_SOCs_1m_kgm2_expoential_relationship_SOC_BD_with_originaldata_20250107_GBM.tif'
)
writeRaster(
  socs_ml_1m,
  '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/canada_SOCs_1m_kgm2_expoential_relationship_SOC_BD_with_originaldata_20250107_XGB.tif'
)
writeRaster(
  socs_ml_1m,
  '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/canada_SOCs_1m_kgm2_expoential_relationship_SOC_BD_with_originaldata_20250107_RF.tif'
)


for (i in 2:4) {
  socs <- list.files(
    pattern = 'tif',
    '/mnt/DataSpace/Projects/Canada_C/SOC_maps/',
    full.names = T
  )
  
  all_maps <- grep('v3.tif', socs, value = T)[c(1, 7, 3, 6, 8, 2, 4, 5)]
  modname <- c('Stack', 'GBM', 'XGBoost', 'RF')[i]
  all_maps2 <- lapply(all_maps, function(x)
    rast(x)[[i]]) # Stack,        GBM,     XGBoost,          RF,
  all_maps3 <- do.call(c, all_maps2)
  all_maps4 <- (all_maps3[[1:7]] + all_maps3[[2:8]]) / 2
  # ca_f <- f5[[1:5]] * all_maps4[[1:5]] * c(5,10,15,30,40)
  f4 <- all_maps4
  BD0 <- 0.08403 + 1.34364 * exp(-0.00874 * f4)
  depth00 <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/ca_bed_depth_of_each_of5layer.tif')
  socs_ml <- f4 * BD0 * depth00
  socs_ml_1m <- sum(socs_ml[[1:5]]) / 100
  
  writeRaster(
    socs_ml_1m,
    paste0(
      '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/canada_SOCs_1m_kgm2_expoential_relationship_SOC_BD_with_originaldata_20250107_',
      modname,
      '_considerbedrock_for_sd.tif'
    )
  )
}

### 20250115
depth00 <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/ca_bed_depth_of_each_of5layer.tif')
socs_ml <- f4 * BD0 * depth00
socs_ml_1m <- sum(socs_ml[[1:5]]) / 100
plot(socs_ml_1m)#
writeRaster(
  socs_ml_1m,
  '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/canada_SOCs_1m_kgm2_expoential_relationship_SOC_BD_with_originaldata_20250115_considerBEdROCK.tif'
)
#169.8796 Pg + landcover 149.1977
### 20250115 end 30cm 95.05967



# then we got

ca_f1 <- sum(ca_f[[1:5]])
ca_pure_area84 <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/ca_cmip6_area_WGS84.tif') #new created
ca_f2 <- ca_f1 * ca_pure_area84
global(ca_f2, function(x)
  sum(x, na.rm = T)) / 10^14
plot(ca_f1 / 100, breaks = c(0, 3, 5, 10, 20, 30, 40, 50, 70, 100, 150, 200, 400))
writeRaster(
  ca_f1 / 100,
  '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/canada_SOCs_1m_kgm2_expoential_relationship_SOC_BD_with_originaldata_20241230.tif',
  overwrite = T
) #both ml for soc and bd

socs <- list.files(
  pattern = 'tif',
  '/mnt/DataSpace/Projects/Canada_C/SOC_maps/',
  full.names = T
)
socs <- list.files(
  pattern = 'tif',
  '/mnt/Fastrun/Canada_C/supplementary_plant_pet_climate_sm_1230',
  full.names = T
)
all_maps <- grep('v3.tif', socs, value = T)[c(1, 7, 3, 6, 8, 2, 4, 5)]
all_maps2 <- lapply(all_maps, function(x)
  rast(x)[[6]]) ##
all_maps3 <- do.call(c, all_maps2)
all_maps4 <- (all_maps3[[1:7]] + all_maps3[[2:8]]) / 2
ca_f <- f5[[1:5]] * all_maps4[[1:5]] * c(5, 10, 15, 30, 40)
ca_f1 <- sum(ca_f[[1:5]])

writeRaster(
  ca_f1 / 100,
  '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/canada_SOCs_1m_kgm2_mean_expoential_relationship_SOC_BD_with_originaldata_20241230.tif',
  overwrite = T
)
socs <- list.files(
  pattern = 'tif',
  '/mnt/DataSpace/Projects/Canada_C/SOC_maps/',
  full.names = T
)
all_maps <- grep('v3.tif', socs, value = T)[c(1, 7, 3, 6, 8, 2, 4, 5)]
all_maps2 <- lapply(all_maps, function(x)
  rast(x)[[5]])
all_maps3 <- do.call(c, all_maps2)
all_maps4_sd <- (all_maps3[[1:7]] + all_maps3[[2:8]]) / 2
upper <- all_maps4 + 1.64 * all_maps4_sd
lower <- all_maps4 - 1.64 * all_maps4_sd
ca_f2 <- f5[[1:5]] * upper[[1:5]] * c(5, 10, 15, 30, 40)
ca_f_up <- sum(ca_f2[[1:5]])
ca_f3 <- f5[[1:5]] * lower[[1:5]] * c(5, 10, 15, 30, 40)
ca_f_low <- sum(ca_f3[[1:5]])

ca_pure_area84 <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/ca_cmip6_area_WGS84.tif') #new created
ca_u <- ca_f_up * ca_pure_area84
global(ca_u, function(x)
  sum(x, na.rm = T)) / 10^14
ca_l <- ca_f_low * ca_pure_area84
global(ca_l, function(x)
  sum(x, na.rm = T)) / 10^14


writeRaster(
  ca_f_up / 100,
  '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/canada_SOCs_1m_kgm2_upper95percent_expoential_relationship_SOC_BD_with_originaldata_20241230.tif',
  overwrite = T
)
writeRaster(
  ca_f_low / 100,
  '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/canada_SOCs_1m_kgm2_lower5percent_expoential_relationship_SOC_BD_with_originaldata_20241230.tif',
  overwrite = T
)


writeRaster(
  BD0,
  '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/BD_based_on_SOC.tif'
) #old

plot(ca_f1 / 100, breaks = c(0, 3, 5, 10, 20, 30, 40, 50, 70, 200, 400))
plot(rast(f2)[[1]])
plot(all_maps4, breaks = c(0, 450, 550, 600))
BD0 <- rast(
  '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/BD_based_on_SOC.tif'
)
plot(
  BD0,
  breaks = c(0, 0.1, 0.3, 0.4, 0.5, 0.6, 1.5, 2),
  col = c('red', 'red', 'blue', 'black', 'black', 'black', 'black')
)
BD0 <- rast(
  '/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/BD_basedon_ml20241214.tif'
)
