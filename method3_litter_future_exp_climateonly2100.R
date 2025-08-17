# 20040801 update 10 year 

#
#
############################################################################
############################################################################
############################################################################
############################################################################

#   outline 
#> 1.  plant ++++  model_gpp  ==  to 2100 and to 2300 
#> 2. all ++++ plant + pet + climate + sm  == to 2100 as, pet just to 2100 
#> 3. water related ++++ pet + precipitation(8) + sm  == 2100 , as pet just to 2100, 
#> #>  some of which we could just precipitation to 2300 for that five models, 
#> #> but from the results of SOC change, we could deduced that temperature is the most important factor
#> 4. temperature related ++++ solar + temperature  == to 2100 first. 9 esm models 
#> 4-1. temperature related ++++ solar + temperature (11) == to 2300 then, 4 esm models, we could know that
#> as mojority change of SOC comes from T, so we use the 4 models alone could fullfill the request. 

############################################################################
############################################################################
############################################################################
############################################################################
#from /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/005_python_EC/0006_GPP_EC_concise_and_used.R
model_gpp_ec <- c('ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','CanESM5-1','CESM2-WACCM','EC-Earth3-CC','EC-Earth3-Veg','EC-Earth3-Veg-LR','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','TaiESM1')

model_gpp <- c("EC-Earth3-Veg-LR", "IPSL-CM6A-LR", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CanESM5-1", "CAS-ESM2-0", "CESM2-WACCM", "CMCC-CM2-SR5", "EC-Earth3-CC", "INM-CM4-8", "INM-CM5-0", "TaiESM1", "EC-Earth3-Veg", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR")
#2300 'CanESM5','IPSL-CM6A-LR','ACCESS-ESM1-5','EC-Earth3-Veg'(both 245 and 585)
# model_pet <- c("ACCESS-CM2", "ACCESS-ESM1-5", "CanESM5", "CanESM5-1", "CAS-ESM2-0", "FGOALS-g3", "GFDL-ESM4", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "EC-Earth3-CC", "MPI-ESM1-2-LR")
model_pet <- c('ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CanESM5-1','CAS-ESM2-0','FGOALS-g3','GFDL-ESM4','IPSL-CM6A-LR','MIROC6','MPI-ESM1-2-HR','MRI-ESM2-0','EC-Earth3-CC','MPI-ESM1-2-LR','EC-Earth3-Veg-LR','INM-CM4-8','INM-CM5-0')

model_2100_new <- c("ACCESS-CM2","CanESM5","CMCC-ESM2","EC-Earth3-Veg","FIO-ESM-2-0","INM-CM4-8","IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0",  "ACCESS-ESM1-5", "EC-Earth3-Veg-LR", "INM-CM5-0" ,"BCC-CSM2-MR","CanESM5-1","CAS-ESM2-0", "EC-Earth3-CC")
model_2300_clim <- c("CanESM5", "EC-Earth3-Veg", "IPSL-CM6A-LR", "MRI-ESM2-0", "ACCESS-ESM1-5")
model_2100_only <- setdiff(model_2100_new,model_2300_clim)
model_soil_water <- c("ACCESS-CM2","ACCESS-ESM1-5", "CanESM5"  , "CMCC-ESM2"  ,  "EC-Earth3-Veg" ,"EC-Earth3-Veg-LR" ,"IPSL-CM6A-LR", "MIROC6","MPI-ESM1-2-HR" ,   "MPI-ESM1-2-LR",    "MRI-ESM2-0" )
model_solar_2300 <- c('IPSL-CM6A-LR','ACCESS-ESM1-5','CanESM5','CESM2-WACCM','EC-Earth3-Veg') # all ssp585, just EC-Earth3-Veg has ssp245 to 2300
model_climate <- intersect(intersect(model_pet,model_2100_new),model_soil_water)
model_climate_2300 <- intersect(intersect(model_pet,model_2300_clim),model_soil_water) # waters all to 2100 as pet is to 2100 , and it is not influence that much so ignor it
model_all <- intersect(intersect(intersect(model_gpp_ec,model_pet),model_2100_new),model_soil_water) #5 to 2100 
model_water_2100 <- intersect(intersect(model_pet,model_2100_new),model_soil_water) # 8 to 2100 without solar, if with solar, only 5 
model_temp_2100 <- intersect(setdiff(model_gpp_ec,model_solar_2300),model_2100_new) # gpp is equal with solar + temperature or just temperature 
model_temp_2300 <- intersect(model_solar_2300,model_2100_new)
# model_gpp <- c("EC-Earth3-Veg-LR", "IPSL-CM6A-LR", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CanESM5-1", "CAS-ESM2-0", "CESM2-WACCM", "CMCC-CM2-SR5", "EC-Earth3-CC", "INM-CM4-8", "INM-CM5-0", "TaiESM1", "EC-Earth3-Veg", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR")
# #2300 'CanESM5','IPSL-CM6A-LR','ACCESS-ESM1-5','EC-Earth3-Veg'(both 245 and 585)
# model_pet <- c("ACCESS-CM2", "ACCESS-ESM1-5", "CanESM5", "CanESM5-1", "CAS-ESM2-0", "FGOALS-g3", "GFDL-ESM4", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "EC-Earth3-CC", "MPI-ESM1-2-LR")
# model_2100 <- c("ACCESS-CM2", "CanESM5", "CMCC-ESM2", "EC-Earth3-Veg", "FIO-ESM-2-0", "GISS-E2-1-G", "GISS-E2-1-H", "HadGEM3-GC31-LL", "INM-CM4-8", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "UKESM1-0-LL", "ACCESS-ESM1-5", "CanESM5-CanOE", "CNRM-CM6-1", "CNRM-CM6-1-HR", "CNRM-ESM2-1", "EC-Earth3-Veg-LR", "INM-CM5-0")
# model_2300 <- c('GISS-E2-1-G','ACCESS-CM2','ACCESS-ESM1-5','IPSL-CM6A-LR',
#                 'MRI-ESM2-0','CanESM5')
# model_all <- unique(c(model_2100,model_2300))
# 
# plant_pet_clim <- intersect(model_gpp,intersect(model_pet,model_2100))
# # we use solar from gpp npp solar, but pet also has solar, so we use pet solar to produce solar
# solars_map <- intersect(model_gpp,model_2100)
# pet_map_solar <- intersect(model_pet,model_2100) # all need 
# pet_solor_map <- intersect(model_gpp,pet_map_solar) # currently we have 
# setdiff(pet_map_solar,pet_solor_map) 
# #> wee still needed see using rsds from /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/001code/00_1_01complementary_PET.R
# #> but function from /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/001code/00_1complementary_experiments.R
# #> see L 157 in /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/001code/00_1complementary_experiments.R line that typed # for /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/001code/01_2_h2o_versoin2_complementary_PET.R
# 
# to2100_all <-  c("IPSL-CM6A-LR" , "ACCESS-ESM1-5", "CanESM5" ,      "MPI-ESM1-2-HR", "MPI-ESM1-2-LR")
# to2300_all <- c("IPSL-CM6A-LR" , "ACCESS-ESM1-5", "CanESM5") # with out pet, all pet to 2100. 
# so the workflow should be 
#> step1 do all gpp npp ndvi evi
#> step2 do all solar pet ai 
#> step3 do all throughout 5 to 2100 
#> step4 do all gpp + climate 3 without pet to 2300 
##########################################################################################
#
#one model for each pixel 
h2_arg <- commandArgs(trailingOnly = T)
stopifnot(length(h2_arg) >0 )
begin <- as.numeric(h2_arg[1])
end <- as.numeric(h2_arg[2])
port0 <- as.numeric(h2_arg[3])
mk <- h2_arg[4]
#2300
#model_ind <- h2_arg[4]

##############################################################################################

#######
unixtools::set.tempdir('/mnt/Fastrun/temp4r')

library(terra)
library(dplyr)
library(tidyr)
library(h2o)
#library(readr)
library(data.table) #fread fwrite
#
#ca_covas <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global_re_align.tif')

# as new 2300 future climate has different extend of the fomer canada_all_covs_global_re_align one ,we make the realign to new future data for future climate 2300 ,
#see /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/001code/00CMIP2300_Preprocesssing2.R

# ca_covas <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global_re_align_2300_d1.tif')

#writeRaster(ca_covas1,'/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global_re_align_2300_d1.tif') #digit =1 version

# names(ca_covas) <- c('NDVI1','EVI2','NPP3','GPP4','PET5','AI6',
#                      'BIO1','BIO10','BIO11','BIO12','BIO13','BIO14','BIO15','BIO16','BIO17','BIO18','BIO19','BIO2','BIO3','BIO4','BIO5','BIO6','BIO7','BIO8','BIO9',
#                      'Solar26','So_seasonal27','Pop_den28','HFP29','Vrm30','Tcurv31','Pcurv32','Ele33','Slope34','Asp_c35','Asp_s36','Est37','Nor38','Roug39','TPI40','TRI41','Dx42','Dxx43','Dy44','Dyy45',
#                      'Maj48','lith.8','lith.24','WTD85','smap','SUMAP','Ref_band1','ref_band2','ref_band3','ref_band4','ref_band5','ref_band6','ref_band7')
#


#cmip6_ca <- list.files(pattern = 'tif',path = '/mnt/DataSpace/Data_pool/CMIP6_1km_ca',recursive = TRUE,full.names = TRUE)

#2300 
# cmip6_ca1 <- list.files(pattern = 'tif$',path = '/mnt/DataSpace/Projects/CMIP2300/Cmip6/clim_bio19_5m',recursive = F,full.names = TRUE)
# #2300 new 4models like access and giss
# cmip6_ca2 <- list.files(pattern = 'tif$',path = '/mnt/DataSpace/Projects/CMIP2300/Cmip6/clim_bio19_5m0',recursive = F,full.names = TRUE)
# cmip6_ca <- c(cmip6_ca1,cmip6_ca2)[c(1:4)] #parallel when 4 cores, it is 1:4, when for the left 3 cores , 5:7 
#
# cmip6_ca_raster <- purrr::map(cmip6_ca,rast)
###


#
dep00 <- c(0,5,15,30,60,100,200,300)
#ca_base_new <- preproce(ca_covas,dep0 = 0)
#fwrite(ca_base_new[ind,],'/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ca_base_new_df.csv')
# hard drive
#ca_base_new <- fread('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ca_base_new_df.csv')#ssd 
#20240329 changed to 
#system.time(ca_base_new <- fread('/mnt/Fastrun/Canada_C/rawdata/004_sensitive_analysis/ca_base_new_df2024.csv'))
#ca_base_new <- fread('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ca_base_new_df20240110.csv')
#20240801 change to 
# system.time(ca_base_new <- fread('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ca_base_new_df20240730.csv')) # from /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/006_Update_after_EC_CSSS/001_h2o_automl_build_SOC_model.R
system.time(ca_base_new <- fread('/mnt/Fastrun/Canada_C/ca_base_new_df20240730.csv'))
# as ca_base_new has already preproce so already loged now just need factor and log when in bio4 12 18 
ca_base_new <- ca_base_new %>% mutate(
  Maj48 = factor(Maj48, levels = 1:10),
  lith.8 = factor(lith.8, levels = 1:7),
  lith.24 = factor(lith.24, levels = 1:15)
) #%>% mutate(across(where(is.numeric), ~ round(., 1)))
#ind <- rowSums(is.na(ca_base_new)) == 0
#saveRDS(ind,'/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ind.rds')
#hand drive
#ind <- readRDS('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ind.rds')
#ssd
# ind <- readRDS('/mnt/Fastrun/Canada_C/rawdata/004_sensitive_analysis/ind2024.rds')
#20240801 change to 
#ind <- readRDS('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_soc_git_repo/004_sensitive_analysis/ind20240730.rds')
ind <- readRDS('/mnt/Fastrun/Canada_C/ind20240730.rds')
#length(which(ind))

#2 load models


h2o.init(nthreads = 12,max_mem_size = '100G',port = port0)
#how ensemble model produced
#https://github.com/h2oai/h2o-tutorials/blob/master/tutorials/ensembles-stacking/README.md
#20230704_121940 is for depth as categorical
# stack_m <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/bestfamiliy_ensemble/StackedEnsemble_BestOfFamily_6_AutoML_1_20230704_151137')
# gbm1 <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/gbm/GBM_lr_annealing_selection_AutoML_1_20230704_151137_select_model')
# xgb2 <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/xgb/XGBoost_lr_search_selection_AutoML_1_20230704_151137_select_grid_model_4')
# deepl4 <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/deeplearn/DeepLearning_grid_1_AutoML_1_20230704_151137_model_4')
# rf3 <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/rf/DRF_1_AutoML_1_20230704_151137')
#
# models <- list(stack_m,gbm1,xgb2,rf3, deepl4) #acturalll only should use gbm xgb deep, we here add rf


#
#test <- as.h2o(x = ca_base_new[ind, ])
#modified based on new version, 20230804
pred_fun0 <- function(model,t_data) {
  maps0 <- h2o.predict(model,t_data)
  var_use_df[ind,'pred'] <- as.data.frame(maps0)['predict']
  plant_r$pred = var_use_df$pred
  return(plant_r$pred)
}
#
#which is larger than 3, which means larger than 1000
lt1000 <- function(x) {
  x[which(x > 2.7634)] <- 2.7634  # som to soc log10(580) ; 1000/1.724
  return(x)
}

#

########


################
#### 2300 new like access-em2 giss 245 etc.this is just for three model to 2300 


# periods0 <- data.frame(a=seq(2021,2300,20),b=seq(2040,2300,20))
# for (i in 1) {
#   cmip6_df0 <- as.data.frame(as.matrix(cmip6_ca_raster[[esm0]][[1:19 + (i-1)*19]]))
#   #system.time(cmip6_df1 <- as(raster(cmip6_ca_raster[[j]]),'SpatialPixelsDataFrame'))
#   for (k in 18:19) {
#   ca_base_new[, paste0('BIO', k)] <-
#     cmip6_df0[, k]  #note the order
#   ca_base_new <- ca_base_new %>% mutate(
#     BIO4 = round(log10(BIO4),1),
#     BIO12 = round(log10(BIO12),1),
#     BIO18 = round(log10(BIO18),1)) %>% mutate(across(where(is.numeric), ~ round(., 1)))
#   for (j in dep00) {
#     stack_m <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/bestfamiliy_ensemble/StackedEnsemble_BestOfFamily_6_AutoML_1_20230704_151137')
#     # gbm1 <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/gbm/GBM_lr_annealing_selection_AutoML_1_20230704_151137_select_model')
#     # xgb2 <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/xgb/XGBoost_lr_search_selection_AutoML_1_20230704_151137_select_grid_model_4')
#     # deepl4 <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/deeplearn/DeepLearning_grid_1_AutoML_1_20230704_151137_model_4')
#     # rf3 <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/rf/DRF_1_AutoML_1_20230704_151137')
#     # models <- list(stack_m,gbm1,xgb2,rf3, deepl4) 
#     ca_base_new$Depth = j #8 depth
#     test <- as.h2o(x = ca_base_new[ind,]) #to h2o frame
#     #maps <- lapply(models, pred_fun0, t_data = test) #predict
#     maps <- pred_fun0(model=stack_m,t_data=test)
#     #maps0 <- do.call(c, maps)
#     maps1 <- app(maps, lt1000, cores = 20)
#     # maps1 <- 10 ^ maps1
#     # maps1$sd <- app(maps1[[2:5]], cores = 24, fun = 'sd')
#     # maps1$mean <- app(maps1[[2:5]], cores = 24, fun = 'mean')
#     # maps1$sd <-
#     #   subst(maps1$sd, 0, minmax(maps1$sd)[2]) #sd is 0 refers to the condition that all layers is 1000, which is larger than normal. we substitute them.
#     # maps1$cv <- maps1$sd / maps1$mean
#     # names(maps1) <- c('Stack', 'GBM', 'XGBoost', 'RF', 'DeepL', 'sd', 'mean', 'cv')
#     names(maps1) <- 'Stack'
#     terra::writeRaster(10^maps1,
#                        paste0('/mnt/DataSpace/Projects/Canada_C/SOC_maps/sensitive_analysis/','BIO',k,'_',j,'cm_','ca_wc2.1_30s_bioc_',
#                               gsub('2020_2300_CA.tif','',basename(cmip6_ca[esm0])),periods0[i,1],'-',periods0[i,2],'.tif'),overwrite = TRUE)
#     h2o.removeAll()
#     h2o.removeAll()
#     h2o.removeAll()
#     }
#   h2o.shutdown(prompt = FALSE)
#   ca_base_new <- preproce(ca_covas,dep0 = 0)#just need time to reinitial h2o so this and the following line
#   ind = rowSums(is.na(ca_base_new)) == 0 #just time to waste
#   h2o.init(nthreads = 24,max_mem_size = '250G', port = port0)
#   stack_m <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/bestfamiliy_ensemble/StackedEnsemble_BestOfFamily_6_AutoML_1_20230704_151137')
#   # gbm1 <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/gbm/GBM_lr_annealing_selection_AutoML_1_20230704_151137_select_model')
#   # xgb2 <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/xgb/XGBoost_lr_search_selection_AutoML_1_20230704_151137_select_grid_model_4')
#   # deepl4 <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/deeplearn/DeepLearning_grid_1_AutoML_1_20230704_151137_model_4')
#   # rf3 <- h2o.loadModel('/mnt/DataSpace/Projects/Canada_C/models/rf/DRF_1_AutoML_1_20230704_151137')
#   #
#   # models <- list(stack_m,gbm1,xgb2,rf3, deepl4)
# }}
# h2o.shutdown(prompt = FALSE)
#39724s


############# new version for each 
#var_p <- c('NDVI1','EVI2','NPP3','GPP4')
#var_p <- c('PET5','AI6','Solar26',paste0('BIO',12:19),'smap','SUMAP')
var_p <- paste0('BIO',1:19)

# to get gpp np 
#20240329 changed to 
#f1 <- list.files(pattern = 'tif',path = '/mnt/Fastrun/Data_pool/01cmip6_processed_output',full.names = T)
#f1 <- list.files(pattern = 'tif',path = '/mnt/DataSpace/Data_pool/01cmip6_processed_output_update_npp',full.names = T)
f1 <- list.files(pattern = 'tif',path = '/mnt/DataSpace/Data_pool/01cmip6_processed_output_EC_GPP',full.names = T)

# to get solar
#f100 <- list.files(pattern = 'tif',path = '/mnt/DataSpace/Data_pool/cmip6_complementary/processed_output',full.names = T)
# f100 <- list.files(pattern = 'tif',path = '/mnt/Fastrun/Data_pool/01cmip6_processed_output',full.names = T)
f100 <- list.files(pattern = 'tif',path = '/mnt/DataSpace/Data_pool/01_rsds_output_10',full.names = T)

# to get pet
#f200 <- list.files(pattern = 'tif',path = '/mnt/DataSpace/Data_pool/PET_cmip/processed_pet_each',full.names = T)
# f200 <- list.files(pattern = 'tif',path = '/mnt/Fastrun/Data_pool/03pet_use_individual',full.names = T)
f200 <- list.files(pattern = 'tif',path = '/mnt/DataSpace/Data_pool/03pet_use_individual_10',full.names = T)

# to get precipitation 
#f300 <- list.files(pattern = 'tif',path = '/mnt/DataSpace/Data_pool/CMIP6_1km_ca',recursive = TRUE,full.names = TRUE)
# f300 <- list.files(pattern = 'tif',path = '/mnt/Fastrun/Data_pool/04_bios_variable',recursive = TRUE,full.names = TRUE)
f300 <- list.files(pattern = 'tif',path = '/mnt/DataSpace/Data_pool/03cmip6_processed_bio19_5min_10',recursive = TRUE,full.names = TRUE)

#f400 <- list.files(pattern = 'tif',path = '/mnt/DataSpace/Data_pool/05_smap')
# f400 <- list.files(pattern = 'tif',path = '/mnt/Fastrun/Data_pool/05_smap',full.names = T)
f400 <- list.files(pattern = 'tif',path = '/mnt/DataSpace/Data_pool/05_smap_10',full.names = T)

#cmip6_ca_raster <- purrr::map(f300,rast)
#cmip6_df0 <- as.data.frame(as.matrix(cmip6_ca_raster[[i]]))  #40s
#system.time(cmip6_df1 <- as(raster(cmip6_ca_raster[[j]]),'SpatialPixelsDataFrame'))
#ca_base_new[, paste0('BIO', 1:19)] <-  cmip6_df0[, 1:19]  #note the order


####

models <- model_2100_new[begin:end] # 18 
#models <- 'EC-Earth3-Veg-LR'
for (i in models) {
  for(scenario in c('ssp245','ssp585')) {
    #f_solar <- grep(paste0('rsds_',i,'_',scenario),f100,value = T)
    #f_pet <- grep(paste0('pet_',i,'_',scenario),f200,value = T)
    f_map <- grep(paste0(i,'_',scenario),f300,value = T) # we changed names in this folder like 
    # f_plant <- grep(paste0(i,'_',scenario),f1,value = T)
    # f_smap <- grep(paste0(i,'_',scenario),f400,value = T)
    
    ##> /mnt/DataSpace/Data_pool/CMIP6_1km_ca/CMIP6_1km_ca_idm/ , files that with ca_disaggregated put that in the end , for example
    ##> from
    ##> [1] "/mnt/DataSpace/Data_pool/CMIP6_1km_ca/CMIP6_1km_ca_idm/ca_disaggregatedwc2.1_30s_bioc_CanESM5_ssp585_2081-2100.tif"
    ##> [2] "/mnt/DataSpace/Data_pool/CMIP6_1km_ca/CMIP6_1km_ca_idm/ca_wc2.1_30s_bioc_CanESM5_ssp585_2021-2040.tif"             
    ##> [3] "/mnt/DataSpace/Data_pool/CMIP6_1km_ca/CMIP6_1km_ca_idm/ca_wc2.1_30s_bioc_CanESM5_ssp585_2041-2060.tif"             
    ##> [4] "/mnt/DataSpace/Data_pool/CMIP6_1km_ca/CMIP6_1km_ca_idm/ca_wc2.1_30s_bioc_CanESM5_ssp585_2061-2080.tif"     
    ##> to
    ##>  
    ##> [1] "/mnt/DataSpace/Data_pool/CMIP6_1km_ca/CMIP6_1km_ca_idm/ca_wc2.1_30s_bioc_CanESM5_ssp585_2021-2040.tif"          
    ##> [2] "/mnt/DataSpace/Data_pool/CMIP6_1km_ca/CMIP6_1km_ca_idm/ca_wc2.1_30s_bioc_CanESM5_ssp585_2041-2060.tif"          
    ##> [3] "/mnt/DataSpace/Data_pool/CMIP6_1km_ca/CMIP6_1km_ca_idm/ca_wc2.1_30s_bioc_CanESM5_ssp585_2061-2080.tif"          
    ##> [4] "/mnt/DataSpace/Data_pool/CMIP6_1km_ca/CMIP6_1km_ca_idm/wc2.1_30s_bioc_CanESM5_ssp585_2081-2100disaggregated.tif"
    
    #if(identical(length(f_solar),length(f_pet),length(f_map))) {
    for (k in 1:10) {  # as pet always 2100
      #    for (k in 1:(length(f3)/3)) {
      # f_ind <- f_plant[seq(k, length(f_plant), length(f_plant) / 3)] #time period
      # gpp <- rast(f_ind[1]) / 10000
      # npp <- rast(f_ind[2]) / 10000
      # #ndvi_f <- readRDS('/mnt/DataSpace/Data_pool/PET_cmip/ndvi_gpp_npp_formula.rds');print(nvdi_f)
      # #evi_f <- readRDS('/mnt/DataSpace/Data_pool/PET_cmip/evi_gpp_npp_formula.rds');print(evi_f)
      # ndvi1 <- npp * 4.379 - gpp * 1.574
      # evi1 <- npp * 1.16274 - gpp * 0.08355
      # #original
      # npp0 <- rast('/mnt/Fastrun/Canada_C/rawdata/forfinalmap_1km/NPP3_ca.tif')
      # gpp0 <- rast('/mnt/Fastrun/Canada_C/rawdata/forfinalmap_1km/GPP4_ca.tif')
      # #old generated ndvi and evi compare to new generated ones
      # ndvi0 <- npp0 * (4.379/10000) - gpp0 * (1.574/10000)
      # evi0 <- npp0 * (1.16274/10000) - gpp0 * (0.08355/10000)
      # #
      # # t_nvdi <- extend(ndvi0,ndvi1)
      # # t_ndvi0 <- crop(ndvi0,ext(ndvi1))
      # # t_evi <- extend(evi0,evi1)
      # # t_evi0 <- crop(evi0,ext(evi1))
      # #
      # ndvi_multip <- ndvi1/ndvi0 # t_nvdi0
      # evi_multip <- evi1/evi0 #t_evi0
      # ndvi_multip <- clamp(ndvi_multip, upper=10, lower= -10, values=T)
      # evi_multip <- clamp(evi_multip, upper=10, lower= -10, values=T)
      # # original vndvi 
      # ndvi00 <- rast('/mnt/Fastrun/Canada_C/rawdata/forfinalmap_1km/NDVI1_ca.tif')
      # evi00 <- rast('/mnt/Fastrun/Canada_C/rawdata/forfinalmap_1km/EVI2_ca.tif')
      # #
      # # t_nvdi <- extend(ndvi00,ndvi1)
      # # t_ndvi00 <- crop(ndvi00,ext(ndvi1))
      # # t_evi <- extend(evi00,evi1)
      # # t_evi00 <- crop(evi00,ext(evi1))
      # #
      # ndvi_new <- ndvi00/10000 * ndvi_multip  #t_ndvi00
      # evi_new <- evi00/10000 * evi_multip  #t_evi00
      #
      # solar26 <- rast(f_solar[k]) / 10000 # based on function preproce
      # pet5 <- round(log10(rast(f_pet[k])),1)
      #   clim_ind <- seq(k,length(f_map),length(f_map)/19)[c(1,12:19,2:11)]
      # map12 <- rast(f_map[clim_ind])[[12]] # original extend, unfortunityly base data csv using current ext with ca_shp2 (not from world) + 1
      # # t10 <- terra::extend(map12,pet5)
      # # t20 <- terra::crop(t10,ext(pet5))
      # # ai6 <- t20/rast(f_pet[k])
      # ai6 <- map12/rast(f_pet[k])
      # ai60 <- clamp(ai6,upper=7,values=TRUE)
      clim19 <- rast(f_map[k])
      names(clim19) <- paste0('BIO',1:19)
      # t19 <- terra::extend(clim19,pet5)
      # t270 <- terra::crop(t19,ext(pet5))
      # t270[[c('BIO4','BIO12','BIO18')]] <- round(log10(t270[[c(4,12,18)]]),1)
      # clim19[[c('BIO4','BIO12','BIO18')]] <- round(log10(clim19[[c(4,12,18)]]),1) #note now 12 and 18 == 1 and 7
      clim19[[12:19]] <- log(clim19[[12:19]]) #note now 12 and 18 == 1 and 7
      clim19[[4]] <- clim19[[4]] / 100 
      #
      # smap <- rast(grep('SMAP',f_smap,value=T)[k])
      # sumap <- rast(grep('SUMAP',f_smap,value=T)[k])
      # par(mfrow=c(2,2))
      # plot(map12)
      # plot(rast(f_pet[k]))
      # plot(ai60)
      # plot(ca_covas[[6]]/10000)
      
      # plant_r <- c(ndvi_new, evi_new, npp, gpp,pet5,ai60,solar26,t270,smap,sumap)
      # plant_r <- c(pet5,ai60,solar26,clim19,smap,sumap)
      plant_r <- clim19
      var_use_df <- as.data.frame(as.matrix(plant_r))
      names(var_use_df) <- var_p
      #
      
      # fls <- list.files(pattern = 'tif',path = '/mnt/DataSpace/Projects/CMIP2300/Cmip6/clim_bio19_each',full.names = T)
      # can_fls <- grep('CanESM',fls,value=TRUE)
      # r_fls <- setdiff(fls,can_fls)
      # bio <- esm0
      # bio_fls <- grep(paste0(bio,"_"),r_fls,value=TRUE)
      
      #
      #creat global variable for function pred00
      #bio_raster1 <- rast(fls[1])
      #cmip6_df1 <- as.data.frame(as.matrix(bio_raster1))
      #
      
      #
      periods0 <- data.frame(a=seq(2001,2300,10),b=seq(2010,2300,10))
      time00 <- paste0(seq(2001,2291,10),'-',seq(2010,2300,10))
      # for (i in bio_fls[begin:end]) {
      #   if (bio %in% c('BIO4', 'BIO12', 'BIO18')) {
      #     bio_raster <- log10(rast(i))
      #     bio_raster <- round(bio_raster,3)
      #   } else {
      #     bio_raster <- rast(i)
      #   }
      #   cmip6_df0 <- as.data.frame(as.matrix(bio_raster))
      #   
      # dep_df <- data.frame(x1=rep(0,length(which(ind))),#56637000
      #                      x2=rep(5,length(which(ind))),
      #                      x3=rep(15,length(which(ind))),
      #                      x4=rep(30,length(which(ind))),
      #                      x5=rep(60,length(which(ind))),
      #                      x6=rep(100,length(which(ind))),
      #                      x7=rep(200,length(which(ind))),
      #                      x8=rep(300,length(which(ind)))
      # )
      # ca_base_new[, bio] <- cmip6_df0[ind, 1]  #note the order
      ca_base_new[, var_p] <- var_use_df[ind, var_p]  #note the order
      
      h2o.init(nthreads = 12,max_mem_size = '100G', port = port0)
      # dep_hex <- as.h2o(dep_df)
      test <- as.h2o(x = ca_base_new) #to h2o frame
      #gc();gc();gc();gc();
      # for (j in 1:8) {
      #   test$Depth <- dep_hex[[paste0('x',j)]]
      # stack_m <- h2o.loadModel('/mnt/Fastrun/Canada_C/rawdata/model/StackedEnsemble_AllModels_6_AutoML_1_20231128_125909')
      #model000 <- grep(pattern = "GBM", x = list.files(path = paste0('/mnt/DataSpace/Projects/Canada_C/turnover_constrain/traind_model/',mk),full.names = T),value=T)
      #20240329 changed to   
      # model000 <- grep(pattern = glob2rx("*GBM*20240731*"), x = list.files(path = paste0('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/traind_model/',mk),full.names = T),value=T) #lu he cardamom 排行第二2 是ensemble 因此就选这个吧
      model000 <- "/mnt/DataSpace/Projects/Canada_C/Canada_C_final/traind_model/litterinput/GBM_grid_1_AutoML_1_20240824_111040_model_10"
      GBM <- h2o.loadModel(model000)
      maps <- pred_fun0(model=GBM,t_data=test)
      #maps0 <- do.call(c, maps)
      #maps1 <- terra::clamp(maps, upper = 2.7634, values = T)
      maps2 <- round(exp(maps),4)
      names(maps2) <- 'GBM'
      terra::writeRaster(maps2,
                         paste0('/mnt/Fastrun/Canada_C/litter_climate_only2100/','climateonly2100_',i,'_',scenario,'_',time00[k],'.tif'),overwrite = TRUE)
      # }
      h2o.shutdown(prompt = FALSE)
    }
  }}

#h2o.shutdown(prompt = FALSE)
#39724s



# for(i in 1:9){
#   print(
#     paste0("sleep ", i * 100,"; Rscript --vanilla '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/006_Update_after_EC_CSSS/010_02_method3_litter_future_exp_climateonly2100.R' ",2*i -1," ",2*i," ",11791 + (i -1)*1000)
#   )
# }

