#20240826 
#----------------------------- 
#from #/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/003_code_lehigh/09_01_turnover_after_chat_add_sasha_peat_BD.R
#one model for each pixel 
h2_arg <- commandArgs(trailingOnly = T)
stopifnot(length(h2_arg) >0 )
model_no1 <- as.numeric(h2_arg[1])
model_no2 <- as.numeric(h2_arg[2])
###############################################
unixtools::set.tempdir('/mnt/temp4r/temp')
library(terra)
library(dplyr)
library(tidyr)
library(h2o)
#library(readr)
library(data.table) #fread fwrite

# future turnover rate after 
#Rscript --vanilla /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/002_code/09_turnover_k_Climateonly2300.R 1 2 21217 stell
landcover <- rast('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/LAND_cover/land_cover_mask_lt15_keep_origin_1km_w84_1214.tif') #> 15  #1214 with new extend as same as bd

# ca_soc_stock <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/base_soc_kgm2_WGS84_1214.tif')
# ca_soc1m <- sum(ca_soc_stock[[1:5]])
ca_soc1m <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/canada_SOCs_1m_kgm2_expoential_relationship_SOC_BD_with_originaldata_20240730.tif')


#ca_stock84_0 <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/base_socs_WGS84_1214.tif') # unit = 1/10 gram,every ten gram, so transfer to every 1000g equals to divided by 100 , transfer to gram need *10
ca_pure_area84 <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/ca_cmip6_area_WGS84.tif') #20231214 new created

peat_ca <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/peatland_area.tif')
peat_area <- ca_pure_area84 * peat_ca  # peat pixel area, and then multiple ca_stock84_0, then sum together

#all
# models <- c('EC-Earth3-Veg-LR','IPSL-CM6A-LR','ACCESS-ESM1-5','CanESM5','MPI-ESM1-2-HR','MPI-ESM1-2-LR')
#expt ='all'
# #climate 2100
# expt = 'plant_climateonly2100'
 expt = 'climateonly2100'
 models <- c("ACCESS-CM2","CanESM5","CMCC-ESM2","EC-Earth3-Veg","FIO-ESM-2-0","INM-CM4-8","IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0",  "ACCESS-ESM1-5", "EC-Earth3-Veg-LR", "INM-CM5-0" ,"BCC-CSM2-MR","CanESM5-1","CAS-ESM2-0", "EC-Earth3-CC")
 # EC is relation, flittersoil_model is actrually models we have , from flittersoil_model_use in  /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/006_Update_after_EC_CSSS/002_future_fveglitter_after_EC.R
 # fveg_mods <- c('BCC-CSM2-MR','CanESM5','CanESM5-1','CESM2-WACCM','CMCC-CM2-SR5','CMCC-ESM2','EC-Earth3-CC','EC-Earth3-Veg','EC-Earth3-Veg-LR','IPSL-CM6A-LR','KIOST-ESM','NorESM2-LM','NorESM2-MM','TaiESM1')
 # models <- intersect(models,fveg_mods)
 # models_p <- c("EC-Earth3-Veg-LR", "IPSL-CM6A-LR", "ACCESS-ESM1-5", "BCC-CSM2-MR",
#               "CanESM5", "CanESM5-1", "CAS-ESM2-0", "CESM2-WACCM", 
#               "CMCC-CM2-SR5", "EC-Earth3-CC", "INM-CM4-8", "INM-CM5-0",
#               "TaiESM1", "EC-Earth3-Veg", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR")
# models <- intersect(models,models_p)
# #climate 2300 
# expt = 'climateonly2300_ave_npp'
# expt = 'climateonly2300_clitterSoil'
# models <- c("CanESM5", "EC-Earth3-Veg", "IPSL-CM6A-LR", "MRI-ESM2-0", "ACCESS-ESM1-5")
# #all with no plants
# expt = 'all_noplant'
# models <- c('EC-Earth3-Veg-LR','IPSL-CM6A-LR','ACCESS-ESM1-5','CanESM5','MPI-ESM1-2-HR','MPI-ESM1-2-LR')
# #
# expt = 'plant_climate2300'
# models <- c("CanESM5", "IPSL-CM6A-LR", "ACCESS-ESM1-5")
# 
# f1 <- list.files(pattern = 'tif',path = '/mnt/Fastrun/Data_pool/01cmip6_processed_output',full.names = T)

#(f_k <- list.files('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/all_2100',full.names = T))

# f_k <- list.files('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/plant_climateonly_2100',full.names = T)
# 
# f_k <- list.files('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/Climateonly_2300',full.names = T)
 f_k <- list.files('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/Climateonly_2100',full.names = T)
# f_k <- list.files('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/plant_climateonly_585_2300',full.names = T)
# f_k <- list.files('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/all2100_noplant',full.names = T)

#
#f_npp0 <- list.files('/mnt/DataSpace/Data_pool/01cmip6_processed_output/turnover_npp/',full.names = T)
#f_npp0 <- list.files(pattern = 'tif',path = '/mnt/DataSpace/Data_pool/01cmip6_processed_output_update_npp',full.names = T)
#(f_npp0 <- list.files(pattern = 'tif',path = '/mnt/DataSpace/Data_pool/01cmip6_processed_output_update_GPP_NPP_from_fluxcom_20240423',full.names = T))
#20240423


#20240826
 #bias correct for the beginning litter0 = rh0, and keep the ratio constant afterward.
 litter_base <- rast('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/litterinput_base_litterinput_1981to2000_GBM_20240824.tif')
 #
 
 lu_k <- rast('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/lu_base_turnover_rate1981to2000_GBM_20240730.tif')
 cardamom_k <- rast('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/cardamom_base_turnover_rate1981to2000_GBM_20240730.tif')
 lu_rh0 <- ca_soc1m * lu_k/10000
 cardamom_rh0 <- ca_soc1m * cardamom_k/10000
 # input_output_base <- c(litter_base,lu_rh0,cardamom_rh0)
 # names(input_output_base) <- c('litter0','lu_rh0','card_rh0')
 # plot(input_output_base,nr=1,pax = list(cex.axis=1.5),plg = list( x= 'bottom'))
 # title(main = 'Base input and outputs')
 # title(sub = 'TMD')
 lu_ratio <- clamp(lu_rh0 *1000/litter_base,upper = 5)
 #plot(lu_ratio)
 #title(main='lu_to_litter')
 card_ratio <- clamp(cardamom_rh0 *1000/litter_base,upper=5)
 #plot(card_ratio)
 #title(main = 'card_to_litter')
 
 
#(f_veglit <- list.files(pattern = 'tif', pat = '/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/litter_all_pet_and_Temp_P',full.names = T))
(f_veglit <- list.files(pattern = 'tif', pat = '/mnt/Fastrun/Canada_C/litter_climate_only2100/',full.names = T))


#20240731
# (f_veglit <- list.files(pattern = 'tif', pat = '/mnt/DataSpace/Data_pool/01cmip6_processed_output_EC_fVegLitter',full.names = T))
# as we assume fveglitter == Rh, there are two Rh = c(lu + cardamom)
#
#ca_covas <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/canada_all_covs_global1214.tif')
#npp1 <- ca_covas[[3]]


#####

#ratio_file <- list.files(pattern = 'tif','/mnt/DataSpace/Data_pool/NPP_litter_ratio/',full.names = T)
# ratio_file <- list.files(pattern = 'tif','/mnt/DataSpace/Data_pool/fveglitter_gpp_ratio_from_fluxcom/',full.names = T)


# model_clitter_soil <-  c("CanESM5", "CanESM5-1","CESM2-WACCM","CMCC-CM2-SR5","EC-Earth3-CC" ,"EC-Earth3-Veg", "EC-Earth3-Veg-LR","IPSL-CM6A-LR",  "TaiESM1")
# models <- intersect(models,model_clitter_soil)
# models <- c("ACCESS-ESM1-5","CanESM5","EC-Earth3-Veg-LR", "IPSL-CM6A-LR")#20240731
# #EC is relation, flittersoil_model is actrually models we have , from flittersoil_model_use in  /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/006_Update_after_EC_CSSS/002_future_fveglitter_after_EC.R
# fveg_mods <- c('BCC-CSM2-MR','CanESM5','CanESM5-1','CESM2-WACCM','CMCC-CM2-SR5','CMCC-ESM2','EC-Earth3-CC','EC-Earth3-Veg','EC-Earth3-Veg-LR','IPSL-CM6A-LR','KIOST-ESM','NorESM2-LM','NorESM2-MM','TaiESM1')
# models <- intersect(models,fveg_mods)

k_type <- c('lu','cardamom')
SOC_ca <- list()
input_ca <- list()
output_ca <- list()
for (km in k_type) {
  k_base <- rast(paste0('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/',km,'_base_turnover_rate1981to2000_GBM_20240730.tif'))
  soc_s <- list()
  # c_input <- list()
  # c_output <- list()
  socs0 <- vector()
  total_peat_storage84 <- vector()
  total_storage84_lc <- list()
  #SOC_s[[1]] <- ca_soc1m + (npp0 / 10000) * 20 - (ca_soc1m * k_base / 10000) * 20 # 20 represent 20 year and 10000 to get back to original units
  #
  for (i in models[model_no1:model_no2]) {
    for (scenario in c('ssp245','ssp585')) { #change c('ssp245','ssp585') to c('ssp585) when 2300
      name0 <- paste0(km,'_',i,'_',scenario)
      cat(name0)
      #ca_soc1m <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/canada_SOCs_1m_kgm2_expoential_relationship_SOC_BD_with_originaldata_20240329.tif')
      ca_soc1m <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/canada_SOCs_1m_kgm2_expoential_relationship_SOC_BD_with_originaldata_20240730.tif')
      # from /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/006_Update_after_EC_CSSS/001_1_h2o_automl_build_turnoverrate_model.R
      
      (C_inp <- grep(pattern = paste0('climateonly2100_',i,'_',scenario),f_veglit,value=T))
      if(km == 'cardamom'){
        ratio00 <- card_ratio
      } else {ratio00 <- lu_ratio}
      
      (f_k0 <- grep(glob2rx(paste0('*',km,'*',i, '_', scenario,'*')), f_k, value = T))
      
      for(k in 1:10){  # all 10 period to 2100
        C_inp_use <- rast(C_inp[k])
        C_out_use <- rast(f_k0[k])
        for(j in 1:10){ # each period 10 years 
          ca_soc1m <- ca_soc1m * (1- C_out_use / 10000) + C_inp_use/1000 # g to kg 
          ca_som2m <- mask(ca_soc1m,landcover) * ca_pure_area84
          ca_som3m <- clamp(ca_som2m,lower=0,value=TRUE)
          socs0[(k-1)*10 + j] <- global(ca_som3m,function(x) sum(x,na.rm = T))/10^12
          #landcover
          #> land cover /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/006_Update_after_EC_CSSS/11_02_method2_pureML_calculation.R
          total_C84_lc <- terra::zonal(ca_som3m,landcover,fun=function(x) sum(x,na.rm=TRUE))
          total_storage84_lc[[(k-1)*10 + j]] <- total_C84_lc/10^14
          #peatland from L68 /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/006_Update_after_EC_CSSS/11_02_method2_pureML_calculation.R
          ca_peat_stock84 <- ca_soc1m * peat_area
          ca_peat_stock84_lc <- mask(ca_peat_stock84,landcover)# macstudio notes 
          #TOTAL STORAGE
          total_peat_C84 <- global(ca_peat_stock84_lc,fun= function(x) sum(x,na.rm=TRUE))
          ####
          total_peat_storage84[(k-1)*10 + j] <- total_peat_C84 * 10^7/10^15 /10^6
          #
          triger0 <- 10 * (k-1) + j
          if(triger0 %in% c(60,100)){  # if you want to output 20 40 60 80 then replace the values. 
            ind1 <- which(c(60,100) == triger0)
            names(ca_soc1m) <- 2000 + triger0
            soc_s[[ind1]] <- ca_soc1m
          }
        }
      }
      # npp0 <- rast(grep(paste0('gpp_',i, '_', scenario), f_npp0, value = T)[1])
      # npp0 <- clamp(npp0,lower = 0 )
      # #npp0 <- rast(grep(paste0(i, '_', scenario), f_npp0, value = T))
      # ratio0 <- rast(grep(paste0(i, '_', scenario), ratio_file, value = T))
      # ratio0 <- clamp(ratio0,lower = 0 , upper = 1, values =T )
      # #SOC_s[[1]] <- ca_soc1m + (npp0 * ratio0[[1]] / 10000) * 20 - (ca_soc1m * k_base / 10000) * 20 
      #  name0 <- paste0(km,'_',i,'_',scenario)
      #  for( w in 1:20){
      #    ca_soc1m <- ca_soc1m * (1 - k_base/10000) + (npp0 * ratio0[[1]] / 1000) # g to kg 
      #    ca_som2m <- mask(ca_soc1m,landcover) * ca_pure_area84
      #    ca_som3m <- clamp(ca_som2m,lower=0,value=TRUE)
      #    socs0[w] <- global(ca_som3m,function(x) sum(x,na.rm = T))/10^12
      #  }
      # # writeRaster(ca_soc1m,paste0('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/SOC_stock_k_cLitterSoil/',name0,'_',expt,'_constraintbulkdensity2020.tif'),overwrite=TRUE)
      #  #c_input[[1]] <- npp0 * ratio0[[1]] / 10000
      #  #c_output[[1]] <- ca_soc1m * k_base / 10000 # because we need output spatial pattern , here we outputs each 20 year output 
      #  # f_plant <- grep(paste0(i, '_', scenario), f1, value = T)
      #  f_plant <- grep(paste0('gpp_',i, '_', scenario), f_npp0, value = T)[2:15]
      #  f_k0 <- grep(glob2rx(paste0('*',km,'*',i, '_', scenario,'*')), f_k, value = T)
      #  for (k in 1:4) {  #change 4 to 14 when 2300
      #    # f_ind <- f_plant[seq(k, length(f_plant), length(f_plant) / 3)] #time period
      #    # npp <- rast(f_ind[2]) 
      #    if(k>=2){
      #      #writeRaster(ca_soc1m,paste0('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/SOC_stock_k_cLitterSoil/',name0,'_',expt,'_constraintbulkdensity',seq(2020,2100,20)[k],'.tif'),overwrite=TRUE)
      #    }
      #    npp <- rast(f_plant[k])
      #    npp <- clamp(npp,lower = 0)
      #    k_current <- rast(f_k0[k])
      #    k_current0 <- clamp(k_current,lower=0,values=T)
      #    #
      #    for(j in 1:20){
      #    #SOC_s[[k + 1]] <- SOC_s[[k]] + (npp * ratio0[[k+1]] / 10000) * 20 - (SOC_s[[k]] * k_current0 / 10000) * 20
      #    ca_soc1m <- ca_soc1m * (1- k_current0 / 10000) + (npp * ratio0[[k + 1]] / 1000) # g to kg 
      #    ca_som2m <- mask(ca_soc1m,landcover) * ca_pure_area84
      #    ca_som3m <- clamp(ca_som2m,lower=0,value=TRUE)
      #    socs0[k*20 + j] <- global(ca_som3m,function(x) sum(x,na.rm = T))/10^12
      #    }
      #    #c_input[[k+1]] <- npp * ratio0[[k+1]]/10000
      #    #c_output[[k+1]] <- ca_soc1m * k_current0 / 10000 #
      #  }
      
      
      #writeRaster(ca_soc1m,paste0('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/SOC_stock_k_cLitterSoil/',name0,'_',expt,'_constraintbulkdensity2100.tif'),overwrite=TRUE)
      
      
      #input_all <- do.call(c,c_input) * ca_pure_area84
      #input_all <- mask(input_all,landcover) # >15 not 
      #input_all0 <- clamp(input_all,lower=0,values=TRUE)
      #input_ca[[name0]] <- global(input_all0,function(x) sum(x,na.rm = TRUE))/10^12 # because of the units is kg/m2
      #writeRaster(do.call(c,c_input),filename = paste0('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/SOC_stock_k_cLitterSoil/',name0,'_',expt,'_constraintbulkdensity_input.tif'),overwrite=TRUE)
      
      #output_all <- do.call(c,c_output) * ca_pure_area84
      #output_all <- mask(output_all,landcover) # >15 not 
      #output_all0 <- clamp(output_all,lower=0,values=T)
      #output_ca[[name0]] <- global(output_all0,function(x) sum(x,na.rm = TRUE))/10^12 # because of the units is kg/m2
      #writeRaster(do.call(c,c_output),filename = paste0('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/SOC_stock_k_cLitterSoil/',name0,'_',expt,'_constraintbulkdensity_output.tif'),overwrite=TRUE)
      
      #SOC_all <- do.call(c,SOC_s) * ca_pure_area84
      #SOC_all <- mask(SOC_all,landcover) # >15 not 
      #SOC_all0 <- clamp(SOC_all,lower=0,values=TRUE)
      #SOC_ca[[name0]] <- global(SOC_all0,function(x) sum(x,na.rm = TRUE))/10^12 # because of the units is kg/m2
      #writeRaster(do.call(c,SOC_s),filename = paste0('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/SOC_stock_k_cLitterSoil/',name0,'_',expt,'_constraintbulkdensity.tif'),overwrite=TRUE)
      #writeRaster(do.call(c,SOC_s),filename = paste0('/mnt/DataSpace/Projects/Canada_C/turnover_constrain/model_outputs/SOC_stock_k_cLitterSoil/',name0,'_',expt,'_constraintbulkdensity.tif'),overwrite=TRUE)
      writeRaster(do.call(c,soc_s),filename = paste0('/mnt/DataSpace/Projects/Canada_C/turnover_constrain/model_outputs/SOC_stock_sasha_2box_balance/',name0,'_',expt,'_based_on_sasa_two_box_balance_2060_2100status.tif'),overwrite=TRUE)
      saveRDS(list(socs0,total_peat_storage84,total_storage84_lc),paste0('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs/',name0,'_SOC_dynamic_fluxcom_gpp_climateonly2100_annual20240826_newlitter.rds'))
    }
  }
}


# f1 <- list.files(pattern = glob2rx('*dynamic*rds'),'/mnt/DataSpace/Projects/Canada_C/Canada_C_final/',full.names = T)
# plot(unlist(readRDS(f1[1])))

