#### ca_bed creating
bedrock <- rast('/mnt/File0/DAAATAAA/bedrock_depth_from2017_tomhengl/BDTICM_M_250m_ll.tif')
ca_soc1m <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/EC_after_CSSS/canada_SOCs_1m_kgm2_expoential_relationship_SOC_BD_with_originaldata_20240730.tif')
nidaye <- rast('/mnt/DataSpace/Projects/Canada_C/bedrock_1km.tif')
ca_bed <- crop(nidaye,ext(ca_soc1m))
ca_bed1 <- resample(ca_bed,ca_soc1m)
ca_bed2 <- ifel(ca_bed1 > 300,300,ca_bed1)
ca_bed3 <- focal(ca_bed2,w=3,na.policy='only',function(x) mean(x,na.rm=T)) 
ca_bed30 <- ifel(ca_bed3 == 0,10^-5,ca_bed3) # make 0 as small as it could, because some ==0 will be NA
ca_bed4 <- mask(ca_bed30,mask = ca_soc1m)
plot(ca_bed4)
writeRaster(ca_bed4,'/mnt/File0/DAAATAAA/forfinalmap_1km/Canada_bedrock_depth0926.tif',overwrite=T)
plot(ifel(ca_bed4 == 0 ,NA,1))
t1 <- ifel(ca_bed4 == 100,1,0)
plot(t1)

ca_bed <- rast('/mnt/File0/DAAATAAA/forfinalmap_1km/Canada_bedrock_depth0926.tif')
# depth of each layer 
lyr5 <- ifel(ca_bed-5 >=0, 5, ca_bed)
lyr15 <- ifel(ca_bed-15 >=0, 15 - 5, ca_bed-5)
lyr30 <- ifel(ca_bed-30 >=0, 30 -15, ca_bed-15)
lyr60 <- ifel(ca_bed-60 >=0, 60 - 30, ca_bed-30)
lyr100 <- ifel(ca_bed-100 >=0, 100 - 60 , ca_bed-60)
lyr200 <- ifel(ca_bed-200 >=0, 200 - 100 , ca_bed-100)
lyr300 <- ifel(ca_bed-300 >=0, 300 - 200 , ca_bed-200)
#plot(lyr100)
dep0 <- c(lyr5,lyr15,lyr30,lyr60,lyr100,lyr200,lyr300)
dep1 <- ifel(dep0 <=0, 0,dep0)
#plot(dep1[[5]])
ca_bed_depths <- writeRaster(dep1,
                             '/mnt/File0/DAAATAAA/forfinalmap_1km/ca_bed_depth_of_each_of5layer.tif',overwrite=T)

