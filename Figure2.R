



###########################################        Figure 2  ###########################
library(terra)
library(sf)
library(scales)
library(fields)
library(spData)
library(dplyr)

data("world")
ca_outline <- st_read('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/2021digitallpr_000b21a_e/lpr_000b21a_e.shp')
ca_shap2 <- st_simplify(ca_outline,dTolerance = 3e3)
ca_shap3 <- project(vect(ca_shap2),'epsg:4326')
######################## M1
######################## M1
######################## M1
######################## M1


############################################   plot  ###################
############################################   plot  ################### 
############################################   plot  ###################
############################################   plot  ###################
############################################   plot  ###################

M1_soc <- rast('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/Processed_output/Visual_materials/M1_soc_two_senario_exp_climonly_2025_new.tif')
M2_soc <- rast('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/Processed_output/Visual_materials/M2_soc_two_senario_exp_climonly_2025_new.tif')
M3_soc <- rast('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/Processed_output/Visual_materials/M3_ML_initial_two_senario_exp_climonly_card_2025_new.tif')
M3_soc1 <- rast('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/Processed_output/Visual_materials/M3_q10_two_senario_exp_climonly_2025_new.tif')
M3_soc2 <- rast('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/Processed_output/Visual_materials/M3_hashi_two_senario_exp_climonly_2025_new.tif')

M1_cv <- rast('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/Processed_output/Visual_materials/M1_cv_eqarea_2025_new.tif')
M2_cv <- rast('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/Processed_output/Visual_materials/M2_cv_eqarea_2025_new.tif')
M3_cv <- rast('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/Processed_output/Visual_materials/M3_cv_eqarea_2025_new.tif')
landcover <- rast('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/LAND_cover/land_cover_equalarea_4visual.tif')
change_percent_mean <- rast('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/Processed_output/Visual_materials/three_methods_change_percent_to_hist_climonly_2025_new.tif')
quant2 <- quantile(na.omit(values(change_percent_mean[[1]])),seq(0,1,0.125),na.rm = T)
quant2 <- c(-28.1,-11.8,-6.89,-2.64,-1.77,0,1.38,3.9,7.8,40,88.22)
quant2 <- c(-28.1,-11.8,-6.89,-2.64,-1.77,0,1.38,3.9,7.8,40,88.22)
quant3 <- c(-28.1,-11.8,-6.89,0,3.9,7.8,88.22) # set less than 5% not change 
quant3 <- c(-28.1,-15,-5,0,5,15,88.22)

breaks_at3 <- seq(from = -28,to = 88,length.out = length(quant3)  )

######> color 
value0 <- na.omit(c(values(M1_soc[[1]]),values(M2_soc[[1]]),values(M3_soc[[1]])))


#2100 
#-58.7439690  -6.0573335  -2.9941118  -1.2042932  -0.1315464   0.6569971   1.5235373   3.6016523  33.9307632 
quant1 <- c(-60,-6.1,-3,-1.2,-0.13,0,0.66,1.52,3.6,34)
breaks_at <- seq(from = -60,to = 34,length.out = length(quant1)  )
colorTable <- designer.colors(length(quant1)-1  , c('#b2182b','#ef8a62','#fddbc7','#f7f7f7','#d1e5f0','#67a9cf','#2166ac') )
colorTable1 <- designer.colors(length(quant2)-1  , c('#b2182b','#ef8a62','#fddbc7','#f7f7f7','#d1e5f0','#67a9cf','#2166ac') )
colorTable2 <- designer.colors(length(quant3)-1  , c('#b2182b','#ef8a62','#fddbc7','#f7f7f7','#d1e5f0','#67a9cf','#2166ac') )

scales::show_col(colorTable)



M2_soc_0.5 <- resample(M2_soc,M1_soc)
M3_soc_0.5 <- resample(M3_soc,M1_soc)


M_all_0.5 <- c(M1_soc,M2_soc_0.5,M3_soc_0.5)
M_mean <- tapp(M_all_0.5,index = rep(1:2,3),function(x) mean(x,na.rm = T))
#
M_agreement <- ifel(M_all_0.5 > 0,1,ifel(M_all_0.5 < 0 ,-1,NA)) 
M_agreement1 <- tapp(M_agreement,index = rep(1:2,3),function(x) sum(x)) 
M_agreement2 <- ifel(M_agreement1 == 3,1,ifel(M_agreement1 == -3,-1,0))
landcover_0.5 <- rast('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/LAND_cover/land_cover_0.5.tif')
landcover_0.5_qa <- project(landcover_0.5,M_agreement2)
M_agreement3 <- mask(M_agreement2,landcover_0.5_qa)
par(mfrow = c(1,2))
plot(M_agreement2[[1]],xlim=c(3.5e6,9.2e6))
plot(M_agreement3[[1]],xlim=c(3.5e6,9.2e6))
C_agreement <- as.polygons(M_agreement3[[1]])

plot(C_agreement,col=c('red','green','blue'))
#

###############cv shadow 
M1_soc_cv_0.5 <- project(clamp(M1_cv,upper = 200,lower = -200),ca_outline)
M2_soc_cv_0.5 <- resample(x = clamp(M2_cv,upper = 100,lower = -100),y = M1_soc[[1]],method = 'med',threads = T)#
M3_soc_cv_0.5 <- resample(x = clamp(M3_cv,upper = 100,lower = -100),y = M1_soc[[1]],method = 'med',threads = T)#



q1 <- quantile(na.omit(abs(values(M1_soc_cv_0.5))),0.5)
q2 <- quantile(na.omit(abs(values(M2_soc_cv_0.5))),0.5)
q3 <- quantile(na.omit(abs(values(M3_soc_cv_0.5))),0.5)


M1_soc_cv_2 <- ifel(abs(M1_soc_cv_0.5) < q1,1,NA)
M2_soc_cv_2 <- ifel(abs(M2_soc_cv_0.5) < q2,1,NA)
M3_soc_cv_2 <- ifel(abs(M3_soc_cv_0.5) < q3,1,NA)

#1
ca_shap3 <- st_simplify(ca_shap2,dTolerance = 2e4)
ca2 <- world %>% filter(name_long =='Canada')
ca3 <- st_transform(x = ca2,  crs = crs(ca_shap2))


remove_small <- function(r,n=40){
  # 3. 
  r_patches <- patches(r, directions = 4)
  

  
  # 4.
  # 4a. 
  freq_table <- freq(r_patches)
  head(freq_table)
  # 
  # 
  small_ids <- freq_table$value[freq_table$count < n]
  
  # 
  r_patches[r_patches %in% small_ids] <- NA
  
  # 4b. 
  # 5. 
  p <- as.polygons(r_patches, dissolve = TRUE)
  
  # 
  p_sf <- st_as_sf(p)
  return(p_sf)
}
par(mfrow=c(1,3))

t1 <- remove_small(M1_soc_cv_2[[1]],n=60);plot(st_geometry(t1),col='red',border = 'white');plot(st_geometry(ca3),add = T)
t2 <- remove_small(M2_soc_cv_2[[1]],n=60);plot(st_geometry(t2),col='red',border = 'white');plot(st_geometry(ca3),add = T)
t3 <- remove_small(M3_soc_cv_2[[1]],n=60);plot(st_geometry(t3),col='red',border = 'white');plot(st_geometry(ca3),add = T)


#
#
#
#
#
library(cartography)
source('/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/006_Update_after_EC_CSSS/sf_patternfun.R')
#

######>plot 
path_fig <- '/mnt/DataSpace/Projects/Canada_C/Canada_C_final/Processed_output/Figures/'
pdf(file = paste0(path_fig,"Fig2_spatial_change_2025new.pdf"),  
    width = 10, 
    height = 7.5) 
par(mfrow = c(2,3))
tiff(filename = '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/figures/Figures_new/fig2_245_climonly.tiff',
     width = 10,
     height = 10,
     units = 'in',
     bg= 'white',
     compression = 'lzw',
     res = 600)
tiff(filename = '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/figures/Figures_new/fig2_585_climonly.pdf',
     width = 10,
     height = 10,
     bg= 'white'
)
#bottom, left, top, right
# 
dev.off()
par(mfrow=c(2,2))
ori_par <- par(no.readonly = T)
par(mar = c(5,4,3,2))
plot(mask(M2_soc[[1]],landcover),mar = c(0,0,0,0),
     breaks = quant1,col = colorTable,xlim=c(3.5e6,9.2e6),
     axes = 0,legend = F)
#if sf
plot(st_geometry(ca_shap2),add=T,lwd = 1,border = scales::alpha(colour = 'black',alpha = 0.5))

mtext(paste0(letters[1],'','\n     Ensemble machine \n     learning'),
      side = 3, 
      font=2, cex = 1.2,
      line = -1.5,
      adj = 0.01 ) 

####
par(mar = c(5,4,3,2))
plot(mask(M3_soc[[1]],landcover),mar = c(0,0,0,0),
     breaks = quant1,col = colorTable,xlim=c(3.5e6,9.2e6),
     axes = 0,legend = F)
#if sf
plot(st_geometry(ca_shap2),add=T,lwd = 1,border = scales::alpha(colour = 'black',alpha = 0.5))

mtext(paste0(letters[2],'','\n     Carbon balance \n     approach'), 
      side = 3, 
      font=2, cex = 1.2,
      line = -1.5,
      adj = 0.01 ) 


#### 3
par(mar = c(5,4,3,2))
plot(M1_soc[[1]],mar = c(0,0,0,0),
     breaks = quant1,col = colorTable,xlim=c(3.5e6,9.2e6),
     axes = 0,legend = F)
#if sf
plot(st_geometry(ca_shap2),add=T,lwd = 1,border = scales::alpha(colour = 'black',alpha = 0.5))

mtext(
  paste0(letters[3],'','\n     Earth system \n     model'), 
  side = 3, 
  font=2, cex = 1.2,
  line = -1.5,
  adj = 0.01 ) 
###legend 
par(mar=c(1,1,1,2))
imagePlot(zlim = c(-60,40),
          #zlim = range(value0),
          legend.only = TRUE, 
          legend.shrink = 0.7,
          smallplot = c(0.10, 0.58, 0.1, 0.13),horizontal = T,
          legend.width = 0.8, legend.mar = 2, 
          col = colorTable,
          breaks = breaks_at, 
          axis.args = list(at= breaks_at,
                           tck=-0.1,
                           mgp=c(0,0.4,0),
                           labels = round(quant1,1), cex.axis=0.9),
          legend.args=list( text=expression(~x~10^-2~"kg C m"^"-2"~yr^-1),col="black", 
                            cex=0.6, side=3, line=0.08, adj=0.35)
)
## for change percent 
par(mar = c(0,0,0,0)) 
plot(change_percent_mean[[1]],mar = c(0,0,0,0),
     breaks = quant3,col = colorTable2,xlim=c(3.5e6,9.2e6), 
     axes = 0,legend = F)
#if sf
plot(st_geometry(ca_shap2),add=T,lwd = 1,border = scales::alpha(colour = 'black',alpha = 0.5))

mtext(paste0(letters[4],'','\n     Mean change \n     percentage'), side = 3, 
      font=2, cex = 1.2,
      line = -4.5,
      adj = 0.01 ) 
patternLayer(subset(st_as_sf(C_agreement),X1==-1),xlim=c(3.5e6,9.2e6), "right2left", density=4, lwd =1.6, col= scales::alpha('black',0.9), add=T)
patternLayer(subset(st_as_sf(C_agreement),X1==1), xlim=c(3.5e6,9.2e6), "vertical", density=4, lwd =1.6, col= scales::alpha('black',0.9), add=T)
legendPattern(
  title.txt = "",
  pos = c(4e6,6e5),
  categ = c("SOC loss agreement","SOC gain agreement"),
  patterns = c("right2left","vertical"),
  pch = 22,
  cex = 2,
  values.cex= 1.,
  frame = F,
  lwd=0.5,
  add=T
)
par(mar=c(1,1,1,2))
imagePlot(zlim = c(-28,88),
          legend.only = TRUE, 
          legend.shrink = 0.7,
          smallplot = c(0.82, 0.87, 0.46, 0.83),horizontal = F, 
          legend.width = 0.8, legend.mar = 2, 
          col = colorTable2,
          breaks = breaks_at3, 
          axis.args = list(at= breaks_at3,
                           tck=-0.1,
                           mgp=c(0,0.4,0),
                           labels = round(quant3,1), cex.axis=0.9),
          legend.args=list( text = '  %',
                            cex=0.8, side=3, line=0.08, adj=0.35)
)
# # # # # # 
# # # # # # 
# # # # # # insert 

par(mar=c(0,0,0,0.5),fig = c(x1 = 0.32, x2 = 0.5,x3 = 0.75,x4 = 1), bg = NA,new = TRUE) ##note the order, bg must be set before new = T. 

plot(st_geometry(t2),col='blue',#'#01b4f2',#'#2988b8' '#17e4db'
     border = 'white');plot(st_geometry(ca3),add=T,border='grey')
par(mar=c(0,0,0,0.5),fig = c(x1 = 0.32 + 0.5, x2 = 0.5 + 0.5,x3 = 0.75,x4 = 1),bg = NA, new = TRUE)

plot(st_geometry(t3),col='blue',#'#01b4f2',#'#2988b8' '#17e4db'
     border = 'white');plot(st_geometry(ca3),add=T,border='grey')

par(mar=c(0,0,0,0.5),fig = c(x1 = 0.32 , x2 = 0.5 ,x3 = 0.75 -0.5,x4 = 1-0.5),bg = NA, new = TRUE)

plot(st_geometry(t1),col='blue',#'#01b4f2',#'#2988b8' '#17e4db'
     border = 'white');plot(st_geometry(ca3),add=T,border='grey')

dev.off()
par(ori_par)


