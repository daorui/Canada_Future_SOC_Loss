




##############################     Figure 1       ##############################
# currently used for submit
# 
# #。      scenario     time upper lower upper1 lower1 anomaly
# 1 Cardamom SSP245    2100 -0.397 -8.83   -1.41  -9.13    -4.61
# 2 Cardamom SSP585    2100  1.72  -9.57   -1.19 -10.1     -3.93
# 
# #ML 
# 
# SSP245    2100  2.71 -15.9  0.592  -14.2   -6.61
# SSP585    2100 -5.62 -17.8 -8.39   -17.5  -11.7 
# 
# #
# 1 SSP245    2100  8.27 -6.68     218.   8.23  -7.10   0.793
# 2 SSP585    2100  8.18 -9.02     217.   6.73  -8.38  -0.421
# 
# ##
# M2 peat 
# scenario  time  upper lower upper1 lower1 anomaly
# 1 SSP245    2100 -0.347 -4.87 -0.887  -4.28   -2.61
# 2 SSP585    2100 -2.50  -6.46 -3.30   -6.07   -4.48
# 
# M3 peat 
# Rh       scenario  time upper lower upper1 lower1 anomaly
# 
# 1 Cardamom SSP245    2100 -5.05 -8.02  -5.41  -8.10   -6.53
# 2 Cardamom SSP585    2100 -5.07 -8.30  -5.95  -8.62   -6.69
# ca_pure_area84 <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/ca_cmip6_area_WGS84.tif') #20231214 new created
# 
# peat_ca <- rast('/mnt/DataSpace/Projects/Canada_C/processed_output/peatland_area.tif')
# peat_area <- ca_pure_area84 * peat_ca  # peat pixel area, and then multiple ca_stock84_0, then sum together
# p_m2 <- global(peat_area,function(x) sum(x,na.rm = T))$global/100
# 2.61*10^15/p_m2  =2.305 kg m-2
# 4.48*10^15/p_m2  =3.956 kg m-2
# 6.53*10^15/p_m2  =5.766 kg m-2
# 6.69*10^15/p_m2  =5.907 kg m-2
# CMIP6 
library(dplyr)
library(ggplot2)
library(terra)
# national 
cmip6_soc <- readRDS('/mnt/DataSpace/Data_pool/CMIP6_cSoil_trend/anomaly_in_future_litter_csoil_totalsoc.rds')
# peatland 
# cmip6_soc <- readRDS('/mnt/DataSpace/Data_pool/CMIP6_cSoil_trend/anomaly_in_future_litter_csoil_totalsoc_peatland.rds')

df_commom <- cmip6_soc[[5]]
df_commom_ave <- cmip6_soc[[6]]
# model mean
df_commom_family <- df_commom %>% mutate(family = substr(model,1,4)) %>% 
  filter(time >= 2000 & time <= 2100) %>% 
  group_by(family,scenario,time) %>% 
  summarize(Totalsoc=mean(Totalsoc),anomaly=mean(anomaly)) 
#
df_common_box <- df_commom_family %>% filter(time == 2100) %>% mutate(scen = case_when(scenario == 'SSP245' ~2110,
                                                                                       scenario == 'SSP585' ~2120))

# family weighted mean 
df_commom_ave_family_ave <- df_commom_family %>% 
  group_by(scenario,time) %>% 
  summarize(upper = mean(anomaly) + 1.645 * sd(anomaly), 
            lower = mean(anomaly) - 1.645 * sd(anomaly),
            Totalsoc=mean(Totalsoc),upper1 = max(anomaly),
            lower1 = min(anomaly),anomaly=mean(anomaly))
#
model_f <- unique(df_commom_family$family)
model0 <- c('ACCESS-ESM1-5','CanESM5(2)','CESM2-WACCM','CMCC(2)','EC-Earth3(3)','IPSL-CM6A-LR','NorESM2-MM','TaiESM1')
model1 <- unlist(lapply(model_f,function(x) grep(x,model0,value=T)))
base_2000 <- df_commom %>% mutate(family = substr(model,1,4)) %>% 
  group_by(family,scenario) %>% 
  summarize(base = base[1]) %>% group_by(scenario) %>% reframe(base0 = mean(base),sd = sd(base))

m3_models <- c('ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5(2)','EC-Earth3(3)','INM-CM(2)','IPSL-CM6A-LR')
m2_models <- c("ACCESS-ESM1-5", "BCC-CSM2-MR"  , "CanESM5(2)"   , "EC-Earth3(3)" , "INM-CM(2)" ,    "IPSL-CM6A-LR" )
m1_models <-  c('ACCESS-ESM1-5','CanESM5(2)','CESM2-WACCM','CMCC(2)','EC-Earth3(3)','IPSL-CM6A-LR','NorESM2-MM','TaiESM1')
a_models <- unique(c(m1_models,m2_models,m3_models))
#
shape_value0 <- c("ACCESS-ESM1-5", "CanESM5(2)", "CESM2-WACCM","CMCC(2)","EC-Earth3(3)","IPSL-CM6A-LR",
                  "NorESM2-MM","TaiESM1","BCC-CSM2-MR","INM-CM(2)")

shape_value1 <- c("ACCE"=1, "CanE"=2 , "CESM"=3,"CMCC"=4,"EC-E"=5,  "IPSL"=6,
                  "NorE"=7,"TaiE"=8,"BCC-"=9,"INM-" =10)

library(ggplot2)
library(extrafont)
library(RColorBrewer)
esm_cols <- colorRampPalette(brewer.pal(12,'Paired'))(length(model0))
names(esm_cols) <- model_f
library(scales)
# show_col(esm_cols)

theme0 <- theme(axis.text = element_text(face = 'bold',family = 'Times',size = 6 * 2.5),
                axis.text.x.bottom = element_text(face = 'bold',margin=margin(t = 10,0,0,0)),
                axis.title.y.left = element_text(family = 'Times',size = 8 * 2.5,margin=margin(0,10,0,0)),
                axis.title.x.bottom = element_text(family = 'Times',size = 8 * 2.5,margin=margin(10,0,0,0)),
                axis.line.x.bottom = element_line(size=1,color = "black"),
                axis.line.y.left = element_line(size=1,color = "black"),
                axis.line.x.top = element_line(size=1,color = 'black'),
                axis.line.y.right = element_line(size=1,color = 'black'),
                axis.ticks.x.top = element_blank(),
                axis.ticks.y.right = element_blank(),
                axis.text.x.top = element_blank(),
                axis.text.y.right = element_blank(),
                plot.background = element_blank(),
                rect = element_rect(fill = "transparent"),
                
                legend.position = c(0.05,0.15),
                legend.justification = 'left',
                
                legend.box = "horizontal", 
                legend.box.just = "top", 
                legend.text = element_text(family = 'Times',size = 4*2.5),
                legend.title = element_text(family = 'Times',size = 5 * 2.5),
                legend.spacing.y = unit(0, "cm"),
                title = element_text(family = 'Times',size = 18)
)


ylim_min1 <- -18
ylim_max1 <- 15

# Calculate positions as a ratio of the y-range

y_pos_01 <- ylim_min1 + 0.96 * (ylim_max1 - ylim_min1)  # 50% height
y_pos_11 <- ylim_min1 + 0.88 * (ylim_max1 - ylim_min1)  # 70% height
y_pos_21 <- ylim_min1 + 0.80 * (ylim_max1 - ylim_min1)  # 50% height



library(patchwork)


p1 <-  
  ggplot() +
  geom_ribbon(data=df_commom_ave_family_ave,
              aes(x= time, ymin = lower, ymax = upper,
                  group= scenario,fill=scenario), alpha = 0.6) +
  geom_segment(aes(x = 1995, xend = 2100, y = 0, yend = 0),
               linetype = "dashed", color = "#454545", lwd = 0.7) +
  annotate(geom = 'text',
           x=rep(2000,3)[1],
           y=10,
           label = ' Earth system model',
            family = 'Times',hjust = 0,
           size = 6,na.rm = T) + 
  scale_y_continuous(breaks=c(-20,-10,0,10),limits = c(-18,12)) +
  scale_x_continuous(breaks=c(seq(2000,2100,25)),
                     limits = c(1995, 2130), expand = c(0, 0)) + 
  labs(y= 'SOC storage change (Pg)',x= 'Year') +
  geom_line(data = df_commom_ave_family_ave,
            aes(x=time,y=anomaly,group=scenario,col=scenario),
            lwd=2) +
  scale_fill_manual(name='Scenario',values = c('SSP245'= '#92c5de','SSP585'= '#fcd1c5'),
                    labels = c('SSP245' = 'SSP2-4.5','SSP585' ='SSP5-8.5'),
                      ) +
  scale_color_manual(name = 'Scenario',values = c('SSP245' = '#70a0cd','SSP585'='#990002'),
                     labels = c('SSP245' = 'SSP2-4.5','SSP585' ='SSP5-8.5')) +
  scale_shape_manual(name = 'ESMs',
                     values = shape_value1,
                     labels = shape_value0,
                     breaks = names(shape_value1),
                     
                     drop= FALSE
  ) + 
  scale_linetype_discrete(guide = guide_legend(order=1,title='Scenario',
                                               show.legend = F,
                                               theme = theme(
                                                 legend.direction = "vertical",
                                                 legend.text = element_text(hjust = 0, vjust = 0.5, angle = 0)
                                               ))) +
  guides(shape = guide_legend(ncol = 2,
                              override.aes = list(size = 2), # shape size 
                              theme = theme(legend.title.position = "top",
                                            legend.byrow = TRUE))) + 
  guides(col = guide_legend(order=1)) +
  guides(fill = guide_legend(order= 1, 
                             theme = theme(
                               legend.direction = "vertical",
                               legend.title.position = "top",
                               legend.text.position = "right",
                               legend.text = element_text(hjust = 0, vjust = 0.5, angle = 0)
                             ))) +
  theme_classic() + 
  ggtitle(label = 'M1: Earth System Models') +
  theme0 +
  theme(axis.line.x = element_blank())


p10 <- 
  p1 + 
  geom_boxplot(data = df_common_box,aes(x = scen, y = anomaly, fill = scenario),
               outlier.shape = NA, alpha = 0.7,size = 0.3 ) +
  geom_jitter(data= df_common_box,  
              aes(x = scen, y = anomaly, 
                  shape =factor(family,levels = names(shape_value1))), 
              width = 0.5,   
              alpha = 1, 
              stroke = 0.5,
              size = 3
  ) 
p10


#### M2

soc_stock_use <- readRDS('/mnt/DataSpace/Projects/Canada_C/Canada_C_final/Processed_output/method2_pureml_change_pattern_landcover/temporal_trend_2025_national.rds')


#before 
# peat = 0
# basev = ifelse(peat==1,28.78664,152.8529)
#now 
peat = 0
basev = ifelse(peat==1,36.37054,143.988)

# for NWT northest terroteries
# NT = 0
# baseev = ifelse(NT==1 & peat == 1,5.291845,23.98221)

exp0 <- 'All'
exp0 <- 'Clim2100'

m2_df_family <- soc_stock_use %>% filter(exp == exp0 ) %>%
  mutate(family = substr(model,1,4),anomaly = stock_1m - basev) %>% #160.753 before bedrock 
  filter(time >= 2000 & time <= 2100) %>% 
  group_by(family,scenario,time) %>% 
  summarize(anomaly = mean(anomaly),time = time_end[1]) %>% 
  group_modify(~ add_row(.,time = 2000,anomaly = 0,.before = 1)) %>% # add base 2000
  mutate(scenario = toupper(scenario))

#
m2_df_box <- m2_df_family %>% filter(time == 2100) %>% mutate(scen = case_when(scenario == 'SSP245' ~2110,
                                                                               scenario == 'SSP585' ~2120))


m2_df_ave_family_ave <- m2_df_family %>% 
  group_by(scenario,time) %>% 
  summarize(upper = mean(anomaly) + 1.645 * sd(anomaly),
            lower = mean(anomaly) - 1.645 * sd(anomaly),
            upper1 = max(anomaly),
            lower1 = min(anomaly),anomaly=mean(anomaly)) %>% mutate(scenario = toupper(scenario))

model_f <- unique(m2_df_family$family)
model_all <- c('ACCESS-ESM1-5','CanESM5','EC-Earth3-Veg-LR','IPSL-CM6A-LR')
model_clim2100 <- c('ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5(2)','EC-Earth3(3)','INM-CM(2)','IPSL-CM6A-LR')
shape_value =c()
get_model <- function(y) {unlist(lapply(model_f,function(x) grep(x,y,value=T)))}
if(exp0 == 'All') {model1 <- get_model(model_all)} else {model1 <- get_model(model_clim2100)}
if(exp0 == 'All') {delta245 = -1.08;delta585 =-1.4 } else {delta245 = -0.15;delta585 = -0.293}
# all 245 : -10.8   +   585: -13.97  # peat ssp245 : -3.4873 + ssp585: - 4.966 

#
esm_cols <- colorRampPalette(brewer.pal(12,'Paired'))(length(model1))
names(esm_cols) <- model_f
library(scales)
# show_col(esm_cols)
#
ylim_min2 <- -20
ylim_max2 <- 18

# Calculate positions as a ratio of the y-range
y_pos_12 <- ylim_min2 + 0.88 * (ylim_max2 - ylim_min2)  
y_pos_22 <- ylim_min2 + 0.96 * (ylim_max2 - ylim_min2)  
p2 <-  
  ggplot() +
   geom_ribbon(data=m2_df_ave_family_ave,
              aes(x= time, ymin = lower, ymax = upper,
                  group= scenario,fill=scenario), alpha = 0.6) +
  geom_segment(aes(x = 1995, xend = 2100, y = 0, yend = 0),
               linetype = "dashed", color = "#454545", lwd = 0.7) +
   annotate(geom = 'text',
           x=rep(2000,2)[1],
           y=10,
           label = ' Ensemble machine learning',
           family = 'Times',hjust = 0,
           size=6,na.rm = T) +
  scale_y_continuous(breaks=c(-20,-10,0,10),limits = c(-18,12)) +
  scale_x_continuous(breaks=c(seq(2000,2100,25)),
                  limits = c(1995, 2130), expand = c(0, 0)) + # not expand， set the xlim exact at 1995 投2130 for segment dash
  labs(y= 'SOC storage change (Pg)',x= 'Year') +
   geom_line(data = m2_df_ave_family_ave,aes(x=time,y=anomaly,group=scenario,col=scenario),lwd=2) +
  scale_fill_manual(name='Scenario',values = c('SSP245'= '#92c5de','SSP585'= '#fcd1c5'),
                    labels = c('SSP245' = 'SSP2-4.5','SSP585' ='SSP5-8.5')) +
  scale_color_manual(name = 'Scenario',values = c('SSP245' = '#70a0cd','SSP585'='#990002'),
                     labels = c('SSP245' = 'SSP2-4.5','SSP585' ='SSP5-8.5')) +
  scale_shape_manual(name = 'ESMs',
                     values = shape_value1,
                     labels = shape_value0,
                   
                     drop=FALSE 
                     
  ) + 
  scale_linetype_discrete(guide = guide_legend(order=1,title='Scenario',show.legend = F,theme = theme(
    legend.direction = "vertical",
    legend.text = element_text(hjust = 0, vjust = 0.5, angle = 0)
  ))) +
  guides(shape = guide_legend(ncol = 2,
                              override.aes = list(size = 3), # shape size 
                              theme = theme(legend.title.position = "top",
                                            legend.byrow = TRUE))) +
  guides(col = guide_legend(order=1)) +
  guides(fill = guide_legend(order= 1,
                             theme = theme(
                               legend.direction = "vertical",
                               legend.title.position = "top",
                               legend.text.position = "right",
                               legend.text = element_text(hjust = 0, vjust = 0.5, angle = 0)
                             ))) +
  theme_classic() + 
  ggtitle(label = 'M2: Ensembled Machine Learning') +
  theme0

p20 <- 
  p2 + 
  geom_boxplot(data = m2_df_box,aes(x = scen, y = anomaly, fill = scenario),
               outlier.shape = NA, alpha = 0.7,size = 0.3) +
  geom_jitter(data= m2_df_box,  
              aes(x = scen, y = anomaly, 
                  shape =factor(family,levels = names(shape_value1))), 
              width = 0.5,   
              alpha = 1,    
              size = 3
  ) 
p20

# method 3 carbon balance 

forpeat = 0

get_exp <- function(x,ssp1,ssp2){
  x1 <- expand.grid(c('lu','cardamom'),x,c(ssp1,ssp2)) %>% 
    mutate(exp = paste0(Var1,"_",Var2,"_",Var3)) %>% pull(exp)
  return(x1)
}
get_data <- function(ye,model,data_file,resp,scen) {
  if(ye ==2100){
    exp0 <- get_exp(model,ssp1 = 'ssp245',ssp2 = 'ssp585')
  } else {exp0 <- get_exp(model,ssp1 = NULL,ssp2 = 'ssp585')}
  if(forpeat){
    data0 <- data.frame(do.call(cbind,lapply(exp0,function(x) unlist(readRDS(grep(x,data_file,value = T))[[2]])))[,grep(glob2rx(paste0(resp,'*',scen)),exp0)])
  } else {
    data0 <- data.frame(do.call(cbind,lapply(exp0,function(x) unlist(readRDS(grep(x,data_file,value = T))[[1]])))[,grep(glob2rx(paste0(resp,'*',scen)),exp0)])
  }
  names(data0) <- grep(glob2rx(paste0(resp,'*',scen)),exp0,value = T)
  return(data0)
}
f1 <- list.files(glob2rx('*annual20240826_newlitter_bedrock*0121.rds$'),path = '/mnt/DataSpace/Projects/Canada_C/Canada_C_final/model_outputs',full.names = T)
model_plant_clim2100 <- c('CanESM5','EC-Earth3-Veg','INM-CM4-8','IPSL-CM6A-LR','ACCESS-ESM1-5','EC-Earth3-Veg-LR','INM-CM5-0','BCC-CSM2-MR','CanESM5-1','EC-Earth3-CC')
m3_clim2100_plant <- grep('kfromplantclimateonly2100_ann',f1,value = T)

m3_clim2100_plant_ca_ssp245_lu <- get_data(ye=2100,model = model_plant_clim2100,data_file = m3_clim2100_plant,resp = 'lu',scen = 'ssp245') %>% mutate(Rh = 'Lu',scenario = 'SSP245',time = 2001:2100)
m3_clim2100_plant_ca_ssp585_lu <- get_data(ye=2100,model = model_plant_clim2100,data_file = m3_clim2100_plant,resp = 'lu',scen = 'ssp585') %>% mutate(Rh = 'Lu',scenario = 'SSP585',time = 2001:2100)
m3_clim2100_plant_ca_ssp245_card <- get_data(ye=2100,model = model_plant_clim2100,data_file = m3_clim2100_plant,resp = 'cardamom',scen = 'ssp245') %>% mutate(Rh = 'Cardamom',scenario = 'SSP245',time = 2001:2100)
m3_clim2100_plant_ca_ssp585_card <- get_data(ye=2100,model = model_plant_clim2100,data_file = m3_clim2100_plant,resp = 'cardamom',scen = 'ssp585') %>% mutate(Rh = 'Cardamom',scenario = 'SSP585',time = 2001:2100)

m3_df_245_1 <- m3_clim2100_plant_ca_ssp245_lu %>% mutate(exp = 'Clim2100') %>% tidyr::pivot_longer(col = starts_with('lu'),names_to = 'experiment',values_to = 'SOC_stock')
m3_df_585_2 <- m3_clim2100_plant_ca_ssp585_lu %>% mutate(exp = 'Clim2100') %>% tidyr::pivot_longer(col = starts_with('lu'),names_to = 'experiment',values_to = 'SOC_stock')
m3_df_245_3 <- m3_clim2100_plant_ca_ssp245_card %>% mutate(exp = 'Clim2100') %>% tidyr::pivot_longer(col = starts_with('card'),names_to = 'experiment',values_to = 'SOC_stock')
m3_df_585_4 <- m3_clim2100_plant_ca_ssp585_card %>% mutate(exp = 'Clim2100') %>% tidyr::pivot_longer(col = starts_with('card'),names_to = 'experiment',values_to = 'SOC_stock')

#
exp0 = 'All'
exp0 = 'NOTALL'
peat000 = FALSE # if for peat
m3_df_family <- rbind(m3_df_245_1,m3_df_585_2,m3_df_245_3,m3_df_585_4) %>% 
  mutate(model = gsub('lu_|cardamom_|_ssp245|_ssp585',replacement = "",experiment)) %>% 
  #mutate(family = substr(model,1,4),anomaly = SOC_stock - 152.8529) %>% #before bedrock 160.753
  {if (peat000) {mutate(.,family = substr(model,1,4),anomaly = SOC_stock - 41.93237) # M3 base from /mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/007_1_for_updateSOC2ndtimes/12_02_visual_combine_spatial2_2025new.R
  } else {mutate(.,family = substr(model,1,4),anomaly = SOC_stock - 149.1977)} } %>%
  #peat 
  #mutate(family = substr(model,1,4),anomaly = SOC_stock - 28.78664) %>% #for peat 
  filter(time >= 2000 & time <= 2100) %>% 
  group_by(family,Rh,scenario,time) %>% 
  summarize(anomaly = mean(anomaly),time = time[1]) %>% mutate(scenario = toupper(scenario))
#
m3_df_box <- m3_df_family %>% filter(time == 2100,Rh=='Cardamom') %>% mutate(scen = case_when(scenario == 'SSP245' ~2110,
                                                                                              scenario == 'SSP585' ~2120))

m3_df_ave_family_ave <- m3_df_family %>% 
  group_by(Rh,scenario,time) %>% 
  summarize(upper = mean(anomaly) + 1.645 * sd(anomaly),
            lower = mean(anomaly) - 1.645 * sd(anomaly),
            upper1 = max(anomaly),
            lower1 = min(anomaly),
            anomaly=mean(anomaly)) %>% mutate(scenario = toupper(scenario))


model_f <- unique(m3_df_family$family)

model_all <- c("ACCESS-ESM1-5","CanESM5","EC-Earth3-Veg-LR","IPSL-CM6A-LR" )

model_clim2100 <- c('ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5(2)','EC-Earth3(3)','INM-CM(2)','IPSL-CM6A-LR')
get_model <- function(y) {unlist(lapply(model_f,function(x) grep(x,y,value=T)))}
if(exp0 == 'All') {model1 <- get_model(model_all)} else {model1 <- get_model(model_clim2100)}
if(exp0 == 'All') {delta245_lu = 4.22;delta585_lu =2.02;delta245_card = -6.81;delta585_card =-5.42 
} else {delta245_lu = 4.22;delta585_lu =2.02;delta245_card = -6.81;delta585_card =-5.42}
# two Rh average ssp585:-25.2 ssp245:-26.2
#
esm_cols <- colorRampPalette(brewer.pal(12,'Paired'))(length(model1))
names(esm_cols) <- model_f
library(scales)

m3_df_ave_family_ave %>% filter(time==2100 & time ==2001)
#
# 
ylim_min3 <- -18
ylim_max3 <- 12

# Calculate positions as a ratio of the y-range
y_pos_13 <- ylim_min3 + 0.88 * (ylim_max3 - ylim_min3)  # 70% height
y_pos_23 <- ylim_min3 + 0.96 * (ylim_max3 - ylim_min3)  # 50% height
p3 <- 
  ggplot() +

  geom_ribbon(data=m3_df_ave_family_ave %>% filter(Rh =='Cardamom'),
              aes(x= time, ymin = lower, ymax = upper,
                  group= interaction(Rh,scenario),fill=scenario), alpha = 0.6) +
  geom_segment(aes(x = 1995, xend = 2100, y = 0, yend = 0),
               linetype = "dashed", color = "#454545", lwd = 0.7) +
  
  annotate(geom = 'text',
           x=rep(2000,2)[1],
            y=10,
           label = ' Carbon balance approach',
           family = 'Times',hjust = 0,
           size=6,na.rm = T) +
  
  scale_y_continuous(breaks=c(-20,-10,0,10),limits = c(-18,12)) +
  scale_x_continuous(breaks=c(seq(2000,2100,25)),
                     limits = c(1995, 2130), expand = c(0, 0)) + 
  labs(y= 'SOC storage change (Pg)',x= 'Year') +
   geom_line(data = m3_df_ave_family_ave %>% filter(Rh =='Cardamom'),
            aes(x=time,y=anomaly,group=interaction(Rh,scenario),col=scenario),lwd=2) +
   scale_fill_manual(name='Scenario',values = c('SSP245'= '#92c5de','SSP585'= '#fcd1c5'),
                    labels = c('SSP245' = 'SSP2-4.5','SSP585' ='SSP5-8.5')) +
  scale_color_manual(name = 'Scenario',values = c('SSP245' = '#70a0cd','SSP585'='#990002'),
                     labels = c('SSP245' = 'SSP2-4.5','SSP585' ='SSP5-8.5')) +
  scale_shape_manual(name = 'ESMs',
                     values = shape_value1,
                     
                     labels= shape_value0,
                   
                     drop= FALSE 
  ) + 
  scale_linetype_discrete(guide = guide_legend(order=1,title='Scenario',show.legend = F,theme = theme(
    legend.direction = "vertical",
    legend.text = element_text(hjust = 0, vjust = 0.5, angle = 0)
  ))) +
  guides(shape = guide_legend(ncol = 2,
                              override.aes = list(size = 3), 
                              theme = theme(legend.title.position = "top",
                                            legend.byrow = TRUE))) + 
  guides(col = guide_legend(order=1)) +
  guides(fill = guide_legend(order= 1, 
                             theme = theme(
                               legend.direction = "vertical",
                               legend.title.position = "top",
                               legend.text.position = "right",
                               legend.text = element_text(hjust = 0, vjust = 0.5, angle = 0)
                             ))) +
  theme_classic() + 
  ggtitle(label = 'M3: Carbon Transient Dynamic Method') +
  theme0
p3

p30 <- 
  p3 + 
  geom_boxplot(data = m3_df_box,aes(x = scen, y = anomaly, fill = scenario),
               outlier.shape = NA, alpha = 0.7,size = 0.3) +
  geom_jitter(data= m3_df_box,  
              aes(x = scen, y = anomaly, 
                  shape =factor(family,levels = names(shape_value1))), 
              width = 0.5,   
              alpha = 1,    
              size = 3
  ) 
p30
#


theme1 <- theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.length.y = unit(0, "pt"),
  axis.text.x.top = element_blank(),   
  axis.ticks.x.top = element_blank(), 
  axis.text.y.right = element_blank(),  
  axis.ticks.y.right = element_blank(),
  axis.ticks = element_line(linewidth = 1, color = "black")  
)
fig1 <- 
  (p20 + labs(x="") + 
     guides(shape = guide_none(),fill = guide_none(),col = guide_none()) +
     theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
           axis.text.x.top = element_blank(),   
           axis.ticks = element_line(linewidth = 1, color = "black"), 
           axis.ticks.x.top = element_blank(),  
           axis.text.y.right = element_blank(),  
           axis.ticks.y.right = element_blank(),
           plot.tag = element_text(hjust = -20,
                                   vjust = -2 )) + 
     scale_x_continuous(breaks=c(seq(2000,2100,25)),labels = c(2000,'',2050,'',2100),
                        limits = c(1995, 2130), expand = c(0, 0),sec.axis = dup_axis(name = NULL, labels = NULL)) +  # 顶部X轴，无标签
     scale_y_continuous(breaks=c(-20,-10,0,10),limits = c(-18,12),sec.axis = dup_axis(name = NULL, labels = NULL)) +  # 右侧Y轴，无标签
     ggtitle("")) + theme(axis.text.x.bottom = element_text(margin = margin(t=4,unit = 'pt'))) +
  (p30 + labs(y="") + 
     guides(shape = guide_none(),fill = guide_none(),col = guide_none())+ 
     theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
           plot.tag = element_text(
             vjust = -2
             
           )  ) +  theme1 + 
     scale_x_continuous(breaks=c(seq(2000,2100,25)),labels = c(2000,'',2050,'',2100),
                        limits = c(1995, 2130), expand = c(0, 0),sec.axis = dup_axis(name = NULL, labels = NULL)) +  # 顶部X轴，无标签
     scale_y_continuous(breaks=c(-20,-10,0,10),limits = c(-18,12),sec.axis = dup_axis(name = NULL, labels = NULL)) +  # 右侧Y轴，无标签
     ggtitle("")) + theme(axis.text.x.bottom = element_text(margin = margin(t=4,unit = 'pt'))) +
  p10+ labs(y="",x="")  + ggtitle("") + 
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(
          vjust = -2
        )) + theme1 +
  theme(
    legend.key.size = unit(0.6, "cm"), 
    legend.margin = margin(0, 0, 0, 0), 
    legend.box.margin = margin(0, 0, 0, 0),
    axis.text.x.bottom = element_text(margin = margin(t=4,unit = 'pt'))
  )+  
  scale_x_continuous(breaks=c(seq(2000,2100,25)),labels = c(2000,'',2050,'',2100),
                     limits = c(1995, 2130), expand = c(0, 0),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +  
  scale_y_continuous(breaks=c(-20,-10,0,10),limits = c(-18,12),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  plot_annotation(tag_levels = 'a')
fig1

path0 <- '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/figures/Figures_new/'
ggsave(plot=fig1,filename = paste0(path0,'Fig1_temporal_trend.tiff'),
       width = 15,height = 6,units = 'in',dpi = 600,compression = 'lzw')
pdf(file = paste0(path0,"fig1_temporal_trend.pdf"),   
    width = 15, 
    height = 6) 
dev.off()


