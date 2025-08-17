



############################           Figure 3         ########################
library(terra)
library(lubridate)
library(sf)
library(dplyr)
library(spData)
library(zoo) 
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
library(RColorBrewer)
quadrotic_data <- readRDS('/mnt/DataSpace/Projects/Canada_C/processed_output/CMIP6_LAI_to_GDD5_and_Observation_relationship_removeland15_ssp245.rds')[[3]]#1 quadrotic 2 gaussian 3 sigmoid 4 exp
quadratic_df <- do.call(rbind,quadrotic_data)
quadratic_df$model <- factor(quadratic_df$model,
                             levels = c('ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','CanESM5-1','CESM2-WACCM',
                                        'EC-Earth3-CC','EC-Earth3-Veg','EC-Earth3-Veg-LR','INM-CM4-8','INM-CM5-0',
                                        'IPSL-CM6A-LR','TaiESM1','Observation'))

cols <- c(colorRampPalette(brewer.pal(8,'Set1'))(length(unique(quadratic_df$model))-1),
          'black')
names(cols) <- levels(quadratic_df$model)

theme0 <- theme(axis.text = element_text(face = 'bold',family = 'Times',size=16),
                axis.text.x.bottom = element_text(face = 'bold',margin=margin(10,00,0,0)),
                axis.title.y.left = element_text(family = 'Times',size=18,margin=margin(0,10,0,0)),
                axis.title.x.bottom = element_text(family = 'Times',size=18,margin=margin(10,10,0,0)),
                panel.border = element_rect(          
                  colour = "black",
                  size   = 2,                            
                  fill   = NA ),
                axis.ticks.x.top = element_blank(),
                axis.ticks.y.right = element_blank(),
                axis.title.x.top = element_blank(),
                axis.title.y.right = element_blank(),
                axis.text.x.top = element_blank(),
                axis.text.y.right = element_blank(),
                plot.background = element_blank(),
                rect = element_rect(fill = "transparent"),
                legend.position = 'right',
                legend.direction = "verticle",
                legend.box = "verticle", 
                legend.box.just = "top", 
                legend.text = element_text(family = 'Times'),
                legend.title = element_text(family = 'Times'),
                legend.spacing.y = unit(0, "cm"),
                title = element_text(family = 'Times',size = 20)
)
#
#
mod_u <- c('ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','CanESM5-1','CAS-ESM2-0','CESM2-WACCM',
           'EC-Earth3-CC','EC-Earth3-Veg','EC-Earth3-Veg-LR','INM-CM4-8','INM-CM5-0',
           'IPSL-CM6A-LR','TaiESM1')
mod_u <- mod_u[-5] 
#
cmip6_esm <- quadratic_df %>% filter(model %in% mod_u)
observation <- quadratic_df %>% filter(model %in% 'Observation')

cols <- c(colorRampPalette(brewer.pal(8,'Set1'))(length(unique(quadratic_df$model))-1),
          'black')
names(cols) <- levels(quadratic_df$model)


### in order to share legend key 
##
csoil00 <- c('ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','CESM2-WACCM','CMCC-ESM2','EC-Earth3-CC','EC-Earth3-Veg','EC-Earth3-Veg-LR','GFDL-ESM4','IPSL-CM6A-LR','TaiESM1')

model_all <- c(unique(c(mod_u,csoil00)),'Observation')
model_all <- c('ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','CanESM5-1','CESM2-WACCM','EC-Earth3-CC','EC-Earth3-Veg','EC-Earth3-Veg-LR','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','TaiESM1','CMCC-ESM2','GFDL-ESM4','Observation')
cols <- setNames(c(colorRampPalette(brewer.pal(10,'Paired'))(length(model_all)-1),
                   'black'),model_all)
shapes <- setNames(c(4:15,17,18,16),model_all)
linetypes <- setNames(c(2:15,1),model_all)

pgdd <- 
  ggplot() +
  geom_point(data = cmip6_esm, aes(
    x = GDD,
    y = LAI_sigmoid,
    group = model,
    shape = model,
    col = model
  ),alpha = 0.4) +
  geom_line(data = cmip6_esm, aes(
    x = GDD,
    y = LAI_sigmoid,
    group = model,
    linetype=model,
    col = model
  ),lwd =1.2 ) +
  geom_line(data = observation, aes(
    x = GDD,
    y = LAI_sigmoid,
    linetype='Observation',
    col = 'Observation',
  ),lwd =2 ) +
  geom_point(data = observation, aes( 
    x = GDD,
    y = LAI_sigmoid,
    shape='Observation',
    col = 'Observation',
  ) ) + 
  scale_color_manual(values = cols  ) + 
  scale_linetype_manual(values = linetypes) + 
  scale_shape_manual(values = shapes) +
  
  labs(x = 'GDD5 (\u00B0C)',y=expression('LAImax ('~m^2~m^-2~')'),
       col ='ESM',linetype='ESM',shape = '') +
  guides(
    shape = guide_legend(
      title = 'ESM',
      label = FALSE,                
      override.aes = list(size = 2.5),  
      order = 1
    ),
    linetype = guide_legend(
      title = '',
      label = TRUE,
      override.aes = list(size = 1.2),
      order = 2
    ),col = guide_legend(
      title = '',
      label = TRUE,
      order = 2
    )) + 
  scale_x_continuous(
    sec.axis = dup_axis()
  ) +
  scale_y_continuous(
    limits = c(0,4.5),
    sec.axis = dup_axis()
  ) +
  theme_bw() + theme(
    legend.position = "right",
    legend.box      = "horizontal",
    legend.spacing.x = unit(0.01, "pt"),    
    legend.title = element_text(size = 10,family = 'Times'),
    legend.text  = element_text(size = 9,family ='Times'),
    legend.box.margin = margin(t=0,r=-10,b=0,l=-10,unit = 'pt'),
    axis.text = element_text(face = 'bold',family = 'Times',size=16),
    axis.text.x.bottom = element_text(face = 'bold',margin=margin(10,00,0,0)), 
    axis.title.y.left = element_text(family = 'Times',size=18,margin=margin(0,10,0,0)),
    axis.title.x.bottom = element_text(family = 'Times',size=18,margin=margin(10,10,0,0)),
    panel.border = element_rect(        
      colour = "black",
      size   = 2,                          
      fill   = NA ),                           
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.title.x.top = element_blank(),
    axis.title.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.text.y.right = element_blank(),
    plot.background = element_blank(),
    rect = element_rect(fill = "transparent"),
    title = element_text(family = 'Times',size = 20)
  ) 
pgdd + theme(plot.margin = margin(t = 5.5, r = 8, b = 8, l = 18,unit = 'pt')) + theme(legend.position = 'none')

##### fig 3b
slope00 <- readRDS('/mnt/DataSpace/Projects/Canada_C/processed_output/CMIP6_LAI_to_GDD5_and_Observation_relationship_slopemax_removeland15_ssp245.rds')[[3]]
slope01 <- as.data.frame(matrix(unlist(slope00),ncol=2,byrow = T))

names(slope01) <- c('Slopemax','model')
slope01$source <- factor(c(rep('ESM',nrow(slope01)-1),'Observation'))
slope01$Slopemax <- round(as.numeric(slope01$Slopemax),4)
slope01$model <- factor(slope01$model,levels = unique(slope01$model))
slope01$ID <- as.numeric(slope01$source)
median(slope01$Slopemax)
mean(c(slope01$Slopemax[c(1,2,5,11,12)],mean(slope01$Slopemax[c(3,4)]),mean(slope01$Slopemax[c(6:8)]),mean(slope01$Slopemax[c(9,10)])))

median(c(slope01$Slopemax[c(1,2,5,11,12)],mean(slope01$Slopemax[c(3,4)]),mean(slope01$Slopemax[c(6:8)]),mean(slope01$Slopemax[c(9,10)])))
slope_fig3b <- 
  ggplot() +
  geom_boxplot(data = slope01, aes(x = source, y = Slopemax * 1000),
               col=c('maroon','black'),cex=1) +
  geom_point(
    data = slope01,
    aes(x = source, y = Slopemax *1000, shape = model,col = model),
    position = position_jitter(width = 0.15, seed = 5),
    show.legend = F,
    alpha = 1,
    stroke = 0.5,
    size = 3.5
  ) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = cols) +
  labs(x = ' ',y=expression(Max~omega~" ("*"×"*10^-3~m^2~m^-2~degree*C^-1*")"),
       col ='ESM',shape='ESM') +
  guides(col = guide_legend(ncol = 3), 
         linetype = guide_legend(ncol = 3)) +  
  scale_y_continuous(
    limits = c(0,18),
    breaks = c(0,4,8,12,16)
  ) +
  theme_bw()  + 
  theme(legend.position = "right",
        legend.box      = "horizontal",
        legend.spacing.x = unit(0.01, "pt"),    
        legend.title = element_text(size = 10,family = 'Times'),
        legend.text  = element_text(size = 9,family ='Times'),
        legend.box.margin = margin(t=0,r=-10,b=0,l=-10,unit = 'pt'),
        axis.text = element_text(face = 'bold',family = 'Times',size=16),
        axis.text.x.bottom = element_text(face = 'bold',margin=margin(10,00,0,0)), 
        axis.title.y.left = element_text(family = 'Times',size=18,margin=margin(0,0,0,0)),
        axis.title.x.bottom = element_text(family = 'Times',size=18,margin=margin(10,10,0,0)),
        panel.border = element_rect(          
          colour = "black",
          size   = 2,                          
          fill   = NA ),                          
        axis.ticks.x.top = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.text.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        plot.background = element_blank(),
        rect = element_rect(fill = "transparent"),
        title = element_text(family = 'Times',size = 20))
slope_fig3b
  

#######fig 3c
df_a <- readRDS('/mnt/DataSpace/Data_pool/CMIP6_cSoil_trend/csoil_rh_turnover_dataframe_withslopes_20250601.rds')
csoil00 <- c('ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','CESM2-WACCM','CMCC-ESM2',
             'EC-Earth3-CC','EC-Earth3-Veg','EC-Earth3-Veg-LR','GFDL-ESM4','IPSL-CM6A-LR',
             'TaiESM1')

two_slopes <- df_a[[2]]
# 
calculate_Q10 <- function(slope, T_low_C = 0, T_high_C = 10,cond = 'Warmer') {
  if(cond == 'Warmer'){
    
    T1 <- T_low_C + 273.15
    T2 <- T_high_C + 273.15
  } else {
    T1 <- T_low_C + 273.15 -10
    T2 <- T_high_C + 273.15 -10
    }
  
  exp(abs(slope) * 10 / (T1 * T2))
}

slopes <- two_slopes$Warmer
slopes <- two_slopes$Cooler

Q10_values_w <- sapply(two_slopes$Warmer, calculate_Q10,cond='Warmer')
Q10_values_c <- sapply(two_slopes$Cooler, calculate_Q10,cond='Cooler')
Q10s <- data.frame(Q10_W = Q10_values_w,Q10_C= Q10_values_c,model = c(csoil00[c(1,3:11)],'Observation')) %>%
  pivot_longer(data = .,cols = starts_with('Q10'),values_to = 'Q10',names_to = 'T')




####
df_int <- df_a[[1]] %>% mutate(model = factor(model,levels = c(csoil00[c(1,3:11)],'Observation')))
df_pred <- df_a[[3]] %>% mutate(model = factor(model,levels = c(csoil00[c(1,3:11)],'Observation')))
head(df_int,2)
head(Q10s)
esm_model <- csoil00[c(1,3:11)]
logKK3 <- 
  ggplot() + 
  geom_point(aes(x = invTK, y = log_k,
                 group = model,
                 shape = model,
                 col=model),data = df_int) +
  geom_line(aes(x = invTK, y = predicted_k,
                group = model,
                col=model),data = df_pred, linewidth = 1.5) +
  geom_vline(xintercept=1/273.15,linetype=2,linewidth = 1.2,col='maroon') +
  annotate('text',label = '0 (\u00B0C)',x = 1/271,y = -5.9,size = 6, family = 'Times') + 
  annotate('text',label = c('Warmer','Colder'),x = c(1/278.15,1/268.15),
           y = -6.25,size =8,family = 'Times') + 
  annotate('segment',x = 1/273.15, xend = 1/271.6,y = -6.2,yend = -6.05,col='maroon',
           size = 1,arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +  
  labs( y = expression("ln k (yr"^{-1}*")"),
        x = expression("Inverse mean annual temperature (1/T, K"^{-1}*")"),
        col ='ESM',shape = 'ESM',linetype='ESM') +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = shapes) + 
  scale_linetype_manual(values = shapes) +
  theme_classic() +
  scale_x_continuous(
     sec.axis = dup_axis()
  ) +
  scale_y_continuous(
    limits = c(-6.5,-2),
    breaks = c(-2,-4,-6),
    sec.axis = dup_axis()
  ) +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(face = 'bold',family = 'Times',size=16),
        axis.text.x.bottom = element_text(face = 'bold',margin=margin(10,00,0,0)), 
        axis.title.y.left = element_text(family = 'Times',size=18,margin=margin(0,0,0,0)),
        axis.title.x.bottom = element_text(family = 'Times',size=18,margin=margin(10,10,0,0)),
        axis.line.x.bottom = element_line(size=1,color = "black"),
        axis.line.y.left = element_line(size=1,color = "black"),
        axis.line.x.top = element_line(size=1,color = 'black'),
        axis.line.y.right = element_line(size=1,color = 'black'),
        axis.ticks.x.top = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.text.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        plot.background = element_blank(),
        rect = element_rect(fill = "transparent"),
        legend.position = c(0.85,0.7),
        legend.box = "horizontal",
        legend.box.just = "top", 
        legend.text = element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'),
        legend.spacing.y = unit(0, "cm"),
        title = element_text(family = 'Times',size = 20)
  ) 
logKK3
##### fig. 3d
head(Q10s)
Q10s <- Q10s %>% mutate(Source = c(rep('ESM',20),rep('Observation',2)),
                        xnum = as.numeric(factor(T)),
                        model = factor(model,levels = unique(Q10s$model))) 
Q10_fig3d <- 
  ggplot() + 
  geom_boxplot(aes(x=T,y=Q10),data = Q10s %>% filter(Source =='ESM'),
                col = scales::alpha(c('#048eca','#f55c39'),alpha = 0.7),
               size = 1) + 
  geom_jitter(aes(x=xnum,
                  y=Q10,
                  shape = model,
                  col= model
                  ),data = Q10s %>% filter(Source =='ESM'),
              width = 0.06,   
              alpha = 0.5, 
              stroke = 0.5,
              size = 3.5
              ) +
  geom_hline(yintercept = 2,linetype=4,linewidth = 0.8,col = 'grey') +
 
  geom_line(aes(x = c(1,2),y = Q10 ),data = Q10s %>% 
              filter(Source =='ESM') %>% group_by(T) %>% 
              summarise(Q10 = median(Q10),.groups = 'drop') %>% as.data.frame(),
            size = 2,col = 'maroon') +
  geom_point(aes(x = xnum,y = Q10,shape = model,col = model),size =3.5, 
             data = Q10s %>% filter(Source =='Observation')) + 
  geom_line(aes(x = xnum,y = Q10),data = Q10s %>% filter(Source =='Observation'),
            linewidth = 2,col = 'black') + 
  labs( y = expression(Q[10]), 
        x = 'Regions',
        col ='ESM',shape = 'ESM',linetype='ESM') +
  scale_shape_manual(values = shapes) + 
  scale_color_manual(values = cols) + 
  scale_linetype_manual(values = linetypes) +
  scale_x_discrete(
    breaks = c('Q10_C','Q10_W'),
    labels = c(expression(paste(Q[10],' colder')),expression(paste(Q[10],' warmer')))
  ) +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(face = 'bold',family = 'Times',size=16),
        axis.text.x.bottom = element_text(face = 'bold',margin=margin(10,00,0,0)), 
        axis.title.y.left = element_text(family = 'Times',size=18,margin=margin(0,10,0,0)),
        axis.title.x.bottom = element_text(family = 'Times',size=18,margin=margin(10,10,0,0)),
        panel.border = element_rect(        
          colour = "black",
          size   = 2,                             
          fill   = NA ),                           
        axis.ticks.x.top = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.text.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        plot.background = element_blank(),
        rect = element_rect(fill = "transparent"),
        legend.position = c(0.85,0.7),
        legend.box = "horizontal", 
        legend.box.just = "top", 
        legend.text = element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'),
        legend.spacing.y = unit(0, "cm"),
        title = element_text(family = 'Times',size = 20)
  ) 
  
logKK3 
Q10_fig3d + theme(axis.text = element_text(face = 'bold',family = 'Times',size=16))

######## for all legend 
legend_df <- data.frame(x = rnorm(n= length(model_all) * 10),
                        y= rnorm(n= length(model_all) * 10),
                        model= factor(model_all,levels = model_all))
legend_fig3 <- 
  ggplot() +
  geom_point(data=legend_df,aes(x=x,y=y,shape=model,linetype = model)) + 
  geom_line(data=legend_df,aes(x=x,y=y,col=model,shape=model,linetype = model)) +
  scale_color_manual(values = cols) +
  scale_shape_manual(values =shapes) + 
  scale_linetype_manual(values = linetypes) + 
  guides(shape = guide_legend(
    title = 'Source',
    label = T,                
    override.aes = list(size = 3),  
    order = 2),
    linetype = guide_legend(
      title = 'Source',
      label = T,
      override.aes = list(size = 2),
      order = 1),
    col = guide_legend(
      title = 'Source',
      label = T,
      order = 1
    )) + theme_bw() + 
  theme(legend.position = "right",
        legend.box      = "verticle",
        legend.box.spacing = unit(0.01,'cm'),
        legend.spacing.x = unit(0.01, "cm"),  
        legend.title = element_text(size = 16,family = 'Times'),
        legend.text  = element_text(size = 12,family ='Times'),
        legend.key.height = unit(0.5,'cm'),
        legend.box.margin = margin(t=0,r=0,b=0,l=0,unit = 'pt'),
        axis.text = element_text(face = 'bold',family = 'Times',size=16),
        axis.text.x.bottom = element_text(face = 'bold',margin=margin(10,00,0,0)), 
        axis.title.y.left = element_text(family = 'Times',size=18,margin=margin(0,10,0,0)),
        axis.title.x.bottom = element_text(family = 'Times',size=18,margin=margin(10,10,0,0)),
        panel.border = element_rect(          
          colour = "black",
          size   = 2,                            
          fill   = NA ),                           
        axis.ticks.x.top = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        axis.text.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        plot.background = element_blank(),
        rect = element_rect(fill = "transparent"),
        title = element_text(family = 'Times',size = 20))
legend_fig3
#####################################

legend_f3 <- cowplot::get_legend(legend_fig3)

#############align
tt1 <- align_plots(slope_fig3b + theme(panel.grid = element_blank(),plot.margin = margin(r=10,l=10,t=5.5,b=5.5,unit='pt')),
                   Q10_fig3d + 
                     theme(legend.position = 'none'), align = "v", axis = "l")

combined_cow3 <- plot_grid(
  plot_grid(
    combined_cow2,
    combined_cow,
    ncol = 2,rel_widths = c(2,1)
  ),legend_f3,
  ncol = 2,rel_widths = c(6,1),
  align      = "h",       
  axis       = "tblr"      
);combined_cow3
ggsave(filename = '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/figures/Figures_new/fig3_new.tiff',
       width = 12.5,height = 8.5,
       units = 'in',
       dpi = 600,bg = 'white',
       compression = 'lzw')
ggsave(filename = '/mnt/File0/DAAATAAA/Data_collection/Canada_soil_archive/Canada_SOC_rebuild/figures/Figures_new/fig3_new.pdf',
       width = 12.5,height = 8.5,units = 'in',
       bg = 'white')