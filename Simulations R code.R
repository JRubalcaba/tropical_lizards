require(NicheMapR)
require(tidyr)
require(dplyr)
require(raster)
require(Rcpp)
library(ncdf4) 
library(rgdal)
require(ggplot2)

############## Create map ###########

dir <- ".../Sources"

shapeTotal <- raster(paste0(dir,"/bio1.bil"))

resolution = 0.5 # Set resolution (deg)

world <- aggregate(shapeTotal, fact=resolution/res(shapeTotal)[1]) 

tropical_forest0 <- shapefile(paste0(dir,"/Shapes/amazon_sensulatissimo_gmm_v1.shp"))
Aameiva_range <- shapefile(paste0("/Shapes/Aameiva_range.shp"))

Amazonia_map <- mask(world, tropical_forest0)
Amazonia_map <- mask(Amazonia_map, Aameiva_range)
Amazonia_map <- crop(Amazonia_map, c(-100, -30, -40, 30)) # Extent

plot(Amazonia_map)

xy <- xyFromCell(Amazonia_map, 1:ncell(Amazonia_map))
# cells <- cellFromPolygon(Amazonia_map, tropical_forest0)[[10]]
vals <- raster::extract(Amazonia_map, xy)
cells <- which(!is.na(vals))
xy.values <- xy[cells,]
nrow(xy.values)

# save(xy.values, file="C:/Users/juanv/Dropbox/Marie Curie/Modeling impacts of deforestation/ProcB/Code/Sources/xy.values_Aameiva.RData")
# save(cells, file="C:/Users/juanv/Dropbox/Marie Curie/Modeling impacts of deforestation/ProcB/Code/Sources/cells_Aameiva.RData")

############## Model parameters and functions ###########

require(Rcpp)
dir <- ".../Sources" # set working directory

# 1. Load model functions
load(paste0(dir,"/Te_model.RData")) # Operative temperature model function
load(paste0(dir,"/Tb_t.RData")) # Transitent body temperature model
sourceCpp(paste0(dir,"/Cpp_functions.cpp")) # Probabilistic behavioral thermoregulation model (in Cpp)

# 2. Load the Leaf Area Index data and geographic coordinates of the tropical forest biome
# load(paste0(dir,"/LAI/LAI_matrix_nov_Aameiva_HD.RData"))
# load(paste0(dir,"/xy.values_HD.RData"))

# Only for sensitivity analysis (no HD)
load(paste0(dir,"/LAI/LAI_matrix_nov_Aameiva_HD.RData"))
load(paste0(dir,"/xy.values_Aameiva_HD.RData"))

# 3. Load the function compute microclimates (requires NicheMapR; Kearney and Porter 2017 Ecography) 
load(paste0(dir, "/compute_microclimates.RData"))

############## Run NicheMapR with simulated increments in air temperature ###########

M = 21.19        # Body mass (g) 
a = 0.88        # skin absorbance 
height = 0.01   # height from ground (cm)
month = 11       # month for simulations
sigma_shade = 1

Tpref_mean = 36  # prefered body temperature (ºC) 
Tpref_sd = 1    # body temperature standard deviation
max_iter = 10   # number of iterations

deltaT <- seq(0, 8, length.out = 10) # simulated temperature increase (between 0 and 4 ºC)

meanTb_matrix_df <- SDTb_matrix_df <- Te_shade_matrix <- Te_sun_matrix <- array(NA, dim=c(1440, length(deltaT)))
meanTb_matrix_list <- SDTb_matrix_list <- Te_shade_matrix_list <- Te_sun_matrix_list <- list()
for(i in 1:nrow(xy.values)){
  t0 <- Sys.time()
  LAI <- LAI_matrix_nov[i]

  # if(is.na(LAI)) LAI <- 1 # only 1 cell ###
  
  for(j in 1:length(deltaT)){
    ## Run NicheMapR at each location (grid cell)
    microclimate <- compute_microclimates(run.gads = 0, lat = xy.values[i,2], long = xy.values[i,1], month = month, height = height, 
                                          LAI = NA, maxshade=90, delta_T = deltaT[j])
    
    microclimate_sun <- microclimate$microclimate_sun_min
    microclimate_shade <- microclimate$microclimate_shade_min
    
    S_sun = microclimate_sun$S
    Ta_sun = microclimate_sun$Ta
    Tsk_sun = microclimate_sun$Tsk
    Tg_sun = microclimate_sun$Tg
    V_sun = microclimate_sun$V
    
    S_shade = microclimate_shade$S
    Ta_shade = microclimate_shade$Ta
    Tsk_shade = microclimate_shade$Tsk
    Tg_shade = microclimate_shade$Tg
    V_shade = microclimate_shade$V
    
    Te_sun <- Te_model(M=M, a=a, ground=T, S=S_sun, Tsk=Tsk_sun, Ta=Ta_sun, Tg=Tg_sun, v=V_sun)
    Te_shade <- Te_model(M=M, a=a, ground=T, S=S_shade, Tsk=Tsk_shade, Ta=Ta_shade, Tg=Tg_shade, v=V_shade)

    ## Run thermoregulation model
    output <- behav_therm(M=M, a=a, ground=1, S_sun=S_sun, Ta_sun=Ta_sun, Tsk_sun=Tsk_sun, Tg_sun=Tg_sun, V_sun=V_sun,
                          S_shade=S_shade, Ta_shade=Ta_shade, Tsk_shade=Tsk_shade, Tg_shade=Tg_shade, V_shade=V_shade,
                          Tpref_mean=Tpref_mean, Tpref_sd=Tpref_sd, sigma_sun=1, sigma_shade=sigma_shade, max_iter=max_iter)

    meanTb_matrix_df[,j] <- output[,1]
    SDTb_matrix_df[,j] <- output[,3]
    Te_shade_matrix[,j] <- Te_shade
    Te_sun_matrix[,j] <- Te_sun
  }
  
  meanTb_matrix_list[[i]] <- meanTb_matrix_df
  SDTb_matrix_list[[i]] <- SDTb_matrix_df
  Te_shade_matrix_list[[i]] <- Te_shade_matrix
  Te_sun_matrix_list[[i]] <- Te_sun_matrix
  
  t1 <- Sys.time()
  elapsed <- t1 - t0
  remaining <- (nrow(xy.values) - i) * elapsed

  print(paste(round(100 * i / nrow(xy.values),2), "%", "remaining:", round(remaining / 60 , 2), "min"))
}

############## Run NicheMapR across shade levels ###########

M = 21.19            # Body mass
a = 0.88         # skin absorbance
height = 0.01
month = 11

Tpref_mean = 34.34343 
Tpref_sd = 1
max_iter = 10

shade_reduction <- seq(0.1, 1, length.out = 10)

meanTb_matrix_df <- SDTb_matrix_df <- Te_shade_matrix <- array(NA, dim=c(1440, length(deltaT)))
meanTb_matrix_list <- SDTb_matrix_list <- Te_shade_matrix_list <- list()
for(i in 1:nrow(xy.values)){
  t0 <- Sys.time()
  LAI <- LAI_matrix_nov[i] 
  
  # if(is.na(LAI)) LAI <- 1 ##
  
  for(j in 1:length(shade_reduction)){
    ## Run NicheMapR
    
    LAI_reduction <- LAI * shade_reduction[j]
  
    microclimate <- compute_microclimates(run.gads = 0, lat = xy.values[i,2], long = xy.values[i,1], month = month, height = height, 
                                          LAI = LAI_reduction, delta_T = 0)
    
    microclimate_sun <- microclimate$microclimate_sun_min
    microclimate_shade <- microclimate$microclimate_shade_min
    
    S_sun = microclimate_sun$S
    Ta_sun = microclimate_sun$Ta
    Tsk_sun = microclimate_sun$Tsk
    Tg_sun = microclimate_sun$Tg
    V_sun = microclimate_sun$V
    
    S_shade = microclimate_shade$S
    Ta_shade = microclimate_shade$Ta
    Tsk_shade = microclimate_shade$Tsk
    Tg_shade = microclimate_shade$Tg
    V_shade = microclimate_shade$V
    
    Te_sun <- Te_model(M=M, a=a, ground=T, S=S_sun, Tsk=Tsk_sun, Ta=Ta_sun, Tg=Tg_sun, v=V_sun)
    Te_shade <- Te_model(M=M, a=a, ground=T, S=S_shade, Tsk=Tsk_shade, Ta=Ta_shade, Tg=Tg_shade, v=V_shade)
    
    ## Run thermoregulation model
    output <- behav_therm(M=M, a=a, ground=1, S_sun=S_sun, Ta_sun=Ta_sun, Tsk_sun=Tsk_sun, Tg_sun=Tg_sun, V_sun=V_sun,
                          S_shade=S_shade, Ta_shade=Ta_shade, Tsk_shade=Tsk_shade, Tg_shade=Tg_shade, V_shade=V_shade,
                          Tpref_mean=Tpref_mean, Tpref_sd=Tpref_sd, sigma_sun=1, sigma_shade=1, max_iter=max_iter)
    
    meanTb_matrix_df[,j] <- output[,1]
    SDTb_matrix_df[,j] <- output[,3]
    Te_shade_matrix[,j] <- Te_shade
  }
  
  meanTb_matrix_list[[i]] <- meanTb_matrix_df
  SDTb_matrix_list[[i]] <- SDTb_matrix_df
  Te_shade_matrix_list[[i]] <- Te_shade_matrix
  
  t1 <- Sys.time()
  elapsed <- t1 - t0
  remaining <- (nrow(xy.values) - i) * elapsed
  
  print(paste(round(100 * i / nrow(xy.values),2), "%", "remaining:", round(remaining / 60 , 2), "min"))
}

############## Run NicheMapR across shade availabilities ###########

M = 21.19            # Body mass
a = 0.88         # skin absorbance
height = 0.01
month = 11

Tpref_mean = 34.34343 
Tpref_sd = 1
max_iter = 10

shade_availability <- seq(0, 1, length.out = 30)

meanTb_matrix_df <- SDTb_matrix_df <- Te_shade_matrix <- array(NA, dim=c(1440, length(shade_availability)))
meanTb_matrix_list <- SDTb_matrix_list <- Te_shade_matrix_list <- list()
for(i in 1:nrow(xy.values)){
  t0 <- Sys.time()
  
  ## Run NicheMapR
  LAI <- LAI_matrix_nov[i]
  # if(is.na(LAI)) LAI <- 1 # only 1 cell ###
  
  microclimate <- compute_microclimates(run.gads = 0, lat = xy.values[i,2], long = xy.values[i,1], month = month, height = height, 
                                        LAI = NA, maxshade = 90, delta_T = 0)
  
  microclimate_sun <- microclimate$microclimate_sun_min
  microclimate_shade <- microclimate$microclimate_shade_min
  
  S_sun = microclimate_sun$S
  Ta_sun = microclimate_sun$Ta
  Tsk_sun = microclimate_sun$Tsk
  Tg_sun = microclimate_sun$Tg
  V_sun = microclimate_sun$V
  
  S_shade = microclimate_shade$S
  Ta_shade = microclimate_shade$Ta
  Tsk_shade = microclimate_shade$Tsk
  Tg_shade = microclimate_shade$Tg
  V_shade = microclimate_shade$V
  
  Te_sun <- Te_model(M=M, a=a, ground=T, S=S_sun, Tsk=Tsk_sun, Ta=Ta_sun, Tg=Tg_sun, v=V_sun)
  Te_shade <- Te_model(M=M, a=a, ground=T, S=S_shade, Tsk=Tsk_shade, Ta=Ta_shade, Tg=Tg_shade, v=V_shade)
  
  for(j in 1:length(shade_availability)){
    sigma_shade = shade_availability[j]
    
    ## Run thermoregulation model
    output <- behav_therm(M=M, a=a, ground=1, S_sun=S_sun, Ta_sun=Ta_sun, Tsk_sun=Tsk_sun, Tg_sun=Tg_sun, V_sun=V_sun,
                          S_shade=S_shade, Ta_shade=Ta_shade, Tsk_shade=Tsk_shade, Tg_shade=Tg_shade, V_shade=V_shade,
                          Tpref_mean=Tpref_mean, Tpref_sd=Tpref_sd, sigma_sun=1, sigma_shade=sigma_shade, max_iter=max_iter)
    
    
    meanTb_matrix_df[,j] <- output[,1]
    SDTb_matrix_df[,j] <- output[,3]
    Te_shade_matrix[,j] <- Te_shade
  }
  
  meanTb_matrix_list[[i]] <- meanTb_matrix_df
  SDTb_matrix_list[[i]] <- SDTb_matrix_df
  Te_shade_matrix_list[[i]] <- Te_shade_matrix
  
  t1 <- Sys.time()
  elapsed <- t1 - t0
  remaining <- (nrow(xy.values) - i) * elapsed
  
  print(paste(round(100 * i / nrow(xy.values),2), "%", "remaining:", round(remaining / 60 , 2), "min"))
}

############## Maps for sensitivity and buffering capacity ###########

dir <- ".../Sources" # set working directory

##### Air temperature increase #####

load(paste0(dir,"/Results/meanTb_matrix_list_deltaT_Aameiva_",resolution,"x",resolution,".RData"))
load(paste0(dir,"/Results/SDTb_matrix_list_deltaT_Aameiva_",resolution,"x",resolution,".RData"))
load(paste0(dir,"/Results/Te_shade_matrix_list_deltaT_Aameiva_",resolution,"x",resolution,".RData"))
load(paste0(dir,"/Results/Te_sun_matrix_list_deltaT_Aameiva_",resolution,"x",resolution,".RData"))

# no LAI

# load(paste0(dir,"/Results/meanTb_matrix_list_deltaT_noLAI_Aameiva_",resolution,"x",resolution,".RData"))
# load(paste0(dir,"/Results/SDTb_matrix_list_deltaT_noLAI_Aameiva_",resolution,"x",resolution,".RData"))
# load(paste0(dir,"/Results/Te_shade_matrix_list_deltaT_noLAI_Aameiva_",resolution,"x",resolution,".RData"))
# load(paste0(dir,"/Results/Te_sun_matrix_list_deltaT_noLAI_Aameiva_",resolution,"x",resolution,".RData"))

deltaT <- seq(0, 8, length.out = 10)

sensitivityT_dat <-Teshade_dat <- Tesun_dat <- numeric(length(meanTb_matrix_list))
for(i in 1:length(sensitivityT_dat)){
  meanTb_max <- colMeans(meanTb_matrix_list[[i]]) - mean(meanTb_matrix_list[[i]][,1])
  
  linear_model <- lm(meanTb_max ~ deltaT)
  
  sensitivityT_dat[i] <- coef(linear_model)[2]
  Teshade_dat[i] <- mean(Te_shade_matrix_list[[i]][,1])
  Tesun_dat[i] <- mean(Te_sun_matrix_list[[i]][,1])
}

sensitivityT_map <- Amazonia_map
j=0
for(i in cells){
  j=j+1
  sensitivityT_map[i] <- sensitivityT_dat[j]
}
plot(sensitivityT_map)

##### Buffering capacity #####

buffSD <- 0.2 # select the deltaT value for which TbSD is lower than this
buffer_dat <- numeric(length(meanTb_matrix_list))
for(i in 1:length(buffer_dat)){
  maxTshade <- which.max(Te_shade_matrix_list[[i]][,1])
  
  # sdTb_max <- SDTb_matrix_list[[i]][maxTshade,]
  sdTb_max <- colMeans(SDTb_matrix_list[[i]])
  if(any(sdTb_max < buffSD)){
    maxTbuff <- min(deltaT[which.min(abs(sdTb_max - buffSD))])
    buffer_dat[i] <- maxTbuff
  }else{
    buffer_dat[i] <- 8
  }
}

buffer_map <- Amazonia_map
j=0
for(i in cells){
  j=j+1
  buffer_map[i] <- buffer_dat[j]
}
plot(buffer_map)

##### Shade availability ##### 

load(paste0(dir,"/Results/meanTb_matrix_list_shadeavail_Aameiva_",resolution,"x",resolution,".RData"))
load(paste0(dir,"/Results/SDTb_matrix_list_shadeavail_Aameiva_",resolution,"x",resolution,".RData"))
load(paste0(dir,"/Results/Te_shade_matrix_list_shadeavail_Aameiva_",resolution,"x",resolution,".RData"))

# no LAI

# load(paste0(dir,"/Results/meanTb_matrix_list_shadeavail_noLAI_Aameiva_",resolution,"x",resolution,".RData"))
# load(paste0(dir,"/Results/SDTb_matrix_list_shadeavail_noLAI_Aameiva_",resolution,"x",resolution,".RData"))
# load(paste0(dir,"/Results/Te_shade_matrix_list_shadeavail_noLAI_Aameiva_",resolution,"x",resolution,".RData"))
# load(paste0(dir,"/Results/Te_sun_matrix_list_shadeavail_noLAI_Aameiva_",resolution,"x",resolution,".RData"))


shade_availability <- seq(0, 1, length.out = 30)

sensitivitySA_dat <- Teshade_dat <- numeric(length(meanTb_matrix_list))
for(i in 1:length(sensitivitySA_dat)){
  meanTb_max <- colMeans(meanTb_matrix_list[[i]]) - mean(meanTb_matrix_list[[i]][,30])
  linear_model <- lm(meanTb_max[2:30] ~ log(1/shade_availability[2:30]))
  sensitivitySA_dat[i] <- coef(linear_model)[2]
  Teshade_dat[i] <- Te_shade_matrix_list[[i]][maxTshade,30]
}

sensitivitySA_map <- Amazonia_map
j=0
for(i in cells){
  j=j+1
  sensitivitySA_map[i] <- sensitivitySA_dat[j]
}
plot(sensitivitySA_map)

##### Shade level #####

load(paste0(dir,"/Results/meanTb_matrix_list_shadelev_Aameiva_",resolution,"x",resolution,".RData"))
load(paste0(dir,"/Results/SDTb_matrix_list_shadelev_Aameiva_",resolution,"x",resolution,".RData"))
load(paste0(dir,"/Results/Te_shade_matrix_list_shadelev_Aameiva_",resolution,"x",resolution,".RData"))

shade_reduction <- seq(0.1, 1, length.out = 10)

sensitivitySL_dat <- Teshade_dat <- numeric(length(meanTb_matrix_list))
for(i in 1:length(sensitivitySL_dat)){
  meanTb_max <- colMeans(meanTb_matrix_list[[i]]) - mean(meanTb_matrix_list[[i]][,1])
  
  linear_model <- lm(meanTb_max~   shade_reduction)
  
  sensitivitySL_dat[i] <- coef(linear_model)[2]
  Teshade_dat[i] <- Te_shade_matrix_list[[i]][maxTshade,10]
}

sensitivitySL_map <- Amazonia_map
j=0
for(i in cells){
  j=j+1
  sensitivitySL_map[i] <- sensitivitySL_dat[j]
}
plot(sensitivitySL_map)

##### Climate change projections #####

# current 

folder <- "C:/globalclimate"
global_climate <- raster::brick(paste0(folder, "/global_climate.nc"))
TMAXX_nov <- global_climate$X61 # T max in nov New et al.
TMAXX_nov <- crop(TMAXX_nov, c(-100, -30, -40, 30)) # Extent
TMAXX_nov <- aggregate(TMAXX_nov, fact=resolution/res(TMAXX_nov)[1]) 
TMAXX_nov <- mask(TMAXX_nov, tropical_forest0)
TMAXX_nov <- mask(TMAXX_nov, Aameiva_range)

plot(TMAXX_nov*0.1)

# future

dir <- ".../Sources" # set working directory
ssp245 <- stack(paste0(dir, "/Climate change/wc2.1_10m_tmax_tmin_MIROC6_ssp245_2021-2040_HD.grd"))

ssp245_nov <- ssp245$wc2.1_10m_tmax_MIROC6_ssp245_2021.2040.11
ssp245_nov <- crop(ssp245_nov, c(-100, -30, -40, 30)) # Extent
# ssp245_nov <- aggregate(ssp245_nov, fact=resolution/res(ssp245_nov)[1]) 
ssp245_nov <- mask(ssp245_nov, tropical_forest0)
ssp245_nov <- mask(ssp245_nov, Aameiva_range)

plot(ssp245_nov)

# Temperature change (tmax in nov)

deltaT_map <- ssp245_nov - TMAXX_nov*0.1

plot(deltaT_map)


############## Figures ###########

require(ggplot2)
require(RColorBrewer)
require(ggpubr)

##### Fig. 4 ##### 
# OPEN THE CODE functions_multiple_scales.R AND LOAD ALL FUNCTIONS

south_america <- crop(world, c(-100, -30, -40, 30)) # Extent

south_america <- south_america/south_america
xy_southamerica <- xyFromCell(south_america, 1:ncell(south_america))
vals<-raster::extract(south_america, xy_southamerica)
vals[which(!is.na(vals))]<-0
vals[which(is.na(vals))]<-1

# sensitivityT_dat[1] <- 1

p1 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=sensitivityT_dat)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())


# sensitivitySA_dat[190] <- 0.15
# sensitivitySA_dat[151] <- 1.3

p2 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=sensitivitySA_dat)) +
  scale_fill_viridis_c()+
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())


p3 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=-1*sensitivitySL_dat)) +
  scale_fill_viridis_c()+
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

ggarrange(p1,p2,p3, nrow=1) # 1000 x 250

p1 <- ggplot() + geom_point(mapping=aes(y=sensitivityT_dat, x=LAI_matrix_nov)) +
  theme_classic() +  ylab("Sensitivity to air temperature") + xlab("Leaf area index") +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=14, colour="black"))

p2 <- ggplot() + geom_point(mapping=aes(y=sensitivitySA_dat, x=LAI_matrix_nov)) +
  theme_classic() + ylab("Sensitivity to shade availability") + xlab("Leaf area index") +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=14, colour="black"))

p3 <- ggplot() + geom_point(mapping=aes(y=-sensitivitySL_dat, x=LAI_matrix_nov)) +
  theme_classic() +  ylab("Sensitivity to shade quality") + xlab("Leaf area index") +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=14, colour="black"))

ggarrange(p1, p2, p3, nrow=1) # 1200 x 400


load(paste0(dir,"/Results/Te_shade_matrix_list_shadeavail_Aameiva_",resolution,"x",resolution,".RData"))
load(paste0(dir,"/Results/Te_sun_matrix_list_shadeavail_Aameiva_",resolution,"x",resolution,".RData"))

folder <- "C:/globalclimate"
global_climate <- raster::brick(paste0(folder, "/global_climate.nc"))
TMAXX_nov <- global_climate$X61 # T max in nov New et al.
TMAXX_nov <- crop(TMAXX_nov, c(-100, -30, -40, 30)) # Extent
TMAXX_nov <- aggregate(TMAXX_nov, fact=resolution/res(TMAXX_nov)[1]) 
TMAXX_nov <- mask(TMAXX_nov, tropical_forest0)
TMAXX_nov <- mask(TMAXX_nov, Aameiva_range)
TMAXX_nov_dat <- raster::extract(TMAXX_nov, xy.values)/10

p1 <- ggplot() + geom_point(mapping=aes(y=sensitivityT_dat, x=TMAXX_nov_dat)) +
  theme_classic() + ylab("Sensitivity to air temperature") + xlab("Temperature") +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=14, colour="black")) 
  # geom_smooth(mapping=aes(y=sensitivityT_dat, x=Te_sun_mean_dat))

p2 <- ggplot() + geom_point(mapping=aes(y=sensitivitySA_dat, x=TMAXX_nov_dat)) +
  theme_classic() + ylab("Sensitivity to shade availability") + xlab("Leaf area index") +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=14, colour="black"))
  # geom_smooth(mapping=aes(y=sensitivitySA_dat, x=Te_sun_mean_dat))

p3 <- ggplot() + geom_point(mapping=aes(y=-sensitivitySL_dat, x=TMAXX_nov_dat)) +
  theme_classic() +  ylab("Sensitivity to shade quality") + xlab("Leaf area index") +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=14, colour="black"))
  # geom_smooth(mapping=aes(y=sensitivitySL_dat, x=Te_sun_mean_dat))


ggarrange(p1, p2, p3, nrow=1) # 1200 x 400


##### Fig. 5 ##### 

deltaT_dat <- raster::extract(deltaT_map, xy.values)

# buffer_dat[212] <- 0

background <- data.frame(x=xy_southamerica[,1], y=xy_southamerica[,2],vals=vals)
ggplot() + geom_tile(background, mapping=aes(x=x, y=y,fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(mapping=aes(x=xy.values[,1], y=xy.values[,2], fill=buffer_dat)) +
  scale_fill_viridis_c()+
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     strip.text.x = element_blank())  # 600x550


ggplot() + geom_tile(background, mapping=aes(x=x, y=y,fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(mapping=aes(x=xy.values[,1], y=xy.values[,2], fill=deltaT_dat)) +
  scale_fill_viridis_c()+
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     strip.text.x = element_blank())  # 600x550

# buffer_dat[212] <- 4
deltaT_diff <- buffer_dat - deltaT_dat
positive <- which(deltaT_diff>0)

ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  #Positive values
  new_scale_fill() +
  geom_tile(aes(x=xy.values[positive,1], y=xy.values[positive,2],fill=deltaT_diff[positive])) +
  scale_fill_gradient(low="#F7FBFF", high="#084594")+
  #Negative values
  new_scale_fill() +
  geom_tile(aes(x=xy.values[-positive,1], y=xy.values[-positive,2],fill=deltaT_diff[-positive])) +
  scale_fill_gradient(low="#CB181D", high="#FFF5F0")+
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank()) # 600x550

















