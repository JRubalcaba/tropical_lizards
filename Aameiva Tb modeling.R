require(NicheMapR)
require(tidyr)
require(dplyr)
require(raster)
require(Rcpp)
library(ncdf4) 
library(rgdal)
require(ggplot2)

dir <- ".../Sources" # set working directory

# 1. Load model functions
load(paste0(dir,"/Te_model.RData")) # Operative temperature model function
load(paste0(dir,"/Tb_t.RData")) # Transitent body temperature model
sourceCpp(paste0(dir,"/Cpp_functions.cpp")) # Probabilistic behavioral thermoregulation model (in Cpp)

# 2. Load the Leaf Area Index data and geographic coordinates of the tropical forest biome
load(paste0(dir,"/LAI/LAI_matrix_july_tropical_forest_HD.RData"))
load(paste0(dir,"/xy.values_HD.RData"))

# Only for sensitivity analysis (no HD)
load(paste0(dir,"/LAI/LAI_matrix_july_tropical_forest.RData"))
load(paste0(dir,"/xy.values.RData"))

# 3. Load the function compute microclimates (requires NicheMapR; Kearney and Porter 2017 Ecography) 
load(paste0(dir, "/compute_microclimates.RData"))

# 4. Load Ameiva ameiva Tb data (Cruz-Neto & Gordo 1996)
Tbdata <- read.csv(paste0(dir_data,"validation_Tb_data.csv"), sep=";")

Tbdata_Aameiva <- subset(Tbdata, Reference == "Cruz-Neto & Gordo 1996")

###### Location and traits ###### 

Tbdata_Aameiva <- subset(Tbdata, Reference == "Cruz-Neto & Gordo 1996")

location = c(x=-59.87527778, y=-2.42222222)
month = 11
height = 0.01
M = 21.19
a = 0.88
Tpref_mean = 34.34343 # 
Tpref_sd = 1 # 

sigma_sun_forest = 0.4444444
sigma_shade_forest = 0.3333333

sigma_sun_edge = 0.8989899 
sigma_shade_edge = 0.2727273

# LAI value

nc_data <- nc_open(paste0(dir,"/LAI/LAI_mean_monthly_1981-2015.nc4"))
lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)
t <- ncvar_get(nc_data, "time")
lai.array <- ncvar_get(nc_data, "LAI")
dim(lai.array)

if(month>=8){m=month-7}else{m=month+5} # 1: aug, 2:sept....
lai.slice <- lai.array[, , m]
r <- raster(t(lai.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
r <- flip(r, direction='y')

LAI_forest = raster::extract(r, cellFromXY(r, location))
LAI_edge = LAI_forest/2

microclimate_forest <- compute_microclimates(run.gads = 0, lat = location[2], long = location[1], month = month, height = height, LAI = LAI_forest, delta_T = 0)
microclimate_edge <- compute_microclimates(run.gads = 0, lat = location[2], long = location[1], month = month, height = height, LAI = LAI_edge, delta_T = 0)

S_sun_forest = microclimate_forest$microclimate_sun_min$S
Ta_sun_forest = microclimate_forest$microclimate_sun_min$Ta
Tsk_sun_forest = microclimate_forest$microclimate_sun_min$Tsky
Tg_sun_forest = microclimate_forest$microclimate_sun_min$Tg
V_sun_forest = microclimate_forest$microclimate_sun_min$V

S_shade_forest = microclimate_forest$microclimate_shade_min$S
Ta_shade_forest = microclimate_forest$microclimate_shade_min$Ta
Tsk_shade_forest = microclimate_forest$microclimate_shade_min$Tsky
Tg_shade_forest = microclimate_forest$microclimate_shade_min$Tg
V_shade_forest = microclimate_forest$microclimate_shade_min$V

S_sun_edge = microclimate_edge$microclimate_sun_min$S
Ta_sun_edge = microclimate_edge$microclimate_sun_min$Ta
Tsk_sun_edge = microclimate_edge$microclimate_sun_min$Tsky
Tg_sun_edge = microclimate_edge$microclimate_sun_min$Tg
V_sun_edge = microclimate_edge$microclimate_sun_min$V

S_shade_edge = microclimate_edge$microclimate_shade_min$S
Ta_shade_edge = microclimate_edge$microclimate_shade_min$Ta
Tsk_shade_edge = microclimate_edge$microclimate_shade_min$Tsky
Tg_shade_edge = microclimate_edge$microclimate_shade_min$Tg
V_shade_edge = microclimate_edge$microclimate_shade_min$V 

Te_sun_forest <- Te_model(M=M, a=a, ground=T, S=S_sun_forest, Tsk=Tsk_sun_forest, Ta=Ta_sun_forest, Tg=Tg_sun_forest, v=V_sun_forest)
Te_shade_forest <- Te_model(M=M, a=a, ground=T, S=S_shade_forest, Tsk=Tsk_shade_forest, Ta=Ta_shade_forest, Tg=Tg_shade_forest, v=V_shade_forest)

Te_sun_edge <- Te_model(M=M, a=a, ground=T, S=S_sun_edge, Tsk=Tsk_sun_edge, Ta=Ta_sun_edge, Tg=Tg_sun_edge, v=V_sun_edge)
Te_shade_edge <- Te_model(M=M, a=a, ground=T, S=S_shade_edge, Tsk=Tsk_shade_edge, Ta=Ta_shade_edge, Tg=Tg_shade_edge, v=V_shade_edge)

output_forest <- behav_therm(M=M, a=a, ground=1, S_sun = S_sun_forest, Ta_sun = Ta_sun_forest, Tsk_sun = Tsk_sun_forest, Tg_sun = Tg_sun_forest, V_sun = V_sun_forest,
                             S_shade = S_shade_forest, Ta_shade = Ta_shade_forest, Tsk_shade = Tsk_shade_forest, Tg_shade = Tg_shade_forest, V_shade = V_shade_forest,
                             Tpref_mean = Tpref_mean, Tpref_sd = Tpref_sd, sigma_sun = sigma_sun_forest, sigma_shade = sigma_shade_forest, max_iter = 100)

output_edge <- behav_therm(M=M, a=a, ground=1, S_sun = S_sun_edge, Ta_sun = Ta_sun_edge, Tsk_sun = Tsk_sun_edge, Tg_sun = Tg_sun_edge, V_sun = V_sun_edge,
                           S_shade = S_shade_edge, Ta_shade = Ta_shade_edge, Tsk_shade = Tsk_shade_edge, Tg_shade = Tg_shade_edge, V_shade = V_shade_edge,
                           Tpref_mean = Tpref_mean, Tpref_sd = Tpref_sd, sigma_sun = sigma_sun_edge, sigma_shade = sigma_shade_edge, max_iter = 100)

Tb_mean_forest <- output_forest[,1]
Tb_sd_forest <- output_forest[,3]

Tb_mean_edge <- output_edge[,1]
Tb_sd_edge <- output_edge[,3]

env <- Tbdata_Aameiva$Notes
meanTb_forest <- Tbdata_Aameiva$Mean[which(env == "FOREST")]
uppTb_forest <- Tbdata_Aameiva$uppSD[which(env == "FOREST")]
lowTb_forest <- Tbdata_Aameiva$lowSD[which(env == "FOREST")]
time_forest <- Tbdata_Aameiva$Hour[which(env == "FOREST")]

env <- Tbdata_Aameiva$Notes
meanTb_edge <- Tbdata_Aameiva$Mean[which(env == "EDGE")]
uppTb_edge <- Tbdata_Aameiva$uppSD[which(env == "EDGE")]
lowTb_edge <- Tbdata_Aameiva$lowSD[which(env == "EDGE")]
time_edge <- Tbdata_Aameiva$Hour[which(env == "EDGE")] 

Tb_plot_forest <- data.frame(Tb=Tb_mean_forest, Tb_sd=Tb_sd_forest, time=1:1440)
Tb_plot_edge <- data.frame(Tb=Tb_mean_edge, Tb_sd=Tb_sd_edge, time=1:1440)

# 800 x 600
ggplot() + theme_classic() + ylab("Temperature (ºC)") + xlab("Time (hour)") + 
  scale_x_continuous(breaks=seq(1,24,by=2), limits = c(8,18)) + 
  ylim(32, 41) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=12, colour="black"), legend.position = "none") +
  
  # FOREST - predicted
  geom_point(aes(y=Tb_mean_forest, x=1:1440/60), size=3, col=grey.colors(1, 0.5, alpha = 0.7)) +
  geom_ribbon(mapping = aes(ymin=Tb_mean_forest-Tb_sd_forest, ymax=Tb_mean_forest+Tb_sd_forest, x=1:1440/60), fill="grey50", alpha=0.2) +
  
  # EDGE - predicted 
  geom_point(aes(y=Tb_mean_edge, x=1:1440/60), size=3, col="goldenrod") +
  geom_ribbon(mapping = aes(ymin=Tb_mean_edge-Tb_sd_edge, ymax=Tb_mean_edge+Tb_sd_edge, x=1:1440/60), fill="goldenrod", alpha=0.2) +
  
  # FOREST - observed
  geom_errorbar(aes(x=time_forest-0.1, ymin=lowTb_forest, ymax=uppTb_forest), width = 0.4) +
  geom_point(aes(y=meanTb_forest, x=time_forest-0.1), fill="black", shape=21, size=5) +
  
  # EDGE - observed
  geom_errorbar(aes(x=time_edge+0.1, ymin=lowTb_edge, ymax=uppTb_edge), width = 0.4) +
  geom_point(aes(y=meanTb_edge, x=time_edge+0.1), fill="white", shape=21, size=5)


###### Simulating warming ###### 

delta_T <- 6

microclimate <- compute_microclimates(run.gads = 0, elev=elev, lat = location[2], long = location[1], month = month, 
                                      height = height, LAI = LAI_forest, delta_T = delta_T)

S_sun = microclimate$microclimate_sun_min$S
Ta_sun = microclimate$microclimate_sun_min$Ta
Tsk_sun = microclimate$microclimate_sun_min$Tsky
Tg_sun = microclimate$microclimate_sun_min$Tg
V_sun = microclimate$microclimate_sun_min$V

S_shade = microclimate$microclimate_shade_min$S
Ta_shade = microclimate$microclimate_shade_min$Ta
Tsk_shade = microclimate$microclimate_shade_min$Tsky
Tg_shade = microclimate$microclimate_shade_min$Tg
V_shade = microclimate$microclimate_shade_min$V

Te_sun <- Te_model(M=M, a=a, ground=T, S=S_sun, Tsk=Tsk_sun, Ta=Ta_sun, Tg=Tg_sun, v=V_sun)
Te_shade <- Te_model(M=M, a=a, ground=T, S=S_shade, Tsk=Tsk_shade, Ta=Ta_shade, Tg=Tg_shade, v=V_shade)

output_forest_warming <- behav_therm(M=M, a=a, ground=1, S_sun = S_sun, Ta_sun = Ta_sun, Tsk_sun = Tsk_sun, Tg_sun = Tg_sun, V_sun = V_sun,
                                     S_shade = S_shade, Ta_shade = Ta_shade, Tsk_shade = Tsk_shade, Tg_shade = Tg_shade, V_shade = V_shade,
                                     Tpref_mean = Tpref_mean, Tpref_sd = Tpref_sd, sigma_shade = sigma_shade_forest, sigma_sun = sigma_sun_forest, max_iter = 100)

Tb_mean_forest_warming <- output_forest_warming[,1]
Tb_sd_forest_warming <- output_forest_warming[,3]

Tb_mean_forest <- output_forest[,1]
Tb_sd_forest <- output_forest[,3]

ggplot() + theme_classic() + ylab("Temperature (ºC)") + xlab("Time (hour)") + scale_x_continuous(breaks=seq(1,24,by=2), limits = c(5,20)) + 
  ylim(20, 50) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=12, colour="black"), legend.position = "none") +
  # geom_point(mapping = aes(y=Te_sun, x=1:1440/60), size=2, colour="gold") +
  # geom_point(mapping = aes(y=Te_shade, x=1:1440/60), size=2, colour="grey30") +
  # geom_point(Tb_plot, mapping = aes(y=Tb, x=time/60), size=2, colour=grey.colors(1, 0.5, alpha = 0.05)) +
  geom_point(aes(y=Tb_mean_forest_warming, x=1:1440/60), col="goldenrod") +
  geom_ribbon(mapping = aes(ymin=Tb_mean_forest_warming-Tb_sd_forest_warming, ymax=Tb_mean_forest_warming+Tb_sd_forest_warming, x=1:1440/60), fill="goldenrod", alpha=0.3) +
  
  geom_point(aes(y=Tb_mean_forest, x=1:1440/60), col=grey.colors(1, 0.3, alpha = 0.1)) +
  geom_ribbon(mapping = aes(ymin=Tb_mean_forest-Tb_sd_forest, ymax=Tb_mean_forest+Tb_sd_forest, x=1:1440/60), fill="black", alpha=0.3)

###### Increasing air temp ###### 

CTmax=45

delta_T_seq <- seq(0, 10, length.out = 100)

Tb_mean_warming <- Tb_sd_warming <- pCTmax_warming <- numeric(length(delta_T_seq))
for(i in 1:length(delta_T_seq)){
  delta_T <- delta_T_seq[i]
  
  microclimate <- compute_microclimates(run.gads = 0, elev=elev, lat = location[2], long = location[1], month = month, height = height, LAI = LAI_forest, delta_T = delta_T)
  
  S_sun = microclimate$microclimate_sun_min$S
  Ta_sun = microclimate$microclimate_sun_min$Ta
  Tsk_sun = microclimate$microclimate_sun_min$Tsky
  Tg_sun = microclimate$microclimate_sun_min$Tg
  V_sun = microclimate$microclimate_sun_min$V
  
  S_shade = microclimate$microclimate_shade_min$S
  Ta_shade = microclimate$microclimate_shade_min$Ta
  Tsk_shade = microclimate$microclimate_shade_min$Tsky
  Tg_shade = microclimate$microclimate_shade_min$Tg
  V_shade = microclimate$microclimate_shade_min$V
  
  Te_sun <- Te_model(M=M, a=a, ground=T, S=S_sun, Tsk=Tsk_sun, Ta=Ta_sun, Tg=Tg_sun, v=V_sun)
  Te_shade <- Te_model(M=M, a=a, ground=T, S=S_shade, Tsk=Tsk_shade, Ta=Ta_shade, Tg=Tg_shade, v=V_shade)
  
  output_forest1 <- behav_therm(M=M, a=a, ground=1, S_sun = S_sun, Ta_sun = Ta_sun, Tsk_sun = Tsk_sun, Tg_sun = Tg_sun, V_sun = V_sun,
                               S_shade = S_shade, Ta_shade = Ta_shade, Tsk_shade = Tsk_shade, Tg_shade = Tg_shade, V_shade = V_shade,
                               Tpref_mean = Tpref_mean, Tpref_sd = Tpref_sd, sigma_shade = sigma_shade_forest, sigma_sun = sigma_sun_forest, max_iter = 100)
  
  daytimeTb_forest <- output_forest1[which(Te_sun > Te_shade),1]
  daytimeTb_SD_forest <- output_forest1[which(Te_sun > Te_shade),3]
  
  Tb_exc_CTmax <- output_forest1[,1][which(Te_sun > CTmax)]
  densityTb <- density(Tb_exc_CTmax)
  densityTb_funct <- approxfun(densityTb$y ~ densityTb$x, rule = 2)
  pCTmax <- integrate(densityTb_funct, lower = CTmax, upper = 100)[[1]]
  
  # focus on the time where Tb_sun is higher: higher risk, better model's performance
  Tb_mean_warming[i] <- output_forest1[which.max(Te_sun),1]
  Tb_sd_warming[i] <- output_forest1[which.max(Te_sun),3]
  pCTmax_warming[i] <- pCTmax
  
  # use the mean value during the daytime
  # Tb_mean_warming[i] <- mean(daytimeTb_forest)
  # Tb_sd_warming[i] <- mean(daytimeTb_SD_forest)
  # pCTmax_warming[i] <- pCTmax
  
  print(i / length(delta_T_seq))
}

# 400  x 350
ggplot(mapping=aes(y=Tb_mean_warming, x=delta_T_seq)) + theme_classic() +
  ylim(30, 50) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=12, colour="black"), legend.position = "none") +
  ylab("Mean body temperature (ºC)") + xlab("Increase in air temperature (ºC)") +
  geom_abline(intercept = 29.2, slope = 1) +
  geom_point(col=grey.colors(1, 0.3, alpha = 0.3)) +
  geom_ribbon(mapping = aes(ymin=Tb_mean_warming-Tb_sd_warming, ymax=Tb_mean_warming+Tb_sd_warming, x=delta_T_seq), fill="black", alpha=0.3)

ggplot(mapping=aes(y=pCTmax_warming, x=delta_T_seq)) + theme_classic() +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=12, colour="black"), legend.position = "none") +
  ylab("Prob (Tb > CTmax)") + xlab("Increase in air temperature (?C)") +
  geom_line(size=1.2) 

###### Lower shade availability ###### 

delta_T <- 0

microclimate <- compute_microclimates(run.gads = 0, elev=elev, lat = location[2], long = location[1], month = month, height = height, LAI = LAI_forest, delta_T = delta_T)

S_sun = microclimate$microclimate_sun_min$S
Ta_sun = microclimate$microclimate_sun_min$Ta
Tsk_sun = microclimate$microclimate_sun_min$Tsky
Tg_sun = microclimate$microclimate_sun_min$Tg
V_sun = microclimate$microclimate_sun_min$V

S_shade = microclimate$microclimate_shade_min$S
Ta_shade = microclimate$microclimate_shade_min$Ta
Tsk_shade = microclimate$microclimate_shade_min$Tsky
Tg_shade = microclimate$microclimate_shade_min$Tg
V_shade = microclimate$microclimate_shade_min$V

Te_sun <- Te_model(M=M, a=a, ground=T, S=S_sun, Tsk=Tsk_sun, Ta=Ta_sun, Tg=Tg_sun, v=V_sun)
Te_shade <- Te_model(M=M, a=a, ground=T, S=S_shade, Tsk=Tsk_shade, Ta=Ta_shade, Tg=Tg_shade, v=V_shade)

output_forest_lowAvail <- behav_therm(M=M, a=a, ground=1, S_sun = S_sun, Ta_sun = Ta_sun, Tsk_sun = Tsk_sun, Tg_sun = Tg_sun, V_sun = V_sun,
                                      S_shade = S_shade, Ta_shade = Ta_shade, Tsk_shade = Tsk_shade, Tg_shade = Tg_shade, V_shade = V_shade,
                                      Tpref_mean = Tpref_mean, Tpref_sd = Tpref_sd, sigma_shade = 0.1, sigma_sun = sigma_sun_forest, max_iter = 100)

Tb_mean_sun <- output_forest_lowAvail[,1]
Tb_sd_sun <- output_forest_lowAvail[,3]

Tb_mean_forest <- output_forest[,1]
Tb_sd_forest <- output_forest[,3]

ggplot() + theme_classic() + ylab("Temperature (ºC)") + xlab("Time (hour)") + scale_x_continuous(breaks=seq(1,24,by=2), limits = c(5,20)) + 
  ylim(20, 50) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=12, colour="black"), legend.position = "none") +
  # geom_point(mapping = aes(y=Te_sun, x=1:1440/60), size=2, colour="gold") +
  # geom_point(mapping = aes(y=Te_shade, x=1:1440/60), size=2, colour="grey30") +
  # geom_point(Tb_plot, mapping = aes(y=Tb, x=time/60), size=2, colour=grey.colors(1, 0.5, alpha = 0.05)) +
  geom_point(aes(y=Tb_mean_sun, x=1:1440/60), col="goldenrod") +
  geom_ribbon(mapping = aes(ymin=Tb_mean_sun-Tb_sd_sun, ymax=Tb_mean_sun+Tb_sd_sun, x=1:1440/60), fill="goldenrod", alpha=0.3) +
  
  geom_point(aes(y=Tb_mean_forest, x=1:1440/60), col=grey.colors(1, 0.3, alpha = 0.1)) +
  geom_ribbon(mapping = aes(ymin=Tb_mean_forest-Tb_sd_forest, ymax=Tb_mean_forest+Tb_sd_forest, x=1:1440/60), fill="black", alpha=0.3)

###### Reducing shade availability ###### 

sigma_sequence <- seq(0, 1, length.out = 100)

Tb_mean_lambda <- Tb_sd_lambda <- pCTmax_lambda <- numeric(length(sigma_sequence))
for(i in 1:length(sigma_sequence)){
  
  sigma <- sigma_sequence[i]
  
  microclimate <- compute_microclimates(run.gads = 0, elev=elev, lat = location[2], long = location[1], month = month, height = height, LAI = LAI_forest, delta_T=0)
                                        
  
  S_sun = microclimate$microclimate_sun_min$S
  Ta_sun = microclimate$microclimate_sun_min$Ta
  Tsk_sun = microclimate$microclimate_sun_min$Tsky
  Tg_sun = microclimate$microclimate_sun_min$Tg
  V_sun = microclimate$microclimate_sun_min$V
  
  S_shade = microclimate$microclimate_shade_min$S
  Ta_shade = microclimate$microclimate_shade_min$Ta
  Tsk_shade = microclimate$microclimate_shade_min$Tsky
  Tg_shade = microclimate$microclimate_shade_min$Tg
  V_shade = microclimate$microclimate_shade_min$V
  
  Te_sun <- Te_model(M=M, a=a, ground=T, S=S_sun, Tsk=Tsk_sun, Ta=Ta_sun, Tg=Tg_sun, v=V_sun)
  Te_shade <- Te_model(M=M, a=a, ground=T, S=S_shade, Tsk=Tsk_shade, Ta=Ta_shade, Tg=Tg_shade, v=V_shade)
  
  output_forest1 <- behav_therm(M=M, a=a, ground=1, S_sun = S_sun, Ta_sun = Ta_sun, Tsk_sun = Tsk_sun, Tg_sun = Tg_sun, V_sun = V_sun,
                               S_shade = S_shade, Ta_shade = Ta_shade, Tsk_shade = Tsk_shade, Tg_shade = Tg_shade, V_shade = V_shade,
                               Tpref_mean = Tpref_mean, Tpref_sd = Tpref_sd, sigma_shade = sigma, sigma_sun = sigma_sun_forest, max_iter = 100)
  
  daytimeTb_forest <- output_forest1[which(Te_sun > Te_shade),1]
  daytimeTb_SD_forest <- output_forest1[which(Te_sun > Te_shade),3]
  
  Tb_exc_CTmax <- output_forest1[,1][which(Te_sun > CTmax)]
  densityTb <- density(Tb_exc_CTmax)
  densityTb_funct <- approxfun(densityTb$y ~ densityTb$x, rule = 2)
  pCTmax <- integrate(densityTb_funct, lower = CTmax, upper = 100)[[1]]
  
  Tb_mean_lambda[i] <- output_forest1[which.max(Te_sun),1] 
  Tb_sd_lambda[i] <- output_forest1[which.max(Te_sun),3]
  pCTmax_lambda[i] <- pCTmax
  
  print(i / length(sigma_sequence))
}

ggplot(mapping=aes(y=Tb_mean_lambda, x=sigma_sequence)) + theme_classic() +
  ylim(30,50) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=12, colour="black"), legend.position = "none") +
  ylab("Mean body temperature (ºC)") + xlab("Shade availability") +
  geom_point(col=grey.colors(1, 0.3, alpha = 0.3)) +
  geom_ribbon(mapping = aes(ymin=Tb_mean_lambda-Tb_sd_lambda, ymax=Tb_mean_lambda+Tb_sd_lambda, x=sigma_sequence), fill="black", alpha=0.3)

ggplot(mapping=aes(y=pCTmax_lambda, x=sigma_sequence)) + theme_classic() +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=12, colour="black"), legend.position = "none") +
  ylab("Prob (Tb > CTmax)") + xlab("Shade availability") +
  geom_line(size=1.2) 

###### Lower shade quality ###### 

delta_T <- 0

LAI_correct <- LAI_forest * 0.3

microclimate <- compute_microclimates(run.gads = 0, elev=elev, lat = location[2], long = location[1], month = month, height = height, LAI = LAI_correct, delta_T=delta_T)
                                      

S_sun = microclimate$microclimate_sun_min$S
Ta_sun = microclimate$microclimate_sun_min$Ta
Tsk_sun = microclimate$microclimate_sun_min$Tsky
Tg_sun = microclimate$microclimate_sun_min$Tg
V_sun = microclimate$microclimate_sun_min$V

S_shade = microclimate$microclimate_shade_min$S
Ta_shade = microclimate$microclimate_shade_min$Ta
Tsk_shade = microclimate$microclimate_shade_min$Tsky
Tg_shade = microclimate$microclimate_shade_min$Tg
V_shade = microclimate$microclimate_shade_min$V

Te_sun <- Te_model(M=M, a=a, ground=T, S=S_sun, Tsk=Tsk_sun, Ta=Ta_sun, Tg=Tg_sun, v=V_sun)
Te_shade <- Te_model(M=M, a=a, ground=T, S=S_shade, Tsk=Tsk_shade, Ta=Ta_shade, Tg=Tg_shade, v=V_shade)

output_forest_lowLAI <- behav_therm(M=M, a=a, ground=1, S_sun = S_sun, Ta_sun = Ta_sun, Tsk_sun = Tsk_sun, Tg_sun = Tg_sun, V_sun = V_sun,
                                    S_shade = S_shade, Ta_shade = Ta_shade, Tsk_shade = Tsk_shade, Tg_shade = Tg_shade, V_shade = V_shade,
                                    Tpref_mean = Tpref_mean, Tpref_sd = Tpref_sd, sigma_sun = sigma_sun_forest, sigma_shade = sigma_shade_forest, max_iter = 100)

Tb_mean_sun <-output_forest_lowLAI[,1]
Tb_sd_sun <- output_forest_lowLAI[,3]

Tb_mean_forest <- output_forest[,1]
Tb_sd_forest <- output_forest[,3]

ggplot() + theme_classic() + ylab("Temperature (ºC)") + xlab("Time (hour)") + scale_x_continuous(breaks=seq(1,24,by=2), limits = c(5,20)) + 
  ylim(20, 50) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=12, colour="black"), legend.position = "none") +
  # geom_point(mapping = aes(y=Te_sun, x=1:1440/60), size=2, colour="gold") +
  # geom_point(mapping = aes(y=Te_shade, x=1:1440/60), size=2, colour="grey30") +
  # geom_point(Tb_plot, mapping = aes(y=Tb, x=time/60), size=2, colour=grey.colors(1, 0.5, alpha = 0.05)) +
  geom_point(aes(y=Tb_mean_sun, x=1:1440/60), col="goldenrod") +
  geom_ribbon(mapping = aes(ymin=Tb_mean_sun-Tb_sd_sun, ymax=Tb_mean_sun+Tb_sd_sun, x=1:1440/60), fill="goldenrod", alpha=0.3) +
  
  geom_point(aes(y=Tb_mean_forest, x=1:1440/60), col=grey.colors(1, 0.3, alpha = 0.1)) +
  geom_ribbon(mapping = aes(ymin=Tb_mean_forest-Tb_sd_forest, ymax=Tb_mean_forest+Tb_sd_forest, x=1:1440/60), fill="black", alpha=0.3)

###### Reducing shade quality ###### 

shadelev_sequence <- seq(0.2, 1, length.out = 100)

Tb_mean_shadelev <- Tb_sd_shadelev <- pCTmax_shadelev <- numeric(length(shadelev_sequence))
for(i in 1:length(shadelev_sequence)){
  
  LAI_correct <- LAI_forest * shadelev_sequence[i]
  
  microclimate <- compute_microclimates(run.gads = 0, elev=elev, lat = location[2], long = location[1], month = month, height = height, LAI = LAI_correct, delta_T = 0)
  
  S_sun = microclimate$microclimate_sun_min$S
  Ta_sun = microclimate$microclimate_sun_min$Ta
  Tsk_sun = microclimate$microclimate_sun_min$Tsky
  Tg_sun = microclimate$microclimate_sun_min$Tg
  V_sun = microclimate$microclimate_sun_min$V
  
  S_shade = microclimate$microclimate_shade_min$S
  Ta_shade = microclimate$microclimate_shade_min$Ta
  Tsk_shade = microclimate$microclimate_shade_min$Tsky
  Tg_shade = microclimate$microclimate_shade_min$Tg
  V_shade = microclimate$microclimate_shade_min$V
  
  Te_sun <- Te_model(M=M, a=a, ground=T, S=S_sun, Tsk=Tsk_sun, Ta=Ta_sun, Tg=Tg_sun, v=V_sun)
  Te_shade <- Te_model(M=M, a=a, ground=T, S=S_shade, Tsk=Tsk_shade, Ta=Ta_shade, Tg=Tg_shade, v=V_shade)
  
  output_forest1 <- behav_therm(M=M, a=a, ground=1, S_sun = S_sun, Ta_sun = Ta_sun, Tsk_sun = Tsk_sun, Tg_sun = Tg_sun, V_sun = V_sun,
                               S_shade = S_shade, Ta_shade = Ta_shade, Tsk_shade = Tsk_shade, Tg_shade = Tg_shade, V_shade = V_shade,
                               Tpref_mean = Tpref_mean, Tpref_sd = Tpref_sd, sigma_sun = sigma_sun_forest, sigma_shade = sigma_shade_forest, max_iter = 100)
  
  daytimeTb_forest <- output_forest1[which(Te_sun > Te_shade),1]
  daytimeTb_SD_forest <- output_forest1[which(Te_sun > Te_shade),3]
  
  Tb_exc_CTmax <- output_forest1[,1][which(Te_sun > CTmax)]
  densityTb <- density(Tb_exc_CTmax)
  densityTb_funct <- approxfun(densityTb$y ~ densityTb$x, rule = 2)
  pCTmax <- integrate(densityTb_funct, lower = CTmax, upper = 100)[[1]]
  
  Tb_mean_shadelev[i] <- output_forest1[which.max(Te_sun),1] 
  Tb_sd_shadelev[i] <- output_forest1[which.max(Te_sun),3]
  pCTmax_shadelev[i] <- pCTmax
  
  print(i / length(shadelev_sequence))
}

ggplot(mapping=aes(y=Tb_mean_shadelev, x=shadelev_sequence*100)) + theme_classic() +
  ylim(30,50) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=12, colour="black"), legend.position = "none") +
  ylab("Mean body temperature (ºC)") + xlab("Shade quality (%)") +
  geom_point(col=grey.colors(1, 0.3, alpha = 0.3)) +
  geom_ribbon(mapping = aes(ymin=Tb_mean_shadelev-Tb_sd_shadelev, ymax=Tb_mean_shadelev+Tb_sd_shadelev, x=shadelev_sequence*100), fill="black", alpha=0.3)

ggplot(mapping=aes(y=pCTmax_shadelev, x=shadelev_sequence*100)) + theme_classic() +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=12, colour="black"), legend.position = "none") +
  ylab("Prob (Tb > CTmax)") + xlab("Shade quality (%)") +
  geom_line(size=1.2) 


###### Exploring patterns #####

delta_T <- 0

microclimate <- compute_microclimates(run.gads = 0, lat = location[2], long = location[1], month = month, height = height, LAI = LAI_forest, delta_T = delta_T)
                                      

S_sun = microclimate$microclimate_sun_min$S
Ta_sun = microclimate$microclimate_sun_min$Ta
Tsk_sun = microclimate$microclimate_sun_min$Tsky
Tg_sun = microclimate$microclimate_sun_min$Tg
V_sun = microclimate$microclimate_sun_min$V

S_shade = microclimate$microclimate_shade_min$S
Ta_shade = microclimate$microclimate_shade_min$Ta
Tsk_shade = microclimate$microclimate_shade_min$Tsky
Tg_shade = microclimate$microclimate_shade_min$Tg
V_shade = microclimate$microclimate_shade_min$V

Te_sun <- Te_model(M=M, a=a, ground=T, S=S_sun, Tsk=Tsk_sun, Ta=Ta_sun, Tg=Tg_sun, v=V_sun)
Te_shade <- Te_model(M=M, a=a, ground=T, S=S_shade, Tsk=Tsk_shade, Ta=Ta_shade, Tg=Tg_shade, v=V_shade)

#

sigma = 0.3
time_exposed_forest_total <- 0
for(i in 1:500){
  output_forest <- behav_therm(M=M, a=a, ground=1, S_sun = S_sun, Ta_sun = Ta_sun, Tsk_sun = Tsk_sun, Tg_sun = Tg_sun, V_sun = V_sun,
                               S_shade = S_shade, Ta_shade = Ta_shade, Tsk_shade = Tsk_shade, Tg_shade = Tg_shade, V_shade = V_shade,
                               Tpref_mean = Tpref_mean, Tpref_sd = Tpref_sd, sigma_shade = sigma, sigma_sun = 1, max_iter = 1)
  
  Tb_forest <- output_forest[,1]
  sigma_forest <- output_forest[,2]
  
  aboveTpref_forest <- sigma_forest[which(Tb_forest > Tpref_mean)]
  
  time_exposed_forest <- 0
  for(i in 1:length(aboveTpref_forest)){
    if(aboveTpref_forest[i] == 1){
      time_exposed_forest[length(time_exposed_forest)] <- time_exposed_forest[length(time_exposed_forest)] + 1
    }
    if(aboveTpref_forest[i] == 0){
      time_exposed_forest <- c(time_exposed_forest, 0)
    }
  }
  time_exposed_forest_total <- c(time_exposed_forest_total,time_exposed_forest[-which(time_exposed_forest == 0)])
  
}

mean(time_exposed_forest_total)
1 / sigma

x <- seq(0, 20, length.out = 500)
y <- sigma * exp(-sigma * x)
ggplot() + geom_histogram(mapping=aes(y=stat(count) / sum(count), x=time_exposed_forest_total), bins=20) + theme_classic() +
  ylab("Density") + xlab("Time exposed (min)") + xlim(0, 20) + ylim(0,1) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=12, colour="black")) +
  geom_vline(xintercept = 1/sigma) +
  geom_line(mapping=aes(y=y, x=x), size=1) # 450 x 250



