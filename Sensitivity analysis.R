require(NicheMapR)
require(tidyr)
require(dplyr)
require(raster)
require(RColorBrewer)
require(ggplot2)
require(ggpubr)

dir <- ".../Sources" # set working directory
shapeTotal <- raster(paste0(dir,"/bio1.bil"))

resolution = 1 # Set resolution (deg) 0.5 or 1 (only sensitivity analyses)

world <- aggregate(shapeTotal, fact=resolution/res(shapeTotal)[1]) 

tropical_forest0 <- shapefile(paste0(dir,"/Shapes/htrop.shp"))
tropical_forest <- aggregate(rbind(tropical_forest0))
Aameiva_range <- shapefile(paste0("/Shapes/Aameiva_range.shp"))

Amazonia_map <- mask(world, Aameiva_range)
Amazonia_map <- crop(Amazonia_map, c(-100, -30, -40, 30)) # Extent

south_america <- crop(world, c(-100, -30, -40, 30)) # Extent

plot(Amazonia_map)

xy <- xyFromCell(Amazonia_map, 1:ncell(Amazonia_map))
cells <- cellFromPolygon(Amazonia_map, Aameiva_range)[[1]]
xy.values <- xy[cells,]
nrow(xy.values)

load(paste0(dir,"/xy.values_Aameiva.RData")) # only sensitivity analysis

######### Temperature increase #########  

parameter = "T36_"

load(paste0(dir,"/Results/meanTb_matrix_list_deltaT_",parameter,resolution,"x",resolution,".RData"))
load(paste0(dir,"/Results/SDTb_matrix_list_deltaT_",parameter,resolution,"x",resolution,".RData"))
load(paste0(dir,"/Results/Te_shade_matrix_list_deltaT_",parameter,resolution,"x",resolution,".RData"))
load(paste0(dir,"/Results/Te_sun_matrix_list_deltaT_",parameter,resolution,"x",resolution,".RData"))

deltaT <- seq(0, 8, length.out = 10)

sensitivityT_dat <- Teshade_dat <- Tesun_dat <- numeric(length(meanTb_matrix_list))
for(i in 1:length(sensitivityT_dat)){
  maxTshade <- which.max(Te_shade_matrix_list[[i]][,1])
  
  # meanTb_max <- meanTb_matrix_list[[i]][maxTshade,] - meanTb_matrix_list[[i]][maxTshade,1]
  meanTb_max <- colMeans(meanTb_matrix_list[[i]]) - mean(meanTb_matrix_list[[i]][,1])
  
  linear_model <- lm(meanTb_max ~ deltaT)
  
  sensitivityT_dat[i] <- coef(linear_model)[2]
  Teshade_dat[i] <- mean(Te_shade_matrix_list[[i]][,1])
  Tesun_dat[i] <- mean(Te_sun_matrix_list[[i]][,1])
}

######### Buffering capacity #########  

buffSD <- 0.2 # select the deltaT value for which TbSD is lower than this
buffer_dat <- numeric(length(meanTb_matrix_list))
for(i in 1:length(buffer_dat)){
  maxTshade <- which.max(Te_shade_matrix_list[[i]][,1])
  
  sdTb_max <- SDTb_matrix_list[[i]][maxTshade,]
  if(any(sdTb_max < buffSD)){
    maxTbuff <- min(deltaT[which.min(abs(sdTb_max - buffSD))])
    buffer_dat[i] <- maxTbuff
  }else{
    buffer_dat[i] <- 8
  }
}

require(ggplot2)

south_america <- south_america/south_america
xy_southamerica <- xyFromCell(south_america, 1:ncell(south_america))
vals<-raster::extract(south_america, xy_southamerica)
vals[which(!is.na(vals))]<-0
vals[which(is.na(vals))]<-1

ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
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

# 

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
                     strip.text.x = element_blank())  # 400x300


if(parameter == "M10_"){
  sensitivityT_dat_M10 <- sensitivityT_dat
  buffer_dat_M10 <- buffer_dat
}
if(parameter == "M30_"){
  sensitivityT_dat_M30 <- sensitivityT_dat
  buffer_dat_M30 <- buffer_dat
}
if(parameter == "a086_"){
  sensitivityT_dat_a086 <- sensitivityT_dat
  buffer_dat_a086 <- buffer_dat
}
if(parameter == "a090_"){
  sensitivityT_dat_a090 <- sensitivityT_dat
  buffer_dat_a090 <- buffer_dat
}
if(parameter == "s05_"){
  sensitivityT_dat_s05 <- sensitivityT_dat
  buffer_dat_s05 <- buffer_dat
}
if(parameter == "s025_"){
  sensitivityT_dat_s025 <- sensitivityT_dat
  buffer_dat_s025 <- buffer_dat
}
if(parameter == "T36_"){
  sensitivityT_dat_T36 <- sensitivityT_dat
  buffer_dat_T36 <- buffer_dat
}
if(parameter == "T32_"){
  sensitivityT_dat_T32 <- sensitivityT_dat
  buffer_dat_T32 <- buffer_dat
}
if(parameter == "M20_"){
  sensitivityT_dat_M20 <- sensitivityT_dat
  sensitivityT_dat_a088 <- sensitivityT_dat
  sensitivityT_dat_s1 <- sensitivityT_dat
  sensitivityT_dat_T34 <- sensitivityT_dat
  
  buffer_dat_M20 <- buffer_dat
  buffer_dat_a088 <- buffer_dat
  buffer_dat_s1 <- buffer_dat
  buffer_dat_T34 <- buffer_dat
}


##### Mass sensitivity #####

matrix_sensit_M <- data.frame(sensit_value = c(sensitivityT_dat_M10, sensitivityT_dat_M20, sensitivityT_dat_M30), 
                            Mass = sort(rep(c(10,20,30), 558)))
matrix_sensit_M <- matrix_sensit_M[which(matrix_sensit_M$sensit_value>0.4),]

ggplot(matrix_sensit_M, aes(y=sensit_value, x=as.factor((Mass)))) + geom_boxplot()

matrix_buff_M <- data.frame(buff_value = c(buffer_dat_M10, buffer_dat_M20, buffer_dat_M30), 
                              Mass = sort(rep(c(10,20,30), 558)))
ggplot(matrix_buff_M, aes(y=buff_value, x=as.factor((Mass)))) + geom_boxplot()

##### Absorbance sensitivity ##### 

matrix_sensit_a <- data.frame(sensit_value = c(sensitivityT_dat_a088, sensitivityT_dat_a086, sensitivityT_dat_a090), 
                            Abs = sort(rep(c(0.86,0.88,0.90), 558)))
matrix_sensit_a <- matrix_sensit_a[which(matrix_sensit_a$sensit_value>0.4),]

ggplot(matrix_sensit_a, aes(y=sensit_value, x=as.factor(Abs))) + geom_boxplot()


matrix_buff_a <- data.frame(buff_value = c(buffer_dat_a088, buffer_dat_a086, buffer_dat_a090), 
                              Abs = sort(rep(c(0.86,0.88,0.90), 558)))

ggplot(matrix_buff_a, aes(y=buff_value, x=as.factor(Abs))) + geom_boxplot()


##### Shade availability sensitivity ##### 

matrix_sensit_s <- data.frame(sensit_value = c(sensitivityT_dat_s1, sensitivityT_dat_s05, sensitivityT_dat_s025), 
                              Shade = sort(rep(c(1,05,025), 558)))
matrix_sensit_s <- matrix_sensit_s[which(matrix_sensit_s$sensit_value>0.4),]

ggplot(matrix_sensit_s, aes(y=sensit_value, x=as.factor(Shade))) + geom_boxplot()

matrix_buff_s <- data.frame(buff_value = c(buffer_dat_s1, buffer_dat_s05, buffer_dat_s025), 
                              Shade = sort(rep(c(1,05,025), 558)))

ggplot(matrix_buff_s, aes(y=buff_value, x=as.factor(Shade))) + geom_boxplot()


##### Tpref sensitivity ##### 

matrix_sensit_T <- data.frame(sensit_value = c(sensitivityT_dat_T32, sensitivityT_dat_T34, sensitivityT_dat_T36), 
                              Tpref = sort(rep(c(27,30,33), 558)))
matrix_sensit_T <- matrix_sensit_T[which(matrix_sensit_T$sensit_value>0.4),]

ggplot(matrix_sensit_T, aes(y=sensit_value, x=as.factor(Tpref))) + geom_boxplot()

matrix_buff_T <- data.frame(buff_value = c(buffer_dat_T32, buffer_dat_T34, buffer_dat_T36), 
                              Tpref = sort(rep(c(27,30,33), 558)))

ggplot(matrix_buff_T, aes(y=buff_value, x=as.factor(Tpref))) + geom_boxplot()

## maps

require(ggplot2)
require(ggpubr)

p1 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=buffer_dat_M10)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

p2 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=buffer_dat_M20)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

p3 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=buffer_dat_M30)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

diffM <- buffer_dat_M10 - buffer_dat_M30
diffM[1] <- -6; diffM[2] <- 6
diff1 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=diffM)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())
# 

p4 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=buffer_dat_a086)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

p5 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=buffer_dat_a088)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

p6 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=buffer_dat_a090)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

diffa <- buffer_dat_a086 - buffer_dat_a090
diffa[1] <- -6; diffa[2] <- 6
diff2 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=diffa)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

# 

p7 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=buffer_dat_s025)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

p8 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=buffer_dat_s05)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

p9 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=buffer_dat_s1)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

diffs <- buffer_dat_s025 - buffer_dat_s1
diffs[1] <- -6; diffs[2] <- 6
diff3 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=diffs)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

# 

p10 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=buffer_dat_T32)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

p11 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=buffer_dat_T34)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

p12 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=buffer_dat_T36)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

diffT <- buffer_dat_T32 - buffer_dat_T36
diffT[1] <- -6; diffT[2] <- 6
diff4 <- ggplot() + geom_tile(aes(x=xy_southamerica[,1], y=xy_southamerica[,2],fill=vals)) +
  scale_fill_gradient(low="grey90", high="white") +
  new_scale_fill() +
  geom_tile(aes(x=xy.values[,1], y=xy.values[,2],fill=diffT)) +
  scale_fill_viridis_c() +
  # geom_polygon(data = tropical_forest, aes(x = long, y = lat), colour="black", fill=NA) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = "right",
                     panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

ggarrange(p1, p2, p3, diff1,
          p4, p5, p6, diff2,
          p7, p8, p9, diff3,
          p10, p11, p12, diff4,
          nrow=4, ncol=4)



