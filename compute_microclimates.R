compute_microclimates <- function(run.gads = 0, lat, long, elev=NA, month, height, LAI, maxshade=NA, delta_T){
  
  ## Compute shade level (see details in Campbell and Norman 2012 and Algar et al. 2018 GEB)
  
  shade <- function(ac=0.8, z, LAI){ 
    Kbc <- sqrt((1 + tan(z*pi/180)^2) / 2.00132)
    omega_p <- exp(-sqrt(ac) * Kbc * LAI)
    return((1-omega_p)*100)
  }
  zenith <- function(sun, long, lat){ 
    long <- as.numeric(long); lat <- as.numeric(lat)
    rad <- pi/180
    hourAngle <- sun$solarTime + long - 180
    cosZenith <- (sin(rad * lat) * sun$sinSolarDec + cos(rad * lat) * sun$cosSolarDec * cos(rad * hourAngle))
    cosZenith[cosZenith > 1] <- 1
    cosZenith[cosZenith < -1] <- -1
    acos(cosZenith)/rad
  }
  
  location <- c(x=long, y=lat)
  
  if(!is.na(LAI)){
    micro_forZEN <- NicheMapR::micro_global(run.gads = 0, loc = location)
    micro_forZEN <- as.data.frame(micro_forZEN$metout)
    Jday = c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)[month] # transform into Julian day
    micro_forZEN = micro_forZEN %>% filter(DOY == Jday) # filter by Julian day
    ZEN <- min(micro_forZEN$ZEN)
    
    maxshade = min(shade(ac=0.8, z=ZEN, LAI=LAI))   # compute shade level from ZEN, LAI
  }
  # if(!is.na(LAI)){
  #   s <- GeoLight::solar(as.POSIXct(paste0("2000-",month,"-15"," ",1:24,":","00")))
  #   ZEN <- zenith(s, long=long, lat=lat)
  #   ZEN[ZEN>90] <- 90   # compute zenith angle value
  #   
  #   maxshade = min(shade(ac=0.8, z=ZEN, LAI=LAI))   # compute shade level using zenith angle (ZEN) and LAI
  # }
  
  ## Run NicheMapR
  micro <- NicheMapR::micro_global(run.gads = run.gads, loc = location, minshade = 0, maxshade = maxshade, Usrhyt = height, warm=delta_T)

  microclim_sun <- as.data.frame(micro$metout)    # microclimatic conditions in the sun
  microclim_shade <- as.data.frame(micro$shadmet) # and in the shade
  
  soil_sun <- as.data.frame(micro$soil)          # soil temperature in the sun
  soil_shade <- as.data.frame(micro$shadsoil)    # and in the shade
  
  ## Extract variables
  
  Ta_sun = microclim_sun$TALOC # Air temperature (<ba>C)
  Ta_shade = microclim_shade$TALOC
  
  Tg_sun = soil_sun$D0cm      # Soil surface temperature
  Tg_shade = soil_shade$D0cm
  
  S_sun = microclim_sun$SOLR  # Solar radiation (Wm-2)
  S_shade = microclim_shade$SOLR * (1 - maxshade/100) 
  
  V_sun = microclim_sun$VLOC # Wind velocity (ms-2)
  V_shade = microclim_shade$VLOC
  
  Tsky_sun = microclim_sun$TSKYC # Sky temperature (thermal radiation)
  Tsky_shade = microclim_shade$TSKYC
  
  DOY = microclim_sun$DOY # day of the year (we will use it later)
  ZEN = microclim_sun$ZEN # zenith angle
  
  microclimate_sun_df <- data.frame(DOY=DOY, ZEN=ZEN, Ta=Ta_sun, Tg=Tg_sun, S=S_sun, V=V_sun, Tsky=Tsky_sun)
  microclimate_shade_df <- data.frame(DOY=DOY, ZEN=ZEN, Ta=Ta_shade, Tg=Tg_shade, S=S_shade, V=V_shade, Tsky=Tsky_shade)
  
  Jday = c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)[month] # transform into Julian day
  
  microclimate_sun_df = microclimate_sun_df %>% filter(DOY == Jday) # filter by Julian day
  microclimate_shade_df = microclimate_shade_df %>% filter(DOY == Jday)
  
  # Temporal interpolation (increase resolution: hours -> minutes)
  
  microclimate_sun_min <- data.frame(array(NA, dim=c(60 * nrow(microclimate_sun_df), 7)))
  for(i in 2:7){
    var_hour = microclimate_sun_df[,i]
    var_hour_interp_fun = approxfun(var_hour ~ 1:nrow(microclimate_sun_df))
    microclimate_sun_min[,i] <- var_hour_interp_fun(seq(1, nrow(microclimate_sun_df), length.out = 60*nrow(microclimate_sun_df)))
  }
  microclimate_shade_min <- data.frame(array(NA, dim=c(60 * nrow(microclimate_sun_df), 7)))
  for(i in 2:7){
    var_hour = microclimate_shade_df[,i]
    var_hour_interp_fun = approxfun(var_hour ~ 1:nrow(microclimate_shade_df))
    microclimate_shade_min[,i] <- var_hour_interp_fun(seq(1, nrow(microclimate_shade_df), length.out = 60*nrow(microclimate_shade_df)))
  }
  microclimate_sun_min[,1] <- microclimate_shade_min[,1] <- sort(rep(microclimate_sun_df$DOY, 60))
  colnames(microclimate_sun_min) <- colnames(microclimate_shade_min) <- colnames(microclimate_shade_df)
  
  return(list(microclimate_sun_min=microclimate_sun_min, microclimate_shade_min=microclimate_shade_min))
}
