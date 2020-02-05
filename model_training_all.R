setwd(getwd()) 


library(devtools)  
library(randomForest)
library(hydroGOF)
library(chillR)
library(caret)
library(quantregForest)
library(imputeTS)
library(raster)

######################################
#### Input data

all<-read.csv("all_new.csv", header = T)
summary(all)


cali<-all[all$Validation=="Training",]
vali<-all[all$Validation=="Validation",]

summary(cali)
summary(vali)


################ Fit QRF

fitControl = trainControl(## 5-fold CV
  method = "cv", 
  number = 5)


names(cali)


qrf = caret::train(x = cali[, c(9:31)],
                   y = cali[,32],
                   na.action = na.omit,
                   trControl = fitControl, 
                   method="qrf")

summary(qrf$finalModel)

model<-qrf$finalModel

varImpPlot(model)

varImpPlot.qrf(model)

qrf_importance<-qrf$finalMode$importance
write.csv(qrf_importance, "qrf_importance.csv", row.names = T)



################# Model statistics

soil_cali_mean  = predict(qrf$finalModel, newdata = cali[, c(9:31)], what = mean)
soil_vali_mean  = predict(qrf$finalModel, newdata = vali[, c(9:31)], what = mean)




sqrt(mean((soil_cali_mean - cali$SSM)^2))  # RMSE.c 
mean(soil_cali_mean - cali$SSM)   # ME.c
cor(soil_cali_mean, cali$SSM)^2   # R2.c


sqrt(mean((vali[!is.na(soil_vali_mean),]$SSM - soil_vali_mean[!is.na(soil_vali_mean)])^2, na.rm = T))   # RMSE.v 
mean(vali$SSM - soil_vali_mean, na.rm = T)    # ME.v
cor(vali[!is.na(soil_vali_mean),]$SSM , soil_vali_mean[!is.na(soil_vali_mean)])^2   # R2.v



sqrt(mean((cali$SMAP - cali$SSM)^2))  # RMSE.c 
mean(cali$SMAP - cali$SSM)   # ME.c
cor(cali$SMAP, cali$SSM)^2   # R2.c
RPD(cali$SMAP, cali$SSM)
RPIQ(cali$SMAP, cali$SSM)
NSE(cali$SMAP, cali$SSM)





sqrt(mean((vali$SMAP - vali$SSM)^2))  # RMSE.v 
mean(vali$SSM - vali$SMAP)    # ME.v
cor(vali$SSM , vali$SMAP)^2   # R2.v
RPD(vali$SMAP, vali$SSM)
RPIQ(vali$SMAP, vali$SSM)
NSE(vali$SMAP, vali$SSM)


plot( cali$SSM,  soil_cali_mean, xlim = c(0, 0.7), ylim = c(0, 0.7))
abline(0,1)

plot( vali$SSM,  soil_vali_mean, xlim = c(0, 0.7), ylim = c(0, 0.7))
abline(0,1)


#################  Output model predictions at the stations

soil_cali_sd  = predict(qrf$finalModel, newdata = cali[, c(9:31)], what = sd)
soil_vali_sd  = predict(qrf$finalModel, newdata = vali[, c(9:31)], what = sd)

soil_cali_upper  = predict(qrf$finalModel, newdata = cali[, c(9:31)], what = c(0.95)) 
soil_vali_upper  = predict(qrf$finalModel, newdata = vali[, c(9:31)], what = c(0.95)) 

soil_cali_lower  = predict(qrf$finalModel, newdata = cali[, c(9:31)], what = c(0.05)) 
soil_vali_lower  = predict(qrf$finalModel, newdata = vali[, c(9:31)], what = c(0.05)) 



out_cali<-cbind(cali, soil_cali_mean, soil_cali_sd, soil_cali_upper, soil_cali_lower )
out_vali<-cbind(vali, soil_vali_mean, soil_vali_sd, soil_vali_upper, soil_vali_lower )

write.csv(out_cali, "out_cali.csv", row.names = F)

write.csv(out_vali, "out_vali.csv", row.names = F)




################### Mapping onto a farm in Australia


library(raster)
require(rgdal)
library(sf)
library(ggplot2)
library(plyr)
library(stringr)
library(reshape)
library(matrixStats)
library(smapr)
library(sp)
library(raster)
require(rgdal)
library(FedData)


############### SMAP


WD_poly<-c("../Validation/All_stations/")



library(stringr)
require(raster)
require(rgdal)
library(plyr)
library(stringr)
library(reshape)
library(matrixStats)
library(imputeTS)
   
   
   WD_ancillary <- "../Validation/ancillary_validation/"
   WD_LC <- "../Validation/LC_validation/"
   WD_Sentinel_angle <-"../Validation/Sentinel/angle_validation/"
   WD_Sentinel_vv <-"../Validation/Sentinel/vv_validation/"
   WD_Sentinel_vh <-"../Validation/Sentinel/vh_validation/"
   WD_Sentinel_dates <-"../Validation/Sentinel/dates_validation/"
   WD_SMAP <- "../Validation/SMAP_validation/"
   
   out_all_vv<-NULL
   out_all_vh<-NULL
   out_all_angle<-NULL
   out_all_ancillary<-NULL
   out_all_LAI<-NULL
   out_all_SMAP<-NULL
   
   
   
   
   ROI_names<-list.files(WD_poly, pattern = ".shp", full.names = T) 
   sites_subset<- shapefile(ROI_names[1])
   
   
   
   # Make a grid
   # Lambert Conformal Conic
   sites_proj <- spTransform(sites_subset, CRS("+proj=utm +zone=55+datum=WGS84"))
   crs(sites_proj) 
   
   grid <- makegrid(sites_proj, cellsize = 100)
   grid.sp <- SpatialPoints(grid, proj4string = CRS(proj4string(sites_proj)))
   
   grid.sp <- grid.sp[!is.na(over(grid.sp,sites_proj)),]
   
   coordinates<-as.data.frame(grid.sp)
   
   grid.sp <- spTransform(grid.sp, proj4string(sites_subset))

   # SMAP
   smap_names<-list.files(WD_SMAP, pattern = ".tif", full.names = T)
   smap_names_2<-list.files(WD_SMAP, pattern = ".tif", full.names = F)
   smap_raster<-lapply(smap_names, raster)
   smap_dates<-substr(smap_names_2, 6, 13)
   

   
   
   
   # plot(grid.sp) 
   
   
   grid_smap.sp <- spTransform(grid.sp,   crs(smap_raster[[1]]))
   extract.smap<-as.data.frame(do.call("cbind",  lapply( smap_raster, extract,  grid_smap.sp)))
   colnames(extract.smap)<- smap_dates
   
   SMAP_SM_filled<-extract.smap
   
   
   for (i in 1:nrow(extract.smap))
     
   {
     if (sum(!is.na.data.frame(as.numeric(extract.smap[i,])))>2)
     {
       SMAP_SM_filled[i,] <- na_ma(as.numeric(extract.smap[i,]), k = 3, weighting = "simple")
       
     }
     
     else
     {
       SMAP_SM_filled[i,] <- NA
       
     }
   }
   
   ROI_names<-"OZNET_1"
   site_names<-ROI_names
   
   ### Sentinel  
   names_Sentinel_vv<-list.files(WD_Sentinel_vv, pattern = site_names, full.names = T)
   Sentinel_vv_raster <- stack(names_Sentinel_vv)
   
   names_dates_vv<-list.files(WD_Sentinel_dates, pattern = site_names, full.names = T)
   dates_vv_raw<- as.character(read.csv(names_dates_vv, header = T)[,5])
   dates_vv_new<-substr(dates_vv_raw, 2, str_length(dates_vv_raw)-1)
   dates_vv_new_2<- strsplit(dates_vv_new, ", ")[[1]]
   dates_vv_final<-substr(dates_vv_new_2, 1, 10)
   dates_vv_final<-paste(substr(dates_vv_final, 1,4), substr(dates_vv_final, 6,7), substr(dates_vv_final, 9,10), sep = "")
   
   
   extract.Sentinel_vv<-as.data.frame(do.call("cbind", lapply(as.list(Sentinel_vv_raster), extract, grid.sp)))
   colnames(extract.Sentinel_vv)<- dates_vv_final
   
   extract.Sentinel_vv<-extract.Sentinel_vv[,colnames(extract.Sentinel_vv)%in%unique(dates_vv_final)]
   extract.Sentinel_vv<-extract.Sentinel_vv[,colnames(extract.Sentinel_vv)%in%unique(dates_vv_final)]
   
   extract.Sentinel_vv$vv_mean<-rowMeans2(as.matrix(extract.Sentinel_vv), na.rm = T)
   extract.Sentinel_vv$vv_sd<-rowSds(as.matrix(extract.Sentinel_vv), na.rm = T)
   extract.Sentinel_vv$vv_min<-rowMins(as.matrix(extract.Sentinel_vv), na.rm = T)
   extract.Sentinel_vv$vv_max<-rowMaxs(as.matrix(extract.Sentinel_vv), na.rm = T)
   
   names_Sentinel_vh<-list.files(WD_Sentinel_vh, pattern = site_names, full.names = T)
   Sentinel_vh_raster <- stack(names_Sentinel_vh)
   
   names_dates_vh<-list.files(WD_Sentinel_dates, pattern = site_names, full.names = T)
   dates_vh_raw<- as.character(read.csv(names_dates_vh, header = T)[,4])
   dates_vh_new<-substr(dates_vh_raw, 2, str_length(dates_vh_raw)-1)
   dates_vh_new_2<- strsplit(dates_vh_new, ", ")[[1]]
   dates_vh_final<-substr(dates_vh_new_2, 1, 10)
   dates_vh_final<-paste(substr(dates_vh_final, 1,4), substr(dates_vh_final, 6,7), substr(dates_vh_final, 9,10), sep = "")
   
   extract.Sentinel_vh<-as.data.frame(do.call("cbind", lapply(as.list(Sentinel_vh_raster), extract, grid.sp)))
   colnames(extract.Sentinel_vh)<- dates_vh_final
   
   ####### Remove repeated dates
   extract.Sentinel_vh<-extract.Sentinel_vh[,colnames(extract.Sentinel_vh)%in%unique(dates_vh_final)]
   extract.Sentinel_vh<-extract.Sentinel_vh[,colnames(extract.Sentinel_vh)%in%unique(dates_vh_final)]
   
   extract.Sentinel_vh$vh_mean<-rowMeans2(as.matrix(extract.Sentinel_vh), na.rm = T)
   extract.Sentinel_vh$vh_sd<-rowSds(as.matrix(extract.Sentinel_vh), na.rm = T)
   extract.Sentinel_vh$vh_min<-rowMins(as.matrix(extract.Sentinel_vh), na.rm = T)
   extract.Sentinel_vh$vh_max<-rowMaxs(as.matrix(extract.Sentinel_vh), na.rm = T)
   
   
   names_Sentinel_angle<-list.files(WD_Sentinel_angle, pattern = site_names, full.names = T)
   Sentinel_angle_raster <- stack(names_Sentinel_angle)
   
   
   names_dates_angle<-list.files(WD_Sentinel_dates, pattern = site_names, full.names = T)
   dates_angle_raw<- as.character(read.csv(names_dates_angle, header = T)[,3])
   dates_angle_new<-substr(dates_angle_raw, 2, str_length(dates_angle_raw)-1)
   dates_angle_new_2<- strsplit(dates_angle_new, ", ")[[1]]
   dates_angle_final<-substr(dates_angle_new_2, 1, 10)
   dates_angle_final<-paste(substr(dates_angle_final, 1,4), substr(dates_angle_final, 6,7), substr(dates_angle_final, 9,10), sep = "")
   
   
   Sentinel_angle_raster_list<-as.list(Sentinel_angle_raster)
   Sentinel_angle_raster_new<-list(NA)
   for (k in 1:length(dates_angle_final))
   {
     Sentinel_angle_raster_new[[k]] <- Sentinel_angle_raster_list[[k+1]]
     
   }
   
   extract.Sentinel_angle<-as.data.frame(do.call("cbind", lapply(Sentinel_angle_raster_new, extract, grid.sp)))
   colnames(extract.Sentinel_angle)<- dates_angle_final
   
   ####### Remove repeated dates
   extract.Sentinel_angle<-extract.Sentinel_angle[,colnames(extract.Sentinel_angle)%in%unique(dates_angle_final)]
   extract.Sentinel_angle<-extract.Sentinel_angle[,colnames(extract.Sentinel_angle)%in%unique(dates_angle_final)]
   
   
   
   
   ############  extract ancillary
   names_ancillary<-list.files(WD_ancillary, pattern = site_names, full.names = T)
   ancillary_raster <- stack(names_ancillary)
   extract.ancillary<-as.data.frame(do.call("cbind", lapply(as.list(ancillary_raster), extract, grid.sp)))
   colnames(extract.ancillary)<- c("elevation", "slope_0", "aspect_0", "hillshade","slope", "aspect", "profc", "tangc", "planc", "meanc",  "roughness", "clay", "sand", "BD", "SOC", "FC", "PWP")
   
   
   
   ######## buffer zones for soil properties
   if (sum(is.na(extract.ancillary$clay))>0)
   {extract.ancillary_soil<-as.data.frame(do.call("cbind", lapply(as.list(ancillary_raster), extract, grid.sp, buffer = 300, fun = mean)))
   colnames(extract.ancillary_soil)<- c("elevation", "slope_0", "aspect_0", "hillshade","slope", "aspect", "profc", "tangc", "planc", "meanc",  "roughness", "clay", "sand", "BD", "SOC", "FC", "PWP")
   extract.ancillary[is.na(extract.ancillary$clay),c(12:17)]<-extract.ancillary_soil[is.na(extract.ancillary$clay),c(12:17)]
   
   
   }
   
   
   extract.ancillary$BD<-extract.ancillary$BD/100.0
   extract.ancillary$SOC<-extract.ancillary$SOC/5.0  ######## 5 times g /kg #
   extract.ancillary$FC<-extract.ancillary$FC/100.0
   extract.ancillary$PWP<-extract.ancillary$PWP/100.0
   
   
   
   ########## LC
   
   names_LC<-list.files(WD_LC, pattern = site_names, full.names = T)
   LC_raster <- raster(names_LC)
   extract.LC<-as.data.frame(do.call("cbind", lapply(as.list(LC_raster), extract, grid.sp)))
   colnames(extract.LC)<- c("LC")
   
   extract.LC$LC2<-extract.LC$LC
   
   if (sum(extract.LC$LC2== 1|extract.LC$LC2== 2|extract.LC$LC2== 3|extract.LC$LC2== 4|extract.LC$LC2== 5))
   {
     extract.LC[extract.LC$LC2== 1|extract.LC$LC2== 2|extract.LC$LC2== 3|extract.LC$LC2== 4|extract.LC$LC2== 5, ]$LC2 <- "Forest"
   }
   
   if (sum(extract.LC$LC2== 6|extract.LC$LC2== 7))
   {
     extract.LC[extract.LC$LC2== 6|extract.LC$LC2== 7, ]$LC2 <- "Shrubland"
   }
   
   if (sum(extract.LC$LC2== 8|extract.LC$LC2== 9))
   {
     extract.LC[extract.LC$LC2== 8|extract.LC$LC2== 9, ]$LC2 <- "Savanna"
   }
   
   if (sum(extract.LC$LC2== 10))
   {
     extract.LC[extract.LC$LC2== 10, ]$LC2 <- "Grassland"
   }
   
   if (sum(extract.LC$LC2== 12|extract.LC$LC2== 14))
   {
     
     extract.LC[extract.LC$LC2== 12|extract.LC$LC2== 14, ]$LC2 <- "Cropland"
   }
   
   if (sum(extract.LC$LC2== 16))
   {
     extract.LC[extract.LC$LC2== 16, ]$LC2 <- "Barren"
   }
   
   
   
   extract.LC<-as.data.frame(extract.LC$LC2)
   colnames(extract.LC)<-"LC"
   
   
   
   
   ############### DEM
   
   WD_DEM <- "../dem/"
   
   names_DEM<- list.files(WD_DEM, pattern = "tif", full.names = T)
   DEM_raster <- stack(names_DEM)
   
   DEM_raster_list<-as.list(DEM_raster)
   
   extract.DEM<-as.data.frame(do.call("cbind", lapply(DEM_raster_list, extract, grid.sp)))
   colnames(extract.DEM)<- c("aspect","elevation","flowdir","slope","tpi","tri")
   
   
   
   ############## Combine
   
   
   
   
   out_static<-cbind(as.data.frame(grid.sp),extract.LC, extract.DEM, extract.ancillary[,12:17], extract.Sentinel_vv[,(ncol(extract.Sentinel_vv)-3):ncol(extract.Sentinel_vv)], extract.Sentinel_vh[,(ncol(extract.Sentinel_vh)-3):ncol(extract.Sentinel_vh)])
   
   dates_Sentinel<-colnames(extract.Sentinel_vh)[1:(ncol(extract.Sentinel_vh)-4)]
   dates_SMAP <-colnames(SMAP_SM_filled)
   dates_map <- intersect(dates_Sentinel, dates_SMAP )
   
   VWC_mean <- matrix(NA, nrow = nrow(SMAP_SM_filled), ncol = length(dates_map))
   colnames(VWC_mean)<-paste("VWC", dates_map, sep = "_")
   VWC_SD <- matrix(NA, nrow = nrow(SMAP_SM_filled), ncol = length(dates_map))
   colnames(VWC_SD)<-paste("VWC", dates_map, sep = "_")
   VWC_CI_lower <- matrix(NA, nrow = nrow(SMAP_SM_filled), ncol = length(dates_map))
   colnames(VWC_CI_lower)<-paste("VWC", dates_map, sep = "_")
   VWC_CI_upper <- matrix(NA, nrow = nrow(SMAP_SM_filled), ncol = length(dates_map))
   colnames(VWC_CI_upper)<-paste("VWC", dates_map, sep = "_")
   
   
   for ( j in 1:length(dates_map))
     
   {
     
     out_dynamic<-cbind(extract.Sentinel_vv[,colnames(extract.Sentinel_vv)==dates_map[j]], extract.Sentinel_vh[,colnames(extract.Sentinel_vh)==dates_map[j]], extract.Sentinel_angle[,colnames(extract.Sentinel_angle)==dates_map[j]], SMAP_SM_filled[,colnames(SMAP_SM_filled)==dates_map[j]])
     
     colnames(out_dynamic) <- c("vv","vh","angle","SMAP")
     
     out_covariate <- cbind( out_static, out_dynamic)
     
     na_index<- !is.na(out_covariate$aspect)&!is.na(out_covariate$clay)&!is.na(out_covariate$vv)&!is.na(out_covariate$vh)&!is.na(out_covariate$SMAP)
     
     VWC_mean[na_index,j] <- predict(qrf$finalModel, newdata = out_covariate[na_index, c(4, 6:27)], what = mean)
     VWC_SD[na_index,j] <- predict(qrf$finalModel, newdata = out_covariate[na_index, c(4, 6:27)], what = sd)
     VWC_CI_lower[na_index,j] <- predict(qrf$finalModel, newdata = out_covariate[na_index, c(4, 6:27)], what = c(0.05))
     VWC_CI_upper[na_index,j] <- predict(qrf$finalModel, newdata = out_covariate[na_index, c(4, 6:27)], what = c(0.95))
     
     
     VWC_raster <- rasterFromXYZ(cbind(coordinates, VWC_mean[,j]), res=c(100,100), crs =   crs(sites_proj))
     writeRaster(VWC_raster, paste("VWC_Oz_", dates_map[j], ".tif", sep = ""))
     
     
   }
   
   
   
   ##################### Output SSM maps in Australia
   
   out_mean <- cbind(coordinates, out_static,  VWC_mean)
   out_SD <- cbind(coordinates, out_static,  VWC_SD)
   out_CI_lower <- cbind(coordinates, out_static,  VWC_CI_lower)
   out_CI_upper <- cbind(coordinates, out_static,  VWC_CI_upper)
   
   
   
   
   write.table(out_mean, "out_Oz_mean.txt", row.names = F, sep = ",")
   write.table(out_SD, "out_Oz_SD.txt", row.names = F, sep = ",")
   write.table(out_CI_lower, "out_Oz_CI_lower.txt", row.names = F, sep = ",")
   write.table(out_CI_upper, "out_Oz_CI_upper.txt", row.names = F, sep = ",")
   
   Upper_raster_1 <- rasterFromXYZ(cbind(coordinates, VWC_CI_upper[,1]), res=c(100,100), crs =   crs(sites_proj))
   writeRaster(Upper_raster_1 , paste("VWC_upper_Oz_", dates_map[1], ".tif", sep = ""))
   Upper_raster_2 <- rasterFromXYZ(cbind(coordinates, VWC_CI_upper[,2]), res=c(100,100), crs =   crs(sites_proj))
   writeRaster(Upper_raster_2 , paste("VWC_upper_Oz_", dates_map[2], ".tif", sep = ""))
   Upper_raster_3 <- rasterFromXYZ(cbind(coordinates, VWC_CI_upper[,3]), res=c(100,100), crs =   crs(sites_proj))
   writeRaster(Upper_raster_3 , paste("VWC_upper_Oz_", dates_map[3], ".tif", sep = ""))
   Upper_raster_4 <- rasterFromXYZ(cbind(coordinates, VWC_CI_upper[,4]), res=c(100,100), crs =   crs(sites_proj))
   writeRaster(Upper_raster_4 , paste("VWC_upper_Oz_", dates_map[4], ".tif", sep = ""))
   Upper_raster_5 <- rasterFromXYZ(cbind(coordinates, VWC_CI_upper[,5]), res=c(100,100), crs =   crs(sites_proj))
   writeRaster(Upper_raster_5 , paste("VWC_upper_Oz_", dates_map[5], ".tif", sep = ""))
   
   
   Lower_raster_1 <- rasterFromXYZ(cbind(coordinates, VWC_CI_lower[,1]), res=c(100,100), crs =   crs(sites_proj))
   writeRaster(Lower_raster_1 , paste("VWC_lower_Oz_", dates_map[1], ".tif", sep = ""))
   Lower_raster_2 <- rasterFromXYZ(cbind(coordinates, VWC_CI_lower[,2]), res=c(100,100), crs =   crs(sites_proj))
   writeRaster(Lower_raster_2 , paste("VWC_lower_Oz_", dates_map[2], ".tif", sep = ""))
   Lower_raster_3 <- rasterFromXYZ(cbind(coordinates, VWC_CI_lower[,3]), res=c(100,100), crs =   crs(sites_proj))
   writeRaster(Lower_raster_3 , paste("VWC_lower_Oz_", dates_map[3], ".tif", sep = ""))
   Lower_raster_4 <- rasterFromXYZ(cbind(coordinates, VWC_CI_lower[,4]), res=c(100,100), crs =   crs(sites_proj))
   writeRaster(Lower_raster_4 , paste("VWC_lower_Oz_", dates_map[4], ".tif", sep = ""))
   Lower_raster_5 <- rasterFromXYZ(cbind(coordinates, VWC_CI_lower[,5]), res=c(100,100), crs =   crs(sites_proj))
   writeRaster(Lower_raster_5 , paste("VWC_lower_Oz_", dates_map[5], ".tif", sep = ""))
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   