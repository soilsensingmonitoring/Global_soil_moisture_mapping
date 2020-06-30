setwd(getwd()) 


library(plyr)
library(chillR)
library(caret)
library(quantregForest)
library(dynatopmodel)
library(mlbench)
library(imputeTS)
library(raster)
library(viridis)
library(dplyr)
library(magrittr)
library(randomForest)
library(earth)
library(quantregForest)


######################################
#### validation by LC

all<-read.csv("data.csv", header = T)

unique(all$name)

names(all)

all$landcover<-as.factor(all$landcover)

summary(all)

all2<-all[(!is.na(all$VWC))&(!is.na(all$elevation))&(!is.na(all$clay))&(!is.na(all$SMAP))&(!is.na(all$vv_10))&(!is.na(all$vv_30))&(!is.na(all$vv_100))&(!is.na(all$vv_500))&(!is.na(all$vv_1000))&(!is.na(all$vv_5000))&(!is.na(all$vv_10000))&(!is.na(all$vv_36000)),]

summary(all2)

names(all2[,17:126])

# preproc1 <-preProcess(all2[,17:71], method=c("range"))
# all2[,17:71] <- predict(preproc1, all2[,17:71])
#summary(all2)


############ Global model


cali<-all2[all2$Validation==0,]
vali<-all2[all2$Validation==1,]

summary(cali)
summary(vali)


fitControl = trainControl(## 5-fold CV
  method = "cv", 
  number = 5)



############# 100 m


names(cali)


index_100<-names(cali) %in% c("elevation", "slope_9", "aspect_9", "hillshade", "meanc","planc","profc","tangc","roughness", 
                              "clay","sand","SOC","BD", "SMAP", 
                              "vv_100", "vh_100", "angle_100", 
                              "Min.vv_100.", "Min.vh_100." ,"Mean.vv_100." , "Mean.vh_100.", "Max.vv_100.", "Max.vh_100.", "StdDev.vv_100.", "StdDev.vh_100.")


qrf_100 = caret::train(x = cali[,index_100],
                         y = cali$VWC,
                         na.action = na.omit,
                         trControl = fitControl, 
                         metric = "RMSE",
                         ntree = 50,
                         nodesize = 10,
                         method="qrf")

qrf_100

summary(qrf_100$finalModel)

model_100<-qrf_100$finalModel



model <- model_100



varImpPlot(model)
importance(model)

# qrf_importance<-qrf$finalMode$importance
# write.csv(qrf_importance, "qrf_importance.csv", row.names = T)


soil_cali_mean  = predict(model, newdata = cali[, index_100], what = mean)
soil_vali_mean  = predict(model, newdata = vali[, index_100], what = mean)




sqrt(mean((soil_cali_mean - cali$VWC)^2))  # RMSE.c 
mean(soil_cali_mean - cali$VWC)   # ME.c
cor(soil_cali_mean, cali$VWC)^2   # R2.c



sqrt(mean((model$predicted - cali$VWC)^2))  # RMSE.c 
mean(model$predicted - cali$VWC)   # ME.c
cor(model$predicted, cali$VWC)^2   # R2.c


sqrt(mean((vali[!is.na(soil_vali_mean),]$VWC - soil_vali_mean[!is.na(soil_vali_mean)])^2, na.rm = T))   # RMSE.v 
mean(vali$VWC - soil_vali_mean, na.rm = T)    # ME.v
cor(vali[!is.na(soil_vali_mean),]$VWC , soil_vali_mean[!is.na(soil_vali_mean)])^2   # R2.v



plot( cali$VWC,  soil_cali_mean, xlim = c(0, 0.7), ylim = c(0, 0.7))
abline(0,1)

plot( cali$VWC,  model$predicted, xlim = c(0, 0.7), ylim = c(0, 0.7))
abline(0,1)

plot( vali$VWC,  soil_vali_mean, xlim = c(0, 0.7), ylim = c(0, 0.7))
abline(0,1)



sqrt(mean((cali$SMAP - cali$VWC)^2))  # RMSE.c 
mean(cali$SMAP - cali$VWC)   # ME.c
cor(cali$SMAP, cali$VWC)^2   # R2.c
RPD(cali$SMAP, cali$VWC)
RPIQ(cali$SMAP, cali$VWC)
#NSE(cali$SMAP, cali$VWC)





sqrt(mean((vali$SMAP - vali$VWC)^2))  # RMSE.v 
mean(vali$VWC - vali$SMAP)    # ME.v
cor(vali$VWC , vali$SMAP)^2   # R2.v
RPD(vali$SMAP, vali$VWC)
RPIQ(vali$SMAP, vali$VWC)
#NSE(vali$SMAP, vali$VWC)


out_cali<-cali
out_vali<-vali
out_cali$VWC_mean<-soil_cali_mean
out_vali$VWC_mean<-soil_vali_mean
out_cali$VWC_cv<-model$predicted
out_vali$VWC_cv<-soil_vali_mean

out<-rbind(out_cali, out_vali)

write.csv(out, "out_100m_full.csv", row.names = F)



#### Drop covariates
{
  
  ############# Drop SMAP
  
  
  
  index_SMAP<-names(cali) %in% c("elevation", "slope_9", "aspect_9", "hillshade", "meanc","planc","profc","tangc","roughness", 
                                 "clay","sand","SOC","BD",
                                 "vv_100", "vh_100", "angle_100", 
                                 "Min.vv_100.", "Min.vh_100." ,"Mean.vv_100." , "Mean.vh_100.", "Max.vv_100.", "Max.vh_100.", "StdDev.vv_100.", "StdDev.vh_100.")
  
  qrf_100_SMAP = caret::train(x = cali[,index_SMAP],
                                y = cali$VWC,
                                na.action = na.omit,
                                trControl = fitControl, 
                                metric = "RMSE",
                                ntree = 50,
                                nodesize = 10,
                                method="qrf")
  
  qrf_100_SMAP
  
  summary(qrf_100_SMAP$finalModel)
  
  model_100_SMAP<-qrf_100_SMAP$finalModel
  
  model_SMAP<-model_100_SMAP
  
  
  
  varImpPlot(model_SMAP)
  importance(model_SMAP)
  
  # qrf_importance<-qrf$finalMode$importance
  # write.csv(qrf_importance, "qrf_importance.csv", row.names = T)
  
  
  soil_cali_mean_SMAP  = predict(model_SMAP, newdata = cali[, index_SMAP], what = mean)
  soil_vali_mean_SMAP  = predict(model_SMAP, newdata = vali[, index_SMAP], what = mean)
  
  
  
  
  sqrt(mean((soil_cali_mean_SMAP - cali$VWC)^2))  # RMSE.c 
  mean(soil_cali_mean_SMAP - cali$VWC)   # ME.c
  cor(soil_cali_mean_SMAP, cali$VWC)^2   # R2.c
  
  
  
  sqrt(mean((model_SMAP$predicted - cali$VWC)^2))  # RMSE.c 
  mean(model_SMAP$predicted - cali$VWC)   # ME.c
  cor(model_SMAP$predicted, cali$VWC)^2   # R2.c
  
  
  sqrt(mean((vali[!is.na(soil_vali_mean_SMAP),]$VWC - soil_vali_mean_SMAP[!is.na(soil_vali_mean_SMAP)])^2, na.rm = T))   # RMSE.v 
  mean(vali$VWC - soil_vali_mean_SMAP, na.rm = T)    # ME.v
  cor(vali[!is.na(soil_vali_mean_SMAP),]$VWC , soil_vali_mean_SMAP[!is.na(soil_vali_mean_SMAP)])^2   # R2.v
  
  
  
  plot( cali$VWC,  soil_cali_mean_SMAP, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  plot( cali$VWC,  model_SMAP$predicted, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  plot( vali$VWC,  soil_vali_mean_SMAP, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  
  
  out_cali_SMAP<-cali
  out_vali_SMAP<-vali
  out_cali_SMAP$VWC_mean<-soil_cali_mean_SMAP
  out_vali_SMAP$VWC_mean<-soil_vali_mean_SMAP
  out_cali_SMAP$VWC_cv<-model_SMAP$predicted
  out_vali_SMAP$VWC_cv<-soil_vali_mean_SMAP
  
  out_SMAP<-rbind(out_cali_SMAP, out_vali_SMAP)
  
  write.csv(out_SMAP, "out_LC_full_Drop_SMAP.csv", row.names = F)
   
  
  
  ############# Drop Sentinel
  
  
  
  index_S1<-names(cali) %in% c("elevation", "slope_9", "aspect_9", "hillshade", "meanc","planc","profc","tangc","roughness", 
                               "clay","sand","SOC","BD",
                               "SMAP" )
  
  qrf_100_S1 = caret::train(x = cali[,index_S1],
                              y = cali$VWC,
                              na.action = na.omit,
                              trControl = fitControl, 
                              metric = "RMSE",
                              ntree = 50,
                              nodesize = 10,
                              method="qrf")
  
  qrf_100_S1
  
  summary(qrf_100_S1$finalModel)
  
  model_100_S1<-qrf_100_S1$finalModel
  
  model_S1<-model_100_S1
  
  
  
  varImpPlot(model_S1)
  importance(model_S1)
  
  # qrf_importance<-qrf$finalMode$importance
  # write.csv(qrf_importance, "qrf_importance.csv", row.names = T)
  
  
  soil_cali_mean_S1  = predict(model_S1, newdata = cali[, index_S1], what = mean)
  soil_vali_mean_S1  = predict(model_S1, newdata = vali[, index_S1], what = mean)
  
  
  
  
  sqrt(mean((soil_cali_mean_S1 - cali$VWC)^2))  # RMSE.c 
  mean(soil_cali_mean_S1 - cali$VWC)   # ME.c
  cor(soil_cali_mean_S1, cali$VWC)^2   # R2.c
  
  
  
  sqrt(mean((model_S1$predicted - cali$VWC)^2))  # RMSE.c 
  mean(model_S1$predicted - cali$VWC)   # ME.c
  cor(model_S1$predicted, cali$VWC)^2   # R2.c
  
  
  sqrt(mean((vali[!is.na(soil_vali_mean_S1),]$VWC - soil_vali_mean_S1[!is.na(soil_vali_mean_S1)])^2, na.rm = T))   # RMSE.v 
  mean(vali$VWC - soil_vali_mean_S1, na.rm = T)    # ME.v
  cor(vali[!is.na(soil_vali_mean_S1),]$VWC , soil_vali_mean_S1[!is.na(soil_vali_mean_S1)])^2   # R2.v
  
  
  
  plot( cali$VWC,  soil_cali_mean_S1, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  plot( cali$VWC,  model_S1$predicted, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  plot( vali$VWC,  soil_vali_mean_S1, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  
  
  out_cali_S1<-cali
  out_vali_S1<-vali
  out_cali_S1$VWC_mean<-soil_cali_mean_S1
  out_vali_S1$VWC_mean<-soil_vali_mean_S1
  out_cali_S1$VWC_cv<-model_S1$predicted
  out_vali_S1$VWC_cv<-soil_vali_mean_S1
  
  out_S1<-rbind(out_cali_S1, out_vali_S1)
  
  write.csv(out_S1, "out_LC_full_Drop_S1.csv", row.names = F)
  
  
  
  
  
  
  
  
  ############# Drop terrain
  
  
  
  index_terrain<-names(cali) %in% c("clay","sand","SOC","BD",
                                    "SMAP", "vv_100", "vh_100", "angle_100", 
                                    "Min.vv_100.", "Min.vh_100." ,"Mean.vv_100." , "Mean.vh_100.", "Max.vv_100.", "Max.vh_100.", "StdDev.vv_100.", "StdDev.vh_100.")
  
  qrf_100_terrain = caret::train(x = cali[,index_terrain],
                                   y = cali$VWC,
                                   na.action = na.omit,
                                   trControl = fitControl, 
                                   metric = "RMSE",
                                   ntree = 50,
                                   nodesize = 10,
                                   method="qrf")
  
  qrf_100_terrain
  
  summary(qrf_100_terrain$finalModel)
  
  model_100_terrain<-qrf_100_terrain$finalModel
  
  model_terrain<-model_100_terrain
  
  
  
  varImpPlot(model_terrain)
  importance(model_terrain)
  
  # qrf_importance<-qrf$finalMode$importance
  # write.csv(qrf_importance, "qrf_importance.csv", row.names = T)
  
  
  soil_cali_mean_terrain  = predict(model_terrain, newdata = cali[, index_terrain], what = mean)
  soil_vali_mean_terrain  = predict(model_terrain, newdata = vali[, index_terrain], what = mean)
  
  
  
  
  sqrt(mean((soil_cali_mean_terrain - cali$VWC)^2))  # RMSE.c 
  mean(soil_cali_mean_terrain - cali$VWC)   # ME.c
  cor(soil_cali_mean_terrain, cali$VWC)^2   # R2.c
  
  
  
  sqrt(mean((model_terrain$predicted - cali$VWC)^2))  # RMSE.c 
  mean(model_terrain$predicted - cali$VWC)   # ME.c
  cor(model_terrain$predicted, cali$VWC)^2   # R2.c
  
  
  sqrt(mean((vali[!is.na(soil_vali_mean_terrain),]$VWC - soil_vali_mean_terrain[!is.na(soil_vali_mean_terrain)])^2, na.rm = T))   # RMSE.v 
  mean(vali$VWC - soil_vali_mean_terrain, na.rm = T)    # ME.v
  cor(vali[!is.na(soil_vali_mean_terrain),]$VWC , soil_vali_mean_terrain[!is.na(soil_vali_mean_terrain)])^2   # R2.v
  
  
  
  plot( cali$VWC,  soil_cali_mean_terrain, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  plot( cali$VWC,  model_terrain$predicted, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  plot( vali$VWC,  soil_vali_mean_terrain, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  
  
  out_cali_terrain<-cali
  out_vali_terrain<-vali
  out_cali_terrain$VWC_mean<-soil_cali_mean_terrain
  out_vali_terrain$VWC_mean<-soil_vali_mean_terrain
  out_cali_terrain$VWC_cv<-model_terrain$predicted
  out_vali_terrain$VWC_cv<-soil_vali_mean_terrain
  
  out_terrain<-rbind(out_cali_terrain, out_vali_terrain)
  
  write.csv(out_terrain, "out_LC_full_Drop_terrain.csv", row.names = F)
  
  
  
  
  
  
  
  
  ############# Drop soil
  
  
  
  index_soil<-names(cali) %in% c("elevation", "slope_9", "aspect_9", "hillshade", "meanc","planc","profc","tangc","roughness", 
                                 "SMAP", "vv_100", "vh_100", "angle_100", 
                                 "Min.vv_100.", "Min.vh_100." ,"Mean.vv_100." , "Mean.vh_100.", "Max.vv_100.", "Max.vh_100.", "StdDev.vv_100.", "StdDev.vh_100.")
  
  qrf_100_soil = caret::train(x = cali[,index_soil],
                                y = cali$VWC,
                                na.action = na.omit,
                                trControl = fitControl, 
                                metric = "RMSE",
                                ntree = 50,
                                nodesize = 10,
                                method="qrf")
  
  qrf_100_soil
  
  summary(qrf_100_soil$finalModel)
  
  model_100_soil<-qrf_100_soil$finalModel
  
  model_soil<-model_100_soil
  
  
  
  varImpPlot(model_soil)
  importance(model_soil)
  
  # qrf_importance<-qrf$finalMode$importance
  # write.csv(qrf_importance, "qrf_importance.csv", row.names = T)
  
  
  soil_cali_mean_soil  = predict(model_soil, newdata = cali[, index_soil], what = mean)
  soil_vali_mean_soil  = predict(model_soil, newdata = vali[, index_soil], what = mean)
  
  
  
  
  sqrt(mean((soil_cali_mean_soil - cali$VWC)^2))  # RMSE.c 
  mean(soil_cali_mean_soil - cali$VWC)   # ME.c
  cor(soil_cali_mean_soil, cali$VWC)^2   # R2.c
  
  
  
  sqrt(mean((model_soil$predicted - cali$VWC)^2))  # RMSE.c 
  mean(model_soil$predicted - cali$VWC)   # ME.c
  cor(model_soil$predicted, cali$VWC)^2   # R2.c
  
  
  sqrt(mean((vali[!is.na(soil_vali_mean_soil),]$VWC - soil_vali_mean_soil[!is.na(soil_vali_mean_soil)])^2, na.rm = T))   # RMSE.v 
  mean(vali$VWC - soil_vali_mean_soil, na.rm = T)    # ME.v
  cor(vali[!is.na(soil_vali_mean_soil),]$VWC , soil_vali_mean_soil[!is.na(soil_vali_mean_soil)])^2   # R2.v
  
  
  
  plot( cali$VWC,  soil_cali_mean_soil, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  plot( cali$VWC,  model_soil$predicted, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  plot( vali$VWC,  soil_vali_mean_soil, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  
  
  out_cali_soil<-cali
  out_vali_soil<-vali
  out_cali_soil$VWC_mean<-soil_cali_mean_soil
  out_vali_soil$VWC_mean<-soil_vali_mean_soil
  out_cali_soil$VWC_cv<-model_soil$predicted
  out_vali_soil$VWC_cv<-soil_vali_mean_soil
  
  out_soil<-rbind(out_cali_soil, out_vali_soil)
  
  write.csv(out_soil, "out_LC_full_Drop_soil.csv", row.names = F)
  
  
  
  
  
  
  
  
  
  ############# Drop SMAP + S1
  
  
  
  index_SMAP_S1<-names(cali) %in% c("elevation", "slope_9", "aspect_9", "hillshade", "meanc","planc","profc","tangc","roughness", 
                                 "clay","sand","SOC","BD")
  
  qrf_100_SMAP_S1 = caret::train(x = cali[,index_SMAP_S1],
                              y = cali$VWC,
                              na.action = na.omit,
                              trControl = fitControl, 
                              metric = "RMSE",
                              ntree = 50,
                              nodesize = 10,
                              method="qrf")
  
  qrf_100_SMAP_S1
  
  summary(qrf_100_SMAP_S1$finalModel)
  
  model_100_SMAP_S1<-qrf_100_SMAP_S1$finalModel
  
  model_SMAP_S1<-model_100_SMAP_S1
  
  
  
  varImpPlot(model_SMAP_S1)
  importance(model_SMAP_S1)
  
  # qrf_importance<-qrf$finalMode$importance
  # write.csv(qrf_importance, "qrf_importance.csv", row.names = T)
  
  
  soil_cali_mean_SMAP_S1  = predict(model_SMAP_S1, newdata = cali[, index_SMAP_S1], what = mean)
  soil_vali_mean_SMAP_S1  = predict(model_SMAP_S1, newdata = vali[, index_SMAP_S1], what = mean)
  
  
  
  
  sqrt(mean((soil_cali_mean_SMAP_S1 - cali$VWC)^2))  # RMSE.c 
  mean(soil_cali_mean_SMAP_S1 - cali$VWC)   # ME.c
  cor(soil_cali_mean_SMAP_S1, cali$VWC)^2   # R2.c
  
  
  
  sqrt(mean((model_SMAP_S1$predicted - cali$VWC)^2))  # RMSE.c 
  mean(model_SMAP_S1$predicted - cali$VWC)   # ME.c
  cor(model_SMAP_S1$predicted, cali$VWC)^2   # R2.c
  
  
  sqrt(mean((vali[!is.na(soil_vali_mean_SMAP_S1),]$VWC - soil_vali_mean_SMAP_S1[!is.na(soil_vali_mean_SMAP_S1)])^2, na.rm = T))   # RMSE.v 
  mean(vali$VWC - soil_vali_mean_SMAP_S1, na.rm = T)    # ME.v
  cor(vali[!is.na(soil_vali_mean_SMAP_S1),]$VWC , soil_vali_mean_SMAP_S1[!is.na(soil_vali_mean_SMAP_S1)])^2   # R2.v
  
  
  
  plot( cali$VWC,  soil_cali_mean_SMAP_S1, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  plot( cali$VWC,  model_SMAP_S1$predicted, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  plot( vali$VWC,  soil_vali_mean_SMAP_S1, xlim = c(0, 0.7), ylim = c(0, 0.7))
  abline(0,1)
  
  
  
  out_cali_SMAP_S1<-cali
  out_vali_SMAP_S1<-vali
  out_cali_SMAP_S1$VWC_mean<-soil_cali_mean_SMAP_S1
  out_vali_SMAP_S1$VWC_mean<-soil_vali_mean_SMAP_S1
  out_cali_SMAP_S1$VWC_cv<-model_SMAP_S1$predicted
  out_vali_SMAP_S1$VWC_cv<-soil_vali_mean_SMAP_S1
  
  out_SMAP_S1<-rbind(out_cali_SMAP_S1, out_vali_SMAP_S1)
  
  write.csv(out_SMAP_S1, "out_LC_full_Drop_SMAP_S1.csv", row.names = F)
  
  
  
}










################# Scales  10 m
{
  
  
  index_10<-names(cali) %in% c("elevation", "slope_9", "aspect_9", "hillshade", "meanc","planc","profc","tangc","roughness", 
                                  "clay","sand","SOC","BD", "SMAP", 
                                  "vv_10", "vh_10", "angle_10", 
                                  "Min.vv_10.", "Min.vh_10." ,"Mean.vv_10." , "Mean.vh_10.", "Max.vv_10.", "Max.vh_10.", "StdDev.vv_10.", "StdDev.vh_10.")
  
  
  qrf_10 = caret::train(x = cali[,index_10],
                           y = cali$VWC,
                           na.action = na.omit,
                           trControl = fitControl, 
                           metric = "RMSE",
                           ntree = 50,
                           nodesize = 10,
                           method="qrf")
  
  qrf_10
  
  summary(qrf_10$finalModel)
  
  model_10<-qrf_10$finalModel
  
  
  
  
  
  
  model <- model_10
  
  
  
  varImpPlot(model)
  importance(model)
  
  # qrf_importance<-qrf$finalMode$importance
  # write.csv(qrf_importance, "qrf_importance.csv", row.names = T)
  
  
  soil_cali_mean  = predict(model, newdata = cali[, index_10], what = mean)
  soil_vali_mean  = predict(model, newdata = vali[, index_10], what = mean)
  
  
  
  
  sqrt(mean((soil_cali_mean - cali$VWC)^2))  # RMSE.c 
  mean(soil_cali_mean - cali$VWC)   # ME.c
  cor(soil_cali_mean, cali$VWC)^2   # R2.c
  
  
  
  sqrt(mean((model$predicted - cali$VWC)^2))  # RMSE.c 
  mean(model$predicted - cali$VWC)   # ME.c
  cor(model$predicted, cali$VWC)^2   # R2.c
  
  
  sqrt(mean((vali[!is.na(soil_vali_mean),]$VWC - soil_vali_mean[!is.na(soil_vali_mean)])^2, na.rm = T))   # RMSE.v 
  mean(vali$VWC - soil_vali_mean, na.rm = T)    # ME.v
  cor(vali[!is.na(soil_vali_mean),]$VWC , soil_vali_mean[!is.na(soil_vali_mean)])^2   # R2.v
  
  
  
}


################# Scales  30 m
{
  
  
  index_30<-names(cali) %in% c("elevation", "slope_9", "aspect_9", "hillshade", "meanc","planc","profc","tangc","roughness", 
                                "clay","sand","SOC","BD", "SMAP", 
                                "vv_30", "vh_30", "angle_30", 
                                "Min.vv_30.", "Min.vh_30." ,"Mean.vv_30." , "Mean.vh_30.", "Max.vv_30.", "Max.vh_30.", "StdDev.vv_30.", "StdDev.vh_30.")
  
  
  qrf_30 = caret::train(x = cali[,index_30],
                         y = cali$VWC,
                         na.action = na.omit,
                         trControl = fitControl, 
                         metric = "RMSE",
                         ntree = 50,
                         nodesize = 10,
                         method="qrf")
  
  qrf_30
  
  summary(qrf_30$finalModel)
  
  model_30<-qrf_30$finalModel
  
  
  
  
  
  
  model <- model_30
  
  
  
  varImpPlot(model)
  importance(model)
  
  # qrf_importance<-qrf$finalMode$importance
  # write.csv(qrf_importance, "qrf_importance.csv", row.names = T)
  
  
  soil_cali_mean  = predict(model, newdata = cali[, index_30], what = mean)
  soil_vali_mean  = predict(model, newdata = vali[, index_30], what = mean)
  
  
  
  
  sqrt(mean((soil_cali_mean - cali$VWC)^2))  # RMSE.c 
  mean(soil_cali_mean - cali$VWC)   # ME.c
  cor(soil_cali_mean, cali$VWC)^2   # R2.c
  
  
  
  sqrt(mean((model$predicted - cali$VWC)^2))  # RMSE.c 
  mean(model$predicted - cali$VWC)   # ME.c
  cor(model$predicted, cali$VWC)^2   # R2.c
  
  
  sqrt(mean((vali[!is.na(soil_vali_mean),]$VWC - soil_vali_mean[!is.na(soil_vali_mean)])^2, na.rm = T))   # RMSE.v 
  mean(vali$VWC - soil_vali_mean, na.rm = T)    # ME.v
  cor(vali[!is.na(soil_vali_mean),]$VWC , soil_vali_mean[!is.na(soil_vali_mean)])^2   # R2.v
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}









##################  By Countries


############ Global model


cali_2<-all2[all2$Validation2==0,]
vali_2<-all2[all2$Validation2==1,]

summary(cali_2)
summary(vali_2)



############# 100 m


names(cali_2)


index_100_Network<-names(cali_2) %in% c("elevation", "slope_9", "aspect_9", "hillshade", "meanc","planc","profc","tangc","roughness", 
                              "clay","sand","SOC","BD", "SMAP", 
                              "vv_100", "vh_100", "angle_100", 
                              "Min.vv_100.", "Min.vh_100." ,"Mean.vv_100." , "Mean.vh_100.", "Max.vv_100.", "Max.vh_100.", "StdDev.vv_100.", "StdDev.vh_100.")


qrf_100_Network = caret::train(x = cali_2[,index_100_Network],
                       y = cali_2$VWC,
                       na.action = na.omit,
                       trControl = fitControl, 
                       metric = "RMSE",
                       ntree = 50,
                       nodesize = 10,
                       method="qrf")

qrf_100_Network

summary(qrf_100_Network$finalModel)

model_100_Network<-qrf_100_Network$finalModel



model_Network <- model_100_Network



varImpPlot(model_Network)
importance(model_Network)

# qrf_importance<-qrf$finalMode$importance
# write.csv(qrf_importance, "qrf_importance.csv", row.names = T)


soil_cali_mean_Network  = predict(model_Network, newdata = cali_2[, index_100_Network], what = mean)
soil_vali_mean_Network  = predict(model_Network, newdata = vali_2[, index_100_Network], what = mean)




sqrt(mean((soil_cali_mean_Network - cali_2$VWC)^2))  # RMSE.c 
mean(soil_cali_mean_Network - cali_2$VWC)   # ME.c
cor(soil_cali_mean_Network, cali_2$VWC)^2   # R2.c



sqrt(mean((model_Network$predicted - cali_2$VWC)^2))  # RMSE.c 
mean(model_Network$predicted - cali_2$VWC)   # ME.c
cor(model_Network$predicted, cali_2$VWC)^2   # R2.c


sqrt(mean((vali_2[!is.na(soil_vali_mean_Network),]$VWC - soil_vali_mean_Network[!is.na(soil_vali_mean_Network)])^2, na.rm = T))   # RMSE.v 
mean(vali_2$VWC - soil_vali_mean_Network, na.rm = T)    # ME.v
cor(vali_2[!is.na(soil_vali_mean_Network),]$VWC , soil_vali_mean_Network[!is.na(soil_vali_mean_Network)])^2   # R2.v



plot( cali_2$VWC,  soil_cali_mean_Network, xlim = c(0, 0.7), ylim = c(0, 0.7))
abline(0,1)

plot( cali_2$VWC,  model_Network$predicted, xlim = c(0, 0.7), ylim = c(0, 0.7))
abline(0,1)

plot( vali_2$VWC,  soil_vali_mean_Network, xlim = c(0, 0.7), ylim = c(0, 0.7))
abline(0,1)



sqrt(mean((cali_2$SMAP - cali_2$VWC)^2))  # RMSE.c 
mean(cali_2$SMAP - cali_2$VWC)   # ME.c
cor(cali_2$SMAP, cali_2$VWC)^2   # R2.c
RPD(cali_2$SMAP, cali_2$VWC)
RPIQ(cali_2$SMAP, cali_2$VWC)
#NSE(cali$SMAP, cali$VWC)





sqrt(mean((vali_2$SMAP - vali_2$VWC)^2))  # RMSE.v 
mean(vali_2$VWC - vali_2$SMAP)    # ME.v
cor(vali_2$VWC , vali_2$SMAP)^2   # R2.v
RPD(vali_2$SMAP, vali_2$VWC)
RPIQ(vali_2$SMAP, vali_2$VWC)
#NSE(vali$SMAP, vali$VWC)


out_cali_Network<-cali_2
out_vali_Network<-vali_2
out_cali_Network$VWC_mean<-soil_cali_mean_Network
out_vali_Network$VWC_mean<-soil_vali_mean_Network
out_cali_Network$VWC_cv<-model_Network$predicted
out_vali_Network$VWC_cv<-soil_vali_mean_Network

out_Network<-rbind(out_cali_Network, out_vali_Network)

write.csv(out_Network, "out_100m_full_Network.csv", row.names = F)


  
  
  
  




######################################
#### Non-Spiking Wu

{
  
  Wu<-read.csv("data_Wu.csv", header = T)
  
  
  
  
  
  Wu2<-Wu[(!is.na(Wu$VWC))&(!is.na(Wu$elevation))&(!is.na(Wu$clay))&(!is.na(Wu$SMAP))&(!is.na(Wu$vv_10))&(!is.na(Wu$vv_30))&(!is.na(Wu$vv_100)),]
  
  
  spiking_grass_index<-unique(Wu2[Wu2$landcover=="Herbaceous", ]$name)
  spiking_forest_index<-unique(Wu2[Wu2$landcover=="Forest", ]$name)
  
  set.seed(123)
  index_1<-sample(spiking_grass_index, 3)
  
  set.seed(123)
  index_2<-sample(spiking_forest_index, 5)
  
  index<-c(index_1,index_2 )
  
  Wu2_spiking<-Wu2[Wu2$name%in%index,]
  Wu2_vali<-Wu2[!Wu2$name%in%index,]
  
  
  
  index_Wu_10<-names(Wu2_vali) %in% c("elevation", "slope_9", "aspect_9", "hillshade", "meanc","planc","profc","tangc","roughness", "clay","sand","SOC","BD", 
                                      "SMAP", 
                                      "vv_10", "vh_10", "angle_10", 
                                      "Min.vv_10.", "Min.vh_10." ,"Mean.vv_10." , "Mean.vh_10.", "Max.vv_10.", "Max.vh_10.", "StdDev.vv_10.", "StdDev.vh_10.")
  
  
  
  Wu2_10_mean  = predict(model_10, newdata = Wu2_vali[, index_Wu_10], what = mean)
  
  
  
  sqrt(mean((Wu2_10_mean - Wu2_vali$VWC)^2))  # RMSE.c 
  mean(Wu2_10_mean - Wu2_vali$VWC)   # ME.c
  cor(Wu2_10_mean, Wu2_vali$VWC)^2   # R2.c
  
  
  
  
  
  
  
  index_Wu_30<-names(Wu2_vali) %in% c("elevation", "slope_9", "aspect_9", "hillshade", "meanc","planc","profc","tangc","roughness", "clay","sand","SOC","BD", 
                                      "SMAP", 
                                      "vv_30", "vh_30", "angle_30", 
                                      "Min.vv_30.", "Min.vh_30." ,"Mean.vv_30." , "Mean.vh_30.", "Max.vv_30.", "Max.vh_30.", "StdDev.vv_30.", "StdDev.vh_30.")
  
  
  Wu2_30_mean  = predict(model_30, newdata = Wu2_vali[, index_Wu_30], what = mean)
  
  
  
  sqrt(mean((Wu2_30_mean - Wu2_vali$VWC)^2))  # RMSE.c 
  mean(Wu2_30_mean - Wu2_vali$VWC)   # ME.c
  cor(Wu2_30_mean, Wu2_vali$VWC)^2   # R2.c
  
  
  
  
  
  
  
  
  
  
  
  index_Wu_100<-names(Wu2_vali) %in% c("elevation", "slope_9", "aspect_9", "hillshade", "meanc","planc","profc","tangc","roughness", "clay","sand","SOC","BD", 
                                       "SMAP", 
                                       "vv_100", "vh_100", "angle_100", 
                                       "Min.vv_100.", "Min.vh_100." ,"Mean.vv_100." , "Mean.vh_100.", "Max.vv_100.", "Max.vh_100.", "StdDev.vv_100.", "StdDev.vh_100.")
  
  
  Wu2_100_mean  = predict(model_100, newdata = Wu2_vali[, index_Wu_100], what = mean)
  
  
  
  sqrt(mean((Wu2_100_mean - Wu2_vali$VWC)^2))  # RMSE.c 
  mean(Wu2_100_mean - Wu2_vali$VWC)   # ME.c
  cor(Wu2_100_mean, Wu2_vali$VWC)^2   # R2.c
  
  
  
  
  
  
  

  
  Wu2_SMAP_mean  = predict(model_SMAP, newdata = Wu2_vali, what = mean)
  
  
  
  sqrt(mean((Wu2_SMAP_mean - Wu2_vali$VWC)^2))  # RMSE.c 
  mean(Wu2_SMAP_mean - Wu2_vali$VWC)   # ME.c
  cor(Wu2_SMAP_mean, Wu2_vali$VWC)^2   # R2.c
  
  
  
  
  
  Wu2_soil_mean  = predict(model_soil, newdata = Wu2_vali, what = mean)
  
  
  
  sqrt(mean((Wu2_soil_mean - Wu2_vali$VWC)^2))  # RMSE.c 
  mean(Wu2_soil_mean - Wu2_vali$VWC)   # ME.c
  cor(Wu2_soil_mean, Wu2_vali$VWC)^2   # R2.c
  
  
  
  
  
  Wu2_terrain_mean  = predict(model_terrain, newdata = Wu2_vali, what = mean)
  
  
  
  sqrt(mean((Wu2_terrain_mean - Wu2_vali$VWC)^2))  # RMSE.c 
  mean(Wu2_terrain_mean - Wu2_vali$VWC)   # ME.c
  cor(Wu2_terrain_mean, Wu2_vali$VWC)^2   # R2.c
  
  
  
  
  
  Wu2_S1_mean  = predict(model_S1, newdata = Wu2_vali, what = mean)
  
  
  
  sqrt(mean((Wu2_S1_mean - Wu2_vali$VWC)^2))  # RMSE.c 
  mean(Wu2_S1_mean - Wu2_vali$VWC)   # ME.c
  cor(Wu2_S1_mean, Wu2_vali$VWC)^2   # R2.c
  
  
  
  
  Wu_out<-cbind(Wu2_vali, Wu2_10_mean, Wu2_30_mean , Wu2_100_mean,  Wu2_SMAP_mean, Wu2_S1_mean, Wu2_soil_mean, Wu2_terrain_mean)
  
  write.csv(Wu_out, "Wu2_out_non_spiking.csv", row.names = F)
  
  
  
  
  
}



################# Spiking

{
  
  
  index_10_spiking<-names(Wu2_spiking) %in% c("elevation", "slope_9", "aspect_9", "hillshade", "meanc","planc","profc","tangc","roughness", 
                               "clay","sand","SOC","BD", "SMAP", 
                               "vv_10", "vh_10", "angle_10", 
                               "Min.vv_10.", "Min.vh_10." ,"Mean.vv_10." , "Mean.vh_10.", "Max.vv_10.", "Max.vh_10.", "StdDev.vv_10.", "StdDev.vh_10.")
  
  
  cali_10_spiking <-rbind(cali[,index_10], Wu2_spiking[,index_10_spiking])
  
  
  
  qrf_10_spiking = caret::train(x = cali_10_spiking,
                        y = c(cali$VWC,Wu2_spiking$VWC),
                        na.action = na.omit,
                        trControl = fitControl, 
                        metric = "RMSE",
                        ntree = 50,
                        nodesize = 10,
                        method="qrf")
  
  qrf_10_spiking
  
  summary(qrf_10_spiking$finalModel)
  
  model_10_spiking<-qrf_10_spiking$finalModel
  
  
  
  
  
 
  Wu2_10_mean_spiking  = predict(model_10_spiking, newdata = Wu2_vali, what = mean)
   
  
  
  sqrt(mean((Wu2_10_mean_spiking - Wu2_vali$VWC)^2))  # RMSE.c 
  mean(Wu2_10_mean_spiking - Wu2_vali$VWC)   # ME.c
  cor(Wu2_10_mean_spiking, Wu2_vali$VWC)^2   # R2.c
  
  
  
  
  
  
  
  
  index_30_spiking<-names(Wu2_spiking) %in% c("elevation", "slope_9", "aspect_9", "hillshade", "meanc","planc","profc","tangc","roughness", 
                                              "clay","sand","SOC","BD", "SMAP", 
                                              "vv_30", "vh_30", "angle_30", 
                                              "Min.vv_30.", "Min.vh_30." ,"Mean.vv_30." , "Mean.vh_30.", "Max.vv_30.", "Max.vh_30.", "StdDev.vv_30.", "StdDev.vh_30.")
  
  
  cali_30_spiking <-rbind(cali[,index_30], Wu2_spiking[,index_30_spiking])
  
  
  
  qrf_30_spiking = caret::train(x = cali_30_spiking,
                                y = c(cali$VWC,Wu2_spiking$VWC),
                                na.action = na.omit,
                                trControl = fitControl, 
                                metric = "RMSE",
                                ntree = 50,
                                nodesize = 30,
                                method="qrf")
  
  qrf_30_spiking
  
  summary(qrf_30_spiking$finalModel)
  
  model_30_spiking<-qrf_30_spiking$finalModel
  
  
  
  
  
  
  Wu2_30_mean_spiking  = predict(model_30_spiking, newdata = Wu2_vali, what = mean)
  
  
  
  sqrt(mean((Wu2_30_mean_spiking - Wu2_vali$VWC)^2))  # RMSE.c 
  mean(Wu2_30_mean_spiking - Wu2_vali$VWC)   # ME.c
  cor(Wu2_30_mean_spiking, Wu2_vali$VWC)^2   # R2.c
  
  
  
  
  
  
  
  
  
  
  
  index_100_spiking<-names(Wu2_spiking) %in% c("elevation", "slope_9", "aspect_9", "hillshade", "meanc","planc","profc","tangc","roughness", 
                                              "clay","sand","SOC","BD", "SMAP", 
                                              "vv_100", "vh_100", "angle_100", 
                                              "Min.vv_100.", "Min.vh_100." ,"Mean.vv_100." , "Mean.vh_100.", "Max.vv_100.", "Max.vh_100.", "StdDev.vv_100.", "StdDev.vh_100.")
  
  
  cali_100_spiking <-rbind(cali[,index_100], Wu2_spiking[,index_100_spiking])
  
  
  
  qrf_100_spiking = caret::train(x = cali_100_spiking,
                                y = c(cali$VWC,Wu2_spiking$VWC),
                                na.action = na.omit,
                                trControl = fitControl, 
                                metric = "RMSE",
                                ntree = 50,
                                nodesize = 100,
                                method="qrf")
  
  qrf_100_spiking
  
  summary(qrf_100_spiking$finalModel)
  
  model_100_spiking<-qrf_100_spiking$finalModel
  
  

  
  Wu2_100_mean_spiking  = predict(model_100_spiking, newdata = Wu2_vali, what = mean)
  
  
  
  sqrt(mean((Wu2_100_mean_spiking - Wu2_vali$VWC)^2))  # RMSE.c 
  mean(Wu2_100_mean_spiking - Wu2_vali$VWC)   # ME.c
  cor(Wu2_100_mean_spiking, Wu2_vali$VWC)^2   # R2.c
  
  
  
  
  
  
  
  
  Wu_out_spiking<-cbind(Wu2_vali, Wu2_10_mean_spiking, Wu2_30_mean_spiking , Wu2_100_mean_spiking )
  
  write.csv(Wu_out_spiking, "Wu2_out_spiking.csv", row.names = F)
  
  
  
}
   
   
    
   
