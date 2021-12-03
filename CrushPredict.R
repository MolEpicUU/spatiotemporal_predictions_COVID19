## Spatio-temporal prediction of COVID-19 test positivity in Uppsala Län, Sweden ##
## Vera van Zoest (Uppsala University), 2021 ##

# set working directory
setwd("D:/CRUNCH/Data_analysis/R_analysis/Crush_predict")
rm(list=ls())
load("D:/CRUNCH/Data_analysis/R_analysis/Crush_predict/.RData")

# load libraries
library(lubridate)
library(INLA)
library(brinla)
library(rgdal)
library(spdep)

# define RMSE function
rmse <- function(x,y) sqrt(mean((x-y)^2, na.rm=T))

# predictions at service point level IMPORT EACH WEEK
data_sp <- as.data.frame(read.csv("prediction_ServicePoints_Vacc_incl_26_2021_4Aug.csv", head = T, sep =",", dec="."))
names(data_sp)[1] <- "id"
names(data_sp)[11] <- "age"
data_sp$sp_id <- as.numeric(as.factor(data_sp$ServicePoint_postalcode))
data_sp$sp_id2 <- data_sp$sp_id # INLA can't use same variable twice, so make copy of ID
data_sp$sp_id3 <- data_sp$sp_id
data_sp$week2 <- data_sp$week # copy temporal ID too
data_sp <- data_sp[data_sp$week != 26,] # remove week 26, 2020 
summary(data_sp)
# END IMPORT EACH WEEK

# import shapefile service point areas
sp_shp <- readOGR("ServicePoints_shape.shp", layer="ServicePoints_shape")
plot(sp_shp)
sp_ids_df <- data.frame("sp"=unique(data_sp$Service_Point_name), "sp_id"=unique(data_sp$sp_id))
sp_shp <- sp::merge(x=sp_shp, y=sp_ids_df, by.x="srvcPn_", by.y="sp", all.x=T, sort=F) 
table(sp_shp@data$sp_id) # no duplicates
nb_sp <- poly2nb(sp_shp, row.names=sp_shp@data$sp_id)


# run iteratively as new data comes in, predict next week
# includes density plot used in paper
IterPred <- function(data, nb_shp) {
  rmse_list <- c()
  week_ids <- unique(data[["week"]])
  pred_weeks <- week_ids[20:(length(week_ids)-1)] # minimum 20 training weeks + 1 prediction week
  for (w in pred_weeks) { 
    temp_data <- data[data[["week"]] <= w,] 
    pred_ids <- which(temp_data[["week"]]==w)
    temp_pos <- temp_data[["positivity_1w_nextweek"]]
    temp_val <- temp_pos[pred_ids]
    temp_pos[pred_ids] <- NA
    assign("temp_pos", value=temp_pos, envir=globalenv()) # INLA can only take y variable from global environment, not function environment
    formula <- temp_pos ~ 1 + aver_pos_neigh + test_per_capita_1w5 + calls_112_covid_perCapita + pop_density + workplaces_GM + top10_code + f(sp_id, model="bym", graph=nb_shp) + f(week, model="ar1") + f(week2, model="iid") + f(sp_id3, week, model="iid")
    inla_temp <- inla(formula, family="gaussian", data=temp_data, 
                      control.predictor=list(compute=TRUE, link=1), control.compute=list(dic=TRUE, cpo=TRUE))
    temp_rmse <- rmse(inla_temp$summary.fitted.values[[1]][pred_ids], temp_val) 
    rmse_list <- append(rmse_list, temp_rmse)
    if (w==max(pred_weeks)) {
      print(inla_temp$summary.random)
      par(mfrow=c(2,3))
      # plot(inla_temp$marginals.fixed[[1]], type="l", main=names(inla_temp$marginals.fixed[1]))
      # abline(v=0, lwd=1, lty=4, col="darkgrey")
      plot(inla_temp$marginals.fixed[[2]], type="l", main="Avg. positivity adjacent areas (t-1)", xlab=expression(~beta~"k"), ylab="Density")
      abline(v=0, lwd=1, lty=4, col="darkgrey")
      plot(inla_temp$marginals.fixed[[3]], type="l", main="Tests per 100,000 inhabitants (t-1)", xlab=expression(~beta~"k"), ylab="Density")
      abline(v=0, lwd=1, lty=4, col="darkgrey")
      plot(inla_temp$marginals.fixed[[4]], type="l", main="COVID-19 related calls to 112 (t-1)", xlab=expression(~beta~"k"), ylab="Density")
      abline(v=0, lwd=1, lty=4, col="darkgrey")
      plot(inla_temp$marginals.fixed[[5]], type="l", main="Population density", xlab=expression(~beta~"k"), ylab="Density")
      abline(v=0, lwd=1, lty=4, col="darkgrey")
      plot(inla_temp$marginals.fixed[[6]], type="l", main="Mobility to/from workplaces (t-1)", xlab=expression(~beta~"k"), ylab="Density")
      abline(v=0, lwd=1, lty=4, col="darkgrey")
      plot(inla_temp$marginals.fixed[[7]], type="l", main="Top10 highest positivity (t-1)", xlab=expression(~beta~"k"), ylab="Density")
      abline(v=0, lwd=1, lty=4, col="darkgrey")
      # legend(x="topright", legend = c("Intercept", 
      #                                 names(inla_temp$marginals.fixed[2]), 
      #                                 names(inla_temp$marginals.fixed[3]), 
      #                                 names(inla_temp$marginals.fixed[4]), 
      #                                 names(inla_temp$marginals.fixed[5]),
      #                                 names(inla_temp$marginals.fixed[6]),
      #                                 names(inla_temp$marginals.fixed[7])), 
      #        col=c(1,2,3,4,5,6,7), lty=1)
    }
    print(inla_temp$summary.fixed)
  }
  rmse_df <- data.frame("week"=pred_weeks+1, "rmse"=rmse_list)
  return(rmse_df)
}

# run function
rmse_df <- IterPred(data=data_sp, nb_shp=nb_sp)
rm(temp_pos)

# same function, save predictions instead of rmse 
IterPredSave <- function(data, nb_shp) {
  data_ori <- data
  data_sub <- data_ori[data_ori[["week"]] > unique(data_ori[["week"]])[19],]
  data_sub <- data_sub[data_sub[["week"]] < max(unique(data_sub[["week"]])),]
  data_sub$order <- c(1:length(data_sub[["id"]]))
  data_sub <- data_sub[order(data_sub$week),]
  #data <- data[order(data$week),]
  pred_list <- c()
  week_ids <- unique(data[["week"]])
  pred_weeks <- week_ids[20:(length(week_ids)-1)] # minimum 20 training weeks + 1 prediction week
  for (w in pred_weeks) { 
    temp_data <- data[data[["week"]] <= w,] 
    pred_ids <- which(temp_data[["week"]]==w)
    temp_pos <- temp_data[["positivity_1w_nextweek"]]
    temp_val <- temp_pos[pred_ids]
    temp_pos[pred_ids] <- NA
    assign("temp_pos", value=temp_pos, envir=globalenv()) # INLA can only take y variable from global environment, not function environment
    formula <- temp_pos ~ 1 + aver_pos_neigh + test_per_capita_1w5 + calls_112_covid_perCapita + pop_density + workplaces_GM + top10_code + f(sp_id, model="bym", graph=nb_shp) + f(week, model="ar1") + f(week2, model="iid") + f(sp_id3, week, model="iid")
    inla_temp <- inla(formula, family="gaussian", data=temp_data, 
                      control.predictor=list(compute=TRUE, link=1), control.compute=list(dic=TRUE, cpo=TRUE))
    temp_pred <- inla_temp$summary.fitted.values[[1]][pred_ids] 
    pred_list <- append(pred_list, temp_pred)
    if (w==max(pred_weeks)) {
      print(inla_temp$summary.random)
      par(mfrow=c(2,3))
      # plot(inla_temp$marginals.fixed[[1]], type="l", main=names(inla_temp$marginals.fixed[1]))
      # abline(v=0, lwd=1, lty=4, col="darkgrey")
      plot(inla_temp$marginals.fixed[[2]], type="l", main=names(inla_temp$marginals.fixed[2]), xlab=expression(~beta~"k"), ylab="Density")
      abline(v=0, lwd=1, lty=4, col="darkgrey")
      plot(inla_temp$marginals.fixed[[3]], type="l", main=names(inla_temp$marginals.fixed[3]), xlab=expression(~beta~"k"), ylab="Density")
      abline(v=0, lwd=1, lty=4, col="darkgrey")
      plot(inla_temp$marginals.fixed[[4]], type="l", main=names(inla_temp$marginals.fixed[4]), xlab=expression(~beta~"k"), ylab="Density")
      abline(v=0, lwd=1, lty=4, col="darkgrey")
      plot(inla_temp$marginals.fixed[[5]], type="l", main=names(inla_temp$marginals.fixed[5]), xlab=expression(~beta~"k"), ylab="Density")
      abline(v=0, lwd=1, lty=4, col="darkgrey")
      plot(inla_temp$marginals.fixed[[6]], type="l", main=names(inla_temp$marginals.fixed[6]), xlab=expression(~beta~"k"), ylab="Density")
      abline(v=0, lwd=1, lty=4, col="darkgrey")
      plot(inla_temp$marginals.fixed[[7]], type="l", main=names(inla_temp$marginals.fixed[7]), xlab=expression(~beta~"k"), ylab="Density")
      abline(v=0, lwd=1, lty=4, col="darkgrey")
    }
    print(inla_temp$summary.fixed)
  }
  pred_list <- pred_list[order(data_sub$order)]
  data_sub <- data_sub[order(data_sub$order),]
  #data <- data[order(data$order),]
  pred_df <- data.frame("week_now"=data_sub[["week"]], "predicted_week"=data_sub[["week"]]+1, "pred"=pred_list) 
  rm(temp_pos)
  return(pred_df)
}

pred_df <- IterPredSave(data=data_sp, nb_shp=nb_sp) 
data_sp_val <- data_sp[data_sp$week > 45,]
data_sp_val <- data_sp_val[data_sp_val$week < 79,]
pred_df <- cbind(data_sp_val[c(2,30,37,6,27)], pred_df)
names(pred_df)[8] <- "pred_INLA"

# import predictions Gradient Boosting (GB) and Random Forest (RF)
pred_RF_GB <- as.data.frame(read.csv("pred_RF_GB.csv", head = T, sep =",", dec="."))
data_sp$order <- c(1:length(data_sp$id)) # preserve old order dataframe
temp_data <- data_sp[order(data_sp$week),] # change order to match RF and GB file
temp_data <- temp_data[temp_data$week %in% c(31:79),] # keep weeks that match RF and GB file
temp_data <- cbind(temp_data, pred_RF_GB[c(4:5)]) # bind data
temp_data <- temp_data[order(temp_data$order),] # order back
temp_data <- temp_data[temp_data$week > 45,]
temp_data <- temp_data[temp_data$week < 79,]
pred_df <- cbind(pred_df, temp_data[c(42:43)])
rm(temp_data)

# ARIMA code
tsdf <- readRDS("D:/CRUNCH/Projects/Crush_covid/Anders/tsdf.rds")
train <- window(tsdf,start = c(2020,27),end= c(2020,45))
fit <- auto.arima(train)
refit <- Arima(tsdf, model=fit)
fc <- window(fitted(refit), start=c(2020,46))
dffc<-as.data.frame(fc)
write.table(dffc,file = "clipboard",dec = ",",sep = "\t")
summary(fit)
# iterate over all service point areas

# import predictions from ARIMA (all iterations)
pred_ARIMA <- as.data.frame(read.csv("predictions_ARIMA.csv", head = T, sep =",", dec="."))
temp_data <- pred_ARIMA[pred_ARIMA$week > 46,] # uses predicted week instead of actual week, so keep weeks 47-79
pred_df <- cbind(pred_df, temp_data[c(4,1)])
names(pred_df)[11] <- "pred_ARIMA"
rm(temp_data)

# calculate RMSE for all models
CalcRMSE <- function(data) {
  rmse_list_INLA <- c()
  rmse_list_RF <- c()
  rmse_list_GB <- c()
  rmse_list_ARIMA <- c()
  week_ids <- unique(data[["predicted_week"]])
  for (w in week_ids) { 
    temp_data <- data[data[["predicted_week"]] == w,] 
    temp_rmse <- rmse(temp_data$pred_INLA, temp_data$positivity_1w_nextweek) 
    rmse_list_INLA <- append(rmse_list_INLA, temp_rmse)
    temp_rmse <- rmse(temp_data$pred_RF, temp_data$positivity_1w_nextweek) 
    rmse_list_RF <- append(rmse_list_RF, temp_rmse)
    temp_rmse <- rmse(temp_data$pred_GB, temp_data$positivity_1w_nextweek) 
    rmse_list_GB <- append(rmse_list_GB, temp_rmse)
    temp_rmse <- rmse(temp_data$pred_ARIMA, temp_data$positivity_1w_nextweek) 
    rmse_list_ARIMA <- append(rmse_list_ARIMA, temp_rmse)
  }
  rmse_df <- data.frame("pred_week"=week_ids, "rmse_INLA"=rmse_list_INLA, "rmse_RF"=rmse_list_RF,
                        "rmse_GB"=rmse_list_GB, "rmse_ARIMA"=rmse_list_ARIMA)
  return(rmse_df)
}

rmse_df_all <- CalcRMSE(pred_df)

# naive model RMSE (= predict the same as last week)
NaivePred <- function(data) {
  rmse_list <- c()
  week_ids <- unique(data[["week"]])
  pred_weeks <- week_ids[20:(length(week_ids)-1)] # minimum 20 training weeks + 1 prediction week
  for (w in pred_weeks) { 
    temp_data <- data[data[["week"]] == w,] 
    temp_rmse <- rmse(temp_data$positivity_1w, temp_data$positivity_1w_nextweek) 
    rmse_list <- append(rmse_list, temp_rmse)
  }
  rmse_df <- data.frame("week"=pred_weeks+1, "rmse"=rmse_list)
  return(rmse_df)
}

rmse_naive <- NaivePred(data_sp)
rmse_df_all$rmse_naive <- rmse_naive$rmse

# import iso weeks for plot labels
iso_week <- as.data.frame(read.csv("iso_weeks.csv", head = T, sep =",", dec="."))
week_labels <- iso_week[iso_week$prediction_week %in% c(45,50,55,60,65,70,75,80),]
week_list <- week_labels$actual_ISO_week

# plot rmse time series
par(mfrow=c(1,1))
plot(rmse_df_all$pred_week, rmse_df_all$rmse_INLA, type="l", xlab="Week", ylab="RMSE", col="orange", lwd=2, ylim=c(0,0.09), xaxt="n")
axis(side=1, at=c(45,50,55,60,65,70,75,80), labels=week_list)
lines(rmse_df_all$pred_week, rmse_df_all$rmse_GB, col="red", lwd=2)
lines(rmse_df_all$pred_week, rmse_df_all$rmse_RF, col="blue", lwd=2)
lines(rmse_df_all$pred_week, rmse_df_all$rmse_ARIMA, col="purple", lwd=2)
lines(rmse_df_all$pred_week, rmse_df_all$rmse_naive, col="darkgrey", lwd=2, lty=2)
legend(x="topright", legend = c("RF", "GB", "INLA", "ARIMA", "Naive"), col=c("blue", "red", "orange", "purple", "darkgrey"), lwd=2, lty=c(1,1,1,1,2))

# plot observed positivity for some areas (based on expert-selected low/high pandemic effects)
plot(data_sp$week[data_sp$sp_id==42 & data_sp$week >=47], data_sp$positivity_1w[data_sp$sp_id==42 & data_sp$week >=47], type="l", ylab="Positivity", xlab="Week", ylim=c(0,0.3), col="purple") # Gottsunda, id 42
lines(data_sp$week[data_sp$sp_id==44 & data_sp$week >=47], data_sp$positivity_1w[data_sp$sp_id==44 & data_sp$week >=47], col="red") # Skutskär, id 44
lines(data_sp$week[data_sp$sp_id==45 & data_sp$week >=47], data_sp$positivity_1w[data_sp$sp_id==45 & data_sp$week >=47], col="green") # Älvkarlaby, id 45
lines(data_sp$week[data_sp$sp_id==12 & data_sp$week >=47], data_sp$positivity_1w[data_sp$sp_id==12 & data_sp$week >=47], col="orange") # Heby, id 12
# selected low pandemic effect area: Heby (id 12)
# selected high pandemic effect area: Skutskär or Älvkarleby (id 44 or 45)

# compare mean RMSE
summary(rmse_df_all)

# df with observed data (pred_week - 1)
data_sp_obs <- data_sp[data_sp$week > 46,]

# plot observed vs predicted positivity in Heby (low pandemic effect)
plot(data_sp_obs$week[data_sp_obs$sp_id==12], data_sp_obs$positivity_1w[data_sp_obs$sp_id==12], type="l", 
     ylab="PCR test positivity (%)", xlab="Week", ylim=c(0,0.3), col="darkgrey", lty=2, lwd=2, main="Heby", xaxt="n") 
axis(side=1, at=c(45,50,55,60,65,70,75,80), labels=week_list)
lines(pred_df$predicted_week[pred_df$sp_id==12], pred_df$pred_INLA[pred_df$sp_id==12], col="orange", lwd=2)
lines(pred_df$predicted_week[pred_df$sp_id==12], pred_df$pred_GB[pred_df$sp_id==12], col="red", lwd=2)
lines(pred_df$predicted_week[pred_df$sp_id==12], pred_df$pred_RF[pred_df$sp_id==12], col="blue", lwd=2)
lines(pred_df$predicted_week[pred_df$sp_id==12], pred_df$pred_ARIMA[pred_df$sp_id==12], col="purple", lwd=2)
legend(x="topright", legend = c("RF", "GB", "INLA", "ARIMA", "observed"), col=c("blue", "red", "orange", "purple", "darkgrey"), lwd=2, lty=c(1,1,1,1,2))

# plot observed vs predicted positivity in Skutskär (high pandemic effect)
plot(data_sp_obs$week[data_sp_obs$sp_id==44], data_sp_obs$positivity_1w[data_sp_obs$sp_id==44], type="l", 
     ylab="PCR test positivity (%)", xlab="Week", ylim=c(0,0.3), col="darkgrey", lty=2, lwd=2, main="Skutskär", xaxt="n") 
axis(side=1, at=c(45,50,55,60,65,70,75,80), labels=week_list)
lines(pred_df$predicted_week[pred_df$sp_id==44], pred_df$pred_INLA[pred_df$sp_id==44], col="orange", lwd=2)
lines(pred_df$predicted_week[pred_df$sp_id==44], pred_df$pred_GB[pred_df$sp_id==44], col="red", lwd=2)
lines(pred_df$predicted_week[pred_df$sp_id==44], pred_df$pred_RF[pred_df$sp_id==44], col="blue", lwd=2)
lines(pred_df$predicted_week[pred_df$sp_id==44], pred_df$pred_ARIMA[pred_df$sp_id==44], col="purple", lwd=2)
legend(x="topright", legend = c("RF", "GB", "INLA", "ARIMA", "observed"), col=c("blue", "red", "orange", "purple", "darkgrey"), lwd=2, lty=c(1,1,1,1,2))

# plot observed vs predicted positivity in Älvkarleby (high pandemic effect, variability in predictions ARIMA)
plot(data_sp_obs$week[data_sp_obs$sp_id==45], data_sp_obs$positivity_1w[data_sp_obs$sp_id==45], type="l", 
     ylab="PCR test positivity (%)", xlab="Week", ylim=c(0,0.3), col="darkgrey", lty=2, lwd=2, main="Älvkarleby", xaxt="n") 
axis(side=1, at=c(45,50,55,60,65,70,75,80), labels=week_list)
lines(pred_df$predicted_week[pred_df$sp_id==45], pred_df$pred_INLA[pred_df$sp_id==45], col="orange", lwd=2)
lines(pred_df$predicted_week[pred_df$sp_id==45], pred_df$pred_GB[pred_df$sp_id==45], col="red", lwd=2)
lines(pred_df$predicted_week[pred_df$sp_id==45], pred_df$pred_RF[pred_df$sp_id==45], col="blue", lwd=2)
lines(pred_df$predicted_week[pred_df$sp_id==45], pred_df$pred_ARIMA[pred_df$sp_id==45], col="purple", lwd=2)
legend(x="topright", legend = c("RF", "GB", "INLA", "ARIMA", "observed"), col=c("blue", "red", "orange", "purple", "darkgrey"), lwd=2, lty=c(1,1,1,1,2))


# add predictions to shapefile for plotting map
pred_wk66 <- pred_df[pred_df$predicted_week==66,] # peak of third wave: week 13, 2021 (week 66)
sp_shp_pred_wk66 <- merge(sp_shp, pred_wk66, by.x="sp_id", by.y="sp_id", sort=F)
map_pal <- colorRampPalette(c("khaki1", "gold", "orange", "orangered", "darkred"))
# entire region
spplot(sp_shp_pred_wk66, zcol=c("positivity_1w_nextweek", "pred_INLA", "pred_RF", "pred_GB", "pred_ARIMA"), 
       names.attr=c("observed", "INLA", "RF", "GB", "ARIMA"), col.regions=map_pal(16), main="Week 13, 2021 (Uppsala region)")
# zoom in Uppsala city
spplot(sp_shp_pred_wk66, zcol=c("positivity_1w_nextweek", "pred_INLA", "pred_RF", "pred_GB", "pred_ARIMA"), 
       names.attr=c("observed", "INLA", "RF", "GB", "ARIMA"), col.regions=map_pal(16), main="Week 13, 2021 (Uppsala city)", 
       xlim=c(17.5,17.8), ylim=c(59.8,59.97))

# iteratively plot all
for (w in c(47:79)) {
  temp_pred <- pred_df[pred_df$predicted_week==w,] 
  temp_shp <- merge(sp_shp, temp_pred, by.x="sp_id", by.y="sp_id", sort=F)
  print(spplot(temp_shp, zcol=c("positivity_1w_nextweek", "pred_INLA", "pred_RF", "pred_GB", "pred_ARIMA"), 
         names.attr=c("observed", "INLA", "RF", "GB", "ARIMA"), col.regions=map_pal(16), 
         main=as.character(iso_week$actual_ISO_week[iso_week$prediction_week==w])))
  rm(temp_pred, temp_shp)
}




### model ensemble ###

# linear model to estimate shares
summary(lm(positivity_1w_nextweek ~ pred_INLA + pred_RF + pred_GB + pred_ARIMA, data=pred_df))
# suggests 0.4 * INLA, 0.4 * RF, 0.2 * GB, no ARIMA 

EnsembleRMSE <- function(data) {
  rmse_list_INLA <- c()
  rmse_list_RF <- c()
  rmse_list_GB <- c()
  rmse_list_ARIMA <- c()
  rmse_list_combi <- c()
  week_ids <- unique(data[["predicted_week"]])
  for (w in week_ids) { 
    temp_data <- data[data[["predicted_week"]] == w,] 
    temp_rmse <- rmse(temp_data$pred_INLA, temp_data$positivity_1w_nextweek) 
    rmse_list_INLA <- append(rmse_list_INLA, temp_rmse)
    temp_rmse <- rmse(temp_data$pred_RF, temp_data$positivity_1w_nextweek) 
    rmse_list_RF <- append(rmse_list_RF, temp_rmse)
    temp_rmse <- rmse(temp_data$pred_GB, temp_data$positivity_1w_nextweek) 
    rmse_list_GB <- append(rmse_list_GB, temp_rmse)
    temp_rmse <- rmse(temp_data$pred_ARIMA, temp_data$positivity_1w_nextweek) 
    rmse_list_ARIMA <- append(rmse_list_ARIMA, temp_rmse)
    temp_data$pred_combi <- ((0.4)*temp_data$pred_INLA + (0.4)*temp_data$pred_RF + (0.2)*temp_data$pred_GB)
    temp_rmse <- rmse(temp_data$pred_combi, temp_data$positivity_1w_nextweek)
    rmse_list_combi <- append(rmse_list_combi, temp_rmse)
  }
  rmse_df <- data.frame("pred_week"=week_ids, "rmse_INLA"=rmse_list_INLA, "rmse_RF"=rmse_list_RF,
                        "rmse_GB"=rmse_list_GB, "rmse_ARIMA"=rmse_list_ARIMA, "rmse_combi"=rmse_list_combi)
  return(rmse_df)
}

rmse_df_combi <- EnsembleRMSE(pred_df)
rmse_df_combi$rmse_naive <- rmse_naive$rmse

# plot rmse time series incl. ensemble
par(mfrow=c(1,1))
plot(rmse_df_combi$pred_week, rmse_df_combi$rmse_INLA, type="l", xlab="Week", ylab="RMSE", col="orange", lwd=2, ylim=c(0,0.09), xaxt="n")
axis(side=1, at=c(45,50,55,60,65,70,75,80), labels=week_list)
lines(rmse_df_combi$pred_week, rmse_df_combi$rmse_GB, col="red", lwd=2)
lines(rmse_df_combi$pred_week, rmse_df_combi$rmse_RF, col="blue", lwd=2)
lines(rmse_df_combi$pred_week, rmse_df_combi$rmse_ARIMA, col="purple", lwd=2)
lines(rmse_df_combi$pred_week, rmse_df_combi$rmse_combi, col="black", lwd=2, lty=2)
lines(rmse_df_combi$pred_week, rmse_df_combi$rmse_naive, col="darkgrey", lwd=2, lty=2)
legend(x="topright", legend = c("RF", "GB", "INLA", "ARIMA", "Naive", "Combi"), col=c("blue", "red", "orange", "purple", "darkgrey", "black"), lwd=2, lty=c(1,1,1,1,2,2))

mean(rmse_df_combi$rmse_combi)
