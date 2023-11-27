#--- Load libraries
library(geoR)
library(mvtnorm)
library(sp)
library(RColorBrewer)
library(gridExtra)
library(gstat)
library(viridis)
library(proxy)
library(matrixcalc)
library(latticeExtra)
library(sf)
library(pbapply)
library(tigris)
library(viridis)
library(tidyverse)
library(rgeoda)
library(extrafont)
seed = 512

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#------------------------------ REAL DATA --------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# ---------------------------------- 1 -----------------------------------------
# -------------------------------- STDs ----------------------------------------
# ------------------------------------------------------------------------------


# Set working folder
path.folder = "your folder"
setwd(path.folder)

#--- Load helper functions
source('./Multivariate-Poisson-cokriging-utils.R')

# ----- Load spatial data

# Load counties and census tracts and reproject
county = counties(state = 'Pennsylvania') %>%  st_transform(26915)
tract = st_read("./Data/PA-Tracts/tl_2018_42_tract.shp") %>%  st_transform(26915)

# Remove unwanted attributes
county = county[,c('GEOID')]
county$FIPS = as.numeric(county$GEOID)
county = county[,c('FIPS')]
tract = tract[,c('GEOID')]

# Scale coordinates to km
county$geometry = county$geometry / 1000
tract$geometry = tract$geometry / 1000

# Plot geometry
plot(county$geometry)
axis(1);axis(2)

# ----- Load disease data

# List data
files = list.files(path = "./Data", pattern = "\\PA.csv$", full.names = T)

# Define objects names
names = c("HIV.Inc","HIV.Prev","Chlamydia","Gonorrhea")

# Load data as a data frame
STDs = lapply(files, function(x) {d = read.csv(x, sep = ";");                         # Read csv
                                  d = d[,c(1,3,9:10)];                                # Remove unwanted attributes
                                  d = d[!grepl("42101", d$FIPS),]                     # Remove Philadelphia county (it is adding noise)
                                  d$Rate.per.100000 = d$Cases/d$Population * 100000;  # Calculate sample rate
                                  d = merge(county, d, by = "FIPS");                  # Merge with spatial object
                                  return(d)
                                  })
# Assign names to data frames 
names(STDs) = names

# ----- Plot disease data

# Get rate quantile ranges
breaks.plot <- lapply(1:length(STDs), function(x){unique(round(c(0.001,1.001,unname(quantile(STDs[[x]]$Rate.per.100000,na.rm = TRUE, prob = c(0.15,0.30,0.45,0.60,0.75,1)))),1))})

# Plot data
plot.diseases <- lapply(1:length(STDs), function(x){plot.data(STDs[[x]], STDs[[x]]$Rate.per.100000, breaks.plot[[x]], "" , "")})
plot.diseases[[4]]
grid.arrange(plot.diseases[[1]], plot.diseases[[2]], plot.diseases[[3]], plot.diseases[[4]],  ncol = 2)


# ----- Load census tract data

# Load census data
pop.census = read.csv("./Data/Population-Census-Tracts.csv", sep = ";", colClasses = c("GEOID"="character"))
pop.census$FIPS = as.numeric(substr(pop.census$GEOID, 1, 5) )


# Merge with spatial object
pop.census.sp = merge(tract,pop.census, by = "GEOID")
pop.census.sp = pop.census.sp[!grepl("42101", pop.census.sp$FIPS),]

# Plot data
break.pop = unique(c(unname(quantile(pop.census.sp$Pop,na.rm = TRUE, prob = c(0,0.15,0.30,0.45,0.60,0.75,1)))))
plot.data(pop.census.sp, pop.census.sp$Pop, break.pop, "Census tracts population", "CDC")


# ----- Load spatial data

# Check correlation
cor(data.frame(STDs$Chlamydia$Rate.per.100000, STDs$Gonorrhea$Rate.per.100000,STDs$HIV.Inc$Rate.per.100000,STDs$HIV.Prev$Rate.per.100000), use = "complete.obs")

#Calculate Empirical Bayesian standardization
eb.rates = lapply(STDs, function(x) {eb.risk(x)})
breaks.eb.plot <- lapply(1:length(STDs), function(x){unique(round(c(unname(quantile(eb.rates[[x]],na.rm = TRUE, prob = c(0,0.15,0.30,0.45,0.60,0.75,1)))),3))})
plot.diseases.eb <- lapply(1:length(STDs), function(x){ plot.data(STDs[[x]], eb.rates[[x]], breaks.eb.plot[[x]], names[[x]] , "CDC")})
grid.arrange(plot.diseases.eb[[1]], plot.diseases.eb[[2]], plot.diseases.eb[[3]], plot.diseases.eb[[4]], ncol = 2)

# ----- Calculate population-weighted centroids

# Calculate weighted centroids per county
pop.census.sp$x = st_coordinates(st_point_on_surface(pop.census.sp))[, 1]
pop.census.sp$y = st_coordinates(st_point_on_surface(pop.census.sp))[, 2]
weighted.centroids <- pop.census.sp %>% dplyr::group_by(FIPS) %>% dplyr::summarize(weightedpop(x,y,Pop)) %>% st_drop_geometry()

# Create data wih the centroids geometry
weighted.centroids <- st_as_sf(weighted.centroids, coords = c("x", "y"), crs = 26915)

# Plot centroids
plot(county$geometry)
plot(st_point_on_surface(STDs$Chlamydia)$geometry, col = 'blue', add = T)
plot( weighted.centroids$geometry, col = 'red', add = T , pch = 3)
legend( title = 'Centroids', "bottomleft", legend = c('Geographic', 'Population-weighted'), col = c('red', 'blue'), pch = c(1, 3), cex = 0.8)

# Assign geometry 
STDs.points =  lapply(STDs, function(x){x$geometry = weighted.centroids$geometry; return(x)})

# Calculate coordinates
STDs.points =  lapply(STDs.points, function(x){coords = st_coordinates(x); x$x =  coords[,1]; x$y = coords[,2]; return(x)})
STDs =  lapply(STDs, function(x){coords = st_coordinates(st_centroid(x)); x$x =  coords[,1]; x$y = coords[,2]; return(x)})

# ------- Model spatial dependence

# Estimate empirical direct semivariograms
distances = c(260,260,260,260)
direct.vars = lapply(1:length(STDs), function(x){pck.variogram(data = STDs[[x]], cases = "Cases", pop = "Population", 
                                                               pop.rate = 100000, maxd = distances[[x]] ,nbins = 12)})

# Plot direct semivariogram
plots.vars = lapply(1:length(direct.vars), function(x){plot.var.pck(direct.vars[[x]])})
grid.arrange(plots.vars[[1]], plots.vars[[2]], plots.vars[[3]], plots.vars[[4]], ncol = 2)

#Fit direct semivariograms)
fit.1 = fit.variogram(direct.vars[[1]], vgm(25, "Exp", 25.808) , fit.method = 2, fit.ranges = F)
attr(fit.1,"SSErr")
plot(direct.vars[[1]], fit.1)

fit.2 = fit.variogram(direct.vars[[2]], vgm(1500, "Exp", 30), fit.method = 2, fit.ranges = F)
attr(fit.2,"SSErr")
plot(direct.vars[[2]], fit.2)

fit.3 = fit.variogram(direct.vars[[3]], vgm(25000, "Exp", 30), fit.method = 2, fit.ranges = F)
attr(fit.3,"SSErr")
plot(direct.vars[[3]], fit.3)

fit.4 = fit.variogram(direct.vars[[4]], vgm(4000, "Exp", 40), fit.method = 2, fit.ranges = F)
attr(fit.4,"SSErr")
plot(direct.vars[[4]], fit.4)

fit.vars = list(fit.1, fit.2, fit.3, fit.4)

# Plot empirical direct semivariograms and fitted ones
title.vars = c("HIV Inc.", "HIV Prev.", "Chlamydia", "Gonorrhea")
plots.fit.vars = lapply(1:length(direct.vars), function(x){plot.pck.variogram(direct.vars[[x]], fit.vars[[x]], title.vars[[x]])})
grid.arrange(plots.fit.vars[[1]], plots.fit.vars[[2]], plots.fit.vars[[3]], plots.fit.vars[[4]], ncol = 2)

#--------- Joint cases

# Define joint risk
shared.risk = list(hi.hp = 0.3, hi.c = 5.9 / 1000, hi.g = 0.5 / 100, hp.c = 10.5 / 100, hp.g = 2.7 / 100, c.g = 0.089 / 100)
names.shared = list()

# Estimate joint risk
for (i in 1:length(shared.risk)){
  
  # Get the disease
  names.lk = unlist(strsplit(names(shared.risk)[[i]], ".", fixed = TRUE))
  name.l = names.lk[1]
  name.k = names.lk[2]
  
  # Obtain risk name
  risk.name =  paste0("R", name.l, name.k)
  
  # Get the data
  data.l = data.get(name.l)
  data.k = data.get(name.k)
  
  # Compute shared risk
  r.lk = rlk(data.l$Rate.per.100000,
             data.k$Rate.per.100000,
             shared.risk[[i]],
             log = F)
  
  # Assign shared risk
  STDs[[data.name(name.l)]][[risk.name]] <- r.lk
  STDs[[data.name(name.k)]][[risk.name]] <- r.lk
  
  # Get names risk column
  names.shared = append(names.shared, risk.name)

}

# Estimate empirical cross semivariograms
STDs.l = list(STDs $HIV.Inc, STDs $HIV.Inc, STDs $HIV.Inc, STDs $HIV.Prev, STDs $HIV.Prev, STDs $Chlamydia)
STDs.k = list(STDs $HIV.Prev, STDs $Chlamydia , STDs $Gonorrhea,   STDs $Chlamydia , STDs $Gonorrhea   , STDs $Gonorrhea)
cross.vars = lapply(1:length(STDs.l), function(x){pck.crossvariogram(STDs.l[[x]], STDs.k[[x]], "Cases", "Cases", names.shared[[x]], "Population", maxd = 260, nbins = 12)})

# Plot empirical cross semivariograms
grid.arrange(plot(cross.vars[[1]]), plot(cross.vars[[2]]), plot(cross.vars[[3]]), plot(cross.vars[[4]]), plot(cross.vars[[5]]), plot(cross.vars[[6]]), ncol = 3)


#Fit direct semivariograms)
fit.cg = fit.variogram(cross.vars[[1]], vgm(150, "Exp", 30), fit.method = 2, fit.ranges = F)
attr(fit.cg,"SSErr")
plot(cross.vars[[1]], fit.cg)

fit.chp = fit.variogram(cross.vars[[2]], vgm(600, "Exp", 25), fit.method = 2, fit.ranges = F)
attr(fit.chp,"SSErr")
plot(cross.vars[[2]], fit.chp)

fit.chi = fit.variogram(cross.vars[[3]], vgm(200, "Exp", 30), fit.method = 2, fit.ranges = F)
attr(fit.chi,"SSErr")
plot(cross.vars[[3]], fit.chi)

fit.ghp = fit.variogram(cross.vars[[4]], vgm(4000, "Exp", 30), fit.method = 2, fit.ranges = F)
attr(fit.ghp,"SSErr")
plot(cross.vars[[4]], fit.ghp)

fit.ghi = fit.variogram(cross.vars[[5]], vgm(1500, "Exp", 25), fit.method = 2, fit.ranges = F)
attr(fit.ghi,"SSErr")
plot(cross.vars[[5]], fit.ghi)

fit.hihp = fit.variogram(cross.vars[[6]], vgm(6000, "Exp", 30), fit.method = 2, fit.ranges = F)
attr(fit.hihp,"SSErr")
plot(cross.vars[[6]], fit.hihp)


fit.crossvars = list(fit.cg, fit.chp, fit.chi, fit.ghp, fit.ghi, fit.hihp)

# Plot empirical direct semivariograms and fitted ones
title.crossvars = c("HIV Inc.|HIV Prev.", "HIV Inc.|Chlamydia", "HIV Inc.|Gonorrhea", "HIV Prev.|Chlamydia", "HIV Prev.|Gonorrhea", "Chlamydia|Gonorrhea")
plots.fit.crossvars = lapply(1:length(cross.vars), function(x){plot.pck.variogram(cross.vars[[x]], fit.crossvars[[x]], title.crossvars[[x]])})
grid.arrange(plots.fit.crossvars[[1]], plots.fit.crossvars[[2]], plots.fit.crossvars[[3]], plots.fit.crossvars[[4]], plots.fit.crossvars[[5]], plots.fit.crossvars[[6]] ,ncol = 3)


# Plot LMC
grid.arrange(plots.fit.vars[[1]]     , plot.new()               , plot.new()         , plot.new() ,
             plots.fit.crossvars[[1]], plots.fit.vars[[2]]      , plot.new()         , plot.new() ,
             plots.fit.crossvars[[2]], plots.fit.crossvars[[4]] , plots.fit.vars[[3]], plot.new() ,
             plots.fit.crossvars[[3]], plots.fit.crossvars[[5]] , plots.fit.crossvars[[6]], plots.fit.vars[[4]],
             ncol = 4
             )


#--- Define initial values


r = 25.808
lmc.init.stds <- list( v1 = vgm(25, "Exp", r),
                  v2 = vgm(1500, "Exp", r),
                  v3 = vgm(25000, "Exp", r),
                  v4 = vgm(4000, "Exp", r),
                  v12 = vgm(150, "Exp", r),
                  v13 = vgm(600, "Exp", r),
                  v14 = vgm(200, "Exp", r),
                  v23 = vgm(4000, "Exp", r),
                  v24 = vgm(1500, "Exp", r),
                  v34 = vgm(6000, "Exp", r))

#--- Linear model of Coregionalization
fitted.lmc.four <- lmc.poisson.cokrige.four(direct.vars[[1]], direct.vars[[2]], direct.vars[[3]], direct.vars[[4]],
                                       cross.vars[[1]], cross.vars[[2]], cross.vars[[3]], cross.vars[[4]], cross.vars[[5]], cross.vars[[6]],
                                       STDs[[1]], STDs[[2]], STDs[[3]], STDs[[4]],
                                       lmc.init.stds, fit.ranges = F)



#------------------ Prediction PCK

#--- Predict

# Select locations without data
hiv.preds = STDs$HIV.Inc[is.na(STDs$HIV.Inc$Cases),]

# Coordinates of prediction locations
coords.hi = cbind(hiv.preds$x, hiv.preds$y)

# Predicting
predictions.hiv = poisson.cokrige.four(STDs$HIV.Inc, STDs$HIV.Prev, STDs$Chlamydia, STDs$Gonorrhea, coords.pred = coords.hi, fitted.lmc.four, 100000)
predictions.hiv.3 = poisson.cokrige.three(STDs$HIV.Inc, STDs$HIV.Prev, STDs$Chlamydia, coords.pred = coords.hi, fitted.lmc.four, 100000)
predictions.hiv.2 = poisson.cokrige.bi(STDs$HIV.Inc, STDs$HIV.Prev, coords.pred = coords.hi, fitted.lmc.four, 100000)
mean(predictions.hiv$var)
mean(predictions.hiv.3$var)
mean(predictions.hiv.2$var)

# Add FIPS
predictions.hiv$FIPS = hiv.preds$FIPS 

# Add prediciton columns to dataset
STDs$HIV.Inc[STDs$HIV.Inc$FIPS %in% predictions.hiv$FIPS,"Rate.per.100000"] <- predictions.hiv$pred
STDs$HIV.Inc[STDs$HIV.Inc$FIPS %in% predictions.hiv$FIPS,"var"] <- predictions.hiv$var

#--- Smooth

# Select locations to smooth
coords.hi.s = cbind(STDs$HIV.Inc[!(is.na(STDs$HIV.Inc$Cases)),]$x, STDs$HIV.Inc[!(is.na(STDs$HIV.Inc$Cases)),]$y)

#Smoothing
smooths.pck <- poisson.cokrige.four(STDs$HIV.Inc, STDs$HIV.Prev, STDs$Chlamydia, STDs$Gonorrhea, coords.pred = coords.hi.s, fitted.lmc.four, 100000)
smooths.pck.3 <- poisson.cokrige.three(STDs$HIV.Inc, STDs$HIV.Prev, STDs$Chlamydia, coords.pred = coords.hi.s, fitted.lmc.four, 100000)
smooths.pck.2 <- poisson.cokrige.bi(STDs$HIV.Inc, STDs$HIV.Prev, coords.pred = coords.hi.s, fitted.lmc.four, 100000)

#AE
plot(smooths.pck$rate.pred, STDs$HIV.Inc[!(is.na(STDs$HIV.Inc$Cases)),]$Rate.per.100000)
mean(smooths.pck$residual)
mean(smooths.pck.3$residual)
mean(smooths.pck.2$residual)

# MSPE
mean(smooths.pck$residual^2)
mean(smooths.pck.3$residual^2)
mean(smooths.pck.2$residual^2)

#Cor
cor(smooths.pck$rate.pred, STDs$HIV.Inc[!(is.na(STDs$HIV.Inc$Cases)),]$Rate.per.100000)
cor(smooths.pck.3$rate.pred, STDs$HIV.Inc[!(is.na(STDs$HIV.Inc$Cases)),]$Rate.per.100000)
cor(smooths.pck.2$rate.pred, STDs$HIV.Inc[!(is.na(STDs$HIV.Inc$Cases)),]$Rate.per.100000)

# Plot data
break.pred = adjust.breaks(breaks.plot[[1]], max(STDs$HIV.Inc$Rate.per.100000))
plot.pred = plot.data(STDs$HIV.Inc$geometry, STDs$HIV.Inc$Rate.per.100000, break.pred, "HIV Estimated Risk", "CDC")
grid.arrange(plot.diseases[[1]] , plot.pred, ncol = 2)

# Integrate smoothing and prediction
STDs$HIV.Inc[!(is.na(STDs$HIV.Inc$Cases)),"preds"] <- smooths.pck$rate.pred
STDs$HIV.Inc[(is.na(STDs$HIV.Inc$Cases)),"preds"] <- predictions.hiv$pred
STDs$HIV.Inc[!(is.na(STDs$HIV.Inc$Cases)),"var"] <- smooths.pck$rate.var

# Plot data
break.pred = adjust.breaks(breaks.plot[[1]], max(STDs$HIV.Inc$preds))
plot.pred = plot.data(STDs$HIV.Inc$geometry, STDs$HIV.Inc$preds, break.pred, "", "")
grid.arrange(plot.diseases[[1]] , plot.pred, ncol = 2)

#Plot variance
break.var = c(unique(round(c(unname(quantile(STDs$HIV.Inc$var,na.rm = TRUE, prob = c(0,0.15,0.30,0.45,0.60,0.75,1))),max(STDs$HIV.Inc$var)),2)))
plot.variance = plot.data(STDs$HIV.Inc$geometry, STDs$HIV.Inc$var, break.var, "", "")
plot.variance

#------------------ Prediction PK

# Predicting
predictions.hiv.pk = poisson.krige.one(STDs$HIV.Inc, coords.pred = coords.hi, fit.1, 100000)
mean(predictions.hiv.pk$var)

#--- Smooth
smooths.pk <- poisson.krige.one(STDs$HIV.Inc, coords.pred = coords.hi.s, fit.1, 100000)
mean(smooths.pk$residual)
mean(smooths.pk$residual^2)
cor(smooths.pk$rate.pred, STDs$HIV.Inc[!(is.na(STDs$HIV.Inc$Cases)),]$Rate.per.100000)
mean(smooths.pk$rate.var)





