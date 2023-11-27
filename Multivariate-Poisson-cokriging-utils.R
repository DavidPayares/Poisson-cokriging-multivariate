library(lattice)
library(pbapply)

# --- Determine shared risk (proportion)
rlk = function(Rl, Rk, rho, log = T) {
  if (log) {
    shared = log(rho * sqrt(exp(Rl) * exp(Rk)))
  }
  else{
    shared = rho * sqrt(Rl * Rk)
  }
  return(shared)
}

# --- Select best sill and range
initial.params <- function(var, data){
  init.sill = mean(c(max(var$gamma), median(var$gamma)))
  data = na.omit(data)
  bb = st_bbox(data)
  init.range = sqrt(((bb$xmax-bb$xmin)^2)+((bb$ymax-bb$ymin)^2)) * 0.1
  values = c(sill = init.sill, range = init.range)
  return(values)
}


# Calculate direct and cross semivariogram
variogram.pk <- function(data, cases = "y", pop = "n", maxd = NULL, nbins = 15 , pop.rate = 1){
  
  
  # Calculate number of observations
  No = nrow(data)
  # Select counts and population size
  y = as.numeric(as.data.frame(data)[,cases])
  n = as.numeric(as.data.frame(data)[,pop])
  
  # Create indexes to use in building of empirical semivariograms
  idx1 = rep(1:(No - 1), times = (No - 1):1)
  idx2 = unlist(sapply(1:(No - 1), function(i) (1:No)[-(1:i)]))
  
  # Calculate distances
  d = c(dist(coordinates(data), diag = T))
  
  # Calculate difference of unique pairs of points
  diff = (pop.rate * (y[idx1] / n[idx1] - y[idx2] / n[idx2])) ^ 2
  
  # Calculate weighted term
  w = (n[idx1] * n[idx2]) / (n[idx1] + n[idx2])
  # bias term
  bias = sum(y) / sum(n) * pop.rate
  
  # Determine maximum distance of estimator
  if (is.null(maxd)) {
    maxd = max(d) / 2
  }
  #print(max(d))
  # Determine distances of tolerance regions for empirical semivariogram estimator
  bindist = seq(0, maxd, len = nbins + 1)
  # Determine classification cuts for unique pairs of points by distance
  dcuts = cut(d, breaks = bindist)
  
  # Calculate population-weighted empirical semivariogram
  np = unlist(lapply(split(d, dcuts), length), use.names = FALSE)
  middle = unlist(lapply(split(d, dcuts), mean), use.names = FALSE)
  wbinds = split(w , dcuts)
  rbinds = split((w * diff - bias), dcuts)
  semivariance = unlist(lapply(1:nbins, function(x) {
    sum(rbinds[[x]]) / (2 * sum(wbinds[[x]]))
  }), use.names = FALSE)
  
  var <- data.frame(as.numeric(np) , as.numeric(middle) ,semivariance,rep(0, nbins),
                    rep(0, nbins),rep(as.factor("var1"), nbins))
  names(var) <-  c("np", "dist", "gamma", "dir.hor", "dir.ver", "id")
  var <- var[!is.na(var$gamma),]
  rownames(var) <- NULL
  class(var) <- c("gstatVariogram", "data.frame")
  return(var)
  
}

# Calculate cross semivariogram
crossvariogram.pk <- function(data, cases.l = "yl", cases.k = "yk" , cases.lk = "ylk", pop = "n", maxd = NULL, nbins = 15 , pop.rate = 1){
  
  # Calculate number of observations
  No = nrow(data)
  # Select counts and population size
  yl = as.numeric(as.data.frame(data)[,cases.l])
  yk = as.numeric(as.data.frame(data)[,cases.k])
  ylk = as.numeric(as.data.frame(data)[,cases.lk])
  n = as.numeric(as.data.frame(data)[,pop])
  
  # Define total means
  Mlk <- sum(ylk) / sum(n)
  
  # Create indexes to use in building of empirical semivariograms
  idx1 = rep(1:No, times = No:1)
  idx2 = c(1,unlist(sapply(1:No, function(i) (1:No)[-(1:(i - 1))])))
  
  # distances for unique pairs of points
  dist = as.matrix(dist(coordinates(data)))
  d = c(dist[lower.tri(dist, diag = T)])
  
  # Calculate difference of unique pairs of points (centered and standardized)
  diff = (pop.rate * ((yl[idx1] / n[idx1] - yl[idx2] / n[idx2])*(yk[idx1] / n[idx1] - yk[idx2] / n[idx2])))
  
  # Calculate weighted term
  w = (n[idx1] * n[idx2]) / (n[idx1] + n[idx2])
  # bias term
  bias = Mlk * pop.rate
  
  # Determine maximum distance of estimator
  if (is.null(maxd)) {
    maxd = max(d) / 2
  }
  # Determine distances of tolerance regions for empirical semivariogram estimator
  bindist = seq(0, maxd, len = nbins + 1)
  # Determine classification cuts for unique pairs of points by distance
  dcuts = cut(d, breaks = bindist)
  
  
  # Calculate population-weighted empirical semivariogram
  np = unlist(lapply(split(d, dcuts), length), use.names = FALSE)
  middle = unlist(lapply(split(d, dcuts), mean), use.names = FALSE)
  wbinds = split(w , dcuts)
  rbinds = split((w * diff - bias), dcuts)
  semivariance = unlist(lapply(1:nbins, function(x) {
    sum(rbinds[[x]]) / (2 * sum(wbinds[[x]]))
  }), use.names = FALSE)
  
  # if(any(semivariance < 0) ){
  #   stop("semivariance values cannot be negative. Check observations.")
  # }
  
  # Transform varioram to object for gstat (easy manipulation of variogram)
  crossvar <- data.frame(as.numeric(np) , as.numeric(middle),semivariance,rep(0, nbins),
                         rep(0, nbins),rep(as.factor("var1.var2"), nbins))
  names(crossvar) <-  c("np", "dist", "gamma", "dir.hor", "dir.ver", "id")
  crossvar <- crossvar[!is.na(crossvar$gamma),]
  rownames(crossvar) <- NULL
  class(crossvar) <- c("gstatVariogram", "data.frame")
  return(crossvar)
  
}

# Adjust LMC
lmc.poisson.cokrige <- function(var.a, var.b, var.c, crossvar.ab, crossvar.ac, crossvar.bc, data, var.params){
  # Change id auxilliary variable
  var.b$id = as.factor("var2")
  var.c$id = as.factor("var3")
  crossvar.ac$id =  as.factor("var1.var3")
  crossvar.bc$id =  as.factor("var2.var3")
  
  names = c(expression(~ gamma[13]), expression(~ gamma[23]), expression(~ gamma[3]), expression(~ gamma[12]), expression(~ gamma[2]), expression(~ gamma[1]))
  my.settings <- list(
    strip.background=list(col="transparent"),
    strip.border=list(col="transparent")
  )
  
  # Integrate variograms
  variograms <- rbind(crossvar.ac, crossvar.bc, var.c, crossvar.ab, var.b, var.a)
  plot(variograms)
  
  # Gstat object
  g <- gstat(id = "var1", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = data$x ,y = data$y, c = data$Ya, p = data$N), model = var.params$v1)
  g <- gstat(g,id = "var2", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = data$x ,y = data$y, c = data$Yb, p = data$N), model = var.params$v2)
  g <- gstat(g,id = "var3", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = data$x ,y = data$y, c = data$Yc, p = data$N), model = var.params$v3)
  g <- gstat(g, id = c("var1","var2"), model = var.params$v12)
  g <- gstat(g, id = c("var1","var3"), model = var.params$v13)
  g <- gstat(g, id = c("var2","var3"), model = var.params$v23)
  fitted.lmc <- fit.lmc(v = variograms, g = g, fit.method = 2, fit.ranges = T)
  print(plot(variograms, fitted.lmc, pch = 20, col = "black", ylab = 'semivariance', par.settings =  my.settings, strip=strip.custom(factor.levels=names)))
  return(fitted.lmc)
}

# ---- Covariance function
Cexp <-  function(h,a,b){b * exp(-h / a)}
Csph = function(h,a,b){ifelse(h <= a, b * (1-1.5*(h/a)+0.5*(h/a)^3), 0)}

#--- Select cov model
covariance.model <- function(model){
  if(model == "Exp"){
    cov = Cexp
  }else if(model == "Sph"){
    cov = Csph
  }else{
    stop("incorrect covariance model")
  }
  return(cov)
}

#---- Distance function
distFun <- function(xi, x0){sqrt((x0[1]- xi[1])^2 + (x0[2]- xi[2])^2)}

# --- Poisson Cokriging helper function
poisson.cokrige.pred <- function(xao, xta, xtb, xtc, Cinv, cov.fun.1, cov.fun.2, data.a, data.b, data.c, pop.rate, par.a, par.ab, par.ac){
  
  # ----- Predict multiple locations
  ya <- data.a$Ya
  yb <- data.b$Yb
  yc <- data.c$Yc
  na <- data.a$N
  nb <- data.b$N
  nc <- data.c$N
  
  distXoa <- unlist(apply(xta,1, distFun, x0 = xao))
  distXob <- unlist(apply(xtb,1, distFun, x0 = xao))
  distXoc <- unlist(apply(xtc,1, distFun, x0 = xao))
  posa <- which(distXoa == 0)
  posb <- which(distXob == 0)
  posc <- which(distXoc == 0)
  
  
  #--- Covariances to prediction
  covaao <-  (cov.fun.1(distXoa,  par.a$range.1,  par.a$sill.1) + cov.fun.2(distXoa,  par.a$range.2,  par.a$sill.2))
  covbao <-  (cov.fun.1(distXob,  par.ab$range.1, par.ab$sill.1)+ cov.fun.2(distXob,  par.ab$range.2,  par.ab$sill.2))
  covcao <-  (cov.fun.1(distXoc,  par.ac$range.1, par.ac$sill.1)+ cov.fun.2(distXoc,  par.ac$range.2,  par.ac$sill.2))
  
  
  
  if(length(posa) > 0){
    covaao[posa] <- covaao[posa] + par.a$nugget.1 + par.a$nugget.2
    covbao[posb] <- covbao[posb] + par.ab$nugget.1 + par.ab$nugget.2
    covcao[posb] <- covcao[posb] + par.ac$nugget.1 + par.ac$nugget.2
  }
  coTotal <- c(covaao, covbao, covcao, 1, 0, 0)
  
  # --- Calculate weights
  lambda <- Cinv%*%coTotal
  #print(max(lambda))
  
  # predict
  pred.ra  = sum(lambda[1:(length(lambda)-3)] * c(ya/na, yb/nb, yc/nc) * pop.rate)
  if(pred.ra < 0){
    pred.ra = 0
  }
  var.ra = (cov.fun.1(0,  par.a$range.1,  par.a$sill.1) + cov.fun.2(0,  par.a$range.2,  par.a$sill.2)) - sum(lambda*coTotal)
  
  predictions = c(c(xao[1]), c(xao[2]) , pred.ra, var.ra)
  names(predictions) <- c('x','y','pred','var')
  return(predictions)
}

# --- Poisson Cokriging
poisson.cokrige <- function(data.a, data.b, data.c, coords.pred, lmc, pop.rate = 1){
  
  
  
  #-- Get covariance models
  cov.fun.1 = covariance.model(as.character(lmc$model$var1$model[1]))
  cov.fun.2 = covariance.model(as.character(lmc$model$var1$model[2]))
  
  #--- Get covariance parameters
  par.a <- data.frame(sill.1 =  lmc$model$var1$psill[1],       range.1 =  lmc$model$var1$range[1],       nugget.1 = 0,
                      sill.2 =  lmc$model$var1$psill[2],       range.2 =  lmc$model$var1$range[2],       nugget.2 = 0)
  par.b <- data.frame(sill.1 =  0,       range.1 =  1,       nugget.1 = 0,
                      sill.2 =  lmc$model$var2$psill[1],       range.2 =  lmc$model$var2$range[1],       nugget.2 = 0)
  # par.b <- data.frame(sill.1 =  lmc$model$var2$psill[1],       range.1 =  lmc$model$var2$range[1],       nugget.1 = 0,
  #                     sill.2 =  lmc$model$var2$psill[2],       range.2 =  lmc$model$var2$range[2],       nugget.2 = 0)
  par.c <- data.frame(sill.1 =  lmc$model$var3$psill[1],       range.1 =  lmc$model$var3$range[1],       nugget.1 = 0,
                      sill.2 =  lmc$model$var3$psill[2],       range.2 =  lmc$model$var3$range[2],       nugget.2 = 0)
  par.ab <- data.frame(sill.1 =  0, range.1 =  1,  nugget.1 = 0,
                       sill.2 =  lmc$model$var1.var2$psill[1], range.2 =  lmc$model$var1.var2$range[1],  nugget.2 = 0)
  # par.ab <- data.frame(sill.1 =  lmc$model$var1.var2$psill[1], range.1 =  lmc$model$var1.var2$range[1],  nugget.1 = 0,
  #                      sill.2 =  lmc$model$var1.var2$psill[2], range.2 =  lmc$model$var1.var2$range[2],  nugget.2 = 0)
  par.ac <- data.frame(sill.1 =  lmc$model$var1.var3$psill[1], range.1 =  lmc$model$var1.var3$range[1],  nugget.1 = 0,
                       sill.2 =  lmc$model$var1.var3$psill[2], range.2 =  lmc$model$var1.var3$range[2],  nugget.2 = 0)
  par.bc <- data.frame(sill.1 =  lmc$model$var2.var3$psill[1], range.1 =  lmc$model$var2.var3$range[1],  nugget.1 = 0,
                       sill.2 =  lmc$model$var2.var3$psill[2], range.2 =  lmc$model$var2.var3$range[2],  nugget.2 = 0)
  
  coords.a = cbind(data.a$x, data.a$y)
  coords.b = cbind(data.b$x, data.b$y)
  coords.c = cbind(data.c$x, data.c$y)
  
  if(nrow(coords.a) != nrow(coords.pred)){
    
    # distances
    distMa = as.matrix(dist(coords.a))
    distMb = as.matrix(dist(coords.b))
    distMc = as.matrix(dist(coords.c))
    distMab = as.matrix(proxy::dist(coords.a, coords.b))
    distMac = as.matrix(proxy::dist(coords.a, coords.c))
    distMbc = as.matrix(proxy::dist(coords.b, coords.c))
    
    # data sizes
    size.a = nrow(data.a)
    size.b = nrow(data.b)
    size.c = nrow(data.c)
    
    
    # ---- Error term W_{lk} 
    Waa <- diag(sum(data.a$Ya)/sum(data.a$N)  * pop.rate/data.a$N)
    Wbb <- diag(sum(data.b$Yb)/sum(data.a$N)  * pop.rate/data.b$N)
    Wcc <- diag(sum(data.c$Yc)/sum(data.a$N)  * pop.rate/data.c$N)
    
    indexab <- which(distMab == 0, arr.ind = TRUE)[,1]
    Wab <- matrix(0, dim(distMab)[1], dim(distMab)[2])
    Wab[which(distMab == 0, arr.ind = TRUE)] <- (sum(data.a$Yab)/sum(data.a$N)  * pop.rate)/data.a$N[indexab]
    
    indexac <- which(distMac == 0, arr.ind = TRUE)[,1]
    Wac <- matrix(0, dim(distMac)[1], dim(distMac)[2])
    Wac[which(distMac == 0, arr.ind = TRUE)] <- (sum(data.a$Yac)/sum(data.a$N) * pop.rate)/data.a$N[indexac]
    
    indexbc <- which(distMbc == 0, arr.ind = TRUE)[,1]
    Wbc <- matrix(0, dim(distMbc)[1], dim(distMbc)[2])
    Wbc[which(distMbc == 0, arr.ind = TRUE)] <- (sum(data.b$Ybc)/sum(data.a$N)  * pop.rate)/data.b$N[indexbc]
    
    
bloc
    
    pred1 <- poisson.cokrige.pred(coords.pred, xta = coords.a, xtb = coords.b, xtc = coords.c, Cinv = Cinv, 
                                  cov.fun.1 =  cov.fun.1 , cov.fun.2 = cov.fun.2, data.a = data.a, data.b = data.b, data.c = data.c, 
                                  pop.rate = pop.rate, par.a = par.a, par.ab = par.ab, par.ac = par.ac)
    
    
    # predict
    preds <- as.data.frame(t(pbapply(coords.pred, 1, poisson.cokrige.pred, xta = coords.a, xtb = coords.b, xtc = coords.c, Cinv = Cinv, 
                                     cov.fun.1 =  cov.fun.1 , cov.fun.2 = cov.fun.2, data.a, data.b, data.c, pop.rate, par.a, par.ab, par.ac)))
    
  }else{
    
    message("Input locations are the same as prediction locations. Performing smoothing instead of interpolation.")
    
    fold = 1:nrow(data.a)
    preds = data.frame(matrix(as.numeric(NA), nrow(data.a), 6))
    names(preds) = c('rate.pred','rate.var', 'observed', 'residual', 'zscore', 'fold')
    
    folds = sort(unique(fold))
    progress_bar = txtProgressBar(min=0, max=length(folds), style = 3, char  = "+")
    
    #Validation
    for (i in folds){
      
      sel = which(fold == i)
      data.cv = data.a[-sel, ]
      data.target = data.a[sel,]
      
      coords.a = cbind(data.cv$x, data.cv$y)
      coords.xo =  cbind(data.target$x, data.target$y)
      
      distMa = as.matrix(dist(coords.a))
      distMb = as.matrix(dist(coords.b))
      distMc = as.matrix(dist(coords.c))
      distMab = as.matrix(proxy::dist(coords.a, coords.b))
      distMac = as.matrix(proxy::dist(coords.a, coords.c))
      distMbc = as.matrix(proxy::dist(coords.b, coords.c))
      
      # data sizes
      size.a = nrow(data.cv)
      size.b = nrow(data.b)
      size.c = nrow(data.c)
      
      
      # ---- Error term W_{lk} 
      Waa <- diag(sum(data.cv$Ya)/sum(data.cv$N)  * pop.rate/data.cv$N)
      Wbb <- diag(sum(data.b$Yb)/sum(data.cv$N)  * pop.rate/data.b$N)
      Wcc <- diag(sum(data.c$Yc)/sum(data.cv$N)  * pop.rate/data.c$N)
      
      indexab <- which(distMab == 0, arr.ind = TRUE)[,1]
      Wab <- matrix(0, dim(distMab)[1], dim(distMab)[2])
      Wab[which(distMab == 0, arr.ind = TRUE)] <- (sum(data.cv$Yab)/sum(data.cv$N)  * pop.rate)/data.cv$N[indexab]
      
      indexac <- which(distMac == 0, arr.ind = TRUE)[,1]
      Wac <- matrix(0, dim(distMac)[1], dim(distMac)[2])
      Wac[which(distMac == 0, arr.ind = TRUE)] <- (sum(data.cv$Yac)/sum(data.cv$N) * pop.rate)/data.cv$N[indexac]
      
      indexbc <- which(distMbc == 0, arr.ind = TRUE)[,1]
      Wbc <- matrix(0, dim(distMbc)[1], dim(distMbc)[2])
      Wbc[which(distMbc == 0, arr.ind = TRUE)] <- (sum(data.b$Ybc)/sum(data.cv$N)  * pop.rate)/data.b$N[indexbc]
      
      # ---- Covariance matrices with W_{lk}
      Caa <-  (cov.fun.1(distMa,  par.a$range.1,  par.a$sill.1) + cov.fun.2(distMa,  par.a$range.2,  par.a$sill.2)) + Waa
      Cbb <-  (cov.fun.1(distMb,  par.b$range.1,  par.b$sill.1) + cov.fun.2(distMb,  par.b$range.2,  par.b$sill.2)) + Wbb
      Ccc <-  (cov.fun.1(distMc,  par.c$range.1,  par.c$sill.1) + cov.fun.2(distMc,  par.c$range.2,  par.c$sill.2)) + Wcc
      Cab <-  (cov.fun.1(distMab, par.ab$range.1, par.ab$sill.1)+ cov.fun.2(distMab,  par.ab$range.2,  par.ab$sill.2)) + Wab
      Cac <-  (cov.fun.1(distMac, par.ac$range.1, par.ac$sill.1)+ cov.fun.2(distMac,  par.ac$range.2,  par.ac$sill.2)) + Wac
      Cbc <-  (cov.fun.1(distMbc, par.bc$range.1, par.bc$sill.1)+ cov.fun.2(distMbc,  par.bc$range.2,  par.bc$sill.2)) +  Wbc    
      Cba <-  t(Cab)
      Cca <-  t(Cac)
      Ccb <-  t(Cbc)
      
      
      # ----- Unbiasedness constrains
      Ctotal <- rbind(cbind(Caa,Cab,Cac),cbind(Cba,Cbb, Cbc), cbind(Cca,Ccb,Ccc))
      Cbias <- cbind(c(rep(1, size.a), rep(0, size.b), rep(0, size.c)),
                     c(rep(0, size.a), rep(1, size.b), rep(0, size.c)),
                     c(rep(0, size.a), rep(0, size.b), rep(1, size.c)))
      Ctotal <- cbind(Ctotal, Cbias)
      Ctotal <- rbind(Ctotal, cbind(t(Cbias), matrix(0,3,3)))
      
      # ---- Inverse matrix
      Cinv <- solve(Ctotal)
      
      #---- Prediction
      pred.target <- as.data.frame(t(apply(coords.xo, 1, poisson.cokrige.pred, xta = coords.a, xtb = coords.b, xtc = coords.c, Cinv = Cinv, 
                                           cov.fun.1 =  cov.fun.1 , cov.fun.2 = cov.fun.2, data.cv, data.b, data.c, pop.rate, par.a, par.ab, par.ac)))
      
      preds[[1]][sel] = pred.target[[3]]
      preds[[2]][sel] = pred.target[[4]]
      preds[[3]][sel] = data.target$Ya/data.target$N * pop.rate
      preds[[4]][sel] = as.numeric(pred.target[[3]] - preds[[3]][sel])
      preds[[5]][sel] = as.numeric(preds[[4]][sel]/sqrt(pred.target[[4]]))
      preds[[6]][sel] = fold[i]
      
      setTxtProgressBar(progress_bar, value = i)
    }
    
    close(progress_bar)
    
  }
  
  return(preds)
}

# --- Poisson kriging helper function
poisson.krige.pred <- function(xo, xt, cov.inv, cov.fun, data , pop.rate, params){
  
  # ----- Predict multiple locations
  ya <- data$Ya
  n <- data$N
  
  #----- Size
  size = nrow(data)
  
  dist.xo <- apply(xt,1, distFun, x0 = xo)
  pos <- which(dist.xo == 0)
  
  #--- Covariances to prediction
  cov.xo <-  (cov.fun(as.matrix(dist.xo), params$range, params$sill))
  
  if(length(pos) > 0){
    cov.xo[pos] <- cov.xo[pos] + params$nugget
  }
  
  cov.xo <- rbind(cov.xo, 1)
  
  # --- Calculate weights
  lambda <- cov.inv%*%cov.xo
  
  # --- Predict
  pred  = sum(lambda[1:size]*(ya/n) * pop.rate)
  if(pred < 0){
    pred = 0
  }
  var = (cov.fun(0, params$range, params$sill) + params$nugget) - sum(lambda*cov.xo)
  
  predictions = c(c(xo[1]), c(xo[2]) , pred, var)
  names(predictions) <- c('x','y','pred','var')
  return(predictions)
}

# --- Poisson kriging
poisson.krige <- function(data, coords.pred, var.fitted, pop.rate = 10000){
  
  #--- check model
  if (nrow(var.fitted) > 1){
    params <- data.frame(sill = var.fitted$psill[2],  range =  var.fitted$range[2],     nugget = var.fitted$psill[1])
    mod = as.character(var.fitted[2,]$model)
  }else{
    params <- data.frame(sill = var.fitted$psill,  range =  var.fitted$range,     nugget = 0)
    mod = as.character(var.fitted$model)
  }
  
  
  #--- Get covariance parameters
  
  
  #--- Get covariance model
  cov.fun = covariance.model(mod)
  
  coords = cbind(data$x, data$y)
  
  
  if(nrow(coords) != nrow(coords.pred)){
    
    # distances
    dist.m = as.matrix(dist(coords))
    
    # data sizes
    size = nrow(data)
    
    # ---- Error term W_{lk} 
    bias <- sum(data$Ya)/sum(data$N) * pop.rate
    error.var <- diag(c(bias/data$N , 0))
    nug.term <- diag(c(rep(params$nugget, times = size) , 0))
    
    # ---- Covariance matrix
    cov.mat <-  (cov.fun(dist.m,  params$range,  params$sill))
    
    # Add unbiassedness conditions
    cov.mat <- cbind(cov.mat, c(rep(1, size)))
    cov.mat <- rbind(cov.mat, c(rep(1, size + 1)))
    cov.mat[size+1,size+1] <- 0
    
    # ---- Inverse matrix
    c.inv <- solve(cov.mat + error.var + nug.term)
    
    # predict
    preds <- as.data.frame(t(pbapply(coords.pred, 1, poisson.krige.pred, xt = coords, cov.inv = c.inv, 
                                     cov.fun =  cov.fun, data, pop.rate, params)))
  }else{
    
    message("Input locations are the same as prediction locations. Performing smoothing instead of interpolation.")
    
    fold = 1:nrow(data)
    preds = data.frame(matrix(as.numeric(NA), nrow(data), 6))
    names(preds) = c('rate.pred','rate.var', 'observed', 'residual', 'zscore', 'fold')
    
    folds = sort(unique(fold))
    progress_bar = txtProgressBar(min=0, max=length(folds), style = 3, char  = "+")
    
    #Validation
    for (i in folds){
      
      # select observation
      sel = which(fold == i)
      data.cv = data[-sel, ]
      data.target = data[sel,]
      
      # get coordinates
      coords = cbind(data.cv$x, data.cv$y)
      coords.xo =  cbind(data.target$x, data.target$y)
      
      # calculate distances
      dist.m = as.matrix(dist(coords))
      
      # calculate rate parameters
      size = nrow(data.cv)
      n = data.cv$N
      y = data.cv$Ya
      
      # bias term
      error.var = diag(c((sum(y) * pop.rate /sum(n))/n, 0))
      nug.term <- diag(c(rep(params$nugget, times = size) , 0))
      
      # covariance matrix
      cov.mat = (cov.fun(dist.m,  params$range,  params$sill))
      cov.mat = cbind(cov.mat, c(rep(1, size)))
      cov.mat = rbind(cov.mat, c(rep(1, size + 1)))
      cov.mat[size+1,size+1] <- 0
      c.inv = solve((cov.mat+error.var + nug.term))
      
      # predict
      pred.target <- as.data.frame(t(apply(coords.xo, 1, poisson.krige.pred, xt = coords, cov.inv = c.inv, 
                                           cov.fun =  cov.fun, data.cv, pop.rate, params)))
      
      preds[[1]][sel] = pred.target[[3]]
      preds[[2]][sel] = pred.target[[4]]
      preds[[3]][sel] = data.target$Ya/data.target$N * pop.rate
      preds[[4]][sel] = as.numeric(pred.target[[3]] - preds[[3]][sel])
      preds[[5]][sel] = as.numeric(preds[[4]][sel]/sqrt(pred.target[[4]]))
      preds[[6]][sel] = fold[i]
      
      setTxtProgressBar(progress_bar, value = i)
    }
    
    close(progress_bar)
    
  }
  
  return(preds)
}


#--------------------------------------------------------------------------------------------------------------------------------

#------ Plot data
plot.data <- function(data, field, breaks, title, subtitle){
  
  cols = RColorBrewer::brewer.pal(9, "Spectral")
  ggplot(data = data, aes(fill = cut(field, breaks = breaks, include.lowest = TRUE, dig.lab = 5))) + geom_sf(col =  'white', size = 0.3) +
    # scale_fill_brewer(type = "qual",
    #                   palette = "Spectral",
    #                   direction = -1,
    #                   name = "rate/10,000 persons-year") +
    scale_fill_viridis(
      discrete = TRUE,
      direction = -1 ,
      option = 'cividis',
      name = 'Risk per 100,0000/population',
      labels = function(breaks) {breaks[is.na(breaks)] <- "No Data"; breaks},
      na.value = "lightgray",
      guide = guide_legend(
        direction = "horizontal",
        reverse = FALSE,
        label.position = "bottom",
        title.position = 'top',
        title.hjust = 0.9,
        label.hjust = 0.9,
        nrow = 1
      )
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = title,
      subtitle = subtitle
    ) +
    theme(
      line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.background = element_blank(),
      #legend.position = c(0.75, 0.15),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.position = "bottom",
      legend.spacing.x = unit(0.0, 'cm'),
      legend.key.width = unit(2.5,"line"),
      legend.key.height = unit(1,"line"),
      legend.title.align = 0.5,
      legend.text.align = 0.5,
      legend.justification = "center"
    ) +
    coord_sf(datum = NA)
}

# Calculate empirical Bayesian standardized rates
eb.risk <- function(data, pop.rate = 100000){

  
  r = data$Cases
  n = data$Population
  x =  r/n
  N =  nrow(data)
  
  m = sum(r, na.rm = T) / sum(n, na.rm = T)
  s2 = sum(n * (x - m) ^ 2, na.rm = T) / sum(n, na.rm = T)
  a = s2 - (m / (sum(n, na.rm = T) / N))
  if(a < (m / (sum(n,na.rm = T) / N))){
    a = 0
  }
  
  v = a + (m/n)
  
  eb.rate = (x - m) / sqrt(v)
  return(eb.rate)
}

# Calculate Poisson semivariogram
pck.variogram <- function(data, cases = 'Y' , pop = 'N', maxd = NULL, nbins = 15 , pop.rate = 10000){
  
  # Remove NAs
  data <- data[!(is.na(data[[cases]])),]
  #data <- data[!grepl("42077", data$FIPS),]
  #print(nrow(data))
  
  #Get coordinates
  coords <- st_coordinates(st_centroid(data))
  rownames(coords) <- data$FIPS
  
  # Calculate number of observations
  data =  data %>% st_drop_geometry()
  No = nrow(data)
  # Select counts and population size
  y = as.numeric(data[,cases])
  n = as.numeric(data[,pop])
  
  
  # Create indexes to use in building of empirical semivariograms
  # idx1 = rep(1:(No - 1), times = (No - 1):1)
  # idx2 = unlist(sapply(1:(No - 1), function(i) (1:No)[-(1:i)]))
  idx1 =  rep(1:No, times= No)
  idx2 = rep(1:No, each = No)
  
  # Reshape distances for unique pairs of points
  d = c(as.matrix(dist(coords, diag = T, upper = T)))
  
  # Calculate difference of unique pairs of points
  diff = ((y[idx1] / n[idx1] * pop.rate - y[idx2] / n[idx2] * pop.rate)) ^ 2
  
  # Calculate weighted term
  w = (n[idx1] * n[idx2]) / (n[idx1] + n[idx2])
  # bias term
  bias = sum(y) / sum(n) * pop.rate

  
  # Determine maximum distance of estimator
  if (is.null(maxd)) {
    maxd = max(d) / 2
  }
  # Determine distances of tolerance regions for empirical semivariogram estimator
  bindist = seq(0, maxd, len = nbins + 1)
  # Determine classification cuts for unique pairs of points by distance
  dcuts = cut(d, breaks = bindist)
  
  
  # Calculate population-weighted empirical semivariogram
  np = unlist(lapply(split(d, dcuts), length), use.names = FALSE)
  middle = unlist(lapply(split(d, dcuts), mean), use.names = FALSE)
  wbinds = split(w , dcuts)
  rbinds = split((w * diff - bias), dcuts)
  semivariance = unlist(lapply(1:nbins, function(x) {
    sum(rbinds[[x]]) / (2 * sum(wbinds[[x]]))
  }), use.names = FALSE)
  # Plot population-weighted semivariogram
  plot(1:nbins, semivariance)
  
  semivariance[semivariance < 0] <- 0
  
  # Transform varioram to object for gstat (easy manipulation of variogram)
  var <- data.frame(as.numeric(np) , as.numeric(middle),semivariance,rep(0, nbins),
                    rep(0, nbins),rep(as.factor("var1"), nbins))
  names(var) <-  c("np", "dist", "gamma", "dir.hor", "dir.ver", "id")
  var <- var[!is.na(var$gamma),]
  rownames(var) <- NULL
  class(var) <- c("gstatVariogram", "data.frame")
  
  return(var)
}

# Calculate  Poisson corss-semivariogram
pck.crossvariogram <-function(data.a, data.b, cases.a = 'Ya', cases.b = 'Yb' , risk.ab = 'Rab', pop = 'N', maxd = NULL, nbins = 15 , pop.rate = 100000){
  
  data.a = data.a[!(is.na(data.a[[risk.ab]])),]
  data.b = data.b[!(is.na(data.b[[risk.ab]])),]
  
  print(nrow(data.a))
  print(nrow(data.b))
  
  if(nrow(data.a) != nrow(data.b)){
    stop("data frames must have the same number of rows")
  }
  
  if(!(risk.ab %in% names(data.a)) || !(risk.ab %in% names(data.b))){
    stop("The shared risk column must be present in both data frames")
  }
  
  if(!identical(data.a[,pop],data.b[,pop])){
    stop("Population values must be the same for both data frames")
  }
  
  #Get coordinates
  coords <- st_coordinates(st_centroid(data.a))
  rownames(coords) <- data.a$FIPS
  
  # Global parameters
  data.a =  data.a %>% st_drop_geometry()
  data.b =  data.b %>% st_drop_geometry()
  No = nrow(data.a)
  
  # Define population sizes
  ya <- as.numeric(data.a[,cases.a])
  yb <- as.numeric(data.b[,cases.b])
  rab <- as.numeric(data.a[,risk.ab])
  # Define process
  n <- as.numeric(data.a[, pop])
  
  # Define total mean
  Mab <- mean(rab)
  
  # Create indexes to use in building of empirical semivariograms
  idx1 = rep(1:No, times = No:1)
  idx2 = c(1,unlist(sapply(1:No, function(i) (1:No)[-(1:(i - 1))])))
  
  # Reshape distances for unique pairs of points
  dist = as.matrix(dist(coords))
  d = c(dist[lower.tri(dist, diag = T)])
  
  # Calculate difference of unique pairs of points (centered and standardized)
  diff = ( pop.rate * pop.rate *((ya[idx1] / n[idx1]) - (ya[idx2] / n[idx2]))*((yb[idx1] / n[idx1]) - (yb[idx2] / n[idx2])))
  
  # Calculate weighted term
  w = (n[idx1] * n[idx2]) / (n[idx2] + n[idx1])
  # bias term
  bias = Mab * pop.rate
  
  # Determine maximum distance of estimator
  if (is.null(maxd)) {
    maxd = max(d) / 2
  }
  # Determine distances of tolerance regions for empirical semivariogram estimator
  bindist = seq(0, maxd, len = nbins + 1)
  # Determine classification cuts for unique pairs of points by distance
  dcuts = cut(d, breaks = bindist)
  
  
  # Calculate population-weighted empirical semivariogram
  np = unlist(lapply(split(d, dcuts), length), use.names = FALSE)
  middle = unlist(lapply(split(d, dcuts), median), use.names = FALSE)
  wbinds = split(w , dcuts)
  rbinds = split((w * diff - bias), dcuts)
  semivariance = unlist(lapply(1:nbins, function(x) {
    sum(rbinds[[x]]) / (2 * sum(wbinds[[x]]))
  }), use.names = FALSE)
  
  semivariance[semivariance < 0 ] <- NA
  
  
  
  # Transform varioram to object for gstat (easy manipulation of variogram)
  crossvar <- data.frame(as.numeric(np) , as.numeric(middle),semivariance,rep(0, nbins),
                         rep(0, nbins),rep(as.factor("var1.var2"), nbins))
  names(crossvar) <-  c("np", "dist", "gamma", "dir.hor", "dir.ver", "id")
  crossvar <- crossvar[!is.na(crossvar$gamma),]
  rownames(crossvar) <- NULL
  class(crossvar) <- c("gstatVariogram", "data.frame")
  return(crossvar)
}

# Calculate population weighted coordinates
weightedpop <- function(x, y, pop) {
  return(data.frame(x = weighted.mean(x, pop), y = weighted.mean(y, pop)))
}

# Plot empirical variogram
plot.var.pck = function(data){
  ggplot(data = data, aes(x = dist, y = gamma)) +
    geom_point(shape = 20, color = "black", size = 2) +
    #geom_line(color = "black", linetype = "dotdash") +
    labs(x = "distance (km)", y = "semivarince") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}

# Plot empirical vs fitted semivariogram
plot.pck.variogram = function(emp.var, fit.var, title){
  
  gamma.fit = covariance.model(fit.var$model[1])
  dist.fit = seq(0, max(emp.var$dist), length.out = 1000)
  cov.fit = gamma.fit(0, fit.var$range, fit.var$psill) - gamma.fit(dist.fit, fit.var$range, fit.var$psill)
  
  ggplot(data = emp.var, aes(x = dist, y = gamma)) +
    geom_point(shape = 20, color = "black", size = 2) +
    geom_line(color = "black", linetype = "dotdash") +
    labs(x = "distance (km)", y = "semivarince") +
    geom_line(data =  data.frame(dist = dist.fit, gamma = cov.fit)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(size=8),
          axis.title = element_text(size=8),
          plot.title = element_text(hjust = 0.5, size= 10),
          text = element_text(family = "Bookman Old Style")) +
    ggtitle(title)
}

# Get data frame to add column
data.get = function(disease) {
  switch(
    disease,
    "c" = STDs$Chlamydia,
    "g" = STDs$Gonorrhea,
    "hp" =  STDs$HIV.Prev,
    "hi" = STDs$HIV.Inc
  )
}

# Calculate shared risk and return values
data.name = function(disease) {
  switch(
    disease,
    "c" = "Chlamydia",
    "g" = "Gonorrhea",
    "hp" = "HIV.Prev",
    "hi" = "HIV.Inc"
  )
}

# Adjust LMC
lmc.poisson.cokrige.four <- function(var.a, var.b, var.c, var.d, crossvar.ab, crossvar.ac, crossvar.ad, crossvar.bc, crossvar.bd, crossvar.cd, data.a, data.b, data.c, data.d , var.params, fit.ranges = F){
  # Change id auxilliary variable
  var.b$id = as.factor("var2")
  var.c$id = as.factor("var3")
  var.d$id = as.factor("var4")
  crossvar.ac$id =  as.factor("var1.var3")
  crossvar.ad$id =  as.factor("var1.var4")
  crossvar.bc$id =  as.factor("var2.var3")
  crossvar.bd$id =  as.factor("var2.var4")
  crossvar.cd$id =  as.factor("var3.var4")
  
  names = c(expression(~ gamma[14]), expression(~ gamma[24]), expression(~ gamma[34]), expression(~ gamma[4]), expression(~ gamma[13]), expression(~ gamma[23]), expression(~ gamma[3]), expression(~ gamma[12]), expression(~ gamma[2]), expression(~ gamma[1]))
  my.settings <- list(
    strip.background=list(col="transparent"),
    strip.border=list(col="transparent")
  )
  
  # Integrate variograms
  variograms <- rbind(crossvar.ad ,crossvar.bd ,crossvar.cd ,var.d ,crossvar.ac, crossvar.bc, var.c, crossvar.ab, var.b, var.a)
  
  #Remove NA observations
  data.a <- data.a[!(is.na(data.a$Cases)),]
  data.b <- data.b[!(is.na(data.b$Cases)),]
  data.c <- data.c[!(is.na(data.c$Cases)),]
  data.d <- data.d[!(is.na(data.d$Cases)),]

  
  # Gstat object
  g <- gstat(  id = "var1", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = data.a$x ,y = data.a$y, c = data.a$Cases, p = data.a$Population), model = var.params$v1)
  g <- gstat(g,id = "var2", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = data.b$x ,y = data.b$y, c = data.b$Cases, p = data.b$Population), model = var.params$v2)
  g <- gstat(g,id = "var3", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = data.c$x ,y = data.c$y, c = data.c$Cases, p = data.c$Population), model = var.params$v3)
  g <- gstat(g,id = "var4", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = data.d$x ,y = data.d$y, c = data.d$Cases, p = data.d$Population), model = var.params$v4)
  
  g <- gstat(g, id = c("var1","var2"), model = var.params$v12)
  g <- gstat(g, id = c("var1","var3"), model = var.params$v13)
  g <- gstat(g, id = c("var1","var4"), model = var.params$v14)
  g <- gstat(g, id = c("var2","var3"), model = var.params$v23)
  g <- gstat(g, id = c("var2","var4"), model = var.params$v24)
  g <- gstat(g, id = c("var3","var4"), model = var.params$v34)
  
  
  fitted.lmc <- fit.lmc(v = variograms, g = g, fit.method = 2, fit.ranges = fit.ranges)
  if(!fit.ranges){
    fitted.lmc$model$var3$psill <- fit.3$psill
  }
  print(plot(variograms, fitted.lmc, pch = 20, col = "black", ylab = 'semivariance', par.settings =  my.settings, strip=strip.custom(factor.levels=names)))
  return(fitted.lmc)
}


# --- Poisson Cokriging helper function
poisson.cokrige.pred.four <- function(xao, xta, xtb, xtc, xtd, Cinv, cov.fun.1, cov.fun.2, data.a, data.b, data.c, data.d, pop.rate, par.a, par.ab, par.ac, par.ad) {
  
  
  # ----- Predict multiple locations
  ya <- data.a$Cases
  yb <- data.b$Cases
  yc <- data.c$Cases
  yd <- data.d$Cases
  na <- data.a$Population
  nb <- data.b$Population
  nc <- data.c$Population
  nd <- data.d$Population
  
  distXoa <- unlist(apply(xta, 1, distFun, x0 = xao))
  distXob <- unlist(apply(xtb, 1, distFun, x0 = xao))
  distXoc <- unlist(apply(xtc, 1, distFun, x0 = xao))
  distXod <- unlist(apply(xtd, 1, distFun, x0 = xao))
  posa <- which(distXoa == 0)
  posb <- which(distXob == 0)
  posc <- which(distXoc == 0)
  posd <- which(distXod == 0)
  
  #--- Covariances to prediction
  covaao <- (cov.fun.1(distXoa, par.a$range.1,  par.a$sill.1)  + cov.fun.2(distXoa, par.a$range.2,  par.a$sill.2))
  covbao <- (cov.fun.1(distXob, par.ab$range.1, par.ab$sill.1) + cov.fun.2(distXob, par.ab$range.2, par.ab$sill.2))
  covcao <- (cov.fun.1(distXoc, par.ac$range.1, par.ac$sill.1) + cov.fun.2(distXoc, par.ac$range.2, par.ac$sill.2))
  covdao <- (cov.fun.1(distXod, par.ad$range.1, par.ad$sill.1) + cov.fun.2(distXod, par.ad$range.2, par.ad$sill.2))
  
  if (length(posa) > 0) {
    covaao[posa] <- covaao[posa] + par.a$nugget.1  + par.a$nugget.2
    covbao[posb] <- covbao[posb] + par.ab$nugget.1 + par.ab$nugget.2
    covcao[posc] <- covcao[posc] + par.ac$nugget.1 + par.ac$nugget.2
    covdao[posd] <- covdao[posd] + par.ad$nugget.1 + par.ad$nugget.2
  }
  
  coTotal <- c(covaao, covbao, covcao, covdao, 1, 0, 0, 0)
  
  # --- Calculate weights
  lambda <- Cinv %*% coTotal
  
  # lambda.no.mu = lambda[1:(length(lambda) - 4)]
  # coTotal.no.mu = coTotal[1:(length(coTotal) - 4)]
  # 
  # sizes = c(nrow(data.a),nrow(data.b),nrow(data.c),nrow(data.d))
  # lambdas = split(lambda.no.mu, rep(1:length(sizes), sizes))
  # coTotals = split(coTotal.no.mu, rep(1:length(sizes), sizes))
  # 
  # lambda.hat = lapply(lambdas, function(x){mean(abs(x[x < 0]))})
  # covo.hat = lapply(1:length(coTotals), function(x){mean(coTotals[[x]][lambdas[[x]] < 0])})
  # 
  # lambda.no.zeros = lapply(1:length(coTotals), function(x){lambda.no.zero(lambdas[[x]], coTotals[[x]], x,  lambda.hat[[x]], covo.hat[[x]])})
  # tao = unlist(lambda.no.zeros)/sum(lambda.no.zeros[[1]])
  # print(tao)

  # predict
  pred.ra <- sum(lambda[1:(length(lambda)-4)] * c(ya/na, yb/nb, yc/nc, yd/nd) * pop.rate)
  # pred.ra <- sum(tao * c(ya/na, yb/nb, yc/nc, yd/nd) * pop.rate)
  
  if (pred.ra < 0) {
    pred.ra <- 0
  }
  
  #var.ra <- (cov.fun.1(0, par.a$range.1, par.a$sill.1) + cov.fun.2(0, par.a$range.2, par.a$sill.2)) - sum(c(tao, lambda[(length(lambda)-3):length(lambda)]) * coTotal)
  var.ra <- (cov.fun.1(0, par.a$range.1, par.a$sill.1) + cov.fun.2(0, par.a$range.2, par.a$sill.2)) - sum(lambda * coTotal)
  
  predictions <- c(c(xao[1]), c(xao[2]), pred.ra, var.ra)
  names(predictions) <- c('x', 'y', 'pred', 'var')
  return(predictions)
  
}

lambda.no.zero = function(ls, cs, pos, lam.hat, cov.hat) {
  tau = c()
  for(l in 1:length(ls)){
    if (ls[[l]] < 0) {
      tau = append(tau, 0)
    } else if (ls[[l]] > 0 & cs[[l]] < cov.hat & ls[[l]] < lam.hat) {
      tau = append(tau, 0)
    } else{
      tau = append(tau, ls[[l]])
    }
  }
  
  # if(pos == 1){
  #   tau.final = tau/sum(tau)
  # }else{
  #   tau.final = tau
  # }
  return(tau)
}

# --- Poisson Cokriging
poisson.cokrige.four <- function(data.a, data.b, data.c, data.d, coords.pred, lmc, pop.rate = 1) {


  #Get covariance models
  cov.fun.1 <- covariance.model(as.character(lmc$model$var1$model[1]))
  #cov.fun.2 <- covariance.model(as.character(lmc$model$var1$model[2]))
  cov.fun.2 <- function(h,a,b){return(0)}
  
  #Get covariance parameters
  par.a <- data.frame(sill.1 = lmc$model$var1$psill[1], range.1 = lmc$model$var1$range[1], nugget.1 = 0,
                      sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.b <- data.frame(sill.1 = lmc$model$var2$psill[1], range.1 = lmc$model$var2$range[1], nugget.1 = 0,
                      sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.c <- data.frame(sill.1 = lmc$model$var3$psill[1], range.1 = lmc$model$var3$range[1], nugget.1 = 0,
                      sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.d <- data.frame(sill.1 = lmc$model$var4$psill[1], range.1 = lmc$model$var4$range[1], nugget.1 = 0,
                      sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.ab <- data.frame(sill.1 = lmc$model$var1.var2$psill[1], range.1 = lmc$model$var1.var2$range[1], nugget.1 = 0,
                       sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.ac <- data.frame(sill.1 = lmc$model$var1.var3$psill[1], range.1 = lmc$model$var1.var3$range[1], nugget.1 = 0,
                       sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.ad <- data.frame(sill.1 = lmc$model$var1.var4$psill[1], range.1 = lmc$model$var1.var4$range[1], nugget.1 = 0,
                       sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.bc <- data.frame(sill.1 = lmc$model$var2.var3$psill[1], range.1 = lmc$model$var2.var3$range[1], nugget.1 = 0,
                       sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.bd <- data.frame(sill.1 = lmc$model$var2.var4$psill[1], range.1 = lmc$model$var2.var4$range[1], nugget.1 = 0,
                       sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.cd <- data.frame(sill.1 = lmc$model$var3.var4$psill[1], range.1 = lmc$model$var3.var4$range[1], nugget.1 = 0,
                       sill.2 = 0, range.2 = 0, nugget.2 = 0)
  
  #Remove NA observations
  data.a <- data.a[!(is.na(data.a$Cases)),]
  data.b <- data.b[!(is.na(data.b$Cases)),]
  data.c <- data.c[!(is.na(data.c$Cases)),]
  data.d <- data.d[!(is.na(data.d$Cases)),]
  
  #Get coordinates parameters
  coords.a <- cbind(data.a$x, data.a$y)
  coords.b <- cbind(data.b$x, data.b$y)
  coords.c <- cbind(data.c$x, data.c$y)
  coords.d <- cbind(data.d$x, data.d$y)
  
  if (!identical(coords.a,coords.pred)) {
    
    message("----------------- Predicting rates ----------------")
    
    # distances
    distMa <- as.matrix(dist(coords.a))
    distMb <- as.matrix(dist(coords.b))
    distMc <- as.matrix(dist(coords.c))
    distMd <- as.matrix(dist(coords.d))
    distMab <- as.matrix(proxy::dist(coords.a, coords.b))
    distMac <- as.matrix(proxy::dist(coords.a, coords.c))
    distMad <- as.matrix(proxy::dist(coords.a, coords.d))
    distMbc <- as.matrix(proxy::dist(coords.b, coords.c))
    distMbd <- as.matrix(proxy::dist(coords.b, coords.d))
    distMcd <- as.matrix(proxy::dist(coords.c, coords.d))
    
    # data sizes
    size.a <- nrow(data.a)
    size.b <- nrow(data.b)
    size.c <- nrow(data.c)
    size.d <- nrow(data.d)
    
    # ---- Error term W_{lk}
    Waa <- diag(sum(data.a$Cases)/sum(data.a$Population) * pop.rate / data.a$Population)
    Wbb <- diag(sum(data.b$Cases)/sum(data.b$Population) * pop.rate / data.b$Population)
    Wcc <- diag(sum(data.c$Cases)/sum(data.c$Population) * pop.rate / data.c$Population)
    Wdd <- diag(sum(data.d$Cases)/sum(data.d$Population) * pop.rate / data.d$Population)
    
    #--------------------------------------------------------------------------------------
    
    indexab <- which(distMab == 0, arr.ind = TRUE)[, 1]
    Wab <- matrix(0, dim(distMab)[1], dim(distMab)[2])
    Wab[which(distMab == 0, arr.ind = TRUE)] <- mean(data.a$Rhihp, na.rm = T) / data.a$Population[indexab]
    
    indexac <- which(distMac == 0, arr.ind = TRUE)[, 1]
    Wac <- matrix(0, dim(distMac)[1], dim(distMac)[2])
    Wac[which(distMac == 0, arr.ind = TRUE)] <- mean(data.a$Rhic, na.rm = T) / data.a$Population[indexac]
    
    indexad <- which(distMad == 0, arr.ind = TRUE)[, 1]
    Wad <- matrix(0, dim(distMad)[1], dim(distMad)[2])
    Wad[which(distMad == 0, arr.ind = TRUE)] <- mean(data.a$Rhig, na.rm = T) / data.a$Population[indexac]
    
    indexbc <- which(distMbc == 0, arr.ind = TRUE)[, 1]
    Wbc <- matrix(0, dim(distMbc)[1], dim(distMbc)[2])
    Wbc[which(distMbc == 0, arr.ind = TRUE)] <- mean(data.b$Rhpc, na.rm = T)/data.b$Population[indexbc]
    
    indexbd <- which(distMbd == 0, arr.ind = TRUE)[, 1]
    Wbd <- matrix(0, dim(distMbd)[1], dim(distMbd)[2])
    Wbd[which(distMbd == 0, arr.ind = TRUE)] <- mean(data.b$Rhpg, na.rm = T) / data.b$Population[indexbd]
    
    indexcd <- which(distMcd == 0, arr.ind = TRUE)[, 1]
    Wcd <- matrix(0, dim(distMcd)[1], dim(distMcd)[2])
    Wcd[which(distMcd == 0, arr.ind = TRUE)] <- mean(data.c$Rcg, na.rm = T) / data.c$Population[indexcd]
    
    # ---- Covariance matrices with W_{lk}
    Caa <- (cov.fun.1(distMa, par.a$range.1, par.a$sill.1) + cov.fun.2(distMa, par.a$range.2, par.a$sill.2)) + Waa
    Cbb <- (cov.fun.1(distMb, par.b$range.1, par.b$sill.1) + cov.fun.2(distMb, par.b$range.2, par.b$sill.2)) + Wbb
    Ccc <- (cov.fun.1(distMc, par.c$range.1, par.c$sill.1) + cov.fun.2(distMc, par.c$range.2, par.c$sill.2)) + Wcc
    Cdd <- (cov.fun.1(distMd, par.d$range.1, par.d$sill.1) + cov.fun.2(distMd, par.d$range.2, par.d$sill.2)) + Wdd
    Cab <- (cov.fun.1(distMab, par.ab$range.1, par.ab$sill.1) + cov.fun.2(distMab, par.ab$range.2, par.ab$sill.2)) + Wab
    Cac <- (cov.fun.1(distMac, par.ac$range.1, par.ac$sill.1) + cov.fun.2(distMac, par.ac$range.2, par.ac$sill.2)) + Wac
    Cad <- (cov.fun.1(distMad, par.ad$range.1, par.ad$sill.1) + cov.fun.2(distMad, par.ad$range.2, par.ad$sill.2)) + Wad
    Cbc <- (cov.fun.1(distMbc, par.bc$range.1, par.bc$sill.1) + cov.fun.2(distMbc, par.bc$range.2, par.bc$sill.2)) + Wbc
    Cbd <- (cov.fun.1(distMbd, par.bd$range.1, par.bd$sill.1) + cov.fun.2(distMbd, par.bd$range.2, par.bd$sill.2)) + Wbd
    Ccd <- (cov.fun.1(distMcd, par.cd$range.1, par.cd$sill.1) + cov.fun.2(distMcd, par.cd$range.2, par.cd$sill.2)) + Wcd
    Cba <- t(Cab)
    Cca <- t(Cac)
    Ccb <- t(Cbc)
    Cda <- t(Cad)
    Cdb <- t(Cbd)
    Cdc <- t(Ccd)
    
    # ----- Unbiasedness constrains
    Ctotal <- rbind(cbind(Caa, Cab, Cac, Cad),
                    cbind(Cba, Cbb, Cbc, Cbd),
                    cbind(Cca, Ccb, Ccc, Ccd),
                    cbind(Cda, Cdb, Cdc, Cdd))
    Cbias <- cbind(c(rep(1, size.a), rep(0, size.b), rep(0, size.c), rep(0, size.d)),
                   c(rep(0, size.a), rep(1, size.b), rep(0, size.c), rep(0, size.d)),
                   c(rep(0, size.a), rep(0, size.b), rep(1, size.c), rep(0, size.d)),
                   c(rep(0, size.a), rep(0, size.b), rep(0, size.c), rep(1, size.d)))
    Ctotal <- cbind(Ctotal, Cbias)
    Ctotal <- rbind(Ctotal, cbind(t(Cbias), matrix(0, 4, 4)))
    
    # ---- Inverse matrix
    Cinv <- solve(Ctotal)
    
    # predict
    preds <- as.data.frame(t(pbapply(coords.pred, 1, poisson.cokrige.pred.four, xta = coords.a, xtb = coords.b, xtc = coords.c, xtd = coords.d, Cinv = Cinv, 
                                     cov.fun.1 = cov.fun.1, cov.fun.2 = cov.fun.2, data.a, data.b, data.c, data.d, pop.rate, par.a, par.ab, par.ac, par.ad)))
  } else {
    
    message("----------------- Smoothing rates ----------------")
    
    fold = 1:nrow(data.a)
    preds = data.frame(matrix(as.numeric(NA), nrow(data.a), 7))
    names(preds) = c('rate.pred','rate.var', 'observed', 'residual', 'zscore', 'fold', 'lambda')
    
    folds = sort(unique(fold))
    progress_bar = txtProgressBar(min=0, max=length(folds), style = 3, char  = "+")
    
    #Validation
    for (i in folds) {
      
      sel = which(fold == i)
      data.cv = data.a[-sel, ]
      data.target = data.a[sel,]
      
      coords.a = cbind(data.cv$x, data.cv$y)
      coords.xo = cbind(data.target$x, data.target$y)
      
      # distances
      distMa <- as.matrix(dist(coords.a))
      distMb <- as.matrix(dist(coords.b))
      distMc <- as.matrix(dist(coords.c))
      distMd <- as.matrix(dist(coords.d))
      distMab <- as.matrix(proxy::dist(coords.a, coords.b))
      distMac <- as.matrix(proxy::dist(coords.a, coords.c))
      distMad <- as.matrix(proxy::dist(coords.a, coords.d))
      distMbc <- as.matrix(proxy::dist(coords.b, coords.c))
      distMbd <- as.matrix(proxy::dist(coords.b, coords.d))
      distMcd <- as.matrix(proxy::dist(coords.c, coords.d))
      
      # data sizes
      size.a = nrow(data.cv)
      size.b = nrow(data.b)
      size.c = nrow(data.c)
      size.d = nrow(data.d)
      
      # ---- Error term W_{lk}
      Waa <- diag(sum(data.cv$Cases)/sum(data.cv$Population) * pop.rate / data.cv$Population)
      Wbb <- diag(sum(data.b$Cases)/sum(data.b$Population) * pop.rate / data.b$Population)
      Wcc <- diag(sum(data.c$Cases)/sum(data.c$Population) * pop.rate / data.c$Population)
      Wdd <- diag(sum(data.d$Cases)/sum(data.d$Population) * pop.rate / data.d$Population)
      
      #--------------------------------------------------------------------------------------
      
      indexab <- which(distMab == 0, arr.ind = TRUE)[, 1]
      Wab <- matrix(0, dim(distMab)[1], dim(distMab)[2])
      Wab[which(distMab == 0, arr.ind = TRUE)] <- mean(data.cv$Rhihp, na.rm = T) / data.cv$Population[indexab]
      
      indexac <- which(distMac == 0, arr.ind = TRUE)[, 1]
      Wac <- matrix(0, dim(distMac)[1], dim(distMac)[2])
      Wac[which(distMac == 0, arr.ind = TRUE)] <- mean(data.cv$Rhic, na.rm = T) / data.cv$Population[indexac]
      
      indexad <- which(distMad == 0, arr.ind = TRUE)[, 1]
      Wad <- matrix(0, dim(distMad)[1], dim(distMad)[2])
      Wad[which(distMad == 0, arr.ind = TRUE)] <- mean(data.cv$Rhig, na.rm = T) / data.cv$Population[indexac]
      
      indexbc <- which(distMbc == 0, arr.ind = TRUE)[, 1]
      Wbc <- matrix(0, dim(distMbc)[1], dim(distMbc)[2])
      Wbc[which(distMbc == 0, arr.ind = TRUE)] <- mean(data.b$Rhpc, na.rm = T)/data.b$Population[indexbc]
      
      indexbd <- which(distMbd == 0, arr.ind = TRUE)[, 1]
      Wbd <- matrix(0, dim(distMbd)[1], dim(distMbd)[2])
      Wbd[which(distMbd == 0, arr.ind = TRUE)] <- mean(data.b$Rhpg, na.rm = T) / data.b$Population[indexbd]
      
      indexcd <- which(distMcd == 0, arr.ind = TRUE)[, 1]
      Wcd <- matrix(0, dim(distMcd)[1], dim(distMcd)[2])
      Wcd[which(distMcd == 0, arr.ind = TRUE)] <- mean(data.c$Rcg) / data.c$Population[indexcd]
      
      # ---- Covariance matrices with W_{lk}
      Caa <- (cov.fun.1(distMa, par.a$range.1, par.a$sill.1) + cov.fun.2(distMa, par.a$range.2, par.a$sill.2)) + Waa
      Cbb <- (cov.fun.1(distMb, par.b$range.1, par.b$sill.1) + cov.fun.2(distMb, par.b$range.2, par.b$sill.2)) + Wbb
      Ccc <- (cov.fun.1(distMc, par.c$range.1, par.c$sill.1) + cov.fun.2(distMc, par.c$range.2, par.c$sill.2)) + Wcc
      Cdd <- (cov.fun.1(distMd, par.d$range.1, par.d$sill.1) + cov.fun.2(distMd, par.d$range.2, par.d$sill.2)) + Wdd
      Cab <- (cov.fun.1(distMab, par.ab$range.1, par.ab$sill.1) + cov.fun.2(distMab, par.ab$range.2, par.ab$sill.2)) + Wab
      Cac <- (cov.fun.1(distMac, par.ac$range.1, par.ac$sill.1) + cov.fun.2(distMac, par.ac$range.2, par.ac$sill.2)) + Wac
      Cad <- (cov.fun.1(distMad, par.ad$range.1, par.ad$sill.1) + cov.fun.2(distMad, par.ad$range.2, par.ad$sill.2)) + Wad
      Cbc <- (cov.fun.1(distMbc, par.bc$range.1, par.bc$sill.1) + cov.fun.2(distMbc, par.bc$range.2, par.bc$sill.2)) + Wbc
      Cbd <- (cov.fun.1(distMbd, par.bd$range.1, par.bd$sill.1) + cov.fun.2(distMbd, par.bd$range.2, par.bd$sill.2)) + Wbd
      Ccd <- (cov.fun.1(distMcd, par.cd$range.1, par.cd$sill.1) + cov.fun.2(distMcd, par.cd$range.2, par.cd$sill.2)) + Wcd
      Cba <- t(Cab)
      Cca <- t(Cac)
      Ccb <- t(Cbc)
      Cda <- t(Cad)
      Cdb <- t(Cbd)
      Cdc <- t(Ccd)
      
      # ----- Unbiasedness constrains
      Ctotal <- rbind(cbind(Caa, Cab, Cac, Cad),
                      cbind(Cba, Cbb, Cbc, Cbd),
                      cbind(Cca, Ccb, Ccc, Ccd),
                      cbind(Cda, Cdb, Cdc, Cdd))
      Cbias <- cbind(c(rep(1, size.a), rep(0, size.b), rep(0, size.c), rep(0, size.d)),
                     c(rep(0, size.a), rep(1, size.b), rep(0, size.c), rep(0, size.d)),
                     c(rep(0, size.a), rep(0, size.b), rep(1, size.c), rep(0, size.d)),
                     c(rep(0, size.a), rep(0, size.b), rep(0, size.c), rep(1, size.d)))
      Ctotal <- cbind(Ctotal, Cbias)
      Ctotal <- rbind(Ctotal, cbind(t(Cbias), matrix(0, 4, 4)))
      
      
      # ---- Inverse matrix
      Cinv <- solve(Ctotal)
      
      #---- Prediction
      pred.target <- as.data.frame(t(apply(coords.xo, 1, poisson.cokrige.pred.four, xta = coords.a, xtb = coords.b, xtc = coords.c, xtd = coords.d, Cinv = Cinv, 
                                           cov.fun.1 = cov.fun.1, cov.fun.2 = cov.fun.2, data.a = data.cv, data.b = data.b, data.c = data.c, data.d = data.d,
                                           pop.rate = pop.rate, par.a = par.a, par.ab = par.ab, par.ac = par.ac, par.ad = par.ad)))
      
      preds[[1]][sel] = pred.target[[3]]
      preds[[2]][sel] = pred.target[[4]]
      preds[[3]][sel] = data.target$Cases/data.target$Population * pop.rate
      preds[[4]][sel] = as.numeric(pred.target[[3]] - preds[[3]][sel])
      preds[[5]][sel] = as.numeric(preds[[4]][sel]/sqrt(pred.target[[4]]))
      preds[[6]][sel] = fold[i]
      
      setTxtProgressBar(progress_bar, value = i)
      
    }
    close(progress_bar)
  }
  return(preds)
}

adjust.breaks = function(vector, max.val){
  new.vector = vector[1:(length(vector)-1)]
  new.vector = append(new.vector, round(max.val,1)+0.1)
  return(new.vector)
}


# --- Poisson kriging
poisson.krige.one <- function(data, coords.pred, var.fitted, pop.rate = 1){
  
  #Remove NA observations
  data <- data[!(is.na(data$Cases)),]
  
  #--- check model
  if (nrow(var.fitted) > 1){
    params <- data.frame(sill = var.fitted$psill[2],  range =  var.fitted$range[2],     nugget = var.fitted$psill[1])
    mod = as.character(var.fitted[2,]$model)
  }else{
    params <- data.frame(sill = var.fitted$psill,  range =  var.fitted$range,     nugget = 0)
    mod = as.character(var.fitted$model)
  }
  
  
  #--- Get covariance parameters
  
  
  #--- Get covariance model
  cov.fun = covariance.model(mod)
  coords = cbind(data$x, data$y)
  
  
  if(!identical(coords,coords.pred)){
    
    message("----------------- Predicting rates ----------------")
    
    # distances
    dist.m = as.matrix(dist(coords))
    
    # data sizes
    size = nrow(data)
    
    # ---- Error term W_{lk} 
    bias <- sum(data$Cases)/sum(data$Population) * pop.rate
    error.var <- diag(c(bias/data$Population , 0))
    nug.term <- diag(c(rep(params$nugget, times = size) , 0))
    
    # ---- Covariance matrix
    cov.mat <-  (cov.fun(dist.m,  params$range,  params$sill))
    
    # Add unbiassedness conditions
    cov.mat <- cbind(cov.mat, c(rep(1, size)))
    cov.mat <- rbind(cov.mat, c(rep(1, size + 1)))
    cov.mat[size+1,size+1] <- 0
    
    # ---- Inverse matrix
    c.inv <- solve(cov.mat + error.var + nug.term)
    
    # predict
    preds <- as.data.frame(t(pbapply(coords.pred, 1, poisson.krige.pred.one, xt = coords, cov.inv = c.inv, 
                                     cov.fun =  cov.fun, data, pop.rate, params)))
  }else{
    
    message("----------------- Smoothing rates ----------------")
    
    fold = 1:nrow(data)
    preds = data.frame(matrix(as.numeric(NA), nrow(data), 6))
    names(preds) = c('rate.pred','rate.var', 'observed', 'residual', 'zscore', 'fold')
    
    folds = sort(unique(fold))
    progress_bar = txtProgressBar(min=0, max=length(folds), style = 3, char  = "+")
    
    #Validation
    for (i in folds){
      
      # select observation
      sel = which(fold == i)
      data.cv = data[-sel, ]
      data.target = data[sel,]
      
      # get coordinates
      coords = cbind(data.cv$x, data.cv$y)
      coords.xo =  cbind(data.target$x, data.target$y)
      
      # calculate distances
      dist.m = as.matrix(dist(coords))
      
      # calculate rate parameters
      size = nrow(data.cv)
      n = data.cv$Population
      y = data.cv$Cases
      
      # bias term
      error.var = diag(c((sum(y) * pop.rate /sum(n))/n, 0))
      nug.term <- diag(c(rep(params$nugget, times = size) , 0))
      
      # covariance matrix
      cov.mat = (cov.fun(dist.m,  params$range,  params$sill))
      cov.mat = cbind(cov.mat, c(rep(1, size)))
      cov.mat = rbind(cov.mat, c(rep(1, size + 1)))
      cov.mat[size+1,size+1] <- 0
      c.inv = solve((cov.mat+error.var + nug.term))
      
      # predict
      pred.target <- as.data.frame(t(apply(coords.xo, 1, poisson.krige.pred.one, xt = coords, cov.inv = c.inv, 
                                           cov.fun =  cov.fun, data.cv, pop.rate, params)))
      
      preds[[1]][sel] = pred.target[[3]]
      preds[[2]][sel] = pred.target[[4]]
      preds[[3]][sel] = data.target$Cases/data.target$Population * pop.rate
      preds[[4]][sel] = as.numeric(pred.target[[3]] - preds[[3]][sel])
      preds[[5]][sel] = as.numeric(preds[[4]][sel]/sqrt(pred.target[[4]]))
      preds[[6]][sel] = fold[i]
      
      setTxtProgressBar(progress_bar, value = i)
    }
    
    close(progress_bar)
    
  }
  
  return(preds)
}

# --- Poisson kriging helper function
poisson.krige.pred.one <- function(xo, xt, cov.inv, cov.fun, data , pop.rate, params){
  
  # ----- Predict multiple locations
  ya <- data$Cases
  n <- data$Population
  
  #----- Size
  size = nrow(data)
  
  dist.xo <- apply(xt,1, distFun, x0 = xo)
  pos <- which(dist.xo == 0)
  
  #--- Covariances to prediction
  cov.xo <-  (cov.fun(as.matrix(dist.xo), params$range, params$sill))
  
  if(length(pos) > 0){
    cov.xo[pos] <- cov.xo[pos] + params$nugget
  }
  
  cov.xo <- rbind(cov.xo, 1)
  
  # --- Calculate weights
  lambda <- cov.inv%*%cov.xo
  
  # --- Predict
  pred  = sum(lambda[1:size]*(ya/n) * pop.rate)
  if(pred < 0){
    pred = 0
  }
  var = (cov.fun(0, params$range, params$sill) + params$nugget) - sum(lambda*cov.xo)
  
  predictions = c(c(xo[1]), c(xo[2]) , pred, var)
  names(predictions) <- c('x','y','pred','var')
  return(predictions)
}


# --- Poisson Cokriging
poisson.cokrige.three <- function(data.a, data.b, data.c, coords.pred, lmc, pop.rate = 1) {
  
  
  #Get covariance models
  cov.fun.1 <- covariance.model(as.character(lmc$model$var1$model[1]))
  #cov.fun.2 <- covariance.model(as.character(lmc$model$var1$model[2]))
  cov.fun.2 <- function(h,a,b){return(0)}
  
  #Get covariance parameters
  par.a <- data.frame(sill.1 = lmc$model$var1$psill[1], range.1 = lmc$model$var1$range[1], nugget.1 = 0,
                      sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.b <- data.frame(sill.1 = lmc$model$var2$psill[1], range.1 = lmc$model$var2$range[1], nugget.1 = 0,
                      sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.c <- data.frame(sill.1 = lmc$model$var3$psill[1], range.1 = lmc$model$var3$range[1], nugget.1 = 0,
                      sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.ab <- data.frame(sill.1 = lmc$model$var1.var2$psill[1], range.1 = lmc$model$var1.var2$range[1], nugget.1 = 0,
                       sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.ac <- data.frame(sill.1 = lmc$model$var1.var3$psill[1], range.1 = lmc$model$var1.var3$range[1], nugget.1 = 0,
                       sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.bc <- data.frame(sill.1 = lmc$model$var2.var3$psill[1], range.1 = lmc$model$var2.var3$range[1], nugget.1 = 0,
                       sill.2 = 0, range.2 = 0, nugget.2 = 0)
  
  #Remove NA observations
  data.a <- data.a[!(is.na(data.a$Cases)),]
  data.b <- data.b[!(is.na(data.b$Cases)),]
  data.c <- data.c[!(is.na(data.c$Cases)),]
  
  #Get coordinates parameters
  coords.a <- cbind(data.a$x, data.a$y)
  coords.b <- cbind(data.b$x, data.b$y)
  coords.c <- cbind(data.c$x, data.c$y)
  
  if (!identical(coords.a,coords.pred)) {
    
    message("----------------- Predicting rates ----------------")
    
    # distances
    distMa <- as.matrix(dist(coords.a))
    distMb <- as.matrix(dist(coords.b))
    distMc <- as.matrix(dist(coords.c))
    distMab <- as.matrix(proxy::dist(coords.a, coords.b))
    distMac <- as.matrix(proxy::dist(coords.a, coords.c))
    distMbc <- as.matrix(proxy::dist(coords.b, coords.c))
    
    # data sizes
    size.a <- nrow(data.a)
    size.b <- nrow(data.b)
    size.c <- nrow(data.c)
    
    # ---- Error term W_{lk}
    Waa <- diag(sum(data.a$Cases)/sum(data.a$Population) * pop.rate / data.a$Population)
    Wbb <- diag(sum(data.b$Cases)/sum(data.b$Population) * pop.rate / data.b$Population)
    Wcc <- diag(sum(data.c$Cases)/sum(data.c$Population) * pop.rate / data.c$Population)
    
    #--------------------------------------------------------------------------------------
    
    indexab <- which(distMab == 0, arr.ind = TRUE)[, 1]
    Wab <- matrix(0, dim(distMab)[1], dim(distMab)[2])
    Wab[which(distMab == 0, arr.ind = TRUE)] <- mean(data.a$Rhihp, na.rm = T) / data.a$Population[indexab]
    
    indexac <- which(distMac == 0, arr.ind = TRUE)[, 1]
    Wac <- matrix(0, dim(distMac)[1], dim(distMac)[2])
    Wac[which(distMac == 0, arr.ind = TRUE)] <- mean(data.a$Rhic, na.rm = T) / data.a$Population[indexac]
    
    indexbc <- which(distMbc == 0, arr.ind = TRUE)[, 1]
    Wbc <- matrix(0, dim(distMbc)[1], dim(distMbc)[2])
    Wbc[which(distMbc == 0, arr.ind = TRUE)] <- mean(data.b$Rhpc, na.rm = T)/data.b$Population[indexbc]
    
    
    # ---- Covariance matrices with W_{lk}
    Caa <- (cov.fun.1(distMa, par.a$range.1, par.a$sill.1) + cov.fun.2(distMa, par.a$range.2, par.a$sill.2)) + Waa
    Cbb <- (cov.fun.1(distMb, par.b$range.1, par.b$sill.1) + cov.fun.2(distMb, par.b$range.2, par.b$sill.2)) + Wbb
    Ccc <- (cov.fun.1(distMc, par.c$range.1, par.c$sill.1) + cov.fun.2(distMc, par.c$range.2, par.c$sill.2)) + Wcc
    Cab <- (cov.fun.1(distMab, par.ab$range.1, par.ab$sill.1) + cov.fun.2(distMab, par.ab$range.2, par.ab$sill.2)) + Wab
    Cac <- (cov.fun.1(distMac, par.ac$range.1, par.ac$sill.1) + cov.fun.2(distMac, par.ac$range.2, par.ac$sill.2)) + Wac
    Cbc <- (cov.fun.1(distMbc, par.bc$range.1, par.bc$sill.1) + cov.fun.2(distMbc, par.bc$range.2, par.bc$sill.2)) + Wbc
    Cba <- t(Cab)
    Cca <- t(Cac)
    Ccb <- t(Cbc)
    
    # ----- Unbiasedness constrains
    Ctotal <- rbind(cbind(Caa, Cab, Cac),
                    cbind(Cba, Cbb, Cbc),
                    cbind(Cca, Ccb, Ccc))
    Cbias <- cbind(c(rep(1, size.a), rep(0, size.b), rep(0, size.c)),
                   c(rep(0, size.a), rep(1, size.b), rep(0, size.c)),
                   c(rep(0, size.a), rep(0, size.b), rep(1, size.c)))
    Ctotal <- cbind(Ctotal, Cbias)
    Ctotal <- rbind(Ctotal, cbind(t(Cbias), matrix(0, 3, 3)))
    
    # ---- Inverse matrix
    Cinv <- solve(Ctotal)
    
    # predict
    preds <- as.data.frame(t(pbapply(coords.pred, 1, poisson.cokrige.pred.three, xta = coords.a, xtb = coords.b, xtc = coords.c, Cinv = Cinv, 
                                     cov.fun.1 = cov.fun.1, cov.fun.2 = cov.fun.2, data.a, data.b, data.c, pop.rate, par.a, par.ab, par.ac)))
  } else {
    
    message("----------------- Smoothing rates ----------------")
    
    fold = 1:nrow(data.a)
    preds = data.frame(matrix(as.numeric(NA), nrow(data.a), 7))
    names(preds) = c('rate.pred','rate.var', 'observed', 'residual', 'zscore', 'fold', 'lambda')
    
    folds = sort(unique(fold))
    progress_bar = txtProgressBar(min=0, max=length(folds), style = 3, char  = "+")
    
    #Validation
    for (i in folds) {
      
      sel = which(fold == i)
      data.cv = data.a[-sel, ]
      data.target = data.a[sel,]
      
      coords.a = cbind(data.cv$x, data.cv$y)
      coords.xo = cbind(data.target$x, data.target$y)
      
      # distances
      distMa <- as.matrix(dist(coords.a))
      distMb <- as.matrix(dist(coords.b))
      distMc <- as.matrix(dist(coords.c))
      distMab <- as.matrix(proxy::dist(coords.a, coords.b))
      distMac <- as.matrix(proxy::dist(coords.a, coords.c))
      distMbc <- as.matrix(proxy::dist(coords.b, coords.c))
      
      # data sizes
      size.a = nrow(data.cv)
      size.b = nrow(data.b)
      size.c = nrow(data.c)
      
      # ---- Error term W_{lk}
      Waa <- diag(sum(data.cv$Cases)/sum(data.cv$Population) * pop.rate / data.cv$Population)
      Wbb <- diag(sum(data.b$Cases)/sum(data.b$Population) * pop.rate / data.b$Population)
      Wcc <- diag(sum(data.c$Cases)/sum(data.c$Population) * pop.rate / data.c$Population)
      
      #--------------------------------------------------------------------------------------
      
      indexab <- which(distMab == 0, arr.ind = TRUE)[, 1]
      Wab <- matrix(0, dim(distMab)[1], dim(distMab)[2])
      Wab[which(distMab == 0, arr.ind = TRUE)] <- mean(data.cv$Rhihp, na.rm = T) / data.cv$Population[indexab]
      
      indexac <- which(distMac == 0, arr.ind = TRUE)[, 1]
      Wac <- matrix(0, dim(distMac)[1], dim(distMac)[2])
      Wac[which(distMac == 0, arr.ind = TRUE)] <- mean(data.cv$Rhic, na.rm = T) / data.cv$Population[indexac]
      
      indexbc <- which(distMbc == 0, arr.ind = TRUE)[, 1]
      Wbc <- matrix(0, dim(distMbc)[1], dim(distMbc)[2])
      Wbc[which(distMbc == 0, arr.ind = TRUE)] <- mean(data.b$Rhpc, na.rm = T)/data.b$Population[indexbc]
      
      # ---- Covariance matrices with W_{lk}
      Caa <- (cov.fun.1(distMa, par.a$range.1, par.a$sill.1) + cov.fun.2(distMa, par.a$range.2, par.a$sill.2)) + Waa
      Cbb <- (cov.fun.1(distMb, par.b$range.1, par.b$sill.1) + cov.fun.2(distMb, par.b$range.2, par.b$sill.2)) + Wbb
      Ccc <- (cov.fun.1(distMc, par.c$range.1, par.c$sill.1) + cov.fun.2(distMc, par.c$range.2, par.c$sill.2)) + Wcc
      Cab <- (cov.fun.1(distMab, par.ab$range.1, par.ab$sill.1) + cov.fun.2(distMab, par.ab$range.2, par.ab$sill.2)) + Wab
      Cac <- (cov.fun.1(distMac, par.ac$range.1, par.ac$sill.1) + cov.fun.2(distMac, par.ac$range.2, par.ac$sill.2)) + Wac
      Cbc <- (cov.fun.1(distMbc, par.bc$range.1, par.bc$sill.1) + cov.fun.2(distMbc, par.bc$range.2, par.bc$sill.2)) + Wbc
      Cba <- t(Cab)
      Cca <- t(Cac)
      Ccb <- t(Cbc)
      
      # ----- Unbiasedness constrains
      Ctotal <- rbind(cbind(Caa, Cab, Cac),
                      cbind(Cba, Cbb, Cbc),
                      cbind(Cca, Ccb, Ccc))
      Cbias <- cbind(c(rep(1, size.a), rep(0, size.b), rep(0, size.c)),
                     c(rep(0, size.a), rep(1, size.b), rep(0, size.c)),
                     c(rep(0, size.a), rep(0, size.b), rep(1, size.c)))
      Ctotal <- cbind(Ctotal, Cbias)
      Ctotal <- rbind(Ctotal, cbind(t(Cbias), matrix(0, 3, 3)))
      
      # ---- Inverse matrix
      Cinv <- solve(Ctotal)
      
      #---- Prediction
      pred.target <- as.data.frame(t(apply(coords.xo, 1, poisson.cokrige.pred.three, xta = coords.a, xtb = coords.b, xtc = coords.c, Cinv = Cinv, 
                                           cov.fun.1 = cov.fun.1, cov.fun.2 = cov.fun.2, data.a = data.cv, data.b = data.b, data.c = data.c,
                                           pop.rate = pop.rate, par.a = par.a, par.ab = par.ab, par.ac = par.ac)))
      
      preds[[1]][sel] = pred.target[[3]]
      preds[[2]][sel] = pred.target[[4]]
      preds[[3]][sel] = data.target$Cases/data.target$Population * pop.rate
      preds[[4]][sel] = as.numeric(pred.target[[3]] - preds[[3]][sel])
      preds[[5]][sel] = as.numeric(preds[[4]][sel]/sqrt(pred.target[[4]]))
      preds[[6]][sel] = fold[i]
      
      setTxtProgressBar(progress_bar, value = i)
      
    }
    close(progress_bar)
  }
  return(preds)
}

# --- Poisson Cokriging helper function
poisson.cokrige.pred.three <- function(xao, xta, xtb, xtc, Cinv, cov.fun.1, cov.fun.2, data.a, data.b, data.c, pop.rate, par.a, par.ab, par.ac) {
  
  # ----- Predict multiple locations
  ya <- data.a$Cases
  yb <- data.b$Cases
  yc <- data.c$Cases
  na <- data.a$Population
  nb <- data.b$Population
  nc <- data.c$Population
  
  distXoa <- unlist(apply(xta, 1, distFun, x0 = xao))
  distXob <- unlist(apply(xtb, 1, distFun, x0 = xao))
  distXoc <- unlist(apply(xtc, 1, distFun, x0 = xao))
  posa <- which(distXoa == 0)
  posb <- which(distXob == 0)
  posc <- which(distXoc == 0)
  
  #--- Covariances to prediction
  covaao <- (cov.fun.1(distXoa, par.a$range.1,  par.a$sill.1)  + cov.fun.2(distXoa, par.a$range.2,  par.a$sill.2))
  covbao <- (cov.fun.1(distXob, par.ab$range.1, par.ab$sill.1) + cov.fun.2(distXob, par.ab$range.2, par.ab$sill.2))
  covcao <- (cov.fun.1(distXoc, par.ac$range.1, par.ac$sill.1) + cov.fun.2(distXoc, par.ac$range.2, par.ac$sill.2))
  
  if (length(posa) > 0) {
    covaao[posa] <- covaao[posa] + par.a$nugget.1  + par.a$nugget.2
    covbao[posb] <- covbao[posb] + par.ab$nugget.1 + par.ab$nugget.2
    covcao[posc] <- covcao[posc] + par.ac$nugget.1 + par.ac$nugget.2$nugget.2
  }
  
  coTotal <- c(covaao, covbao, covcao, 1, 0, 0)
  
  # --- Calculate weights
  lambda <- Cinv %*% coTotal
  
  # lambda.no.mu = lambda[1:(length(lambda) - 4)]
  # coTotal.no.mu = coTotal[1:(length(coTotal) - 4)]
  # 
  # sizes = c(nrow(data.a),nrow(data.b),nrow(data.c),nrow(data.d))
  # lambdas = split(lambda.no.mu, rep(1:length(sizes), sizes))
  # coTotals = split(coTotal.no.mu, rep(1:length(sizes), sizes))
  # 
  # lambda.hat = lapply(lambdas, function(x){mean(abs(x[x < 0]))})
  # covo.hat = lapply(1:length(coTotals), function(x){mean(coTotals[[x]][lambdas[[x]] < 0])})
  # 
  # lambda.no.zeros = lapply(1:length(coTotals), function(x){lambda.no.zero(lambdas[[x]], coTotals[[x]], x,  lambda.hat[[x]], covo.hat[[x]])})
  # tao = unlist(lambda.no.zeros)/sum(lambda.no.zeros[[1]])
  # print(tao)
  
  # predict
  pred.ra <- sum(lambda[1:(length(lambda)-3)] * c(ya/na, yb/nb, yc/nc) * pop.rate)
  # pred.ra <- sum(tao * c(ya/na, yb/nb, yc/nc, yd/nd) * pop.rate)
  
  if (pred.ra < 0) {
    pred.ra <- 0
  }
  
  #var.ra <- (cov.fun.1(0, par.a$range.1, par.a$sill.1) + cov.fun.2(0, par.a$range.2, par.a$sill.2)) - sum(c(tao, lambda[(length(lambda)-3):length(lambda)]) * coTotal)
  var.ra <- (cov.fun.1(0, par.a$range.1, par.a$sill.1) + cov.fun.2(0, par.a$range.2, par.a$sill.2)) - sum(lambda * coTotal)
  
  predictions <- c(c(xao[1]), c(xao[2]), pred.ra, var.ra)
  names(predictions) <- c('x', 'y', 'pred', 'var')
  return(predictions)
  
}


# --- Poisson Cokriging
poisson.cokrige.bi <- function(data.a, data.b, coords.pred, lmc, pop.rate = 1) {
  
  
  #Get covariance models
  cov.fun.1 <- covariance.model(as.character(lmc$model$var1$model[1]))
  #cov.fun.2 <- covariance.model(as.character(lmc$model$var1$model[2]))
  cov.fun.2 <- function(h,a,b){return(0)}
  
  #Get covariance parameters
  par.a <- data.frame(sill.1 = lmc$model$var1$psill[1], range.1 = lmc$model$var1$range[1], nugget.1 = 0,
                      sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.b <- data.frame(sill.1 = lmc$model$var2$psill[1], range.1 = lmc$model$var2$range[1], nugget.1 = 0,
                      sill.2 = 0, range.2 = 0, nugget.2 = 0)
  par.ab <- data.frame(sill.1 = lmc$model$var1.var2$psill[1], range.1 = lmc$model$var1.var2$range[1], nugget.1 = 0,
                       sill.2 = 0, range.2 = 0, nugget.2 = 0)
  
  #Remove NA observations
  data.a <- data.a[!(is.na(data.a$Cases)),]
  data.b <- data.b[!(is.na(data.b$Cases)),]
  
  #Get coordinates parameters
  coords.a <- cbind(data.a$x, data.a$y)
  coords.b <- cbind(data.b$x, data.b$y)
  
  if (!identical(coords.a,coords.pred)) {
    
    message("----------------- Predicting rates ----------------")
    
    # distances
    distMa <- as.matrix(dist(coords.a))
    distMb <- as.matrix(dist(coords.b))
    distMab <- as.matrix(proxy::dist(coords.a, coords.b))
    
    # data sizes
    size.a <- nrow(data.a)
    size.b <- nrow(data.b)
    
    # ---- Error term W_{lk}
    Waa <- diag(sum(data.a$Cases)/sum(data.a$Population) * pop.rate / data.a$Population)
    Wbb <- diag(sum(data.b$Cases)/sum(data.b$Population) * pop.rate / data.b$Population)
    
    #--------------------------------------------------------------------------------------
    
    indexab <- which(distMab == 0, arr.ind = TRUE)[, 1]
    Wab <- matrix(0, dim(distMab)[1], dim(distMab)[2])
    Wab[which(distMab == 0, arr.ind = TRUE)] <- mean(data.a$Rhihp, na.rm = T) / data.a$Population[indexab]

    
    
    # ---- Covariance matrices with W_{lk}
    Caa <- (cov.fun.1(distMa, par.a$range.1, par.a$sill.1) + cov.fun.2(distMa, par.a$range.2, par.a$sill.2)) + Waa
    Cbb <- (cov.fun.1(distMb, par.b$range.1, par.b$sill.1) + cov.fun.2(distMb, par.b$range.2, par.b$sill.2)) + Wbb
    Cab <- (cov.fun.1(distMab, par.ab$range.1, par.ab$sill.1) + cov.fun.2(distMab, par.ab$range.2, par.ab$sill.2)) + Wab
    Cba <- t(Cab)
    
    # ----- Unbiasedness constrains
    Ctotal <- rbind(cbind(Caa, Cab),
                    cbind(Cba, Cbb))
    Cbias <- cbind(c(rep(1, size.a), rep(0, size.b)),
                   c(rep(0, size.a), rep(1, size.b)))
    Ctotal <- cbind(Ctotal, Cbias)
    Ctotal <- rbind(Ctotal, cbind(t(Cbias), matrix(0, 2, 2)))
    
    # ---- Inverse matrix
    Cinv <- solve(Ctotal)
    
    # predict
    preds <- as.data.frame(t(pbapply(coords.pred, 1, poisson.cokrige.pred.bi, xta = coords.a, xtb = coords.b, Cinv = Cinv, 
                                     cov.fun.1 = cov.fun.1, cov.fun.2 = cov.fun.2, data.a = data.a, data.b = data.b, pop.rate, par.a, par.ab)))
  } else {
    
    message("----------------- Smoothing rates ----------------")
    
    fold = 1:nrow(data.a)
    preds = data.frame(matrix(as.numeric(NA), nrow(data.a), 7))
    names(preds) = c('rate.pred','rate.var', 'observed', 'residual', 'zscore', 'fold', 'lambda')
    
    folds = sort(unique(fold))
    progress_bar = txtProgressBar(min=0, max=length(folds), style = 3, char  = "+")
    
    #Validation
    for (i in folds) {
      
      sel = which(fold == i)
      data.cv = data.a[-sel, ]
      data.target = data.a[sel,]
      
      coords.a = cbind(data.cv$x, data.cv$y)
      coords.xo = cbind(data.target$x, data.target$y)
      
      # distances
      # distances
      distMa <- as.matrix(dist(coords.a))
      distMb <- as.matrix(dist(coords.b))
      distMab <- as.matrix(proxy::dist(coords.a, coords.b))
      
      # data sizes
      size.a = nrow(data.cv)
      size.b = nrow(data.b)
      
      # ---- Error term W_{lk}
      Waa <- diag(sum(data.cv$Cases)/sum(data.cv$Population) * pop.rate / data.cv$Population)
      Wbb <- diag(sum(data.b$Cases)/sum(data.b$Population) * pop.rate / data.b$Population)
      
      #--------------------------------------------------------------------------------------
      
      indexab <- which(distMab == 0, arr.ind = TRUE)[, 1]
      Wab <- matrix(0, dim(distMab)[1], dim(distMab)[2])
      Wab[which(distMab == 0, arr.ind = TRUE)] <- mean(data.cv$Rhihp, na.rm = T) / data.cv$Population[indexab]
      
      
      # ---- Covariance matrices with W_{lk}
      Caa <- (cov.fun.1(distMa, par.a$range.1, par.a$sill.1) + cov.fun.2(distMa, par.a$range.2, par.a$sill.2)) + Waa
      Cbb <- (cov.fun.1(distMb, par.b$range.1, par.b$sill.1) + cov.fun.2(distMb, par.b$range.2, par.b$sill.2)) + Wbb
      Cab <- (cov.fun.1(distMab, par.ab$range.1, par.ab$sill.1) + cov.fun.2(distMab, par.ab$range.2, par.ab$sill.2)) + Wab
      Cba <- t(Cab)
      
      # ----- Unbiasedness constrains
      Ctotal <- rbind(cbind(Caa, Cab),
                      cbind(Cba, Cbb))
      Cbias <- cbind(c(rep(1, size.a), rep(0, size.b)),
                     c(rep(0, size.a), rep(1, size.b)))
      Ctotal <- cbind(Ctotal, Cbias)
      Ctotal <- rbind(Ctotal, cbind(t(Cbias), matrix(0, 2, 2)))
      
      # ---- Inverse matrix
      Cinv <- solve(Ctotal)
      
      #---- Prediction
      pred.target <- as.data.frame(t(apply(coords.xo, 1, poisson.cokrige.pred.bi, xta = coords.a, xtb = coords.b, Cinv = Cinv, 
                                           cov.fun.1 = cov.fun.1, cov.fun.2 = cov.fun.2, data.a = data.cv, data.b = data.b,
                                           pop.rate = pop.rate, par.a = par.a, par.ab = par.ab)))
      
      preds[[1]][sel] = pred.target[[3]]
      preds[[2]][sel] = pred.target[[4]]
      preds[[3]][sel] = data.target$Cases/data.target$Population * pop.rate
      preds[[4]][sel] = as.numeric(pred.target[[3]] - preds[[3]][sel])
      preds[[5]][sel] = as.numeric(preds[[4]][sel]/sqrt(pred.target[[4]]))
      preds[[6]][sel] = fold[i]
      
      setTxtProgressBar(progress_bar, value = i)
      
    }
    close(progress_bar)
  }
  return(preds)
}

# --- Poisson Cokriging helper function
poisson.cokrige.pred.bi <- function(xao, xta, xtb, Cinv, cov.fun.1, cov.fun.2, data.a, data.b, pop.rate, par.a, par.ab) {
  
  # xao = coords.a[1,]
  # xta = coords.a
  # xtb = coords.b
  # xtc = coords.c
  # xtd = coords.d
  
  # ----- Predict multiple locations
  ya <- data.a$Cases
  yb <- data.b$Cases
  na <- data.a$Population
  nb <- data.b$Population
  
  distXoa <- unlist(apply(xta, 1, distFun, x0 = xao))
  distXob <- unlist(apply(xtb, 1, distFun, x0 = xao))
  posa <- which(distXoa == 0)
  posb <- which(distXob == 0)
  
  #--- Covariances to prediction
  covaao <- (cov.fun.1(distXoa, par.a$range.1,  par.a$sill.1)  + cov.fun.2(distXoa, par.a$range.2,  par.a$sill.2))
  covbao <- (cov.fun.1(distXob, par.ab$range.1, par.ab$sill.1) + cov.fun.2(distXob, par.ab$range.2, par.ab$sill.2))
  
  if (length(posa) > 0) {
    covaao[posa] <- covaao[posa] + par.a$nugget.1  + par.a$nugget.2
    covbao[posb] <- covbao[posb] + par.ab$nugget.1 + par.ab$nugget.2
  }
  
  coTotal <- c(covaao, covbao, 1, 0)
  
  # --- Calculate weights
  lambda <- Cinv %*% coTotal
  
  # lambda.no.mu = lambda[1:(length(lambda) - 4)]
  # coTotal.no.mu = coTotal[1:(length(coTotal) - 4)]
  # 
  # sizes = c(nrow(data.a),nrow(data.b),nrow(data.c),nrow(data.d))
  # lambdas = split(lambda.no.mu, rep(1:length(sizes), sizes))
  # coTotals = split(coTotal.no.mu, rep(1:length(sizes), sizes))
  # 
  # lambda.hat = lapply(lambdas, function(x){mean(abs(x[x < 0]))})
  # covo.hat = lapply(1:length(coTotals), function(x){mean(coTotals[[x]][lambdas[[x]] < 0])})
  # 
  # lambda.no.zeros = lapply(1:length(coTotals), function(x){lambda.no.zero(lambdas[[x]], coTotals[[x]], x,  lambda.hat[[x]], covo.hat[[x]])})
  # tao = unlist(lambda.no.zeros)/sum(lambda.no.zeros[[1]])
  # print(tao)
  
  # predict
  pred.ra <- sum(lambda[1:(length(lambda)-2)] * c(ya/na, yb/nb) * pop.rate)
  # pred.ra <- sum(tao * c(ya/na, yb/nb, yc/nc, yd/nd) * pop.rate)
  
  if (pred.ra < 0) {
    pred.ra <- 0
  }
  
  #var.ra <- (cov.fun.1(0, par.a$range.1, par.a$sill.1) + cov.fun.2(0, par.a$range.2, par.a$sill.2)) - sum(c(tao, lambda[(length(lambda)-3):length(lambda)]) * coTotal)
  var.ra <- (cov.fun.1(0, par.a$range.1, par.a$sill.1) + cov.fun.2(0, par.a$range.2, par.a$sill.2)) - sum(lambda * coTotal)
  
  predictions <- c(c(xao[1]), c(xao[2]), pred.ra, var.ra)
  names(predictions) <- c('x', 'y', 'pred', 'var')
  return(predictions)
  
}

