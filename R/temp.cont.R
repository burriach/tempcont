temp.cont <- function (model, driver, random, plot=TRUE){
  dataframe<- model$data
  resp<- as.character(model$call$fixed[[2]])

  ## Use all fixed effects terms if missing
  if (missing(driver)){
    driver<- as.character(model$call$fixed[[3]])
    driver<- sapply(driver, strsplit, split="\\s*\\+\\s*")
    driver<- do.call(c, driver)
    driver<- intersect(driver, names(dataframe))
  }

  # random<- as.character(model$call$random[[2]])
  ## TODO: not exactly the same value as expected by the current input
  ## use update() instead and forget dataframe, resp and random

  dataframe$site <- dataframe[[random]]
  dataframe$seeds <- dataframe[[resp]]
  newdat <- dataframe
  newdat$seeds <- newdat[[resp]]
  newdat$pred.nbp <- as.matrix(predict(model, type="response"))
  pred.mod.nbp <- aggregate(newdat$pred.nbp, list(year=newdat$year), FUN=mean)
  pred.n <- aggregate(newdat$pred.nbp, list(year=newdat$year), FUN=length)
  real.seeds <- aggregate(newdat$seeds, list(year=newdat$year), FUN=mean)

  if (plot){
    plot(x ~ year, data=real.seeds, type="l", lwd=2)
    lines(V1 ~ year, data=pred.mod.nbp, lwd=2, col="blue")
  }

  # ctrl <- lmeControl(opt='optim')
  # formFix<- formula(paste(resp, "~", "year"))
  # formRand<- formula(paste("year", "~", random))

  tt<- try(mod.real<- update(model, fixed=seeds ~ year, random=~ year | site, control=list(maxIter=10000, niterEM=10000), method="REML"), silent=TRUE)
  # tt <- try(mod.real <-lme (seeds ~ year,  random=~ year|site, data=dataframe, correlation = corAR1(form=~year|site), control=list(maxIter=10000, niterEM=10000), method="REML"), silent = T)
  if (class(tt)=="try-error") {mod.real <-lme (seeds ~ year,  random=~ year|site, data=dataframe, correlation = corAR1(form=~year|site), control=list(maxIter=10000, niterEM=10000, opt='optim'), method="REML")}

  slope.real <- summary (mod.real)$tTable[2,1]
  slope.real.se <- summary (mod.real)$tTable[2,2]
  slope.real.p <- summary (mod.real)$tTable[2,5]

  tt <- try(mod.pred.full <- update(model, fixed=pred.nbp ~ year,  random=~ year | site, data=newdat, control=list(maxIter=10000, niterEM=10000), method="REML"), silent=TRUE)
  # tt <- try(mod.pred.full <-lme (pred.nbp ~ year,  random=~ year |site, data=newdat, correlation = corAR1(form=~year|site), control=list(maxIter=10000, niterEM=10000), method="REML"), silent=T)
  if (class(tt)=="try-error") {tt <-try(mod.pred.full <-lme (pred.nbp ~ year,  random=~ year |site, data=newdat, correlation = corAR1(form=~year|site), control=list(maxIter=10000, niterEM=10000, opt='optim'), method="REML"), silent=T)}
  if (class(tt)=="try-error") {tt <-try(mod.pred.full <-lme (pred.nbp ~ year,  random=~ year |site, data=newdat, control=list(maxIter=10000, niterEM=10000), method="ML"), silent=T)}
  if (class(tt)=="try-error") {mod.pred.full <-lme (pred.nbp ~ year,  random=~ year |site, data=newdat, control=list(maxIter=10000, niterEM=10000, opt='optim'), method="ML")}


  mod.contr <- summary (mod.pred.full)$tTable[2,1]
  mod.contr.se <- summary (mod.pred.full)$tTable[2,2]
  mod.contr.t <- summary (mod.pred.full)$tTable[2,4]
  mod.contr.p <- summary (mod.pred.full)$tTable[2,5]

  # tot menys x var contribution
  # newdat <- cbind(as.data.frame(model.matrix(model)), "seeds"=dataframe[["seeds"]], "year"=dataframe[["year"]], "site"=dataframe[["site"]])
  out<- lapply(driver, function(dr){
    newdat <- dataframe

    #calculem mediana de predictora per pixel
    z <- aggregate(newdat[[dr]], by=list("site"=newdat$site), median)
    colnames(z) <- c("site", "varx.pred")
    newdat <- merge(newdat, z, by="site", all.x=T, all.y=F)
    newdat[[dr]] <- newdat$varx.pred

    # newdat[[dr]] <- mean(newdat[[dr]])
    newdat$pred.nbp <- as.numeric(predict(model, newdata=newdat, type="response"))
    pred.mod.nbp.cdioxide <- aggregate(newdat$pred.nbp, list(year=newdat$year), FUN=mean)

    tt <- try(mod.varx <- update(model, fixed=pred.nbp ~ year,  random=~year|site, data=newdat, control=list(maxIter=10000, niterEM=10000), method="REML"), silent=TRUE)
    # tt <- try(mod.varx <-lme (pred.nbp ~ year,  random=~year|site, data=newdat, correlation = corAR1(form=~year|site), control=list(maxIter=10000, niterEM=10000), method="REML"), silent=T)
    if (class(tt)=="try-error") {tt <- try(mod.varx <-lme (pred.nbp ~ year,  random=~year|site, data=newdat, correlation = corAR1(form=~year|site), control=list(maxIter=10000, niterEM=10000, opt='optim'), method="REML"),silent=T)}
    if (class(tt)=="try-error") {tt <- try(mod.varx <-lme (pred.nbp ~ year,  random=~1|site, data=newdat, correlation = corAR1(form=~year|site), control=list(maxIter=10000, niterEM=10000, opt='optim'), method="ML"),silent=T)}
    if (class(tt)=="try-error") {mod.varx <-lme (pred.nbp ~ year,  random=~1|site, data=newdat, correlation = corAR1(form=~year|site), control=list(maxIter=10000, niterEM=10000), method="ML")}

    if (plot){
      lines(x ~ year, data=pred.mod.nbp.cdioxide, lwd=2.5, col="red")
    }
    cdiox.contr <- summary (mod.varx)$tTable[2,1]
    cdiox.contr.se <- summary (mod.varx)$tTable[2,2]

    a <- as.numeric(mod.contr - cdiox.contr) # contr
    b <- as.numeric(sqrt((mod.contr.se^2)+(cdiox.contr.se^2))) # se
    c <- a/b

    newdat <- dataframe
    newdat$pppp <- dataframe[[dr]]
    zz <- try(mod.trendx <- update(model, fixed=pppp ~ year,  random=~year|site, data=newdat, control=list(maxIter=10000, niterEM=10000), method="REML"), silent=TRUE)
    # zz <- try(mod.trendx <-lme (pppp ~ year,  random=~year|site, data=newdat, correlation = corAR1(form=~year|site), control=list(maxIter=10000, niterEM=10000), method="REML"), silent=T)
    if (class(zz)=="try-error") {zz <- try(mod.trendx <-lme (pppp ~ year,  random=~year|site, data=newdat, correlation = corAR1(form=~year|site), control=list(maxIter=10000, niterEM=10000, opt='optim'), method="REML"),silent=T)}
    if (class(zz)=="try-error") {zz <- try(mod.trendx <-lme (pppp ~ year,  random=~year|site, data=newdat, correlation = corAR1(form=~year|site), control=list(maxIter=10000, niterEM=10000, opt='optim'), method="ML"),silent=T)}
    if (class(zz)=="try-error") {zz <- try(mod.trendx <-lme (pppp ~ year,  random=~year|site, data=newdat, correlation = corAR1(form=~year|site), control=list(maxIter=10000, niterEM=10000), method="ML"),silent=T)}
    if (class(zz)=="try-error") {monguer <- aggregate(newdat$pppp, by=list(year=newdat$year), mean)
    colnames(monguer)<- c("year", "x")
    mod.trendx <-gls (x ~ year,  data=monguer, control=list(maxIter=10000, niterEM=10000),correlation = corAR1(form=~1), method="ML")}


    d <- summary (mod.trendx)$tTable[2,1]
    e <- summary (mod.trendx)$tTable[2,2]
    f <- d/e

    h <- (a/d)*sqrt(((e/d)^2)+((b/a)^2))
    res <- data.frame(mod.slope=mod.contr, mod.slope.se=mod.contr.se, mod.slope.t=mod.contr.t, mod.slope.p=mod.contr.p,
                      temp.contr=a, temp.contr.se=b, temp.contr.t=c, temp.contr.p= (1-pt(abs(c),df=2*nrow(pred.mod.nbp.cdioxide)-4)),
                      pred.trend=d, pred.trend.se=e, pred.trend.t=f, pred.trend.p= (1-pt(abs(f),df=2*nrow(pred.mod.nbp.cdioxide)-4)),
                      sensit=a/d, sensit.se=h, sensit.t=(a/d)/h, sensit.p= (1-pt(abs((a/d)/h),df=2*nrow(pred.mod.nbp.cdioxide)-4)),
                      row.names=dr)
    return(res)
  })

  res<- do.call(rbind, out)

  return(res)
} # works well
