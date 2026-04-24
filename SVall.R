SVall <- function(exposure, numerator = NULL, denominator,
                  id, tstart, timevar, data, trunc= NULL){
    
#######################################################################################
#   Filename    :	SVall.R
#
#   Description :   modifed ipwtm function
#                   computes the weights for all time points in the dataset, 
#                   that is event times and time changes in the treatment of interest status
#    
#   Required package  :    survival
#
#   Usage       :   SVall(exposure, numerator, denominator,
#                         id, tstart, timevar, data, trunc)
#                   Same usage as in ipwtm (of the ipw package) with 
#                   family = survival & type = all  
#    
#   Value : same as the ones of ipwtm
########################################################################################
    
    tempcall <- match.call()  
    
    tempdat <- data.frame(id = data[ , as.character(tempcall$id)],
                          tstart = data[ , as.character(tempcall$tstart)],     
                          timevar = data[ , as.character(tempcall$timevar)], 
                          exposure = data[ , eval(tempcall$exposure)])
    
    row.data <- nrow(tempdat)
    # if numerator = NULL (unstabilized weights), then the weight equals 1 
    if (is.null(tempcall$numerator))
        tempdat$w.numerator <- 1
    # Numerator
    else {
        # Andersen-Gill model
        mod1 <- coxph(formula = eval(parse(text = paste("Surv(",
                                                        deparse(tempcall$tstart),
                                                        ",",
                                                        deparse(tempcall$timevar, width.cutoff = 500),  
                                                        ",",
                                                        eval(tempcall$exposure),
                                                        ")",eval(tempcall$numerator),
                                                        #deparse(tempcall$numerator, width.cutoff = 500),
                                                        "+ cluster (",
                                                        deparse(tempcall$id, width.cutoff = 500),
                                                        ")"))),
                        data = data)
        #print(summary(mod1))
        # Relative risks
        tempdat$risk.numerator <- vector("numeric", nrow(tempdat))
        tempdat$risk.numerator <- predict(mod1, data, type = "risk")
        # Baseline hazard
        tauxBasCum <- basehaz(mod1, centered = TRUE)
        # all tstart and timevar
        tempo.times  <- sort(unique(c(tempdat$tstart[tempdat$tstart >= 0], tempdat$timevar)))
        if( ! all(tempo.times %in% tauxBasCum$time)){
            tauxBasCum <- data.frame(time = tempo.times,
                                     hazard = approx(x = tauxBasCum$time, y = tauxBasCum$hazard,
                                                     xout = tempo.times, 
                                                     method = "constant", rule = 2)$y)
        }
        # split by patient 
        tabi   <- split(tempdat, tempdat$id)
        L.tabi <- length(tabi)
        tablist <- lapply(1:L.tabi, function(i){
            lignes.tabi <- nrow(tabi[[i]])
            # for each [t1;t2]: \Lambda(t2) - \Lambda(t1)
            bashaz.numerator <- vector(length = lignes.tabi)
            # when time = t1
            taux.cum.start.numerator <- approx(x = tauxBasCum$time, y = tauxBasCum$hazard,
                                               xout = tabi[[i]][ ,"tstart"])$y 
            # when time = t2
            taux.cum.stop.numerator <- approx(x = tauxBasCum$time, y = tauxBasCum$hazard,
                                              xout = tabi[[i]][ ,"timevar"])$y

            bashaz.numerator[1] <- taux.cum.stop.numerator[1]
            if(lignes.tabi > 1){
                bashaz.numerator[2:lignes.tabi] <- (taux.cum.stop.numerator - taux.cum.start.numerator)[-1]
            }

            indices0 <- which(tabi[[i]]$exposure == 0) # non-exposed patients
            indices1 <- which(tabi[[i]]$exposure == 1) # exposed patients

            tabi[[i]]$p.numerator <- vector("numeric", lignes.tabi)
            # when exposition = 0 => exp( - (\Lambda(t2) - \Lambda(t1)) * RR )
            tabi[[i]]$p.numerator[indices0] <- exp( - bashaz.numerator[indices0] * tabi[[i]]$risk.numerator[indices0])
            # when exposition = 1 => 1 - exp ( - (\Lambda(t2) - \Lambda(t1)) * RR )
            tabi[[i]]$p.numerator[indices1] <- 1 - exp( - bashaz.numerator[indices1] * tabi[[i]]$risk.numerator[indices1])
            # weights obtained by the cumulative product
            tabi[[i]]$w.numerator <- cumprod(tabi[[i]]$p.numerator)
            
            return(tabi[[i]])
        } )
        rm(tauxBasCum)
        
        result.numerator <- do.call(rbind, tablist)
    }
    # Denominator
    # Andersen-Gill model 
    mod2 <- coxph(formula = eval(parse(text = paste("Surv(",
                                                    deparse(tempcall$tstart, width.cutoff = 500),
                                                    ",",
                                                    deparse(tempcall$timevar, width.cutoff = 500),
                                                    ",",
                                                    eval(tempcall$exposure),
                                                    ")",eval(tempcall$denominator),
                                                    #deparse(tempcall$denominator, width.cutoff = 500),
                                                    "+ cluster (",
                                                    deparse(tempcall$id, width.cutoff = 500),
                                                    ")"))),
                   data = data)
    
    # Relative risks
    tempdat$risk.denominator <- vector("numeric", nrow(tempdat))
    tempdat$risk.denominator <- predict(mod2, data, type = "risk") 
    
    # Baseline hazard
    tauxBasCum <- basehaz(mod2, centered = TRUE)
    # all tstart and timevar
    tempo.times  <- sort( unique( c( tempdat$tstart[tempdat$tstart >= 0] , tempdat$timevar ) ) )
    tempo.times  <- sort( unique( c( tempdat$tstart[ tempdat$tstart >= 0 ] , tempdat$timevar ) ) )
    
    if( ! all( tempo.times %in% tauxBasCum$time ) ){        
        tauxBasCum <- data.frame( time = tempo.times,
                                hazard = approx( x = tauxBasCum$time, y = tauxBasCum$hazard,
                                                xout = tempo.times, 
                                                method = "constant", rule = 2 )$y )  
    }
    
    # split by patient
    tabi   <- split(tempdat, tempdat$id)
    L.tabi <- length(tabi)
    tablist <- lapply(1:L.tabi, function(i){
        lignes.tabi <- nrow(tabi[[i]])
        # for each [t1;t2]: \Lambda(t2) - \Lambda(t1)
        bashaz.denominator <- vector(length = lignes.tabi)
        # when time = t1
        taux.cum.start.denominator <- approx(x = tauxBasCum$time, y = tauxBasCum$hazard,
                                             xout = tabi[[i]][ , "tstart"])$y
        # when time = t2
        taux.cum.stop.denominator <- approx(x = tauxBasCum$time, y = tauxBasCum$hazard,
                                            xout = tabi[[i]][ , "timevar"])$y
        if(min(tabi[[i]]$tstart ) < 0){ 
            bashaz.denominator[1] <- tauxBasCum[1, "hazard"]
            if(lignes.tabi > 1){
                bashaz.denominator[2:lignes.tabi] <- (taux.cum.stop.denominator - taux.cum.start.denominator)[-1]
            }
        }else{
            bashaz.denominator[1] <- taux.cum.stop.denominator[1]
            if(lignes.tabi > 1){
                bashaz.denominator[2:lignes.tabi] <- (taux.cum.stop.denominator - taux.cum.start.denominator)[-1]
            }
        }
        indices0 <- which(tabi[[i]]$exposure == 0) # non-exposed patients
        indices1 <- which(tabi[[i]]$exposure == 1) # exposed patients

        tabi[[i]]$p.denominator <- vector("numeric", lignes.tabi)
        # when exposition = 0 => exp( - (\Lambda(t2) - \Lambda(t1)) * RR )
        tabi[[i]]$p.denominator[indices0] <- exp( - bashaz.denominator[indices0] * tabi[[i]]$risk.denominator[indices0])
        # when exposition = 1 => 1 - exp ( - (\Lambda(t2) - \Lambda(t1)) * RR )
        tabi[[i]]$p.denominator[indices1] <- 1 - exp( - bashaz.denominator[indices1] * tabi[[i]]$risk.denominator[indices1])
        # weights obtained by the cumulative product
        tabi[[i]]$w.denominator <- cumprod(tabi[[i]]$p.denominator)

        return( tabi[[i]] )
    } )
    rm(tauxBasCum)
    result.denominator <- do.call( rbind, tablist )
    #print("DONE")
    # Final weights
    if (is.null(tempcall$numerator)){
        tempdat$ipw.weights <- tempdat$w.numerator / result.denominator$w.denominator
    }else{
        tempdat$w.numerator <- result.numerator$w.numerator
        tempdat$w.denominator <- result.denominator$w.denominator
        tempdat$ipw.weights <- result.numerator$w.numerator / result.denominator$w.denominator
    }
    
    # Truncated weights (optional)
    if (!(is.null(tempcall$trunc))) {
        tempdat$weights.trunc <- tempdat$ipw.weights
        tempdat$weights.trunc[tempdat$ipw.weights <= quantile(tempdat$ipw.weights,
                                                              0 + trunc)] <- quantile(tempdat$ipw.weights, 0 +
                                                                                          trunc)
        tempdat$weights.trunc[tempdat$ipw.weights > quantile(tempdat$ipw.weights,
                                                             1 - trunc)] <- quantile(tempdat$ipw.weights, 1 -
                                                                                         trunc)
    }
    
    list(tempdat, summary(mod1), summary(mod2))
}