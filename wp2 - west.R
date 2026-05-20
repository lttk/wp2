library(tidyr)
library(nnet)
library(AICcmodavg)
require(VGAM)
library(caret)
library(tidyverse)
library(survminer)
library(adjustedCurves)
library(ggplot2)
library(dplyr)
library(tableone)
library(openxlsx)
library(boot)
## PS matching
library(Matching)
## Weighted analysis
library(survey)
## Reorganizing data
library(reshape2)
library(fitdistrplus)
library(cobalt)
library(WeightIt)
library(chisquare)
library(pROC)
library(ResourceSelection)
library(survival)
library(ipw)
library(stats)
options(dplyr.summarise.inform = FALSE)
library(gamlss)
source("/work/ttkle/wp2/SVall.R")

shift_left <- function(row) {
  vals <- row[!is.na(row)]
  c(vals, rep(NA, length(row) - length(vals)))
}
shift_right_mat <- function(mat) {
  t(apply(mat, 1, function(row) {
    vals <- row[!is.na(row)]
    c(rep(NA, length(row) - length(vals)), vals)
  }))
}
reverse_diagonal_matrix_blocks <- function(mat, z) {
  n <- nrow(mat)
  m <- ncol(mat)
  num_blocks <- m %/% z
  result <- matrix(NA, nrow = n, ncol = m)
  for (i in 1:n) {
    blocks <- list()
    for (j in 0:(min(i - 1, num_blocks - 1))) {
      start_col <- j * z + 1
      end_col <- start_col + z - 1
      block <- mat[i, start_col:end_col]
      blocks[[length(blocks) + 1]] <- block
    }
    reversed_blocks <- rev(blocks)
    row_vals <- unlist(reversed_blocks)
    len_row <- length(row_vals)
    result[i, 1:len_row] <- row_vals
  }
  return(result)
}
w.treat = function(ori,t.gap,var_X,var_lgnX,sub_var_lgnX,var_y,stab,func.treat) {
    ps.out = data.frame(id_patient = unique(ori$id_patient))
    ps.out.l = data.frame()
    for (t.est in seq(0,max(ori$t),t.gap)) {
        ps.l = data.frame(id_patient = unique(ori$id_patient))
        psnum.l = data.frame(id_patient = unique(ori$id_patient))
        ori_ = ori[ori$mask_cen==0,]
        df.train = ori_[ori_$t<=t.est&ori_$t>=(t.est-len.win),]
        idx = 0
        var.lgn.fml = c()
        y.lgn.fml = c()
        for (t in unique(df.train$t)) {
            if (idx == 0) {
                df_ = df.train[df.train$t==t,c(var_X,var_lgnX,c('id_patient'))]
                if (t<0) {
                    colnames(df_)[which(names(df_) %in% var_lgnX)] = paste(var_lgnX,paste0('neg',abs(t)),sep='_')
                    var.lgn.fml = append(var.lgn.fml,paste(sub_var_lgnX,paste0('neg',abs(t)),sep='_'))
                    if (t!=t.est) {y.lgn.fml = append(y.lgn.fml,paste(var_y,paste0('neg',abs(t)),sep='_'))}
                } else {colnames(df_)[which(names(df_) %in% var_lgnX)] = paste(var_lgnX,t,sep='_')
                        var.lgn.fml = append(var.lgn.fml,paste(sub_var_lgnX,t,sep='_'))
                        if (t!=t.est) {y.lgn.fml = append(y.lgn.fml,paste(var_y,t,sep='_'))}}
                df.t = df_
            } else {
                df_ = df.train[df.train$t==t,c(var_lgnX,c('id_patient'))]
                if (t<0) {
                    colnames(df_)[which(names(df_) %in% var_lgnX)] = paste(var_lgnX,paste0('neg',abs(t)),sep='_')
                    var.lgn.fml = append(var.lgn.fml,paste(sub_var_lgnX,paste0('neg',abs(t)),sep='_'))
                    if (t!=t.est) {y.lgn.fml = append(y.lgn.fml,paste(var_y,paste0('neg',abs(t)),sep='_'))}
                } else {colnames(df_)[which(names(df_) %in% var_lgnX)] = paste(var_lgnX,t,sep='_')
                        var.lgn.fml = append(var.lgn.fml,paste(sub_var_lgnX,t,sep='_'))
                        if (t!=t.est) {y.lgn.fml = append(y.lgn.fml,paste(var_y,t,sep='_'))}}
                df.t = merge(df.t,df_,by='id_patient')
            }
            idx = idx+1
        }
        df = data.frame(df.t)
        for (i in (1:length(var_y))) {
            if (func.treat == 'binomial') {
                if (is.null(y.lgn.fml)) {
                    mod1 = glm(as.formula(paste0(paste(var_y[i],t.est,sep='_'),' ~ ',
                                    paste(var_X[var_X != 'bio_pha_atc'],collapse=' + '),' + ',
                                    paste(var.lgn.fml,collapse=' + '))),
                            data = df, family = binomial)
                    mod1.num = glm(as.formula(paste0(paste(var_y[i],t.est,sep='_'),' ~ ',
                                    paste(var_X[var_X != 'bio_pha_atc'],collapse=' + '))),
                            data = df, family = binomial)
                } else {
                    mod1 = glm(as.formula(paste0(paste(var_y[i],t.est,sep='_'),' ~ ',
                                    paste(var_X,collapse=' + '),' + ',
                                    paste(var.lgn.fml,collapse=' + '),' + ',
                                    paste(y.lgn.fml[startsWith(y.lgn.fml,var_y[i])],collapse=' + '))),
                            data = df, family = binomial)
                    mod1.num = glm(as.formula(paste0(paste(var_y[i],t.est,sep='_'),' ~ ',
                                    paste(var_X,collapse=' + '),' + ',
                                    paste(y.lgn.fml[startsWith(y.lgn.fml,var_y[i])],collapse=' + '))),
                            data = df, family = binomial)
                }
                ps = predict(mod1,df,type="response")
                df[,paste(var_y[i],'ps',sep='_')] = ps
                ps.l = merge(ps.l,df[,c('id_patient',paste(var_y[i],'ps',sep='_'))],by='id_patient',all.x=TRUE)
                ps.l$t = t.est
                df[,paste(var_y[i],t.est,'ps',sep='_')] = ps
                ps.out = merge(ps.out,df[,c('id_patient',paste(var_y[i],t.est,'ps',sep='_'))],by='id_patient',all.x=TRUE)
                if (stab == TRUE) {
                    psnum = predict(mod1.num,df,type="response")
                    df[,paste(var_y[i],'psnum',sep='_')] = psnum
                    ps.l = merge(ps.l,df[,c('id_patient',paste(var_y[i],'psnum',sep='_'))],by='id_patient',all.x=TRUE)
                    df[,paste(var_y[i],t.est,'psnum',sep='_')] = psnum
                    ps.out = merge(ps.out,df[,c('id_patient',paste(var_y[i],t.est,'psnum',sep='_'))],by='id_patient',all.x=TRUE)
                } else {
                    df[,paste(var_y[i],'psnum',sep='_')] = 1
                    ps.l = merge(ps.l,df[,c('id_patient',paste(var_y[i],'psnum',sep='_'))],by='id_patient',all.x=TRUE)
                    df[,paste(var_y[i],t.est,'psnum',sep='_')] = 1
                    ps.out = merge(ps.out,df[,c('id_patient',paste(var_y[i],t.est,'psnum',sep='_'))],by='id_patient',all.x=TRUE)}} 
            else if (func.treat == 'multinomial') {
                df[,paste(var_y[i],t.est,sep='_')] = as.factor(df[,paste(var_y[i],t.est,sep='_')])
                if (is.null(y.lgn.fml)) {
                    mod1 = multinom(as.formula(paste0(paste(var_y[i],t.est,sep='_'),' ~ ',
                                    paste(var_X[var_X != 'bio_pha_atc'],collapse=' + '),' + ',
                                    paste(var.lgn.fml,collapse=' + '))),
                            data = df,maxit = 1000,trace=FALSE)
                    mod1.num = multinom(as.formula(paste0(paste(var_y[i],t.est,sep='_'),' ~ ',
                                    paste(var_X[var_X != 'bio_pha_atc'],collapse=' + '))),
                            data = df,maxit = 1000,trace=FALSE)
                } else {
                    mod1 = multinom(as.formula(paste0(paste(var_y[i],t.est,sep='_'),' ~ ',
                                    paste(var_X,collapse=' + '),' + ',
                                    paste(var.lgn.fml,collapse=' + '),' + ',
                                    paste(y.lgn.fml[startsWith(y.lgn.fml,var_y[i])],collapse=' + '))),
                            data = df,maxit = 1000,trace=FALSE)
                    mod1.num = multinom(as.formula(paste0(paste(var_y[i],t.est,sep='_'),' ~ ',
                                    paste(var_X,collapse=' + '),' + ',
                                    paste(y.lgn.fml[startsWith(y.lgn.fml,var_y[i])],collapse=' + '))),
                            data = df,maxit = 1000,trace=FALSE)
                }
                if (var_y[i]=="lgn_mtx_sld_2.5") {
                  s <- summary(mod1)
                  coefs <- t(s$coefficients)
                  colnames(coefs) = c("Estimate1", "Estimate2")
                  ses   <- t(s$standard.errors)
                  colnames(ses) = c("Std.err1", "Std.err2")
                  summary_table <- as.data.frame(cbind(coefs,ses))
                  print(summary_table)
                  na_vars <- summary_table %>%
                    filter(Estimate1 == 0| Estimate2 == 0 | abs(Estimate1) > 5 | abs(Estimate2) > 5) %>%
                    rownames()
                  na_vars <- ifelse(startsWith(na_vars, "bio_pha_atc"), "bio_pha_atc", na_vars)
                  na_vars = na_vars[na_vars != "(Intercept)"]
                  df_na = df[,c(na_vars,paste(var_y[i],t.est,sep='_'))]
                  if (t.est>=3) {
                    print(table(df[,paste(var_y[i],t.est,sep='_')], df[,paste(var_y[i],t.est-3,sep='_')]))
                  }
                  for (col in na_vars) {
                    print(col)
                    print(table(df_na[,paste(var_y[i],t.est,sep='_')], df_na[,col]))
                  }
                }
                ps = predict(mod1,df,type="probs")
                ref = as.numeric(df[,paste(var_y[i],t.est,sep='_')])
                vec = sapply(1:dim(ps)[1], function(i) ps[i,as.numeric(ref[i])])
                df[,paste(var_y[i],'ps',sep='_')] = vec
                ps.l = merge(ps.l,df[,c('id_patient',paste(var_y[i],'ps',sep='_'))],by='id_patient',all.x=TRUE)
                ps.l$t = t.est
                df[,paste(var_y[i],t.est,'ps',sep='_')] = vec
                ps.out = merge(ps.out,df[,c('id_patient',paste(var_y[i],t.est,'ps',sep='_'))],by='id_patient',all.x=TRUE)
                if (stab == TRUE) {
                    psnum = predict(mod1.num,df,type="probs")
                    vecnum = sapply(1:dim(psnum)[1], function(i) psnum[i,as.numeric(ref[i])])
                    df[,paste(var_y[i],'psnum',sep='_')] = vecnum
                    ps.l = merge(ps.l,df[,c('id_patient',paste(var_y[i],'psnum',sep='_'))],by='id_patient',all.x=TRUE)
                    df[,paste(var_y[i],t.est,'psnum',sep='_')] = vecnum
                    ps.out = merge(ps.out,df[,c('id_patient',paste(var_y[i],t.est,'psnum',sep='_'))],by='id_patient',all.x=TRUE)} 
                else {
                    df[,paste(var_y[i],'psnum',sep='_')] = 1
                    ps.l = merge(ps.l,df[,c('id_patient',paste(var_y[i],'psnum',sep='_'))],by='id_patient',all.x=TRUE)
                    df[,paste(var_y[i],t.est,'psnum',sep='_')] = 1
                    ps.out = merge(ps.out,df[,c('id_patient',paste(var_y[i],t.est,'psnum',sep='_'))],by='id_patient',all.x=TRUE)}
                }
            else if (func.treat == 'lm') {
                if (is.null(y.lgn.fml)) {
                    mod1 = gamlss(as.formula(paste0(paste(var_y[i],t.est,sep='_'),' ~ ',
                                    paste(var_X[var_X != 'bio_pha_atc'],collapse=' + '),' + ',
                                    paste(var.lgn.fml,collapse=' + '))),
                            family = BEINF,
                            data = df,
                            trace = FALSE)                    
                    mod1.num = gamlss(as.formula(paste0(paste(var_y[i],t.est,sep='_'),' ~ ',
                                    paste(var_X[var_X != 'bio_pha_atc'],collapse=' + '))),
                            family = BEINF,
                            data = df,
                            trace = FALSE)                    
                } else {
                    mod1 = gamlss(as.formula(paste0(paste(var_y[i],t.est,sep='_'),' ~ ',
                                    paste(var_X,collapse=' + '),' + ',
                                    paste(var.lgn.fml,collapse=' + '),' + ',
                                    paste(y.lgn.fml[startsWith(y.lgn.fml,var_y[i])],collapse=' + '))),
                            family = BEINF,
                            data = df,
                            trace = FALSE)                    
                    mod1.num = gamlss(as.formula(paste0(paste(var_y[i],t.est,sep='_'),' ~ ',
                                    paste(var_X,collapse=' + '),' + ',
                                    paste(y.lgn.fml[startsWith(y.lgn.fml,var_y[i])],collapse=' + '))),
                            family = BEINF,
                            data = df,
                            trace = FALSE)                    
                }
                

                psnum = dBEINF(df[,paste(var_y[i],t.est,sep='_')], mu = fitted(mod1.num, "mu"), 
                                sigma = fitted(mod1.num, "sigma"), 
                                nu = fitted(mod1.num, "nu"), 
                                tau = fitted(mod1.num, "tau"))
                df[,paste(var_y[i],'psnum',sep='_')] = psnum
                ps.l = merge(ps.l,df[,c('id_patient',paste(var_y[i],'psnum',sep='_'))],by='id_patient',all.x=TRUE)
                df[,paste(var_y[i],t.est,'psnum',sep='_')] = psnum
                ps.out = merge(ps.out,df[,c('id_patient',paste(var_y[i],t.est,'psnum',sep='_'))],by='id_patient',all.x=TRUE)

              
                ps = dBEINF(df[,paste(var_y[i],t.est,sep='_')], mu = fitted(mod1, "mu"), 
                            sigma = fitted(mod1, "sigma"), 
                            nu = fitted(mod1, "nu"), 
                            tau = fitted(mod1, "tau"))
                df[,paste(var_y[i],'ps',sep='_')] = ps
                ps.l = merge(ps.l,df[,c('id_patient',paste(var_y[i],'ps',sep='_'))],by='id_patient',all.x=TRUE)
                ps.l$t = t.est
                df[,paste(var_y[i],t.est,'ps',sep='_')] = ps
                ps.out = merge(ps.out,df[,c('id_patient',paste(var_y[i],t.est,'ps',sep='_'))],by='id_patient',all.x=TRUE)
            }
        }
        ps.out.l = rbind(ps.out.l,ps.l)
    }
    #### calculate weights
    if (FALSE) {
        print(var_y[i])
        s <- summary(mod1)
        coefs <- t(s$coefficients)
        colnames(coefs) = c("Estimate1", "Estimate2")
        ses   <- t(s$standard.errors)
        colnames(ses) = c("Std.err1", "Std.err2")
        summary_table <- as.data.frame(cbind(coefs,ses))
        print(summary_table)
        na_vars <- summary_table %>%
            filter(Estimate1 == 0| Estimate2 == 0 | abs(Estimate1) > 5 | abs(Estimate2) > 5) %>%
            rownames()
        na_vars = na_vars[na_vars != "(Intercept)"]
        df_na = df.fn[,c(na_vars,paste(var_y[i],t.est,sep='_'))]
        for (col in na_vars) {
            print(col)
            print(table(df_na[,paste(var_y[i],t.est,sep='_')], df_na[,col]))
        }
    }
    

    block.t = ori[,c('id_patient','t','comp_cen',var_y)]
    block.t = merge(block.t,ps.out.l,by=c('id_patient','t'))
    df = data.frame(block.t)
    df = df[order(df$id_patient,df$t),]
    row.names(df) = NULL
    for (tar in var_y) {
        if (func.treat == 'lm') {
            df[,paste(tar,'wt',sep='_')] = df[,paste(tar,'psnum',sep='_')]/df[,paste(tar,'ps',sep='_')]
        } else {
            df[,paste(tar,'wt',sep='_')] = ifelse(df[,tar]==0,(1-df[,paste(tar,'psnum',sep='_')])/(1-df[,paste(tar,'ps',sep='_')]),df[,paste(tar,'psnum',sep='_')]/df[,paste(tar,'ps',sep='_')])
        }
    }
    #### calculate cumulative weights
    df_ = data.frame(df)
    cumwt = df_ %>% group_by(id_patient) %>% summarise(lgn_mtx_ovr_2.5_cumwt = cumprod(lgn_mtx_ovr_2.5_wt),
                                                    lgn_mtx_sld_2.5_cumwt = cumprod(lgn_mtx_sld_2.5_wt),
                                                    lgn_mtx_ovr_2.0_cumwt = cumprod(lgn_mtx_ovr_2.0_wt),
                                                    lgn_mtx_sld_2.0_cumwt = cumprod(lgn_mtx_sld_2.0_wt),
                                                    lgn_mtx_ovr_1.4_cumwt = cumprod(lgn_mtx_ovr_1.4_wt),
                                                    lgn_mtx_sld_1.4_cumwt = cumprod(lgn_mtx_sld_1.4_wt))
    cumwt$t = df$t
    df = merge(df,cumwt,by=c('id_patient','t'))
    df = df[order(df$id_patient,df$t),]
    row.names(df) = NULL
    
return(df[,c('id_patient','t',paste(rep(var_y),'cumwt',sep='_'))])}

w.cen = function(ori,t.gap,var_X,var_lgnX,sub_var_lgnX,var_y,stab,func.cen,subgroup) {
    if (func.cen == 'logit') {
        pscen.out = data.frame(id_patient = unique(ori$id_patient))
        pscen.out.l = data.frame()
        for (t.est in seq(t.gap,(max(ori$t)-t.gap),t.gap)) {
            pscen.l = data.frame(id_patient = unique(ori$id_patient))
            ori_ = data.frame(ori)
            df.train = ori_[ori_$t<=t.est&ori_$t>=(t.est-len.win),]
            idx = 0
            cen.lgn.fml = c()
            y.lgn.fml = c()
            for (t in unique(df.train$t)) {
                if (idx == 0) {
                    df_ = df.train[df.train$t==t,c(var_X,var_lgnX,c('id_patient'))]
                    if (t<0) {
                        colnames(df_)[which(names(df_) %in% var_lgnX)] = paste(var_lgnX,paste0('neg',abs(t)),sep='_')
                        if (t!=t.est) {
                            cen.lgn.fml = append(cen.lgn.fml,paste(sub_var_lgnX,paste0('neg',abs(t)),sep='_'))
                            y.lgn.fml = append(y.lgn.fml,paste(var_y,paste0('neg',abs(t)),sep='_'))}
                    } else {colnames(df_)[which(names(df_) %in% var_lgnX)] = paste(var_lgnX,t,sep='_')
                            if (t!=t.est) {
                                cen.lgn.fml = append(cen.lgn.fml,paste(sub_var_lgnX,t,sep='_'))
                                y.lgn.fml = append(y.lgn.fml,paste(var_y,t,sep='_'))}}
                    df.t = df_
                } else {
                    df_ = df.train[df.train$t==t,c(var_lgnX,c('id_patient'))]
                    if (t<0) {
                        colnames(df_)[which(names(df_) %in% var_lgnX)] = paste(var_lgnX,paste0('neg',abs(t)),sep='_')
                        if (t!=t.est) {
                            cen.lgn.fml = append(cen.lgn.fml,paste(sub_var_lgnX,paste0('neg',abs(t)),sep='_'))
                            y.lgn.fml = append(y.lgn.fml,paste(var_y,paste0('neg',abs(t)),sep='_'))}
                    } else {colnames(df_)[which(names(df_) %in% var_lgnX)] = paste(var_lgnX,t,sep='_')
                            if (t!=t.est) {
                                cen.lgn.fml = append(cen.lgn.fml,paste(sub_var_lgnX,t,sep='_'))
                                y.lgn.fml = append(y.lgn.fml,paste(var_y,t,sep='_'))}}
                    df.t = merge(df.t,df_,by='id_patient')
                }
                idx = idx+1
            }
            df = df.t
            df.test = ori[ori$t==t.est,c('id_patient','comp_cen')]
            df = merge(df,df.test,by='id_patient')
            df = df%>%drop_na()
            df$comp_cen = abs(df$comp_cen-1)
            for (treat in var_y) {
                if (is.null(y.lgn.fml)) {
                    mod2 = glm(as.formula(paste0('comp_cen ~ ',
                            paste(var_X[var_X != 'bio_pha_atc'],collapse=' + '))),
                    data = df, family = binomial)
                    mod2.num = glm(as.formula(paste0('comp_cen ~ ',
                            paste(var_X[var_X != 'bio_pha_atc'],collapse=' + '))),
                    data = df, family = binomial)
                } else {
                    mod2 = glm(as.formula(paste0('comp_cen ~ ',
                            paste(var_X,collapse=' + '),' + ',
                            paste(cen.lgn.fml,collapse=' + '),' + ',
                            paste(y.lgn.fml[startsWith(y.lgn.fml,treat)],collapse=' + '))),
                    data = df, family = binomial)
                    mod2.num = glm(as.formula(paste0('comp_cen ~ ',
                                paste(var_X,collapse=' + '),' + ',
                                paste(y.lgn.fml[startsWith(y.lgn.fml,treat)],collapse=' + '))),
                        data = df, family = binomial)
                }
                ps.cen = predict(mod2,df,type="response")
                df[,paste(treat,'pscen',sep='_')] = ps.cen
                pscen.l = merge(pscen.l,df[,c('id_patient',paste(treat,'pscen',sep='_'))],by='id_patient',all.x=TRUE)
                pscen.l$t = t.est
                df[,paste(treat,t.est,'pscen',sep='_')] = ps.cen
                pscen.out = merge(pscen.out,df[,c('id_patient',paste(treat,t.est,'pscen',sep='_'))],by='id_patient',all.x=TRUE)
                if (stab == TRUE) {
                    psnumcen = predict(mod2.num,df,type="response")
                    df[,paste(treat,'psnumcen',sep='_')] = psnumcen
                    pscen.l = merge(pscen.l,df[,c('id_patient',paste(treat,'psnumcen',sep='_'))],by='id_patient',all.x=TRUE)
                    df[,paste(treat,t.est,'psnumcen',sep='_')] = psnumcen
                    pscen.out = merge(pscen.out,df[,c('id_patient',paste(treat,t.est,'psnumcen',sep='_'))],by='id_patient',all.x=TRUE)
                } else {
                    df[,paste(treat,'psnumcen',sep='_')] = 1
                    pscen.l = merge(pscen.l,df[,c('id_patient',paste(treat,'psnumcen',sep='_'))],by='id_patient',all.x=TRUE)
                    df[,paste(treat,t.est,'psnumcen',sep='_')] = 1
                    pscen.out = merge(pscen.out,df[,c('id_patient',paste(treat,t.est,'psnumcen',sep='_'))],by='id_patient',all.x=TRUE)
                }
            }    
            pscen.out.l = rbind(pscen.out.l,pscen.l)
        }
        #### calculate weights
        block.t = ori[,c('id_patient','t','comp_cen',var_y)]
        block.t = merge(block.t,pscen.out.l,by=c('id_patient','t'),all.x=TRUE)
        block.t[is.na(block.t)] = 1
        df = data.frame(block.t)
        df = df[order(df$id_patient,df$t),]
        row.names(df) = NULL
        for (tar in var_y) {
            df[,paste(tar,'wc',sep='_')] = ifelse(df[,'comp_cen']==0,df[,paste(tar,'psnumcen',sep='_')]/df[,paste(tar,'pscen',sep='_')],(1-df[,paste(tar,'psnumcen',sep='_')])/(1-df[,paste(tar,'pscen',sep='_')]))        
        }       
        #### calculate cumulative weights
        df_ = data.frame(df)
        cumwc = df_ %>% group_by(id_patient) %>% summarise(lgn_mtx_ovr_2.5_cumwc = cumprod(lgn_mtx_ovr_2.5_wc),
                                                        lgn_mtx_sld_2.5_cumwc = cumprod(lgn_mtx_sld_2.5_wc),
                                                        lgn_mtx_ovr_2.0_cumwc = cumprod(lgn_mtx_ovr_2.0_wc),
                                                        lgn_mtx_sld_2.0_cumwc = cumprod(lgn_mtx_sld_2.0_wc),
                                                        lgn_mtx_ovr_1.4_cumwc = cumprod(lgn_mtx_ovr_1.4_wc),
                                                        lgn_mtx_sld_1.4_cumwc = cumprod(lgn_mtx_sld_1.4_wc))
        cumwc$t = df$t
        df = merge(df,cumwc,by=c('id_patient','t'),all.x=TRUE)
        df = df[order(df$id_patient,df$t),]
        row.names(df) = NULL
    
    } else if (func.cen=='survival') {
        df.tplus = ori[,c('id_patient','t','time','comp_cen')]
        df.tplus$t = df.tplus$t - t.gap
        df.t = ori[,c('id_patient','t',var_X,var_lgnX)]
        df.t$tstart = df.t$t * 30
        df = merge(df.t,df.tplus,by=c('id_patient','t'))
        df$time[df$tstart==df$time] = df$time[df$tstart==df$time] + 1 
        df = df[order(df$id_patient,df$tstart),]
        row.names(df) = NULL
        for (treat in var_y) {
            df$treat = df[,treat]
            if (subgroup == "on") {
              mod = SVall(exposure = comp_cen,
                          numerator = ~ ben_sex_cod + age + bln_cci + bln_aspsa + bln_dld + bln_hpt + bln_ibd + bln_otinfla + bln_crtc + bln_nsaid + bln_mtx + bln_clpr + bln_acitretin + treat,
                          denominator = ~ ben_sex_cod + age + bln_cci + bln_aspsa + bln_dld + bln_hpt + bln_ibd + bln_otinfla + bln_crtc + bln_nsaid + bln_mtx + bln_clpr + bln_acitretin + lgn_cci + lgn_aspsa + lgn_dyslipidemia + lgn_hyper + lgn_ibd + lgn_otcroinfla + lgn_incrbio + lgn_sys + lgn_vitd + lgn_topsteroid + treat, id = id_patient,
                          tstart = tstart, timevar = time, data = df, trunc = NULL)
            } else {
              mod = SVall(exposure = comp_cen,
                          numerator = ~ ben_sex_cod + age + bln_cci + bln_aspsa + bln_dld + bln_hpt + bln_ibd + bln_otinfla + bln_crtc + bln_nsaid + bln_mtx + bln_clpr + bln_acitretin + treat,
                          denominator = ~ ben_sex_cod + age + bln_cci + bln_aspsa + bln_dld + bln_hpt + bln_ibd + bln_otinfla + bln_crtc + bln_nsaid + bln_mtx + bln_clpr + bln_acitretin + lgn_cci + lgn_aspsa + lgn_dyslipidemia + lgn_hyper + lgn_ibd + lgn_otcroinfla + lgn_incrbio + lgn_sys + lgn_vitd + lgn_topsteroid + treat, id = id_patient,
                          tstart = tstart, timevar = time, data = df, trunc = NULL)
            }
            
            # 
            df[,paste(treat,'cumwc',sep='_')] = mod$ipw.weights
        }
        df$t = df$t + t.gap
        df.bas = data.frame(id_patient = unique(df$id_patient),
                            t = rep(0,length(unique(df$id_patient))))
        df = bind_rows(df,df.bas)
        df[is.na(df)] = 1
        df = df[order(df$id_patient,df$t),]
        row.names(df) = NULL
    }
return(df[,c('id_patient','t',paste(rep(var_y),'cumwc',sep='_'))])}


w.est = function(ori,t.gap,var_X,var_lgnX,sub_var_lgnX,var_y,trunc,func.cen,func.treat,stab,subgroup) {
    wt = w.treat(ori=ori,t.gap=t.gap,var_X=var_X,var_lgnX=var_lgnX,sub_var_lgnX=sub_var_lgnX,var_y=var_y,stab=stab,func.treat=func.treat)
    wt = wt %>% drop_na()
    wc = w.cen(ori=ori,t.gap=t.gap,var_X=var_X,var_lgnX=var_lgnX,sub_var_lgnX=sub_var_lgnX,var_y=var_y,stab=stab,func.cen=func.cen,subgroup = subgroup)
    
    df = merge(wt,wc,by=c('id_patient','t'),all.x=TRUE)
    df[,paste(rep(var_y),'w',sep='_')] = df[,paste(rep(var_y),'cumwc',sep='_')] * df[,paste(rep(var_y),'cumwt',sep='_')]
    df.tplus = ori[ori$t>0,c('id_patient','t','comp_otc')]
    df.tplus$t = df.tplus$t - t.gap
    df.t = ori[,c('id_patient','t',var_y)]
    df.t = merge(df.t,df.tplus,by=c('id_patient','t'))
    df = merge(df.t,df,by = c('id_patient','t'))
    df = df[order(df$id_patient,df$t),]
    row.names(df) = NULL
    #### truncate
    for (col in var_y) {
        for (t in unique(df$t)) {
            w = df[df$t==t,paste(col,'w',sep='_')]
            q3 = quantile(w,(1-trunc))    
            q1 = quantile(w,trunc)
            w[w>q3] = q3
            w[w<q1] = q1
            df[df$t==t,paste(col,'w',sep='_')] = w
            #png(filename=paste0('mtx/graph/weights/w_',col,'_',t,'.png'))
            #boxplot(df[df$t==t,paste(col,'w',sep='_')] ~ df[df$t==t,col])
            #dev.off()
        }
        #png(filename=paste0('mtx/graph/weights/logw_',col,'.png'))
        #boxplot(log(df[,paste(col,'w',sep='_')]) ~ df$t)
        #dev.off()
    }
return(df)}

delta.cal = function(df, mode, func.treat, scale) {
    df_ = data.frame(df)
    if (func.treat=='multinomial') {
        df_[,c('lgn_mtx_ovr_2.5','lgn_mtx_sld_2.5','lgn_mtx_ovr_2.0','lgn_mtx_sld_2.0','lgn_mtx_ovr_1.4','lgn_mtx_sld_1.4')] = 
    ifelse(df_[,c('lgn_mtx_ovr_2.5','lgn_mtx_sld_2.5','lgn_mtx_ovr_2.0','lgn_mtx_sld_2.0','lgn_mtx_ovr_1.4','lgn_mtx_sld_1.4')]==2,1,0)
    } 
    #else if (func.treat=='lm') {
    #    df_[,c('lgn_mtx_ovr_2.5','lgn_mtx_sld_2.5','lgn_mtx_ovr_2.0','lgn_mtx_sld_2.0','lgn_mtx_ovr_1.4','lgn_mtx_sld_1.4')] =
    #ifelse(df_[,c('lgn_mtx_ovr_2.5','lgn_mtx_sld_2.5','lgn_mtx_ovr_2.0','lgn_mtx_sld_2.0','lgn_mtx_ovr_1.4','lgn_mtx_sld_1.4')]>c(75,85,85,100,100,100)/(100/#scale),1,0)
    #}
    
    delta = df_ %>% group_by(id_patient) %>% summarise(lgn_mtx_ovr_2.5_delta = cumsum(lgn_mtx_ovr_2.5),
                                                    lgn_mtx_sld_2.5_delta = cumsum(lgn_mtx_sld_2.5),
                                                    lgn_mtx_ovr_2.0_delta = cumsum(lgn_mtx_ovr_2.0),
                                                    lgn_mtx_sld_2.0_delta = cumsum(lgn_mtx_sld_2.0),
                                                    lgn_mtx_ovr_1.4_delta = cumsum(lgn_mtx_ovr_1.4),
                                                    lgn_mtx_sld_1.4_delta = cumsum(lgn_mtx_sld_1.4))

    delta$t = df$t
    delta$t = delta$t + t.gap
    bas.del = data.frame(id_patient = unique(df$id_patient),
                        t = rep(0,length(unique(df$id_patient))),
                        lgn_mtx_ovr_2.5_delta = rep(0,length(unique(df$id_patient))),
                        lgn_mtx_sld_2.5_delta = rep(0,length(unique(df$id_patient))),
                        lgn_mtx_ovr_2.0_delta = rep(0,length(unique(df$id_patient))),
                        lgn_mtx_sld_2.0_delta = rep(0,length(unique(df$id_patient))),
                        lgn_mtx_ovr_1.4_delta = rep(0,length(unique(df$id_patient))),
                        lgn_mtx_sld_1.4_delta = rep(0,length(unique(df$id_patient))))
    delta = rbind(bas.del,delta)
    df = merge(df,delta,by=c('id_patient','t'))
    df = df[order(df$id_patient,df$t),]
    row.names(df) = NULL
    if (mode =='double' & func.treat=='multinomial') {
        df_ = data.frame(df)
        df_[,c('lgn_mtx_ovr_2.5','lgn_mtx_sld_2.5','lgn_mtx_ovr_2.0','lgn_mtx_sld_2.0','lgn_mtx_ovr_1.4','lgn_mtx_sld_1.4')] = 
        ifelse(df_[,c('lgn_mtx_ovr_2.5','lgn_mtx_sld_2.5','lgn_mtx_ovr_2.0','lgn_mtx_sld_2.0','lgn_mtx_ovr_1.4','lgn_mtx_sld_1.4')]==1,1,0)
        
        delta = df_ %>% group_by(id_patient) %>% summarise(lgn_mtx_ovr_2.5_delta1 = cumsum(lgn_mtx_ovr_2.5),
                                                        lgn_mtx_sld_2.5_delta1 = cumsum(lgn_mtx_sld_2.5),
                                                        lgn_mtx_ovr_2.0_delta1 = cumsum(lgn_mtx_ovr_2.0),
                                                        lgn_mtx_sld_2.0_delta1 = cumsum(lgn_mtx_sld_2.0),
                                                        lgn_mtx_ovr_1.4_delta1 = cumsum(lgn_mtx_ovr_1.4),
                                                        lgn_mtx_sld_1.4_delta1 = cumsum(lgn_mtx_sld_1.4))

        delta$t = df$t
        delta$t = delta$t + t.gap
        bas.del = data.frame(id_patient = unique(df$id_patient),
                            t = rep(0,length(unique(df$id_patient))),
                            lgn_mtx_ovr_2.5_delta1 = rep(0,length(unique(df$id_patient))),
                            lgn_mtx_sld_2.5_delta1 = rep(0,length(unique(df$id_patient))),
                            lgn_mtx_ovr_2.0_delta1 = rep(0,length(unique(df$id_patient))),
                            lgn_mtx_sld_2.0_delta1 = rep(0,length(unique(df$id_patient))),
                            lgn_mtx_ovr_1.4_delta1 = rep(0,length(unique(df$id_patient))),
                            lgn_mtx_sld_1.4_delta1 = rep(0,length(unique(df$id_patient))))
        delta = rbind(bas.del,delta)
        df = merge(df,delta,by=c('id_patient','t'))
        df = df[order(df$id_patient,df$t),]
        row.names(df) = NULL
    }

    bidelta = ifelse(df[,paste(rep(var_y),'delta',sep='_')]>0,1,0)
    colnames(bidelta) = paste(rep(var_y),'bidelta',sep='_')
    df = cbind(df,bidelta)

    if (mode =='double' & func.treat=='multinomial') {
        catdelta = ifelse(df[,paste(rep(var_y),'delta',sep='_')]==df$t/t.gap & df$t>0,2,0)
        catdelta1 = ifelse(df[,paste(rep(var_y),'delta1',sep='_')]>0,1,0)
        res.catdelta = catdelta + catdelta1
        colnames(res.catdelta) = paste(rep(var_y),'catdelta',sep='_')
        df = cbind(df,res.catdelta)}
    return(df)}


