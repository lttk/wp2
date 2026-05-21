library("readxl")
library("nleqslv")
library("MASS")
library("tidyr")
library(nnet)
library(AICcmodavg)
library(VGAM)
library(caret)
library(tidyverse)
library(survminer)
library(adjustedCurves)
library(ggplot2)
library(dplyr)
library(tableone)
library(openxlsx)
library(boot)
library(Matching)
library(survey)
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
library(data.table)
library(tidyr)
library(brant)
library(broom)
source("/work/ttkle/wp2/wp2 - west.R")

df = read.csv("/work/ttkle/wp2/mtx/fin_lgn_tnf_upmed_csDMARDs_AG3.csv")
df1 = read.csv("/work/ttkle/wp2/mtx/fin_lgn_tnf_upmed_syscorti_AG3.csv")
pori = read.csv("/work/ttkle/wp2/mtx/fin_lgn_semicat_comp2_tnf_conti_AG3_sum.csv")
pori = merge(pori,df[,c('id_patient','t','lgn_aza','lgn_mmf','lgn_mer','lgn_hdcq','lgn_lfnm','lgn_sfsl')],by=c('id_patient','t'),all.x = TRUE)
pori = merge(pori,df1[,c('id_patient','t','lgn_syscorti')],by=c('id_patient','t'),all.x = TRUE)
pori = pori%>%drop_na()
pori$t = as.numeric(pori$t)
pori$id_patient = as.numeric(as.factor((pori$id_patient)))
# refine outcomes and censporing
for (i in 1:2) {
  pori[,paste0('comp_cen_',i)] = pori$comp_cen
  pori[,paste0('comp_otc_',i)] = pori$comp_otc
  pori[,paste0('typ_cen_',i)] = pori$typ_cen
  pori[,paste0('typ_otc_',i)] = pori$typ_otc
  if (i==1) {
    pori$comp_cen_1[pori$typ_cen_1<3] = 0
    pori$comp_cen_1[pori$typ_cen_1==5&pori$time==360] = 0
    #pori$mask_cen_1[pori$typ_cen_1==5&pori$time==360] = 0
    pori[pori$typ_cen_1>pori$typ_otc_1 & pori$typ_otc_1>0, c('comp_otc_1','typ_otc_1')] = 0
  } else {
    pori$comp_otc_2[(pori$typ_cen_2==4) & (pori$comp_cen_2==1)] = 1
    pori$typ_otc_2[(pori$typ_cen_2==4) & (pori$comp_cen_2==1)] = 4
    pori[(pori$typ_cen_2==4) & (pori$comp_cen_2==1), c('comp_cen_2','typ_cen_2')] = 0
    
    pori$comp_cen_2[pori$typ_cen_2<3] = 0
    pori$comp_cen_2[pori$typ_cen_2==5&pori$time==360] = 0
    #pori$mask_cen_2[pori$typ_cen_2==5&pori$time==360] = 0
    pori[pori$typ_cen_2>pori$typ_otc_2 & pori$typ_otc_2>0, c('comp_otc_2','typ_otc_2')] = 0
  }
}
pori$lgn_otim = ifelse(pori$lgn_aza + pori$lgn_mmf + pori$lgn_mer + pori$lgn_hdcq + pori$lgn_lfnm + pori$lgn_sfsl>0,1,0)
pori[ , c('lgn_sys','lgn_syscorti','lgn_vitd','lgn_topsteroid')] = lapply(pori[ , c('lgn_sys','lgn_syscorti','lgn_vitd','lgn_topsteroid')], function(x) ifelse(x > 0, 1, x))
pori$lgn_ibdpsa  = ifelse(pori$lgn_ibd + pori$lgn_aspsa + pori$lgn_otim>0,1,0)
# define population without IBD and other csDMARDs
pori = subset(pori,id_patient %in% unique(pori$id_patient[pori$lgn_ibd==0&pori$lgn_otim==0&pori$t==0])) 
#### define subgroup: without infliximab
for (sub in c('psa','pso','wow_ifx','ada')) {
if (sub == 'wow_ifx') {
  ori = pori[pori$bio_pha_atc=="L04AB04"|pori$bio_pha_atc=="L04AB01"|pori$bio_pha_atc=="L04AB05",]
} else if (sub == 'psa') {
  ori = subset(pori,id_patient %in% unique(pori$id_patient[lgn_aspsa==1&pori$t==0]))
} else if (sub == 'pso') {
  ori = subset(pori,id_patient %in% unique(pori$id_patient[lgn_aspsa==0&pori$t==0]))
} else if (sub == 'ada') {
  ori = pori[pori$bio_pha_atc=="L04AB04",]
}
# define variables
var_X = c('ben_sex_cod','age','bln_nsaid','bln_clpr','bln_acitretin') 
var_lgnX = c('lgn_cci','lgn_dyslipidemia','lgn_hyper','lgn_ibdpsa','lgn_incrbio','lgn_sys','lgn_syscorti','lgn_vitd','lgn_topsteroid','lgn_mtx_ovr_2.5','lgn_mtx_sld_2.5','lgn_mtx_ovr_2.0','lgn_mtx_sld_2.0','lgn_mtx_ovr_1.4','lgn_mtx_sld_1.4')
var_y = c('lgn_mtx_ovr_2.5','lgn_mtx_sld_2.5','lgn_mtx_ovr_2.0','lgn_mtx_sld_2.0','lgn_mtx_ovr_1.4','lgn_mtx_sld_1.4')
sub_var_lgnX = setdiff(var_lgnX, var_y)

for (col in var_y) {
  threshold = median(ori[ori[,col] > 0 & ori$t==0,col])
  ori[ori[,col]>=threshold,col] = 20
  ori[ori[,col]<threshold & ori[,col]>0,col] = 10
  ori[ori[,col]==0,col] = 0
}
ori[,var_y] = ori[,var_y]/10
for (t in unique(ori$t)) {
  print(t)
  for (col in var_y) {
    print(col)
    print(as.data.frame(t(table(ori[ori$t==t & ori$mask_cen==0,col]))))
  }
}
t.gap = 3
idx=0
if (FALSE) {
for (t in unique(ori$t)) {
  print(t)
  for (col in var_y) {
    print(col)
    print(as.data.frame(t(table(ori[ori$t==t & ori$mask_cen==0,col]))))
  }
}}

for (t in unique(ori$t)) {
  df = ori[ori$t==t,c(c('id_patient','t'),var_lgnX)]
  colnames(df)[which(names(df) %in% var_lgnX)] = paste(var_lgnX,t,sep='_')
  if (idx == 0) {
    df.t = data.frame(df)
  } else {
    df.t = bind_rows(df.t,df)
  }
  idx = idx+1
}
df.treat.time = df.t[order(df.t$id_patient,df.t$t),]
row.names(df.treat.time) = NULL
setDT(df.treat.time)
all_cols <- apply(expand.grid(var_lgnX, unique(df.treat.time$t)), 1, function(x) paste(x[1],as.numeric(x[2]), sep = "_"))
missing_cols <- setdiff(all_cols, names(df.treat.time))
if(length(missing_cols) > 0) df.treat.time[, (missing_cols) := NA]
for (col in all_cols) {
  setorder(df.treat.time, id_patient, t)
  df.treat.time[, (col) := nafill(get(col), type = "locf"), by = id_patient]
}
value_cols <- names(df.treat.time)[3:ncol(df.treat.time)]
cols <- apply(expand.grid(var_lgnX, c(paste0("neg",c(0,1,2,3,4)))), 1, function(x) paste(x[1], x[2], sep = "_"))
result_list <- df.treat.time[, {
  mat <- as.matrix(.SD)
  shifted <- reverse_diagonal_matrix_blocks(mat, z = length(var_lgnX))
  shifted_df <- as.data.frame(shifted)
  shifted_df <- setNames(shifted_df, cols)
  shifted_df$t <- t
  shifted_df
}, by = id_patient, .SDcols = value_cols]
df.fn <- as.data.frame(result_list)
df.fn[is.na(df.fn)] = 0
# start-stop
start.stop = ori[,c('id_patient','t','time','comp_cen_1','comp_otc_1','comp_cen_2','comp_otc_2',var_X)]
start.stop$t = start.stop$t - t.gap
start.stop$tstart = start.stop$t * 30
colnames(start.stop)[3] = c("tstop")
start.stop = start.stop[,order(colnames(start.stop))]
df.fn = merge(start.stop,df.fn,by=c('id_patient','t'))
df.fn[df.fn$tstart==df.fn$tstop,'tstop'] = df.fn[df.fn$tstart==df.fn$tstop,'tstop'] + 1
df.fn = df.fn[order(df.fn$id_patient,df.fn$tstart),]
row.names(df.fn) = NULL
df = data.frame(df.fn)

for (i in (1:length(var_y))) {
    var.lgn.fml = apply(expand.grid(sub_var_lgnX, c('neg0','neg1','neg2','neg3')), 1, function(x) paste(x[1], x[2], sep = "_"))
    y.lgn.fml = c(paste(var_y[i],c('neg1','neg2','neg3'),sep='_'))
    y.lgn = c(paste(var_y[i],c('neg0','neg1','neg2','neg3'),sep='_'))
    var.lgn.fml.cen = apply(expand.grid(sub_var_lgnX, c('neg1','neg2','neg3')), 1, function(x) paste(x[1], x[2], sep = "_"))
    # delta calculation
    df[,paste(var_y[i],'delta2',sep='_')] = rowSums(df[,y.lgn.fml] == 2)
    df[,paste(var_y[i],'delta1',sep='_')] = rowSums(df[,y.lgn.fml] == 1)
    df[,paste(var_y[i],'deltacat',sep='_')] = apply(df[,c('t',y.lgn.fml)], 1, function(x) {
      last = x[1]/3
      if (last>0) {
        if (all(x[2:(last+1)] == 2)) {2} 
        else if (all(x[2:(last+1)] == 0)) {0} 
        else {1}}
      else {0}}
    )
    
    for (col in y.lgn) {
      df[,col] = factor(df[,col], levels = c("0","1","2"))
    }
    mod1 = multinom(as.formula(paste0(paste(var_y[i],'neg0',sep='_'),' ~ ',
                                      paste(var_X,collapse=' + '),' + ',
                                      paste(var.lgn.fml,collapse=' + '),' + ',
                                      paste(y.lgn.fml,collapse=' + '),' + ',
                                      'tstart')),
                              data = df,maxit = 1000,trace=FALSE)
    mod1.num = multinom(as.formula(paste0(paste(var_y[i],'neg0',sep='_'),' ~ ',
                                  paste(var_X,collapse=' + '),' + ',
                                  paste(y.lgn.fml,collapse=' + '),' + ',
                                  'tstart')),
                          data = df,maxit = 1000,trace=FALSE)
    # treatment weights
    ref = as.numeric(df[,paste(var_y[i],'neg0',sep='_')])
    ps.denominator = predict(mod1,df,type="probs")
    vec.denominator = sapply(1:dim(ps.denominator)[1], function(i) ps.denominator[i,as.numeric(ref[i])])
    ps.numerator = predict(mod1.num,df,type="probs")
    vec.numerator = sapply(1:dim(ps.numerator)[1], function(i) ps.numerator[i,as.numeric(ref[i])])
    df[,paste(var_y[i],'wt',sep='_')] = vec.numerator / vec.denominator
    df[,paste0(var_y[i], "_cumwt")] = ave(df[,paste0(var_y[i], "_wt")], df$id_patient, FUN = cumprod)
    for (j in 1:2) {
      mod = SVall(exposure = paste0('comp_cen_',j),
                numerator = paste0('~ ',paste(var_X,collapse = " + "),' + ',paste(y.lgn.fml,collapse = " + ")), 
                denominator = paste0('~ ',paste(var_X,collapse = " + "),' + ',paste(y.lgn.fml,collapse = " + "), ' + ',paste(var.lgn.fml.cen,collapse = ' + ')),
                id = id_patient, tstart = tstart, timevar = tstop, data = df, trunc = NULL)
      # censoring weights
      df[,paste(var_y[i],'cumwc',j,sep='_')] = mod[[1]]$ipw.weights
      # final weights
      df[,paste(var_y[i],'w',j,sep='_')] = df[,paste(var_y[i],'cumwt',sep='_')] * df[,paste(var_y[i],'cumwc',j,sep='_')]
    }  
}
#### truncate
trunc = 0.01
for (col in var_y) { 
  for (j in 1:2) {
    w = df[,paste(col,'w',j,sep='_')]
    q3 = quantile(w,(1-trunc))
    q1 = quantile(w,trunc)
    w[w>q3] = q3
    w[w<q1] = q1
    df[,paste(col,'w',j,sep='_')] = w
  }
}

hr = list()
ard = list()
for (treat in var_y) {
    df[,paste(treat,'deltacat',sep='_')] = factor(df[,paste(treat,'deltacat',sep='_')], levels = c("0","1","2"))
    for (col in c(paste(treat,c(paste0("neg",c(0,1,2,3))),sep='_'))) {
        df[,col] = factor(df[,col], levels = c("0","1","2"))
    }
    for (j in 1:2) {    
      mod.hr.main = coxph(as.formula(paste0('Surv(tstart, tstop, comp_otc_',j,') ~ ',paste(treat,"neg0",sep='_'),
                                                                              ' + ',paste(treat,"delta2",sep='_'),
                                                                              ' + ',paste(treat,"delta1",sep='_'),
                                                                              ' + cluster(id_patient)')),
                                              data = df, weights = df[,paste(treat,'w',j,sep='_')])
      mod.hr.cat = coxph(as.formula(paste0('Surv(tstart, tstop, comp_otc_',j,') ~ ',paste(treat,"neg0",sep='_'),
                                                                              ' + ',paste(treat,"deltacat",sep='_'),
                                                                              ' + cluster(id_patient)')),
                          data = df, weights = df[,paste(treat,'w',j,sep='_')])
      hr.main = c(summary(mod.hr.main)$coef[,"exp(coef)"])
      hr.cat = c(summary(mod.hr.cat)$coef[,"exp(coef)"])
      hr[[paste(treat,j,sep='_')]] = do.call(rbind,list(hr[[paste(treat,j,sep='_')]],c(hr.main,hr.cat)))
      # ard
      v_neg0   <- paste(treat, "neg0", sep='_')
      v_delta1 <- paste(treat, "delta1", sep='_')
      v_delta2 <- paste(treat, "delta2", sep='_')
    
      df.full <- expand.grid(id_patient = unique(df$id_patient), interval = 1:4) %>%
                    left_join(df %>% group_by(id_patient) %>% slice(1) %>% select(-c(tstart, tstop, all_of(v_neg0), all_of(v_delta1), all_of(v_delta2))), by = "id_patient") %>%
                    mutate(tstart = (interval - 1) * 90,
                          tstop = interval * 90,
                          !!v_neg0 := factor("2", levels=c("0","1","2")),
                          !!v_delta1 := 0,
                          !!v_delta2 := interval - 1)
      df.full$cum_haz <- predict(mod.hr.main, newdata = df.full, type = "expected")
      s2.365 <- df.full %>%
        group_by(id_patient) %>%
        summarise(total_H = sum(cum_haz)) %>%
        summarise(mean_S = mean(exp(-total_H))) %>%
        pull(mean_S)
      
      df.partial <- expand.grid(id_patient = unique(df$id_patient), interval = 1:4) %>%
        left_join(df %>% group_by(id_patient) %>% slice(1) %>% select(-c(tstart, tstop, all_of(v_neg0), all_of(v_delta1), all_of(v_delta2))), by = "id_patient") %>%
        mutate(tstart = (interval - 1) * 90,
               tstop = interval * 90,
               !!v_neg0 := factor("1", levels=c("0","1","2")),
               !!v_delta1 := interval - 1,
               !!v_delta2 := 0)
      df.partial$cum_haz <- predict(mod.hr.main, newdata = df.partial, type = "expected")
      s1.365 <- df.partial %>%
        group_by(id_patient) %>%
        summarise(total_H = sum(cum_haz)) %>%
        summarise(mean_S = mean(exp(-total_H))) %>%
        pull(mean_S)
      
      df.non <- expand.grid(id_patient = unique(df$id_patient), interval = 1:4) %>%
                  left_join(df %>% group_by(id_patient) %>% slice(1) %>% select(-c(tstart, tstop, all_of(v_neg0), all_of(v_delta1), all_of(v_delta2))), by = "id_patient") %>%
                  mutate(tstart = (interval - 1) * 90,
                        tstop = interval * 90,
                        !!v_neg0 := factor("0", levels=c("0","1","2")),
                        !!v_delta1 := 0,
                        !!v_delta2 := 0)
      df.non$cum_haz <- predict(mod.hr.main, newdata = df.non, type = "expected")
      s0.365 <- df.non %>%
        group_by(id_patient) %>%
        summarise(total_H = sum(cum_haz)) %>%
        summarise(mean_S = mean(exp(-total_H))) %>%
        pull(mean_S)
      
      ard.2.365 <- (1-s2.365) - (1-s0.365)
      ard.1.365 <- (1-s1.365) - (1-s0.365)
      ard[[paste(treat,j,sep='_')]] = do.call(rbind,list(ard[[paste(treat,j,sep='_')]],c(ard.2.365,ard.1.365)))
    }
}

n_boot = 999
n_phase = 10
for (phase in 1:n_phase) {
  for (j in split(c(1:n_boot), ceiling(seq_along(c(1:n_boot))/ceiling(n_boot/n_phase)))[[phase]]) {
    set.seed(j)
    cat(j," \r")
    flush.console()
    x = length(unique(df.fn$id_patient))
    bootid = unique(df.fn$id_patient)[sample(x,x,replace=TRUE)]
    num = sapply(unique(bootid),function(i) sum(bootid==i))
    idx = 0
    df = data.frame()
    while (any(num>=1)) {
      b_ = subset(df.fn,id_patient %in% unique(bootid)[num>=1])
      b_$id_patient = ifelse(idx==0,1,-1)*(b_$id_patient*10**ceiling(idx/10)+idx)
      df = rbind(df,b_)
      num = num - 1
      idx = idx + 1
    }
    for (i in (1:length(var_y))) {
      var.lgn.fml = apply(expand.grid(sub_var_lgnX, c('neg0','neg1','neg2','neg3')), 1, function(x) paste(x[1], x[2], sep = "_"))
      y.lgn.fml = c(paste(var_y[i],c('neg1','neg2','neg3'),sep='_'))
      y.lgn = c(paste(var_y[i],c('neg0','neg1','neg2','neg3'),sep='_'))
      var.lgn.fml.cen = apply(expand.grid(sub_var_lgnX, c('neg1','neg2','neg3')), 1, function(x) paste(x[1], x[2], sep = "_"))
      # delta calculation
      df[,paste(var_y[i],'delta2',sep='_')] = rowSums(df[,y.lgn.fml] == 2)
      df[,paste(var_y[i],'delta1',sep='_')] = rowSums(df[,y.lgn.fml] == 1)
      df[,paste(var_y[i],'deltacat',sep='_')] = apply(df[,c('t',y.lgn.fml)], 1, function(x) {
        last = x[1]/3
        if (last>0) {
          if (all(x[2:(last+1)] == 2)) {2} 
          else if (all(x[2:(last+1)] == 0)) {0} 
          else {1}}
        else {0}}
      )
      for (col in y.lgn) {
        df[,col] = factor(df[,col], levels = c("0","1","2"))
      }
      mod1 = multinom(as.formula(paste0(paste(var_y[i],'neg0',sep='_'),' ~ ',
                                        paste(var_X,collapse=' + '),' + ',
                                        paste(var.lgn.fml,collapse=' + '),' + ',
                                        paste(y.lgn.fml,collapse=' + '),' + ',
                                        'tstart')),
                                data = df,maxit = 1000,trace=FALSE)
      mod1.num = multinom(as.formula(paste0(paste(var_y[i],'neg0',sep='_'),' ~ ',
                                    paste(var_X,collapse=' + '),' + ',
                                    paste(y.lgn.fml,collapse=' + '),' + ',
                                    'tstart')),
                            data = df,maxit = 1000,trace=FALSE)
      # treatment weights
      ref = as.numeric(df[,paste(var_y[i],'neg0',sep='_')])
      ps.denominator = predict(mod1,df,type="probs")
      vec.denominator = sapply(1:dim(ps.denominator)[1], function(i) ps.denominator[i,as.numeric(ref[i])])
      ps.numerator = predict(mod1.num,df,type="probs")
      vec.numerator = sapply(1:dim(ps.numerator)[1], function(i) ps.numerator[i,as.numeric(ref[i])])
      df[,paste(var_y[i],'wt',sep='_')] = vec.numerator / vec.denominator
      df[,paste0(var_y[i], "_cumwt")] = ave(df[,paste0(var_y[i], "_wt")], df$id_patient, FUN = cumprod)
      for (k in 1:2) {
        mod = SVall(exposure = paste0('comp_cen_',k),
                    numerator = paste0('~ ',paste(var_X,collapse = " + "),' + ',paste(y.lgn.fml,collapse = " + ")), 
                    denominator = paste0('~ ',paste(var_X,collapse = " + "),' + ',paste(y.lgn.fml,collapse = " + "), ' + ',paste(var.lgn.fml.cen,collapse = ' + ')),
                    id = id_patient, tstart = tstart, timevar = tstop, data = df, trunc = NULL)
        # censoring weights
        df[,paste(var_y[i],'cumwc',k,sep='_')] = mod[[1]]$ipw.weights
        # final weights
        df[,paste(var_y[i],'w',k,sep='_')] = df[,paste(var_y[i],'cumwt',sep='_')] * df[,paste(var_y[i],'cumwc',k,sep='_')]
      }
      
      
      }
    #### truncate
    for (col in var_y) {
      for (k in 1:2) {
        w = df[,paste(col,'w',k,sep='_')]
        q3 = quantile(w,(1-trunc))
        q1 = quantile(w,trunc)
        w[w>q3] = q3
        w[w<q1] = q1
        df[,paste(col,'w',k,sep='_')] = w
      }
    }
      
    for (treat in var_y) {
      df[,paste(treat,'deltacat',sep='_')] = factor(df[,paste(treat,'deltacat',sep='_')], levels = c("0","1","2"))
      for (col in c(paste(treat,c(paste0("neg",c(0,1,2,3))),sep='_'))) {
          df[,col] = factor(df[,col], levels = c("0","1","2"))
      }
      for (k in 1:2) {
        mod.hr.main = coxph(as.formula(paste0('Surv(tstart, tstop, comp_otc_',k,') ~ ',paste(treat,"neg0",sep='_'),
                                                                                ' + ',paste(treat,"delta2",sep='_'),
                                                                                ' + ',paste(treat,"delta1",sep='_'),
                                                                                ' + cluster(id_patient)')),
                                                data = df, weights = df[,paste(treat,'w',k,sep='_')])
        mod.hr.cat = coxph(as.formula(paste0('Surv(tstart, tstop, comp_otc_',k,') ~ ',paste(treat,"neg0",sep='_'),
                                                                                ' + ',paste(treat,"deltacat",sep='_'),
                                                                                ' + cluster(id_patient)')),
                            data = df, weights = df[,paste(treat,'w',k,sep='_')])
        
        hr.main = c(summary(mod.hr.main)$coef[,"exp(coef)"])
        hr.cat = c(summary(mod.hr.cat)$coef[,"exp(coef)"])
        hr[[paste(treat,k,sep='_')]] = do.call(rbind,list(hr[[paste(treat,k,sep='_')]],c(hr.main,hr.cat)))
        # ard
        v_neg0   <- paste(treat, "neg0", sep='_')
        v_delta1 <- paste(treat, "delta1", sep='_')
        v_delta2 <- paste(treat, "delta2", sep='_')
        
        df.full <- expand.grid(id_patient = unique(df$id_patient), interval = 1:4) %>%
                      left_join(df %>% group_by(id_patient) %>% slice(1) %>% select(-c(tstart, tstop, all_of(v_neg0), all_of(v_delta1), all_of(v_delta2))), by = "id_patient") %>%
                      mutate(tstart = (interval - 1) * 90,
                            tstop = interval * 90,
                            !!v_neg0 := factor("2", levels=c("0","1","2")),
                            !!v_delta1 := 0,
                            !!v_delta2 := interval - 1)
        df.full$cum_haz <- predict(mod.hr.main, newdata = df.full, type = "expected")
        s2.365 <- df.full %>%
          group_by(id_patient) %>%
          summarise(total_H = sum(cum_haz)) %>%
          summarise(mean_S = mean(exp(-total_H))) %>%
          pull(mean_S)
        
        df.partial <- expand.grid(id_patient = unique(df$id_patient), interval = 1:4) %>%
          left_join(df %>% group_by(id_patient) %>% slice(1) %>% select(-c(tstart, tstop, all_of(v_neg0), all_of(v_delta1), all_of(v_delta2))), by = "id_patient") %>%
          mutate(tstart = (interval - 1) * 90,
                tstop = interval * 90,
                !!v_neg0 := factor("1", levels=c("0","1","2")),
                !!v_delta1 := interval - 1,
                !!v_delta2 := 0)
        df.partial$cum_haz <- predict(mod.hr.main, newdata = df.partial, type = "expected")
        s1.365 <- df.partial %>%
          group_by(id_patient) %>%
          summarise(total_H = sum(cum_haz)) %>%
          summarise(mean_S = mean(exp(-total_H))) %>%
          pull(mean_S)
        
        df.non <- expand.grid(id_patient = unique(df$id_patient), interval = 1:4) %>%
                    left_join(df %>% group_by(id_patient) %>% slice(1) %>% select(-c(tstart, tstop, all_of(v_neg0), all_of(v_delta1), all_of(v_delta2))), by = "id_patient") %>%
                    mutate(tstart = (interval - 1) * 90,
                          tstop = interval * 90,
                          !!v_neg0 := factor("0", levels=c("0","1","2")),
                          !!v_delta1 := 0,
                          !!v_delta2 := 0)
        df.non$cum_haz <- predict(mod.hr.main, newdata = df.non, type = "expected")
        s0.365 <- df.non %>%
          group_by(id_patient) %>%
          summarise(total_H = sum(cum_haz)) %>%
          summarise(mean_S = mean(exp(-total_H))) %>%
          pull(mean_S)
        
        ard.2.365 <- (1-s2.365) - (1-s0.365)
        ard.1.365 <- (1-s1.365) - (1-s0.365)
        ard[[paste(treat,j,sep='_')]] = do.call(rbind,list(ard[[paste(treat,j,sep='_')]],c(ard.2.365,ard.1.365)))
        }
      }
  }
  exx = apply(expand.grid(var_y, c(1:2)), 1, function(x) paste(x[1], x[2], sep = "_"))
  ci.hr =  sapply(exx, function(i) apply(hr[[i]], 2, quantile, probs = c(0.025,0.50,0.975), na.rm = TRUE ))
  val.hr = sapply(exx, function(i) hr[[i]][1,])
  p.value = sapply(exx, function(i) 2 * exp(pnorm(abs(log(hr[[i]][1,]) / apply(log(hr[[i]]),2,sd)), lower.tail = FALSE, log.p = TRUE) ) )
  
  ci.ard =  sapply(exx, function(i) apply(ard[[i]], 2, quantile, probs = c(0.025,0.50,0.975), na.rm = TRUE ))
  val.ard = sapply(exx, function(i) ard[[i]][1,])
  p.value.ard = sapply(exx, function(i) 2 * exp(pnorm(abs(log(ard[[i]][1,]) / apply(log(ard[[i]]),2,sd)), lower.tail = FALSE, log.p = TRUE) ) )
  
  hr_output_ = list(as.data.frame(val.hr),as.data.frame(ci.hr),as.data.frame(p.value))
  ard_output_ = list(as.data.frame(val.ard),as.data.frame(ci.ard),as.data.frame(p.value.ard))
  wb <- createWorkbook()
  for (r in names(hr)) {
    addWorksheet(wb, paste0(r)) 
    writeData(wb, paste0(r), as.data.frame(hr[[r]]))
  }
  for (r in names(ard)) {
    addWorksheet(wb, paste0(r,"_ard")) 
    writeData(wb, paste0(r,"_ard"), as.data.frame(ard[[r]]))
  }
  saveWorkbook(wb, paste0("/work/ttkle/wp2/raw_",sub,"_out_.xlsx"), overwrite = TRUE)
  wbr <- createWorkbook()
  addWorksheet(wbr, paste("hr_results_"))
  # start row counter
  row_start <- 1
  for(o in seq_along(hr_output_)) {
    if (is.null(hr_output_[[o]])) next
    writeData(wbr, paste("hr_results_"), hr_output_[[o]], startRow = row_start)
    row_start <- row_start + nrow(hr_output_[[o]]) + 2  # leave blank row between
  }
  addWorksheet(wbr, paste("ard_results_"))
  # start row counter
  row_start <- 1
  for(o in seq_along(ard_output_)) {
    if (is.null(ard_output_[[o]])) next
    writeData(wbr, paste("ard_results_"), ard_output_[[o]], startRow = row_start)
    row_start <- row_start + nrow(ard_output_[[o]]) + 2  # leave blank row between
  }
  saveWorkbook(wbr, paste0("/work/ttkle/wp2/",sub,"out_.xlsx"), overwrite = TRUE)
}
}
