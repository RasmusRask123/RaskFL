#######################################
### Simulation 6 mean value changes ###
#######################################
#1000 simulations
#20,20,30,...,220 pts pr site
#20 sites
#follow_up 1:21 year



#################
### Libraries ###
#################---------------------------------------------------------------
library(tidyverse)
library(survival)
library(simsurv)
library(splines)
library(coxed)
library(devtools)
library(riskRegression)
library(Epi)
library(popEpi)
library(mgcv)
library(testpackage)
library(foreach)
library(doParallel)
library(future)
library(future.apply)
library(pec)
#library(mypackage)
#library(evmix)
###########
### CPP ###
###########---------------------------------------------------------------------
library(Rcpp)
#sourceCpp('./C++/AFT_newton_fkt.cpp')

######################
### Specifications ###
######################
ayear <- 365


#################
### Functions ###
#################---------------------------------------------------------------

source("./Functions for simulation/SIM_functions.R")

sim_fkt_sim6 <- function(n_pat=100,
                    True_beta =c(X1=log(1.2),
                                 X2=log(0.8),
                                 X3=log(2.0),
                                 X4=log(0.7)),
                    mean_val=0,
                    SD_Val=1,
                    lambda=0.5,
                    uni_max_time=5,
                    gammas=2,
                    max_time=1000,ayear=365.25 ){
  
  covs <-  data.frame(id = 1:n_pat,
                      X1 = rnorm(n_pat,mean = mean_val,sd = SD_Val),
                      X2 = rnorm(n_pat,mean = mean_val,sd = SD_Val),
                      X3 = rnorm(n_pat,mean = mean_val,sd = SD_Val),
                      X4 = rnorm(n_pat,mean = mean_val,sd = SD_Val)
  )
  
  dat <- simsurv::simsurv(dist = "weibull",#"weibull",
                          lambdas = lambda, #Scale parameter
                          gammas = gammas, #shape parameter
                          betas = True_beta,
                          x = covs,
                          maxt = max_time,
                          interval = c(0,max_time+1),
  )
  
  
  #simulating the censoring mechanism with an exponential dist
  censoring <- runif(n = 1:nrow(dat),min=0,max=uni_max_time )
  
  
  #Making a censoring time
  dat$censoring <- censoring#$eventtime
  #Selecting if eventtime or censoring times occurs first as the status variable
  dat$status <- ifelse(dat$eventtime < dat$censoring, 1, 0)
  dat$eventtime <- pmin(dat$eventtime, dat$censoring)
  dat$censoring <- NULL
  
  #cmbinning into a dataframe
  covs$id <- NULL
  dat$id <- NULL
  dat <- cbind(dat,covs) %>% as.data.frame()
  dat$eventtime_true <- dat$eventtime
  dat$eventtime <- ceiling(dat$eventtime*ayear)
  return(dat) 
}

SIM6_ftk <- function(N_sim,n_pts,N_sites,follow_up,N_cov,Test_set_nr_pts,predtimes,True_beta,lambda=0.5,gammas=2,ayear=365.25,mean_val=0,SD_val=1,
                     envv=ls(globalenv())){
  
  
  #n.cores <- parallel::detectCores() - 1
  n.cores <-50
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)
  
  #loop though each simulation
  PS <- foreach(
    i = 1:N_sim,
    .combine = 'rbind',
    .packages = (.packages()),
    .export = envv
  )%dopar% {#%do%{##
    
    # for (i in 1:N_sim) {
    
    
    #Simulating a list of datasets 
    Sim_data <- lapply(X=1:N_sites,FUN =  function(X){sim_fkt_sim6(n_pat = n_pts,
                                                                   max_time = max(predtimes),
                                                                   uni_max_time = follow_up, #followup
                                                                   True_beta =True_beta,
                                                                   SD_Val = SD_val,
                                                                   mean_val =mean_val[X], 
                                                                   lambda = lambda,
                                                                   gammas = gammas )} )
    
    #Check if there is a df without any events
    Sim_data <- lapply(X = Sim_data,FUN=function(i){
      if (sum(i$status)>0) {return(i)}
      NULL
    })
    
    #removes the dataset without events
    Sim_data <- Sim_data[(which(sapply(Sim_data,negate(is.null) ),arr.ind=TRUE))]
    
    #updating numbers of sites
    if (length(Sim_data)!=N_sites) {
      N_sites <- length(Sim_data)
    }
    
    #Old 
    # Test_set <- sim_fkt(n_pat = Test_set_nr_pts,
    #                     True_beta =True_beta,
    #                     lambda = lambda,
    #                     gammas = gammas )
    
    #Building a testset with equal number of pts from each site with differen mean value for the simulated pats, to mimic if a test set was pooled
    Test_set <- lapply(X=1:N_sites,FUN =  function(X){sim_fkt_sim6(n_pat = 1000/N_sites, #To ensure an equal number of pts from each site, and that the total number of pts are 1000
                                                                   max_time = max(predtimes),
                                                                   uni_max_time = follow_up, #followup
                                                                   True_beta =True_beta,
                                                                   SD_Val = SD_val,
                                                                   mean_val =mean_val[X], 
                                                                   lambda = lambda,
                                                                   gammas = gammas )} )
    Test_set <- do.call("rbind",Test_set)
    
    #Making the full data frame with all data
    full_DF <- do.call("rbind",Sim_data)
    full_DF.s <- Full_Df_fkt(full_DF_input = full_DF )
    
    t.kn_full <- full_DF$eventtime %>% quantile(probs = seq(0,0.99,1/5) ) %>% as.numeric()
    
    
    #Preparing simulated data so it on lexis form. 
    preClean_Simulated_df <- lapply(X=Sim_data,FUN = function(X){Lexis_df_cleaning(X,knots = t.kn_full)})
    
    #Initial beta
    beta_old <- rep(0,ncol(preClean_Simulated_df[[1]]$DF))
    beta_old_cox <- rep(0,4)
    
    #Beta poisson
    beta_fl <- Beta_AFT_FL_fkt(preClean_Simulated_df =preClean_Simulated_df,
                               beta_old = beta_old )
    
    
    #Beta Naive cox 
    beta_COX <- Beta_Cox_FL_fkt(Sim_data = Sim_data,beta_old = beta_old_cox)
    
    FIT <- glm(cbind(lex.Xst == 1, lex.dur)~ Ns(tfe, knots = t.kn_full) + X1+X2+X3+X4,
               family = poisreg,
               data = full_DF.s,
               x = T)
    FIT_COX <- coxph(Surv(eventtime,status)~ X1+X2+X3+X4,data=full_DF,x=T)
    
    
    #Prediction matriices
    Qmatrix_FL_pois <- Pred_FL_Poisson_fkt(df = Test_set,
                                           beta_old = beta_fl,
                                           predtimes = predtimes,
                                           t.kn = t.kn_full)
    
    Qmatrix_Pool_pois <- Pred_FL_Poisson_fkt(df = Test_set,
                                             beta_old =  matrix(coef(FIT)),
                                             predtimes = predtimes,
                                             t.kn = t.kn_full)
    
    
    QMatrix_FL_Cox <- pred_FL_COX_fkt(Sim_data = Sim_data,
                                      test_df = Test_set,
                                      beta_old = beta_COX,
                                      predtimes = predtimes)
    
    QMatrix_Pool_COX <- predictSurvProb(FIT_COX,newdata = Test_set,times = predtimes)
    
    
    #Fitting the true model
    QMatrix_True <- TRUE_model_pred(df = Test_set,
                                    lambda =lambda,
                                    gamma = gammas,
                                    True_beta =True_beta,
                                    ayear = ayear,
                                    pred_time =predtimes,
                                    test = Test_set )
    
    
    
    #Cindex
    C_index_FL_pois <- rcorr.cens(Qmatrix_FL_pois[,max(predtimes)],with(Test_set,Surv(eventtime,status)))[1]
    C_index_Pool_pois <- rcorr.cens(Qmatrix_Pool_pois[,max(predtimes)],with(Test_set,Surv(eventtime,status)))[1]
    C_index_FL_COX <- rcorr.cens(QMatrix_FL_Cox[,max(predtimes)],with(Test_set,Surv(eventtime,status)))[1]
    C_index_Pool_COX <- rcorr.cens(QMatrix_Pool_COX[,max(predtimes)],with(Test_set,Surv(eventtime,status)))[1]
    C_index_TRUE <- rcorr.cens(QMatrix_True[,max(predtimes)],with(Test_set,Surv(eventtime_true,status)))[1]
    
    
    
    
    #IBS 
    IBS_FL_pois <- ibs(pec(Qmatrix_FL_pois,
                           data = Test_set,
                           exact = F,
                           formula = Surv(eventtime,status)~1,
                           times = predtimes),start = 0,times =max(predtimes) )[[2]]
    IBS_pool_pois <- ibs(pec(Qmatrix_Pool_pois,
                             data = Test_set,
                             exact = F,
                             formula = Surv(eventtime,status)~1,
                             times = predtimes),start = 0,times =max(predtimes) )[[2]]
    IBS_FL_Cox <- ibs(pec(QMatrix_FL_Cox,
                          data = Test_set,
                          exact = F,
                          formula = Surv(eventtime,status)~1,
                          times = predtimes),start = 0,times =max(predtimes) )[[2]]
    IBS_Pool_COX <- ibs(pec(QMatrix_Pool_COX,
                            data = Test_set,
                            exact = F,
                            formula = Surv(eventtime,status)~1,
                            times = predtimes),start = 0,times =max(predtimes) )[[2]]
    IBS_true_model <- ibs(pec(QMatrix_True,
                              data = Test_set,
                              exact = F,
                              formula = Surv(eventtime,status)~1,
                              times = predtimes),start = 0,times =max(predtimes) )[[2]]
    
    
    #The average survival difference
    Avg_Surv_diff_FL_POIS <- abs(QMatrix_True-Qmatrix_FL_pois)%>% rowSums() %>% mean()
    Avg_Surv_diff_Pool_Pois <- abs(QMatrix_True-Qmatrix_Pool_pois)%>% rowSums() %>% mean()
    
    Avg_Surv_diff_FL_Cox <- abs(QMatrix_True-QMatrix_FL_Cox)%>% rowSums() %>% mean()
    Avg_Surv_diff_Pool_cox <- abs(QMatrix_True-QMatrix_Pool_COX)%>% rowSums() %>% mean()
    
    
    #Collecting all performance scores
    performance_Scores <- data.frame(
      C_index_FL_pois = C_index_FL_pois,
      C_index_Pool_pois = C_index_Pool_pois,
      C_index_FL_COX = C_index_FL_COX,
      C_index_Pool_COX = C_index_Pool_COX,
      C_index_TRUE = C_index_TRUE,
      
      IBS_FL_pois = IBS_FL_pois,
      IBS_pool_pois = IBS_pool_pois,
      IBS_FL_Cox = IBS_FL_Cox,
      IBS_Pool_COX = IBS_Pool_COX,
      IBS_true_model = IBS_true_model,
      
      Avg_Surv_diff_FL_POIS = Avg_Surv_diff_FL_POIS,
      Avg_Surv_diff_Pool_Pois = Avg_Surv_diff_Pool_Pois,
      Avg_Surv_diff_FL_Cox = Avg_Surv_diff_FL_Cox,
      Avg_Surv_diff_Pool_cox = Avg_Surv_diff_Pool_cox
      
    )
    performance_Scores
  }
  parallel::stopCluster(cl = my.cluster)
  # C_index <- data.frame(C_index_FL=C_index_FL,
  #                       C_index_FULL=C_index_FULL)
  # 
  return(PS)
}




#######################
### Simulating Data ###
#######################
#Specifications---
n_pts <- 20
N_sites <- 20
follow_up <- 5
N_cov <- 4
Test_set_nr_pts <- 1000
predtimes <- 0:1000
mean_val <- seq(0,5,0.2)[1:N_sites]
True_beta =c(X1=log(1.2),
             X2=log(0.8),
             X3=log(2.0),
             X4=log(0.7))
lambda=0.5
gammas=2




system.time({
  sim6_test <- SIM6_ftk(
    N_sim = 2,
    n_pts = n_pts,
    mean_val = mean_val,
    N_sites = N_sites,
    N_cov = N_cov,
    follow_up = follow_up,
    Test_set_nr_pts = Test_set_nr_pts,
    predtimes = predtimes,
    True_beta = True_beta
  )
  
})

saveRDS(sim6_test,"./Simulation 6 - Mean value/sim6_test.rds")




#################
### Test area ###
#################
C_index_sim4 <- C_index_ceaning_fkt(df = sim6_test)
IBS_sim4 <- IBS_df_cleaning_fkt(df = sim6_test)
ABS_Surv_diff_sim4 <- ABS_AVG_SD_fkt(df = sim6_test)

C_index_S4 <- Cindex_box_plot_fkt(C_index_res =C_index_sim4,title="Simulation 4",exclude_y_axis = T,ylim = c(0.70,0.9) )
IBS_S4 <- IBS_box_plot_fkt(IBS_sim4,title="Simulation 4",exclude_y_axis = T ,ylim = c(0.05,0.1) )
ABS_Avg_s4 <- ABS_surv_df_box_plot_fkt(df =ABS_Surv_diff_sim4,title="Simulation 4",exclude_y_axis = T )

sim6_test <- lapply(X=1:N_sites,FUN =  function(X){sim_fkt_sim6(n_pat = n_pts,
                                                                max_time = max(predtimes),
                                                                uni_max_time = follow_up, #followup
                                                                True_beta =True_beta,
                                                                mean_val =mean_val[X], 
                                                                lambda = lambda,
                                                                gammas = gammas )} )


dff <- do.call("rbind",sim6_test)

coxph(Surv(eventtime,status)~X1+X2+X3+X4,data=dff)


cox_l <- lapply(X=1:N_sites,FUN=function(X){
  coef(coxph(Surv(eventtime,status)~X1+X2+X3+X4,data=sim6_test[[X]]))  
  
})

do.call("rbind", cox_l) %>% colMeans()

coxph(Surv(eventtime,status)~X1+X2+X3+X4,data=sim6_test[[2]])
