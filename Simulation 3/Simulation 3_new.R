####################
### Simulation 3 ###
####################
#1000 simulations
#10,20,30,...,200 pts pr site
#20 sites
#follow_up 5 year



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
library(flexsurv)
#library(mypackage)
#library(evmix)
###########
### CPP ###
###########---------------------------------------------------------------------
library(Rcpp)
#sourceCpp('./C++/AFT_newton_fkt.cpp')

#################
### Functions ###
#################---------------------------------------------------------------


source("./Functions for simulation/SIM_functions.R")


SIM3_ftk <- function(N_sim,n_pts,N_sites,follow_up,N_cov,Test_set_nr_pts,predtimes,True_beta,lambda=0.5,gammas=2,ayear=365.25,
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
    Sim_data <- lapply(X=1:N_sites,FUN =  function(X){sim_fkt(n_pat = n_pts[X],
                                                              max_time = max(predtimes),
                                                              uni_max_time = follow_up,
                                                              True_beta =True_beta,
                                                              lambda = lambda,
                                                              gammas = gammas )} )
    
    
    Test_set <- sim_fkt(n_pat = Test_set_nr_pts,
                        True_beta =True_beta,
                        lambda = lambda,
                        gammas = gammas )
    
    
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
    #Avg_Surv_diff_POIS <- abs(Qmatrix-Qmatrix_fit)%>% rowSums() %>% mean()
    #Avg_Surv_diff_COX <- abs(QMatrix_Cox-QMatrix_FIT_COX)%>% rowSums() %>% mean()
    
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
n_pts <- seq(20,220,10)
N_sites <- 20
N_cov <- 4
Test_set_nr_pts <- 1000
predtimes <- 0:1000
follow_up <- 5

True_beta =c(X1=log(1.2),
             X2=log(0.8),
             X3=log(2.0),
             X4=log(0.7))
lambda=0.5
gammas=2

sim3_test <- SIM3_ftk(
  N_sim = 50,
  n_pts = n_pts,
  N_sites = N_sites,
  N_cov = N_cov,
  follow_up = follow_up,
  Test_set_nr_pts = Test_set_nr_pts,
  predtimes = predtimes,
  True_beta = True_beta
)



system.time({
  sim3_test <- SIM3_ftk(
    N_sim = 1000,
    n_pts = n_pts,
    N_sites = N_sites,
    N_cov = N_cov,
    follow_up = follow_up,
    Test_set_nr_pts = Test_set_nr_pts,
    predtimes = predtimes,
    True_beta = True_beta
  ) 
  
})
saveRDS(sim3_test,"./Simulation 3/NewSimulations/sim3_test_new_weightsBZ.rds")
saveRDS(sim3_test,"C:/Users/rasmu/OneDrive/Skrivebord/R ting/Simulationer - zip/Simulationer - zip/Simulation 3/sim3_test.rds")
saveRDS(sim3_test,"./Simulation 3/sim3_test_new.rds")
saveRDS(sim3_test,"./Simulation 3/NewSimulations/sim3_test_new.rds")

