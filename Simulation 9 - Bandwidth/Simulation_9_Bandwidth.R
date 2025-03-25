###############################################
### Simulation 9 Bandwidth changes ###
###############################################
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
#library(devtools)
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
library(survC1)
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

sim_fkt_sim8 <- function(n_pat=100,
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

SIM9_ftk <- function(N_sim,n_pts,N_sites,follow_up,
                     N_cov,Test_set_nr_pts,predtimes,
                     True_beta,
                     lambda=0.5,gammas=2,ayear=365.25,
                     mean_val=0,SD_val=1,bandwidth=10,
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
    
    #Simulating a list of datasets 
    Sim_data <- lapply(X=1:N_sites,FUN =  function(X){sim_fkt_sim8(n_pat = n_pts,
                                                                   max_time = max(predtimes),
                                                                   uni_max_time = follow_up, #followup
                                                                   True_beta =True_beta,
                                                                   SD_Val = SD_val,
                                                                   mean_val =mean_val, 
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
    

    #building a seperate testset
    Test_set <- sim_fkt_sim8(n_pat = 1000, 
                 max_time = max(predtimes),
                 uni_max_time = follow_up, #followup
                 True_beta =True_beta,
                 SD_Val = SD_val,
                 mean_val =mean_val, 
                 lambda = lambda,
                 gammas = gammas )
  
    #Making the full data frame with all data
    full_DF <- do.call("rbind",Sim_data)

    
    #Initial beta
    beta_old_cox <- rep(0,4)
    
    #Beta Naive cox 
    beta_COX <- Beta_Cox_FL_fkt(Sim_data = Sim_data,beta_old = beta_old_cox)
    

    FIT_COX <- coxph(Surv(eventtime,status)~ X1+X2+X3+X4,data=full_DF,x=T)
    
  
    #Prediction matriices
    #FL cox predictions with different bandwidth parameter
    QMatrix_FL_Cox_list <- lapply(X = 1:length(bandwidth),
                                  function(X) {
                                    pred_FL_COX_fkt(
                                      Sim_data = Sim_data,
                                      test_df = Test_set,
                                      beta_old = beta_COX,
                                      predtimes = predtimes,
                                      b = bandwidth[X]
                                    )
    })
    
    
    
    #Normal pooled Cox model
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
    #C-index for the fl cox models with different bandwidths
    C_index_FL_COX_list <- lapply(X=QMatrix_FL_Cox_list, function(X){
      Est.Cval(mydata =cbind(Test_set$eventtime,Test_set$status,-X[,max(predtimes)]),tau=1000 ,nofit = T)$Dhat
    })
   
    C_index_Pool_COX   <- Est.Cval(mydata =cbind(Test_set$eventtime,Test_set$status,-QMatrix_Pool_COX[,max(predtimes)]),tau=1000,nofit = T )$Dhat
    C_index_TRUE       <- Est.Cval(mydata =cbind(Test_set$eventtime,Test_set$status,-QMatrix_True[,max(predtimes)]),tau=1000,nofit = T )$Dhat
    
    
    
    
    #IBS 
    IBS_FL_Cox_list <- lapply(X=QMatrix_FL_Cox_list, function(X){
      ibs(pec(X,
              data = Test_set,
              exact = F,
              formula = Surv(eventtime,status)~1,
              times = predtimes),start = 0,times =max(predtimes) )[[2]]
    })
    
    
    
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
    Avg_Surv_diff_FL_Cox_list <- lapply(X=QMatrix_FL_Cox_list, function(X){
      abs(QMatrix_True-X)%>% rowSums() %>% mean()
    })
    
    
    
    #Combinning perofrmance measures for the list of different bandwidth parameters
    C_index_FL_COX       <- as.data.frame(do.call("cbind", C_index_FL_COX_list))
    IBS_FL_Cox           <- as.data.frame(do.call("cbind", IBS_FL_Cox_list))
    Avg_Surv_diff_FL_Cox <- as.data.frame(do.call("cbind", Avg_Surv_diff_FL_Cox_list))

    #Renaming colloumns to somthing sensible
    names(C_index_FL_COX) <- paste("C-FL-cox (b=",bandwidth,")",sep = "")
    names(IBS_FL_Cox) <- paste("IBS-FL-cox (b=",bandwidth,")",sep = "")
    names(Avg_Surv_diff_FL_Cox) <- paste("MIAD-FL-cox (b=",bandwidth,")",sep = "")
    
    
    
    Avg_Surv_diff_Pool_cox <- abs(QMatrix_True-QMatrix_Pool_COX)%>% rowSums() %>% mean()
    
  
    #Collecting all performance scores
    performance_Scores <- data.frame(
     
      C_index_Pool_COX = C_index_Pool_COX,
      C_index_TRUE = C_index_TRUE,
      
      IBS_Pool_COX = IBS_Pool_COX,
      IBS_true_model = IBS_true_model,
      
      Avg_Surv_diff_Pool_cox = Avg_Surv_diff_Pool_cox
      
    )
    performance_Scores <-cbind(performance_Scores,
                               C_index_FL_COX,
                               IBS_FL_Cox,
                               Avg_Surv_diff_FL_Cox)
    
    performance_Scores
  }
  parallel::stopCluster(cl = my.cluster)
  
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
mean_val <- 0
SD_val <-1
True_beta =c(X1=log(1.2),
             X2=log(0.8),
             X3=log(2.0),
             X4=log(0.7))
bandwidth <- c(4,10,20,40,50)

lambda=0.5
gammas=2


system.time({
  sim9_test <- SIM9_ftk(
    N_sim = 1000,
    n_pts = n_pts,
    mean_val = mean_val,
    bandwidth = bandwidth,
    SD_val = SD_val,
    N_sites = N_sites,
    N_cov = N_cov,
    follow_up = follow_up,
    Test_set_nr_pts = Test_set_nr_pts,
    predtimes = predtimes,
    True_beta = True_beta
  )
  
})

saveRDS(sim9_test,"./Simulation 9 - Bandwidth/sim9_test.rds")


#################
### TEST site ###
#################
Sim_data <- lapply(X=1:N_sites,FUN =  function(X){sim_fkt_sim8(n_pat = n_pts,
                                                               max_time = max(predtimes),
                                                               uni_max_time = follow_up, #followup
                                                               True_beta =True_beta,
                                                               SD_Val = SD_val,
                                                               mean_val =mean_val, 
                                                               lambda = lambda,
                                                               gammas = gammas )} )


qq_mat1 <- pred_FL_COX_fkt(Sim_data = Sim_data,
                test_df = Sim_data[[1]],
                beta_old = True_beta,
                predtimes = predtimes,b=10)
qq_mat2 <- pred_FL_COX_fkt(Sim_data = Sim_data,
                          test_df = Sim_data[[1]],
                          beta_old = True_beta,
                          predtimes = predtimes,b=20)
qq_mat3 <- pred_FL_COX_fkt(Sim_data = Sim_data,
                           test_df = Sim_data[[1]],
                           beta_old = True_beta,
                           predtimes = predtimes,b=30)
qq_mat4 <- pred_FL_COX_fkt(Sim_data = Sim_data,
                           test_df = Sim_data[[1]],
                           beta_old = True_beta,
                           predtimes = predtimes,b=50)
qq_mat5 <- pred_FL_COX_fkt(Sim_data = Sim_data,
                           test_df = Sim_data[[1]],
                           beta_old = True_beta,
                           predtimes = predtimes,b=100)

plot(qq_mat1[4,],type="s")
lines(qq_mat2[4,],col="red")
lines(qq_mat3[4,],col="blue")
lines(qq_mat4[4,],col="green")
lines(qq_mat5[4,],col="purple")
