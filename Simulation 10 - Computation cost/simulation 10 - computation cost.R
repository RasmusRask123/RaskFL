###############################################
### Simulation 10 Iteration cost             ###
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

sim_fkt_sim10 <- function(n_pat=100,
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

SIM10_ftk <- function(N_sim,n_pts,N_sites,follow_up,N_cov,Test_set_nr_pts,predtimes,True_beta,lambda=0.5,gammas=2,ayear=365.25,mean_val=0,SD_val=1,
                     envv=ls(globalenv())){
  
  
  n.cores <- parallel::detectCores() - 1
  #n.cores <-50
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
  )%dopar% {#%do%{##%dopar% {#
    
    # for (i in 1:N_sim) {
    
    
    #Simulating a list of datasets 
    Sim_data <- lapply(X=1:N_sites,FUN =  function(X){sim_fkt_sim10(n_pat = n_pts[X],
                                                                   max_time = max(predtimes),
                                                                   uni_max_time = follow_up[X], #followup
                                                                   True_beta =True_beta,
                                                                   SD_Val = SD_val[X],
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
    
    Test_set <- sim_fkt_sim10(n_pat = 1000, #To ensure an equal number of pts from each site, and that the total number of pts are 1000
                             max_time = max(predtimes),
                             uni_max_time = follow_up, #followup
                             True_beta =True_beta,
                             SD_Val = SD_val,
                             mean_val =mean_val, 
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
    
    #Beta poisson (returns number of computations)
    beta_fl <- Beta_AFT_FL_fkt(preClean_Simulated_df =preClean_Simulated_df,
                               beta_old = beta_old,Computation=TRUE )
    
    
    #Beta Naive cox 
    beta_COX <- Beta_Cox_FL_fkt(Sim_data = Sim_data,beta_old = beta_old_cox,Computation=TRUE)

    Iteration_df <- data.frame(
      FL_pois_iterations=beta_fl,
      FL_cox_iterations=beta_COX
    )
  
    
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
mean_val <- 0
SD_val <-1
True_beta =c(X1=log(1.2),
             X2=log(0.8),
             X3=log(2.0),
             X4=log(0.7))

lambda=0.5
gammas=2


system.time({
  sim10_test <- SIM10_ftk(
    N_sim = 1000,
    n_pts = n_pts,
    mean_val = mean_val,
    SD_val = SD_val,
    N_sites = N_sites,
    N_cov = N_cov,
    follow_up = follow_up,
    Test_set_nr_pts = Test_set_nr_pts,
    predtimes = predtimes,
    True_beta = True_beta
  )
  
})

saveRDS(sim10_test,"./Simulation 10 - Computation cost/sim10_test.rds")

############################
### SIM 1 Iteration cost ###
############################
N_sites <- 30
n_pts <- rep(20,N_sites)
follow_up <- rep(5,N_sites)
N_cov <- 4
Test_set_nr_pts <- 1000
predtimes <- 0:1000
mean_val <- rep(0,N_sites)
SD_val <-rep(1,N_sites)
True_beta =c(X1=log(1.2),
             X2=log(0.8),
             X3=log(2.0),
             X4=log(0.7))

lambda=0.5
gammas=2

system.time({
  sim1_test_it <- SIM10_ftk(
    N_sim = 1000,
    n_pts = n_pts,
    mean_val = mean_val,
    SD_val = SD_val,
    N_sites = N_sites,
    N_cov = N_cov,
    follow_up = follow_up,
    Test_set_nr_pts = Test_set_nr_pts,
    predtimes = predtimes,
    True_beta = True_beta
  )
  
})
saveRDS(sim1_test_it,"./Simulation 10 - Computation cost/sim1_test_it.rds")


############################
### SIM 2 Iteration cost ###
############################
N_sites <- 10
n_pts <- rep(100,N_sites)
follow_up <- 1:N_sites
N_cov <- 4
Test_set_nr_pts <- 1000
predtimes <- 0:1000
mean_val <- rep(0,N_sites)
SD_val <-rep(1,N_sites)
True_beta =c(X1=log(1.2),
             X2=log(0.8),
             X3=log(2.0),
             X4=log(0.7))

lambda=0.5
gammas=2

system.time({
  sim2_test_it <- SIM10_ftk(
    N_sim = 1000,
    n_pts = n_pts,
    mean_val = mean_val,
    SD_val = SD_val,
    N_sites = N_sites,
    N_cov = N_cov,
    follow_up = follow_up,
    Test_set_nr_pts = Test_set_nr_pts,
    predtimes = predtimes,
    True_beta = True_beta
  )
  
})
saveRDS(sim2_test_it,"./Simulation 10 - Computation cost/sim2_test_it.rds")



############################
### SIM 3 Iteration cost ###
############################
N_sites <- 20
n_pts <- seq(20,220,10)
follow_up <- rep(5,N_sites)
N_cov <- 4
Test_set_nr_pts <- 1000
predtimes <- 0:1000
mean_val <- rep(0,N_sites)
SD_val <-rep(1,N_sites)
True_beta =c(X1=log(1.2),
             X2=log(0.8),
             X3=log(2.0),
             X4=log(0.7))

lambda=0.5
gammas=2

system.time({
  sim3_test_it <- SIM10_ftk(
    N_sim = 1000,
    n_pts = n_pts,
    mean_val = mean_val,
    SD_val = SD_val,
    N_sites = N_sites,
    N_cov = N_cov,
    follow_up = follow_up,
    Test_set_nr_pts = Test_set_nr_pts,
    predtimes = predtimes,
    True_beta = True_beta
  )
  
})
saveRDS(sim3_test_it,"./Simulation 10 - Computation cost/sim3_test_it.rds")



############################
### SIM 4 Iteration cost ###
############################
N_sites <- 20
n_pts <- seq(20,220,10)
follow_up <-  1:20
N_cov <- 4
Test_set_nr_pts <- 1000
predtimes <- 0:1000
mean_val <- rep(0,N_sites)
SD_val <-rep(1,N_sites)
True_beta =c(X1=log(1.2),
             X2=log(0.8),
             X3=log(2.0),
             X4=log(0.7))

lambda=0.5
gammas=2

system.time({
  sim4_test_it <- SIM10_ftk(
    N_sim = 1000,
    n_pts = n_pts,
    mean_val = mean_val,
    SD_val = SD_val,
    N_sites = N_sites,
    N_cov = N_cov,
    follow_up = follow_up,
    Test_set_nr_pts = Test_set_nr_pts,
    predtimes = predtimes,
    True_beta = True_beta
  )
  
})
saveRDS(sim4_test_it,"./Simulation 10 - Computation cost/sim4_test_it.rds")


############################
### SIM 5 Iteration cost ###
############################

#We have not done sim 5, ###


############################
### SIM 6 Iteration cost ###
############################
N_sites <- 20
n_pts <- rep(20,N_sites)
follow_up <-  rep(5,N_sites)
N_cov <- 4
Test_set_nr_pts <- 1000
predtimes <- 0:1000
mean_val <- seq(0,5,0.2)[1:N_sites]
SD_val <-rep(1,N_sites)
True_beta =c(X1=log(1.2),
             X2=log(0.8),
             X3=log(2.0),
             X4=log(0.7))

lambda=0.5
gammas=2

system.time({
  sim6_test_it <- SIM10_ftk(
    N_sim = 1000,
    n_pts = n_pts,
    mean_val = mean_val,
    SD_val = SD_val,
    N_sites = N_sites,
    N_cov = N_cov,
    follow_up = follow_up,
    Test_set_nr_pts = Test_set_nr_pts,
    predtimes = predtimes,
    True_beta = True_beta
  )
  
})
saveRDS(sim6_test_it,"./Simulation 10 - Computation cost/sim5_test_it.rds")




############################
### SIM 7 Iteration cost ###
############################
N_sites <- 20
n_pts <- rep(20,N_sites)
follow_up <-  rep(5,N_sites)
N_cov <- 4
Test_set_nr_pts <- 1000
predtimes <- 0:1000
mean_val <- rep(0,N_sites)
SD_val <- seq(1,20,0.2)[1:N_sites]
True_beta =c(X1=log(1.2),
             X2=log(0.8),
             X3=log(2.0),
             X4=log(0.7))

lambda=0.5
gammas=2

system.time({
  sim7_test_it <- SIM10_ftk(
    N_sim = 1000,
    n_pts = n_pts,
    mean_val = mean_val,
    SD_val = SD_val,
    N_sites = N_sites,
    N_cov = N_cov,
    follow_up = follow_up,
    Test_set_nr_pts = Test_set_nr_pts,
    predtimes = predtimes,
    True_beta = True_beta
  )
  
})
saveRDS(sim7_test_it,"./Simulation 10 - Computation cost/sim7_test_it.rds")









############################
### SIM 8 Iteration cost ###
############################
############################
### SIM 9 Iteration cost ###
############################
############################
### SIM 10 Iteration cost ###
############################


















############
### Test ###
############