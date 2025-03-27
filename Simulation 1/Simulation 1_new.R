####################
### Simulation 1 ###
####################
#1000 simulations
#20 pts pr site
#30 sites




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
library(ggsci)
library(gridExtra)
library(ggpubr)

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
#Function for simulating data
source("./Functions for simulation/SIM_functions.R")



SIM1_ftk <- function(N_sim,n_pts,N_sites,N_cov,Test_set_nr_pts,predtimes,True_beta,lambda=0.5,gammas=2,ayear=365.25,
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
  )%do%{##%dopar% {#
    
    # for (i in 1:N_sim) {
    
 
    #Simulating a list of datasets 
    Sim_data <- lapply(X=1:N_sites,FUN =  function(X){sim_fkt(n_pat = n_pts,
                                                              max_time = max(predtimes),
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
# Simulating Data ####
#######################
#Specifications
n_pts <- 20
N_sites <- 30
N_cov <- 4
Test_set_nr_pts <- 1000
predtimes <- 0:1000
True_beta =c(X1=log(1.2),
             X2=log(0.8),
             X3=log(2.0),
             X4=log(0.7))
lambda=0.5
gammas=2


sim1_test <- SIM1_ftk(
  N_sim = 1,
  n_pts = n_pts,
  N_sites = N_sites,
  N_cov,
  Test_set_nr_pts = Test_set_nr_pts,
  predtimes = predtimes,
  True_beta = True_beta
)


system.time({
  sim1_test <- SIM1_ftk(N_sim = 1000,n_pts = n_pts,N_sites = N_sites,N_cov,Test_set_nr_pts = Test_set_nr_pts,predtimes = predtimes,True_beta = True_beta )   
})

saveRDS(sim1_test,file = "./Simulation 1/NewSimulations/sim1_test_results_new_weightsBZ.rds")


saveRDS(sim1_test,"C:/Users/rasmu/OneDrive/Skrivebord/R ting/Simulationer - zip/Simulationer - zip/Simulation 1/sim1_test.rds")
saveRDS(sim1_test,file = "./Simulation 1/sim1_test_new2.rds")
saveRDS(sim1_test,file = "./Simulation 1/sim1_test_new3.rds")

########################
### Kaplanmier plots ###
########################
n_pts <- 20
N_sites <- 30
N_cov <- 4
Test_set_nr_pts <- 1000
predtimes <- 0:1000
True_beta =c(X1=log(1.2),
             X2=log(0.8),
             X3=log(2.0),
             X4=log(0.7))
sd_multi <- 3
lambda <- 0.5
gammas <- 2
ayear <- 365.25

#new function
sim_fkt2 <- function(n_pat=100,
                    True_beta =c(X1=log(1.2),
                                 X2=log(0.8),
                                 X3=log(2.0),
                                 X4=log(0.7)
                    ),
                    lambda=0.5,
                    uni_max_time=5,
                    gammas=2,
                    max_time=1000,ayear=365.25,sd_multi=1 ){
  
  covs <-  data.frame(id = 1:n_pat,
                      X1 = rnorm(n_pat,mean = 0,sd = 1*sd_multi),
                      X2 = rnorm(n_pat,mean = 0,sd = 1*sd_multi),
                      X3 = rnorm(n_pat,mean = 0,sd = 1*sd_multi),
                      X4 = rnorm(n_pat,mean = 0,sd = 1*sd_multi)
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

#Simulating a list of datasets 
Sim_data <- lapply(X=1:N_sites,FUN =  function(X){sim_fkt2(n_pat = n_pts,
                                                          max_time = max(predtimes),
                                                          True_beta =True_beta,
                                                          lambda = lambda,
                                                          gammas = gammas,sd_multi = sd_multi )} )

Validation_data <- sim_fkt2(n_pat = Test_set_nr_pts,
                    True_beta =True_beta,
                    lambda = lambda,
                    gammas = gammas,sd_multi = sd_multi )

full_DF <- do.call("rbind",Sim_data)
full_DF.s <- Full_Df_fkt(full_DF_input = full_DF )

t.kn_full <- full_DF$eventtime %>% quantile(probs = seq(0,0.99,1/5) ) %>% as.numeric()


#Preparing simulated data so it on lexis form. 
preClean_Simulated_df <- lapply(X=Sim_data,
                                FUN = function(X){Lexis_df_cleaning(X,knots = t.kn_full)})

#Initial beta
beta_old <- rep(0,ncol(preClean_Simulated_df[[1]]$DF))
beta_old_cox <- rep(0,4)

#Beta poisson
beta_fl <- Beta_AFT_FL_fkt(preClean_Simulated_df =preClean_Simulated_df,
                           beta_old = beta_old )

Validat <- Lexis_df_cleaning(full_DF,knots = t.kn_full)

beta_Pool_pois <- Beta_AFT_FL_fkt(preClean_Simulated_df =list(Validat),
                                  beta_old = beta_old )

#Beta Naive cox 
beta_COX <- Beta_Cox_FL_fkt(Sim_data = Sim_data,beta_old = beta_old_cox)


FIT <- glm(cbind(lex.Xst == 1, lex.dur)~ Ns(tfe, knots = t.kn_full) + X1+X2+X3+X4,
           family = poisreg,
           data = full_DF.s,
           x = T)
FIT_COX <- coxph(Surv(eventtime,status)~ X1+X2+X3+X4,data=full_DF,x=T)

#makring the prediction matrix
# Qmatrix <- pred_FL_AFT_fkt(df = Test_set,
#                            beta_old = beta_fl,
#                            predtimes = predtimes,
#                            t.kn = t.kn_full)
#Prediction matrices
Qmatrix_FL_pois <- Pred_FL_Poisson_fkt(df = Validation_data,
                                       beta_old = beta_fl,
                                       predtimes = predtimes,
                                       t.kn = t.kn_full)

Qmatrix_Pool_pois <- Pred_FL_Poisson_fkt(df = Validation_data,
                                         beta_old =  matrix(coef(FIT)),
                                         predtimes = predtimes,
                                         t.kn = t.kn_full)
# Qmatrix_Pool_pois <- Pred_FL_Poisson_fkt(df = Validation_data,
#                                          beta_old =  beta_Pool_pois,
#                                          predtimes = predtimes,
#                                          t.kn = t.kn_full)

QMatrix_FL_Cox <- pred_FL_COX_fkt(Sim_data = Sim_data,
                                  test_df = Validation_data %>% as.data.frame(),
                                  beta_old = beta_COX,
                                  predtimes = predtimes)

QMatrix_Pool_COX <- predictSurvProb(FIT_COX,newdata = Validation_data,times = predtimes)
#Fitting the true model
QMatrix_True <- TRUE_model_pred(df = Validation_data,
                                lambda =lambda,
                                gamma = gammas,
                                True_beta =True_beta,
                                ayear = ayear,
                                pred_time =predtimes,
                                test = Validation_data )
library(ggplot2)

pts_nr <- 1
pts_nr <- 4
pred_fl_pois <- data.frame(Model="FL-Poisson",
                           Pred=Qmatrix_FL_pois[pts_nr,],
                           Predtime=predtimes
)
pred_Pool_pois <- data.frame(Model="Poisson-Pool",
                             Pred=Qmatrix_Pool_pois[pts_nr,],
                             Predtime=predtimes
)
pred_fl_cox <- data.frame(Model="FL-Cox",
                          Pred=QMatrix_FL_Cox[pts_nr,],
                          Predtime=predtimes
)
pred_Pool_cox <- data.frame(Model="Cox-Pool",
                            Pred=QMatrix_Pool_COX[pts_nr,],
                            Predtime=predtimes
)
pred_true <- data.frame(Model="TRUE",
                            Pred=QMatrix_True[pts_nr,],
                            Predtime=predtimes
)

pred_df <- rbind(pred_fl_pois,pred_fl_cox,pred_Pool_cox,pred_Pool_pois,pred_true) %>% as.data.frame()


p1 <- ggplot(data = pred_df)+
    geom_line(aes(x=Predtime,y=Pred,col=Model),size=1)  +
    xlim(x=c(0,250))+
    ylim(c(0,1))

p1

#############################
### Plots of hazard scale ###
#############################
#run the section before
######################
### Specifications ###
######################
n_pts <- seq(20,20*10,20)
N_sites <- 10
N_cov <- 4
Test_set_nr_pts <- 1000
predtimes <- 0:1000
True_beta =c(X1=log(1.2),
             X2=log(0.8),
             X3=log(2.0),
             X4=log(0.7))
sd_multi <- 3
lambda <- 0.5
gammas <- 2
ayear <- 365.25


height <-  2160*1.1
width <- 3840
res <- 300

predtimes <- 0:1000
Sim_data <- lapply(X=1:N_sites,FUN =  function(X){sim_fkt2(n_pat = n_pts[X],
                                                           max_time = max(predtimes),
                                                           True_beta =True_beta,
                                                           lambda = lambda,
                                                           gammas = gammas,sd_multi = sd_multi )} )


fit <-  coxph(Surv(eventtime, status) ~ X1+X2+X3+X4, data = Sim_data[[3]], init = c(beta_COX),iter.max=0,x=T,y=T)

bz_true <- basehaz(fit,centered = F)

bz <- predictCox(object = fit,type = "hazard",centered = F,keep.times = T,times = predtimes)
bz_haz_true <- data.frame(time=bz$time,Hazard=bz$hazard)
bz_kernel1 <- BZKernel_fkt(cox_obj =fit,
                          pred_seq = predtimes,
                          df=Sim_data[[3]],b = 25)


bz_kernel2 <- BZKernel_fkt(cox_obj =fit,
                          pred_seq = predtimes,
                          df=Sim_data[[3]],b = 10)

plot(x=bz_kernel1$time,y=bz_kernel1$smoothhazard,type="l")
lines(x=bz$time,y=bz$hazard,lty=2)

BZ_test1 <- BZ_kernel_FL(BZ_list =list(bz_kernel1),
                        pred_seq=predtimes )
BZ_test2 <- BZ_kernel_FL(BZ_list =list(bz_kernel2),
                        pred_seq=predtimes )

plot(x=BZ_test1$time,y=BZ_test1$CumHazard,type="s")
lines(x=bz_true$time,y=bz_true$hazard,type="s")
plot(x=BZ_test2$time,y=BZ_test2$CumHazard,type="s")
lines(x=bz_true$time,y=bz_true$hazard,type="s")


p1 <- ggplot()+
  geom_line(data=bz_kernel1,aes(x=time,y=smoothhazard,color="a"),size=1.2) +
  geom_line(data=bz_haz_true,aes(x=time,y=Hazard/100,color="b"),size=1.2)+
  scale_y_continuous(
    name="Baseline hazard",
    sec.axis=sec_axis(~.*100,name="Baseline hazard")
  )+
  theme_classic() +
  scale_color_manual(values = c("#A73030FF","#3B3B3BFF"),
                     labels = expression( "Kernel Smooth (b=25)","Breslow"),
                     name="")+
  theme_bw()+
  theme(#axis.title = element_blank(),
    #axis.title.x =  element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "top",
    legend.key.size = unit(1,"cm"),
    #legend.key = element_text(size=2),
    legend.text = element_text(size=20),
    axis.title.y.left = element_text(size=20,color="#A73030FF"),
    axis.title.y.right = element_text(size=20,color="#3B3B3BFF"),
    #axis.ticks = element_blank(),
    axis.text.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.x = element_text(size=20),
    panel.background = element_blank(),
    text = element_text(size=20),
    plot.title = element_text(hjust=0.5,size = 20) )+
  labs(x = "",
       y = "Baseline hazard")+
  guides(colour=guide_legend(override.aes = list(size=5)))
p1
  
p11 <- ggplot()+
  geom_line(data=bz_kernel2,aes(x=time,y=smoothhazard,color="a"),size=1.2) +
  geom_line(data=bz_haz_true,aes(x=time,y=Hazard/25,color="b"),size=1.2)+
  scale_y_continuous(
    name="Baseline hazard",
    sec.axis=sec_axis(~.*25,name="Baseline hazard")
  )+
  theme_classic() +
  scale_color_manual(values = c("#A73030FF","#3B3B3BFF"),
                     labels = expression( "Kernel Smooth (b=10)","Breslow"),
                     name="")+
  theme_bw()+
  theme(#axis.title = element_blank(),
    #axis.title.x =  element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "top",
    legend.key.size = unit(1,"cm"),
    #legend.key = element_text(size=2),
    legend.text = element_text(size=20),
    axis.title.y.left = element_text(size=20,color="#A73030FF"),
    axis.title.y.right = element_text(size=20,color="#3B3B3BFF"),
    #axis.ticks = element_blank(),
    axis.text.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.x = element_text(size=20),
    panel.background = element_blank(),
    text = element_text(size=20),
    plot.title = element_text(hjust=0.5,size = 20) )+
  labs(x = "",
       y = "Baseline hazard")+
  guides(colour=guide_legend(override.aes = list(size=5)))
p1
p11

p2 <- ggplot()+
  geom_line(data=BZ_test1,aes(x=time,y=CumHazard,color="a"),color="#A73030FF",size=1.2)+
  geom_step(data=bz_true,aes(x=time,y=hazard,color="b"),color="#3B3B3BFF",size=1.2) +
  theme_classic() +
 # scale_color_manual(values = c("#A73030FF","#3B3B3BFF"),
#                     labels = expression("Kernel Smooth (b=25)", "Breslow"),
#                     name="")+
  theme_bw()+
  theme(#axis.title = element_blank(),
    #axis.title.x =  element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "top",
    legend.key.size = unit(1,"cm"),
    #legend.key = element_text(size=2),
    legend.text = element_text(size=20),
    axis.title.y = element_text(size=20),
    #axis.ticks = element_blank(),
    axis.text.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.x = element_text(size=20),
    panel.background = element_blank(),
    text = element_text(size=20),
    plot.title = element_text(hjust=0.5,size = 20) )+
  labs(x = "Time",
       y = "Cummulated baseline hazard")+
  guides(colour=guide_legend(override.aes = list(size=5)))+
  xlim(c(0,1000))+
  ylim(c(0,2))   
p2
p22 <- ggplot()+
  geom_line(data=BZ_test2,aes(x=time,y=CumHazard,color="a"),color="#A73030FF",size=1.2)+
  geom_step(data=bz_true,aes(x=time,y=hazard,color="b"),color="#3B3B3BFF",size=1.2) +
  theme_classic() +
  theme_bw()+
  theme(#axis.title = element_blank(),
    #axis.title.x =  element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "top",
    legend.key.size = unit(1,"cm"),
    #legend.key = element_text(size=2),
    legend.text = element_text(size=20),
    axis.title.y = element_text(size=20),
    #axis.ticks = element_blank(),
    axis.text.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.x = element_text(size=20),
    panel.background = element_blank(),
    text = element_text(size=20),
    plot.title = element_text(hjust=0.5,size = 20) )+
  labs(x = "Time",
       y = "Cummulated baseline hazard")+
  guides(colour=guide_legend(override.aes = list(size=5)))+
  xlim(c(0,1000))+
  ylim(c(0,2))   

p2
p22


png(filename = "./Results/Figures/haz_thesis_Smooth_B_25.png",
    width =width,
    height =height ,
    res = res)
p1
dev.off()
png(filename = "/Results/Figures/haz_thesis_Smooth_B_10.png",
    width =width,
    height =height ,
    res = res)
p11
dev.off()

png(filename = "./Results/Figures/CUM_haz_thesis_Smooth_B_25.png",
    width =width,
    height =height ,
    res = res)
p2
dev.off()
png(filename = "./Results/Figures/CUM_haz_thesis_Smooth_B_10.png",
    width =width,
    height =height ,
    res = res)
p22
dev.off()

png(filename = "./Results/Figures/CombinedHAZZARDS_Smooth_B.png",
    width =width,
    height =height ,
    res = res)
ggarrange(p1,p11,
          p2,p22,
          ncol=2,nrow=2,align="hv",legend = "top")
dev.off()
ggarrange(p1,p2,
          ncol=1,nrow=2,common.legend = T,align = "hv")

####################
### Working area ###
####################
n_pts <- seq(20,20*10,20)
N_sites <- 10
N_cov <- 4
Test_set_nr_pts <- 1000
predtimes <- 0:1000
True_beta =c(X1=log(1.2),
             X2=log(0.8),
             X3=log(2.0),
             X4=log(0.7))
sd_multi <- 3
lambda <- 0.5
gammas <- 2
ayear <- 365.25


height <-  2160*1.1
width <- 3840
res <- 300

predtimes <- 0:1000
Sim_data <- lapply(X=1:N_sites,FUN =  function(X){sim_fkt2(n_pat = n_pts[X],
                                                           max_time = max(predtimes),
                                                           True_beta =True_beta,
                                                           lambda = lambda,
                                                           gammas = gammas,sd_multi = sd_multi )} )


fit1 <-  coxph(Surv(eventtime, status) ~ X1+X2+X3+X4, data = Sim_data[[1]], init = c(beta_COX),iter.max=0,x=T,y=T)
fit2 <-  coxph(Surv(eventtime, status) ~ X1+X2+X3+X4, data = Sim_data[[10]], init = c(beta_COX),iter.max=0,x=T,y=T)

bz_true <- basehaz(fit,centered = F)

bz <- predictCox(object = fit,type = "hazard",centered = F,keep.times = T,times = predtimes)

bz_haz_true <- data.frame(time=bz$times,Hazard=bz$hazard)

predtimes1 <- 0:1000
predtimes2 <- seq(0,1000,50)

bz_kernel1 <- BZKernel_fkt(cox_obj =fit1,
                           pred_seq = predtimes1,
                           df=as.data.frame(Sim_data[[1]]),b = 25,weight =1 )

bz_kernel2 <- BZKernel_fkt(cox_obj =fit1,
                           pred_seq = predtimes2,
                           df=as.data.frame(Sim_data[[1]]),b = 25,weight =1 )


p1 <- ggplot()+
  geom_point(data=bz_kernel1,aes(x=time,y=smoothhazard,color="a"),size=1.2) +
  geom_point(data=bz_kernel2,aes(x=time,y=smoothhazard,color="b"),size=2) +
  geom_line(data=bz_haz_true,aes(x=time,y=Hazard/100,color="c"),size=1.2)+
  scale_y_continuous(
    name="Baseline hazard",
    sec.axis=sec_axis(~.*100,name="Baseline hazard")
  )+
  theme_classic() +
  scale_color_manual(values = c("#A73030FF","#3B3B3BFF","blue"),
                     labels = expression( "Site 1","Site 2","Breslow"),
                     name="")+
  theme_bw()+
  theme(#axis.title = element_blank(),
    #axis.title.x =  element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "top",
    legend.key.size = unit(1,"cm"),
    #legend.key = element_text(size=2),
    legend.text = element_text(size=20),
    axis.title.y.left = element_text(size=20,color="#A73030FF"),
    axis.title.y.right = element_text(size=20,color="#3B3B3BFF"),
    #axis.ticks = element_blank(),
    axis.text.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.x = element_text(size=20),
    panel.background = element_blank(),
    text = element_text(size=20),
    plot.title = element_text(hjust=0.5,size = 20) )+
  labs(x = "",
       y = "Baseline hazard")+
  guides(colour=guide_legend(override.aes = list(size=5)))
p1





dd <- sim_fkt(n_pat = n_pts,max_time = max(predtimes),uni_max_time = 1)

dd.2 <- Full_Df_fkt(full_DF_input = dd )
knotss <- dd.2$eventtime %>% quantile(probs = seq(0,0.99,1/5) ) %>% as.numeric()
Lexis_df_cleaning(df = dd,knots = knotss)


##################################
### Some fl for a presentation ###
##################################
ggplot()+
  geom_line(data=bz_kernel1,aes(x=time,y=smoothhazard,color="a"),size=1.2) +
  geom_line(data=bz_haz_true,aes(x=time,y=Hazard/25,color="b"),size=1.2)+
  scale_y_continuous(
    name="Baseline hazard",
    sec.axis=sec_axis(~.*25,name="Baseline hazard")
  )+
  theme_classic() +
  scale_color_manual(values = c("#A73030FF",
                                "#3B3B3BFF"),
                     labels = expression( "Kernel Smooth (b=25)",
                                          "Breslow"),
                     name="")+
  theme_bw()+
  theme(#axis.title = element_blank(),
    #axis.title.x =  element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "top",
    legend.key.size = unit(1,"cm"),
    #legend.key = element_text(size=2),
    legend.text = element_text(size=20),
    axis.title.y.left = element_text(size=20,color="#A73030FF"),
    axis.title.y.right = element_text(size=20,color="#3B3B3BFF"),
    #axis.ticks = element_blank(),
    axis.text.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.x = element_text(size=20),
    panel.background = element_blank(),
    text = element_text(size=20),
    plot.title = element_text(hjust=0.5,size = 20) )+
  labs(x = "Time",
       y = "Baseline hazard")+
  guides(colour=guide_legend(override.aes = list(size=5)))



