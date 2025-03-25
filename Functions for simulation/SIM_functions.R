#############################
### Simulation Functions  ###
#############################
#This script contains the functions for simulations that is not simulation specific but used in all simulations



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
#Function for simulating data
sim_fkt <- function(n_pat=100,
                    True_beta =c(X1=log(1.2),
                                 X2=log(0.8),
                                 X3=log(2.0),
                                 X4=log(0.7)
                    ),
                    lambda=0.5,
                    uni_max_time=5,
                    gammas=2,
                    max_time=1000,ayear=365.25 ){
  
  covs <-  data.frame(id = 1:n_pat,
                      X1 = rnorm(n_pat,mean = 0,sd = 1),
                      X2 = rnorm(n_pat,mean = 0,sd = 1),
                      X3 = rnorm(n_pat,mean = 0,sd = 1),
                      X4 = rnorm(n_pat,mean = 0,sd = 1)
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
#Lexis cleaning functions takes a simulated dataset from sim_fkt and prepares it for Bendix carstens model
Lexis_df_cleaning <- function(df,knots){
  
  #Make an Lexis object
  df <- Lexis(exit = list(tfe = eventtime),exit.status = factor(status,labels = unique(df$status)),data = df)
  
  #Finding splits in all event times  
  #tfe_split <- unique(df$eventtime )[sample(1:length(unique(df$eventtime )),size = 100,replace = T )]
  #tfe_split <- seq(10,1825,10)#we choose a secquence of time points/bins acroos all sites
  
  #df <- splitMulti(df, tfe = c(0, sort(tfe_split)))
  tfe_split <- seq(10,1825,20)#we choose a secquence of time points/bins acroos all sites
  #df <- splitMulti(df, tfe = c(0, sort(unique(df$eventtime)))) 
  df <- splitMulti(df, tfe = c(0, sort(tfe_split)))
  
  #t.kn <- df$eventtime %>% quantile(probs = seq(0,0.99,1/5) ) %>% as.numeric()
  #t.kn <- knots
  
  X <- model.matrix(~Ns(tfe, knots = knots)+X1+X2+X3+X4,df)
  
  t_time <- df$lex.dur %>% as.numeric()
  d_status <- df$lex.Xst %>% as.numeric()-1
  list(DF=X,t_time=t_time,d_status=d_status)
}

#Making risksets
GetRiskSet <- function(time_of_interest,entry_vector,time_vector,event_vector) {
  
  ret <- which((time_of_interest >= entry_vector) & 
                 ((time_vector == time_of_interest & event_vector == 1) |
                    (time_vector > time_of_interest)))
  return(ret)
}
COX_grad_fkt <- function(i,df,X_data,entry_vector,Ts,event,beta_old,risk_mat,uniq_event_times){
  #X_D <- df[(df$eventtime==i & df$status==1),3:(ncol(df)]
  X_D <- df[(df$eventtime==i & df$status==1),c("X1","X2","X3","X4")]
  d_i <- X_D %>% nrow()
  
  sum_X_D <- colSums(X_D)
  
  #Making the risk set 
  X_riskSet <- as.matrix(X_data[risk_mat[,which(uniq_event_times==i)],])
  
  EXP_XB <- exp(X_riskSet %*% beta_old)
  t1 <- t(X_riskSet) %*% EXP_XB
  t2 <-  sum(EXP_XB)
  
  return(list(t1=t1,t2=t2,EXP_XB=EXP_XB,X_riskSet=X_riskSet,sum_X_D=sum_X_D,d_i=d_i))
}

cox_update <- function(df,beta_old,times = "eventtime",event = "status"){
  
  #All eventtimes
  Ts <-df[,times]
  
  #Status
  event <- df[,event]
  
  #Design matrix
  #X_data <- df[,3:ncol(df)]
  X_data <-df %>% select(X1,X2,X3,X4)
  
  #All unique eventtimes
  uniq_event_times <- unique(df[event==1,times])
  
  #Entry vector
  entry_vector <- rep(0,nrow(X_data))
  H_old <- matrix(0,ncol=length(beta),nrow=length(beta))
  
  #Making a risk matrix, for all riskSets at each unique eventtime
  RISK_SET_mat <- sapply(uniq_event_times, function(toi) Ts> 0& Ts >= toi)
  
  #Initial Gradient and HEssian
  Grad_old <- rep(0,length(beta_old) )
  H_old <- matrix(0,ncol=length(beta_old),nrow=length(beta_old))
  
  #looping through all uniqe eventtimes (i.e eventtimes were an event happens)
  for (i in uniq_event_times) {
    
    #Calculating gradient
    grad_list <- COX_grad_fkt(
      i = i,
      df = df,
      X_data = X_data,
      entry_vector = entry_vector ,
      Ts = Ts,
      event = event,
      beta_old = beta_old,
      risk_mat = RISK_SET_mat,
      uniq_event_times = uniq_event_times
    )
    
    #Gradient, 
    #Iteratively adding the hessians canculated at every unique eventtime
    Grad_old <- Grad_old+(grad_list$sum_X_D-grad_list$d_i*(grad_list$t1/grad_list$t2))
    
    #Hessian
    #Iteratively adding the hessians canculated at every unique eventtime
    H_old <- H_old+Hessian_fkt(X_riskSet = grad_list$X_riskSet,t2 =grad_list$t2 ,EXP_XB =grad_list$EXP_XB ,n_cov =length(beta_old),d_i = grad_list$d_i )
  }
  
  return(list(Gradient=Grad_old,Hessian=H_old))
}

#Fkt to combine everything nicely
Grad_var_fkt <- function(beta,data,times,event,var_inv=F,return_what=""){
  X <- data %>% select(X1,X2,X3,X4) %>% as.matrix() #,X5,X6,X7,X8
  
  #beta <-   matrix(rep(beta, ncol(X)), nrow = ncol(X))
  times <-  data %>% select(times) %>% as.matrix
  event <-  data %>% select(event) %>% as.matrix
  entry <- rep(0, dim(times)[1]) %>% as.matrix
  
  #Only retunr gradient or variance if specificed
  if (return_what=="Gradient") {
    gradient <- CoxGradient(beta1 =beta,Xs =X,entry = entry,Ts =times ,event = event)  
    return(gradient)
  }
  if (return_what=="Variance") {
    variance <- CoxVariance(beta1 =beta,Xs =X,entry = entry,Ts =times ,event = event,var_inv=var_inv)  
    return(variance)
  }
  
  
  gradient <- CoxGradient(beta1 =beta,Xs =X,entry = entry,Ts =times ,event = event) 
  variance <- CoxVariance(beta1 =beta,Xs =X,entry = entry,Ts =times ,event = event,var_inv=var_inv)  
  ret <- list(gradient=gradient,variance=variance)
  return(ret)
}
Full_Df_fkt <- function(full_DF_input){
  #Making the full data frame with all data
  
  full_DF_input <- Lexis(exit = list(tfe = eventtime),
                         exit.status = factor(status,labels = sort(unique(full_DF_input$status)) ),
                         data = full_DF_input)
  #tfe_split <- seq(10,1825,10)
  #tfe_split_full <- unique(full_DF_input$eventtime )[sample(1:length(unique(full_DF_input$eventtime )),size = 100,replace = T )]
  full_DF_input<- splitMulti(full_DF_input, tfe = c(0, sort(unique(full_DF_input$eventtime ))))
  #full_DF_input<- splitMulti(full_DF_input, tfe = c(0, tfe_split))
  
  
  
  return(full_DF_input)
}
grad_exstract_fkt <- function(updates){
  #Dimmension of covariates
  dimmension <- length(updates[[1]]$Gradient)
  
  #Empty gradients and hessians
  Gradient <- rep(0,dimmension)
  Hessian <- matrix(0,ncol=dimmension,nrow=dimmension)
  
  #adding gradients and hessians
  for (i in 1:length(updates)) {
    Gradient <- Gradient+updates[[i]]$Gradient
    Hessian <- Hessian+updates[[i]]$Hessian
  }
  return(list(Gradient=Gradient,Hessian=Hessian))  
}
Baseline_hazard_fkt <- function(df,beta,times = "eventtime",event = "status"){
  Hazard_fkt <- function(i,df,X_data,beta,risk_mat,uniq_event_times){
    X_D <- df[(df$eventtime==i & df$status==1),3:ncol(df)]
    d_i <- X_D %>% nrow()
    
    #Making the risk set 
    X_riskSet <- as.matrix(X_data[risk_mat[,which(uniq_event_times==i)],])
    
    EXP_XB <- exp(X_riskSet %*% beta)
    return(d_i/sum(EXP_XB))
    
  }
  #beta <- coef(fit)
  #All eventtimes
  Ts <-df[,times]
  
  #Status
  event <- df[,event]
  
  #Design matrix
  X_data <- df[,3:ncol(df)]
  
  #All unique eventtimes
  #uniq_event_times <- sort(unique(df[event==1,times]))
  uniq_event_times <- sort(unique(df[,times]))
  
  #Entry vector
  
  #Making a risk matrix, for all riskSets at each unique eventtime
  RISK_SET_mat <- sapply(uniq_event_times, function(toi) Ts> 0& Ts >= toi)
  
  hazard <- lapply(uniq_event_times,FUN = function(i){
    Hazard_fkt(
      i = i ,
      df = df ,
      X_data = X_data ,
      beta = beta,
      risk_mat = RISK_SET_mat,
      uniq_event_times = uniq_event_times)
  })
  hazard <- data.frame(times=uniq_event_times,Hazard=do.call(what="c",hazard))
  
  return(hazard)
}

Beta_AFT_FL_fkt <- function(preClean_Simulated_df,beta_old,Computation=F){
  
  for (k in 1:100) {
    
    #Calculating gradinet and HESSIAN for each site
    FL_update <- lapply(preClean_Simulated_df,
                        FUN = function(X){AFT_Newton_FL_fkt(beta = beta_old,
                                                            df = X$DF,
                                                            time = X$t_time,
                                                            status = X$d_status)
                        } )
    
    #gradients
    GRAD <- lapply(X=FL_update,
                   FUN = function(X)X$Gradient ) %>% do.call(what="cbind") %>% rowSums()
    #Hessian
    Hessian <- lapply(X=FL_update,
                      FUN = function(X)X$Hessian) %>% Reduce(f = '+')
    
    #Inverting the hessain
    H_inv <- solve(Hessian)
    
    #updating beta
    beta_new <- beta_old - H_inv%*%GRAD
    
    #if no large change then stop algorithm
    if ( norm(beta_new-beta_old,type="2") <0.0000000000001) {
      print(k)
      
      #Returns the numbers of iterations for each site
      if (Computation) {
        return(k)
      }
      
      break
    }
    
    beta_old <- beta_new
    
    
  }
  return(beta_old)
}
Beta_Cox_FL_fkt <- function(Sim_data,beta_old,Computation=F){
  
  Total_nr_patinets <- nrow(do.call("rbind", Sim_data))
  for (k in 1:100) {
    
    #Gradients and hessian from each dataset
    updates <- lapply(X=Sim_data, function(X){
      #Weighting
      Weight <- (nrow(X)/Total_nr_patinets)
      
      #Calculate gradients and hessians
      update <-  cox_update(df = X,beta = beta_old) 
      
      #Weighting the gradients and hessians
      update$Gradient <- update$Gradient*Weight
      update$Hessian <- update$Hessian*Weight
      
      return(update)
    })
    
    #combining all gradients and Hessian
    grad_exstract <- grad_exstract_fkt(updates = updates)
    
    #Combinning the graidnets and variances by the weights
    GRADIENT <- grad_exstract$Gradient#note that weight has already been done at each site so we just sum the gradients
    Hessian <- grad_exstract$Hessian
    
    #Negative inverse hessian for the update
    H_inv <- solve(-Hessian)# The negative sign is correct, we forgot to take the negative when making the hessian. 
    
    #Updating Beta
    beta_new <-  beta_old-H_inv%*%GRADIENT
    #if no large change then stop algorithm
    if ( norm(beta_new-beta_old,type="2") <0.0000000000001) {
      print(k)
      
      #Returns the numbers of iterations for each site
      if (Computation) {
        return(k)
      }
      
      break
    }
    print(beta_new)
    beta_old <- beta_new
  }
  return(beta_old)
}

#Function to make a prediction for the Bendix carstens model after newton raphson FL
pred_FL_AFT_fkt <- function(df,beta_old,predtimes=0:100,t.kn){
  
  #Prediction time points
  tfe <- predtimes
  
  #Numbers of patients
  n_cut <- nrow(df)
  
  #Making the design matrix
  x <- model.matrix(~1+Ns(tfe, knots = t.kn)+X1+X2+X3+X4, merge(data.frame(tfe=tfe),df ))
  
  #Predicting, exp(-exp(lambda)*t) note: t is replicated n_cut to fit the predictions
  pred <- exp(-exp(x%*%beta_old)*rep(tfe,n_cut) )
  
  #A sequence to split up for the predictions in order to align the right predictions in the prediciton vecotr as its done for alle patients at once
  split_length <- 1:length(pred)
  
  #Splitting up predictied survival for
  splits <- split(pred, cut_number(split_length,n = n_cut))
  
  #Making a survival prediction matrix, each row is predictied survival 
  Surv_pred <- do.call("cbind",splits) %>% t()
  
  return(Surv_pred)
}#This is an old one
pred_FL_COX_fkt <- function(Sim_data,test_df,beta_old,predtimes,...){
  
  #Making a weighting argument
  Total_nr_patinets <- nrow(do.call("rbind", Sim_data))

    
  Local_BZ_cox <- lapply(X=Sim_data, function(X){
    fit <-  coxph(Surv(eventtime, status) ~ X1+X2+X3+X4, data = X, init = c(beta_old),iter.max=0,x=T,y=T)
    #Weighting
    Weight <- (nrow(X)/Total_nr_patinets)
    
    BZKernel_fkt(cox_obj =fit,
                 pred_seq = predtimes,df=X,weight=Weight,...)
  })
  
  
  BZ_test <- BZ_kernel_FL(BZ_list =Local_BZ_cox,
                          pred_seq=predtimes )
  
  QMatrix <- Cox_FL_survival_prediction(BZ_kernel_FL = BZ_test,
                                        FL_Beta = beta_old,
                                        testData =test_df )
  
}
Pred_FL_Poisson_fkt <- function (beta_old,df,predtimes = seq(0,1000,1),t.kn){
  #Making time vairable for the whole dataset
  dff <- merge(data.frame(tfe=predtimes),df )
  #Coefficients
  
  
  #Design matrix
  X <- model.matrix(~1+Ns(tfe, knots = t.kn)+X1+X2+X3+X4, dff)
  
  #Calculating X*Beta
  #ct <- X %*% cf
  ct <- X %*% beta_old
  
  #Calculating how many subjects we have for cutting XB
  n_cut <- nrow(df)
  split_length <- 1:(length(predtimes)*n_cut)
  
  #Finding the splits
  subject_index_list <- split(split_length, cut_number(split_length,n = n_cut))
  
  #taking exp to XB
  ct <- exp(ct)
  #lower triangle matrix filled only with 1's. This is for matrix vector multiplication and making a "time varying" prediction
  cum.mat <- lower.tri(diag(predtimes),diag = T)
  #Making predictions for each pts
  pred_list <-  lapply(1:n_cut,FUN=function(i){
    
    pred <- c(exp(- (cum.mat %*% ct[subject_index_list[[i]]])))
    pred <- c(1,pred[1:(length(pred)-1)])
    #pred <- c(exp(-pred))
    #ct <- data.frame(Pred=exp(-ct) ,Time=ctr.mat$tfe)
    return(pred)
  } )
  #Combinning predictions into a matrix
  pred_mat <- do.call(what = "rbind",pred_list)
  #Returning a prediction matrix, rows are predictions for each subject so coloum 1 is predictions for time 1 in predtimes.
  return(pred_mat)
  
}


Hazard_fkt <- function(i,df,X_data,beta,risk_mat,uniq_event_times){
  
  X_D <- df[(df$eventtime==i & df$status==1),3:ncol(df)]
  d_i <- nrow(X_D)
  
  #Making the risk set 
  X_riskSet <- as.matrix(X_data[risk_mat[,which(uniq_event_times==i)],])
  
  EXP_XB <- exp(X_riskSet %*% beta)
  return(d_i/sum(EXP_XB))
  
}
kdgaussian_fkt <- function (x = 0, lambda = 4, bw = NULL){
  z = (x)/lambda
  
  dnorm(z)/lambda
}
K_indput <- function(t,T_i,b){
  (t-T_i)/b
}
Baseline_hazard_fkt <- function(df,beta,times = "eventtime",event = "status"){
  
  #Design matrix
  X_data <- df[,c("X1","X2","X3","X4")]
  
  #All unique eventtimes
  uniq_event_times <- sort(unique(df[,times]))
  
  #Entry vector
  
  #Making a risk matrix, for all riskSets at each unique eventtime
  RISK_SET_mat <- sapply(uniq_event_times, function(toi) df[,times]> 0& df[,times] >= toi)
  
  hazard <- lapply(uniq_event_times,FUN = function(i){
    Hazard_fkt(
      i = i ,
      df = df ,
      X_data = X_data ,
      beta = beta,
      risk_mat = RISK_SET_mat,
      uniq_event_times = uniq_event_times)
  })
  #hazard <- data.table::data.table(times=uniq_event_times,Hazard=do.call(what="c",hazard))
  
  return(list(hazard,uniq_event_times) )
}
BZKernel_fkt <- function(cox_obj,b=4,pred_seq,df,weight=1){
  
  #exstracting the baseline hazard 
  bz <- predictCox(object = cox_obj,type = "hazard",centered = F,keep.times = T)
  BZ_new <- Baseline_hazard_fkt(df = df,beta = coef(cox_obj)) #Note: disse giver det samme output, nemlig at lave risk-settet i kernel estimatoren
  #browser()
  #We only look at up until the last event time, otherwise it becomes 0 for a loong time
  max_time <- max(df[df$status==1,"eventtime"])
  
  #Applying the kernel function
  #Kernel_function <-lapply(pred_seq,FUN=function(x){evmix::kdgaussian(K_indput(t=x,T_i = bz$time,b=b),bw = b)})
  Kernel_function <-lapply(pred_seq,FUN=function(x){kdgaussian_fkt(K_indput(t=x,T_i = bz$time,b=b),bw = b)})
  K <- do.call("cbind", Kernel_function)
  
  
  #Allowing for multi processing
  #K %<-% Kernel_function(pred_seq = pred_seq)
  #K <- do.call("cbind", K)
  
  #calculating the new smooth baseline hazard function
  h0 <- (1/b)*t(K)%*%bz$hazard
  h0_df <- data.frame(smoothhazard=h0*weight,time=pred_seq) %>% filter(smoothhazard>.Machine$double.eps)
  
  return(h0_df)
}

BZ_kernel_FL <- function(BZ_list,pred_seq){
  
  #Combine all baseline hazards and arrange by tuime
  BZdf <- do.call("rbind",BZ_list)%>% arrange(time) #%>% as.data.frame()
  
  #Making sure no value is equal to 0 other wise we are in trouble 
  BZdf[BZdf$smoothhazard==0,"smoothhazard"] <- .Machine$double.eps
  #BZdf_avg <- BZdf %>% group_by(time) %>% summarise(smoothhazard_avg=mean(smoothhazard))
  BZdf_sum <- BZdf %>% group_by(time) %>% summarise(smoothhazard_sum=sum(smoothhazard))
  fit_loess <- loess(log(smoothhazard_sum) ~time,data = BZdf_sum,span = 0.5)
  
  #fit_loess <- loess(log(smoothhazard_avg) ~time,data = BZdf_avg,span = 0.5)
  
  #Smooth out all hazards with loess function, on log-scale (so we dont end up with negative predictions later on)
  #fit_loess <- loess(log(smoothhazard) ~time,data = BZdf,span = 0.5)
  
  h0_smooth <- exp(predict(fit_loess,data.frame(time=pred_seq)))
  
  #Function to integrate over, converted back using the exponent
  # h0_smooth <- function(t){
  #   exp(predict(fit_loess,newdata = data.frame(time=t)))
  # }
  
  h0_smooth <- data.frame(h0_smooth=h0_smooth,time=pred_seq)
  
  H0_smoot_fkt <- function(i){
    h0_smooth$h0_smooth[i+1]
  }
  
  integrate_h0 <- lapply(pred_seq, function(i){integrate(H0_smoot_fkt,lower = 0,upper = i,subdivisions = 1000000,stop.on.error = F)$value})
  FL_cumHazard <- data.frame(CumHazard=do.call("c",integrate_h0),
                             time=pred_seq)
  
  return(FL_cumHazard)
  
}
Cox_FL_survival_prediction <- function(BZ_kernel_FL,FL_Beta,testData){
  
  #Design matrix
  X <- testData %>% select(X1,X2,X3,X4) %>% as.matrix
  
  #X*beta
  FL_XB <- X%*%matrix(FL_Beta,ncol=1)
  
  
  
  #Survival for each patients at each point
  FL_surv %<-% lapply(X=BZ_kernel_FL$CumHazard,
                      FUN = function(X){
                        exp(-X)^exp(FL_XB)
                      })
  
  FL_surv <- do.call(what="cbind",FL_surv)
  return(FL_surv)
}


#The ture model
TRUE_model_pred <- function(df,lambda,gamma,True_beta,ayear,pred_time,test){

  #Fit weibull model
  mod_weib <- flexsurvspline(Surv(eventtime_true, status) ~ X1+X2+X3+X4,
                             data = df,
                             k = 0,
                             scale="hazard")
  
  #Updating parameters to the true parameters
  mod_weib[["res"]][,"est"] <- c(log(lambda),gamma,True_beta )
  
  #Making predictions
  pred <- summary(object = mod_weib,
                  newdata = test,
                  type = "survival",
                  ci = F, 
                  t = pred_time/ayear,
                  tidy = F,
                  start = 0, 
                  conf.int = FALSE, conf.level = 0.95,
                  se.fit = FALSE,se = F,B = 0)
  
  #Exstracting predicted probabilites
  pred <- lapply(X=pred,FUN = function(x){
    x[,2]
  })
  pred <- do.call(what = "rbind",pred)
  rownames(pred) <- NULL
  return(pred)
}

