##########################
### Simulation Results ###
##########################
#This scripts is for analysis the simulation results

################
### Packages ###
################
library(reshape2)
library(ggplot2)
library(tidyverse)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(tableone)
######################
### Specifications ###
######################
ayear <- 365.25
height <-  2160
width <- 3840
res <-  300


#################
### Functions ###
#################---------------------------------------------------------------
C_index_ceaning_fkt <- function(df){
  C_index_df <- df %>%
    select(C_index_FL_pois,C_index_Pool_pois,C_index_FL_COX,C_index_Pool_COX,C_index_TRUE) %>%
    melt(measure.vars=c('C_index_FL_pois',
                        'C_index_Pool_pois',
                        "C_index_FL_COX",
                        "C_index_Pool_COX","C_index_TRUE")) %>% 
    mutate(
      Models=case_when(
        variable=="C_index_FL_pois"~"FL-Poisson",
        variable=="C_index_Pool_pois"~"Poisson-Pool",
        variable=="C_index_FL_COX"~"FL-Cox",
        variable=="C_index_Pool_COX"~"Cox-Pool",
        variable=="C_index_TRUE"~"True"
        ) %>% factor(levels=c("FL-Poisson","Poisson-Pool","FL-Cox","Cox-Pool","True")) )
  return(C_index_df)
}
IBS_df_cleaning_fkt <- function(df){
  IBS_DF <- df %>%
    select(IBS_FL_pois  ,IBS_pool_pois,IBS_FL_Cox,IBS_Pool_COX,IBS_true_model) %>% 
    melt(measure.vars=c('IBS_FL_pois',
                        'IBS_pool_pois',
                        "IBS_FL_Cox",
                        "IBS_Pool_COX","IBS_true_model")) %>% 
    mutate(
      Models=case_when(
        variable=="IBS_FL_pois"~"FL-Poisson",
        variable=="IBS_pool_pois"~"Poisson-Pool",
        variable=="IBS_FL_Cox"~"FL-Cox",
        variable=="IBS_Pool_COX"~"Cox-Pool",
        variable=="IBS_true_model"~"True"
      ) %>% factor(levels=c("FL-Poisson","Poisson-Pool","FL-Cox","Cox-Pool","True"))) 
  return(IBS_DF) 
}
ABS_AVG_SD_fkt <- function(df){
  ABS_AVG_survDIFF_df <- df %>%
    select(Avg_Surv_diff_FL_POIS,Avg_Surv_diff_Pool_Pois,Avg_Surv_diff_FL_Cox,Avg_Surv_diff_Pool_cox ) %>% 
    melt(measure.vars=c('Avg_Surv_diff_FL_POIS',"Avg_Surv_diff_Pool_Pois","Avg_Surv_diff_FL_Cox", 'Avg_Surv_diff_Pool_cox') ) %>% 
    mutate(Models=case_when(
      variable=="Avg_Surv_diff_FL_POIS"~"FL-Poisson",
      variable=="Avg_Surv_diff_Pool_Pois"~"Poisson-Pool",
      variable=="Avg_Surv_diff_FL_Cox"~"FL-Cox",
      variable=="Avg_Surv_diff_Pool_cox"~"Cox-Pool"
      
      ) %>% factor(levels=c("FL-Poisson","Poisson-Pool","FL-Cox","Cox-Pool","True")))
  return(ABS_AVG_survDIFF_df)
}

Coefficient_analyse_cleaning_fkt <- function(sim11_test_list,substract_beta="TRUE_beta"){
  
  pois_fl_beta   <- do.call("cbind", sim11_test_list[,1])[6:9,]
  cox_fl_beta    <- do.call("cbind", sim11_test_list[,2]) 
  pois_pool_beta <- do.call("cbind", sim11_test_list[,3])[6:9,]
  cox_pool_beta  <- do.call("cbind", sim11_test_list[,4]) 
  
  
  
  if (substract_beta=="TRUE_beta") {
    poisson_fl_diff   <- as.data.frame(pois_fl_beta-True_beta)
    poisson_pool_diff <- as.data.frame(pois_pool_beta-True_beta)
    cox_fl_diff       <- as.data.frame(cox_fl_beta-True_beta)
    cox_pool_diff     <- as.data.frame(cox_pool_beta-True_beta)
  }
  if (substract_beta=="Pool") {
    poisson_diff   <- as.data.frame(pois_fl_beta-pois_pool_beta)
    Cox_diff <- as.data.frame(cox_fl_beta-cox_pool_beta)
    
    poisson_diff$Beta    <- c("Beta-1","Beta-2","Beta-3","Beta-4")
    Cox_diff$Beta  <- c("Beta-1","Beta-2","Beta-3","Beta-4")
    
    poisson_diff$Model <- "Poisson FL vs Pool"
    Cox_diff$Model <- "Cox FL vs Pool"
    
    poisson_diff <- poisson_diff %>% melt()
    Cox_diff <- Cox_diff %>% melt()
    
    Beta_diff_list <- list(poisson_diff=poisson_diff,
                           Cox_diff=Cox_diff)
    return(Beta_diff_list)
  }
  
  
  poisson_fl_diff$Beta    <- c("Beta-1","Beta-2","Beta-3","Beta-4")
  poisson_pool_diff$Beta  <- c("Beta-1","Beta-2","Beta-3","Beta-4")
  cox_fl_diff$Beta        <- c("Beta-1","Beta-2","Beta-3","Beta-4")
  cox_pool_diff$Beta      <- c("Beta-1","Beta-2","Beta-3","Beta-4")
  
  poisson_fl_diff$Model <- "FL-Poisson"
  poisson_pool_diff$Model <- "Poisson-Pool"
  cox_fl_diff$Model <- "FL-Cox"
  cox_pool_diff$Model <- "Cox-Pool"
  
  
  poisson_fl_diff <- poisson_fl_diff %>% melt()
  poisson_pool_diff <- poisson_pool_diff %>% melt()
  cox_fl_diff <- cox_fl_diff %>% melt()
  cox_pool_diff <- cox_pool_diff %>% melt()
  
  Beta_diff_list <- list(poisson_fl_diff=poisson_fl_diff,
                         poisson_pool_diff=poisson_pool_diff,
                         cox_fl_diff=cox_fl_diff,
                         cox_pool_diff=cox_pool_diff
  )
  
  return(Beta_diff_list)
}




Cindex_box_plot_fkt <- function(C_index_res,title="",x_title="",exclude_y_axis=F,ylim=c(0.6,0.8)){
p <-   ggplot(C_index_res,aes(y = value,fill=Models))+
    geom_boxplot()+
    ylim(ylim)+
    labs(y="C-index",
         title =title, 
         x=x_title)+
    theme_classic() +
    scale_fill_jco()+
    theme(legend.position = "bottom",
          plot.title = element_text(size=20,vjust = 0,hjust = 0.5),
          legend.text = element_text(size=12),
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=20)
    )
  if (exclude_y_axis) {
    p <- p+theme(axis.title.y =  element_blank())
    return(p)
  }else{
    return(p) 
  }
  
}
IBS_box_plot_fkt <- function(df,title="",x_title="",exclude_y_axis=F,ylim=c(0.1,0.15)){
 p <-  ggplot(df,aes(y = value,fill=Models))+
    geom_boxplot()+
    labs(y="IBS",
         title =title, 
         x=x_title)+
    theme_classic() +
    scale_fill_jco()+
    theme(legend.position = "bottom",
          plot.title = element_text(size=20,vjust = 0,hjust = 0.5),
          legend.text = element_text(size=12),
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=20)
    )+ylim(ylim)   
 if (exclude_y_axis) {
   p <- p+theme(axis.title.y =  element_blank())
   return(p)
 }else{
   return(p) 
 }
}
ABS_surv_df_box_plot_fkt <- function(df,title="",x_title="",exclude_y_axis=F){
 p <-  ggplot(df,aes(y = value,fill=Models))+
    geom_boxplot()+
    ylim(c(0,100))+
    labs(y="MIAD",
         title =title, 
         x=x_title)+
    theme_classic() +
    scale_fill_jco()+
    theme(legend.position = "bottom",
          plot.title = element_text(size=20,vjust = 0,hjust = 0.5),
          legend.text = element_text(size=12),
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=20)
    )   
  if (exclude_y_axis) {
    p <- p+theme(axis.title.y =  element_blank())
    return(p)
  }else{
    return(p) 
  }
  
}
Coef_Box_plot_fkt <- function(Coef_diff_Df,title="",x_title="",exclude_y_axis=F,ylim=c(0,0.2)){

 p <-  ggplot(Coef_diff_Df)+
    geom_boxplot(aes(y = value,fill=Beta))+
    ylim(ylim)+
    labs(y="Difference",
         title =title, 
         x=x_title)+
    theme_classic() +
    scale_fill_jco()+
    theme(legend.position = "bottom",
          plot.title = element_text(size=20,vjust = 0,hjust = 0.5),
          legend.text = element_text(size=12),
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=20)
    )   
  if (exclude_y_axis) {
    p <- p+theme(axis.title.y =  element_blank())
    return(p)
  }else{
    return(p) 
  }
}

####################
### Load results ###
####################------------------------------------------------------------
#sim1_test <- readRDS("./Simulation 1/sim1_test_new3.rds")
#sim2_test <- readRDS("./Simulation 2/sim2_test_new.rds")
#sim3_test <- readRDS("./Simulation 3/sim3_test_new.rds")
#sim4_test <- readRDS("./Simulation 4/sim4_test_new.rds")
#sim5_test <- readRDS("./Simulation 5/sim5_test_new.rds")

sim1_test <- readRDS("./Simulation 1/NewSimulations/sim1_test_results.rds")
sim2_test <- readRDS("./Simulation 2/NewSimulations/sim2_test_results.rds")
sim3_test <- readRDS("./Simulation 3/NewSimulations/sim3_test_new.rds")
sim4_test <- readRDS("./Simulation 4/NewSimulations/sim4_test_new.rds")
#sim5_test <- readRDS("./Simulation 5/sim5_test_new.rds")



#######################
### Result Cleaning ###
#######################
#Cleaning c-index for all models and simulations
C_index_sim1 <- C_index_ceaning_fkt(df = sim1_test)
C_index_sim2 <- C_index_ceaning_fkt(df = sim2_test)
C_index_sim3 <- C_index_ceaning_fkt(df = sim3_test)
C_index_sim4 <- C_index_ceaning_fkt(df = sim4_test)
#C_index_sim5 <- C_index_ceaning_fkt(df = sim5_test)


#Cleaning integrated brier score IBS for all simulations and models
IBS_sim1 <- IBS_df_cleaning_fkt(df = sim1_test)
IBS_sim2 <- IBS_df_cleaning_fkt(df = sim2_test)
IBS_sim3 <- IBS_df_cleaning_fkt(df = sim3_test)
IBS_sim4 <- IBS_df_cleaning_fkt(df = sim4_test)
#IBS_sim5 <- IBS_df_cleaning_fkt(df = sim5_test)

#Cleaning the absolute average survival difference between poisson methods and cox methods for each simulaitons
ABS_Surv_diff_sim1 <- ABS_AVG_SD_fkt(df = sim1_test)
ABS_Surv_diff_sim2 <- ABS_AVG_SD_fkt(df = sim2_test)
ABS_Surv_diff_sim3 <- ABS_AVG_SD_fkt(df = sim3_test)
ABS_Surv_diff_sim4 <- ABS_AVG_SD_fkt(df = sim4_test)
#ABS_Surv_diff_sim5 <- ABS_AVG_SD_fkt(df = sim5_test)


###############
### Results ###
###############-----------------------------------------------------------------

#######################
### C-index results ###
#######################---------------------------------------------------------

#Plotting
C_index_S1 <- Cindex_box_plot_fkt(C_index_res =C_index_sim1,title="Simulation 1")
C_index_S2 <- Cindex_box_plot_fkt(C_index_res =C_index_sim2,title="Simulation 2",exclude_y_axis = T )
C_index_S3 <- Cindex_box_plot_fkt(C_index_res =C_index_sim3,title="Simulation 3",exclude_y_axis = T )
C_index_S4 <- Cindex_box_plot_fkt(C_index_res =C_index_sim4,title="Simulation 4",exclude_y_axis = T )
#C_index_S5 <- Cindex_box_plot_fkt(C_index_res =C_index_sim5,title="Simulation 5",exclude_y_axis = T,ylim=c(0.8,1) )

C_index_plot <- ggarrange(C_index_S1,
                          C_index_S2,
                          C_index_S3,
                          C_index_S4,
                          ncol=4,nrow=1,common.legend = T,legend = "bottom")



###################
### IBS Results ###
###################-------------------------------------------------------------
IBS_S1 <- IBS_box_plot_fkt(IBS_sim1,title="Simulation 1"  )
IBS_S2 <- IBS_box_plot_fkt(IBS_sim2,title="Simulation 2",exclude_y_axis = T )
IBS_S3 <- IBS_box_plot_fkt(IBS_sim3,title="Simulation 3",exclude_y_axis = T  )
IBS_S4 <- IBS_box_plot_fkt(IBS_sim4,title="Simulation 4",exclude_y_axis = T  )
#IBS_S5 <- IBS_box_plot_fkt(IBS_sim5,title="Simulation 5",exclude_y_axis = T ,ylim = c(0,0.15) )
IBS_plot <- ggarrange(
  IBS_S1,
  IBS_S2,
  IBS_S3,
  IBS_S4,
  ncol = 4,
  nrow = 1,
  common.legend = T,
  legend = "bottom"
)


#################################################
### Abs average survival difference - Results ###
#################################################-------------------------------
ABS_Avg_s1 <- ABS_surv_df_box_plot_fkt(df =ABS_Surv_diff_sim1,title="Simulation 1",exclude_y_axis =  F )
ABS_Avg_s2 <- ABS_surv_df_box_plot_fkt(df =ABS_Surv_diff_sim2,title="Simulation 2",exclude_y_axis = T )
ABS_Avg_s3 <- ABS_surv_df_box_plot_fkt(df =ABS_Surv_diff_sim3,title="Simulation 3",exclude_y_axis = T )
ABS_Avg_s4 <- ABS_surv_df_box_plot_fkt(df =ABS_Surv_diff_sim4,title="Simulation 4",exclude_y_axis = T )
#ABS_Avg_s5 <- ABS_surv_df_box_plot_fkt(df =ABS_Surv_diff_sim5,title="Simulation 5",exclude_y_axis = T )

ABS_Avg_plot <- ggarrange(
  ABS_Avg_s1,
  ABS_Avg_s2,
  ABS_Avg_s3,
  ABS_Avg_s4,
  ncol = 4,
  nrow = 1,
  common.legend = T,
  legend = "bottom"
)

#Median values
ABS_Avg_s1$data %>% group_by(Models) %>% summarise(med=median(value) )
ABS_Avg_s2$data %>% group_by(Models) %>% summarise(med=median(value) )
ABS_Avg_s3$data %>% group_by(Models) %>% summarise(med=median(value) )
ABS_Avg_s4$data %>% group_by(Models) %>% summarise(med=median(value) )
#ABS_Avg_s5$data %>% group_by(Models) %>% summarise(med=median(value) )


####################
### Saving Plots ###
####################------------------------------------------------------------
png(
  filename = "Q:/AalbUH-Haema-AnnonymisationOfPredictiveModels/Results/Figures/C_index_results.png",
  width = width,
  height = height,
  res = res
)
C_index_plot
dev.off()

png(
  filename = "Q:/AalbUH-Haema-AnnonymisationOfPredictiveModels/Results/Figures/IBS_results.png",
  width = width,
  height = height,
  res = res
)
IBS_plot
dev.off()

png(
  filename = "Q:/AalbUH-Haema-AnnonymisationOfPredictiveModels/Results/Figures/ABS_AVG_results.png",
  width = width,
  height = height,
  res = res
)
ABS_Avg_plot
dev.off()




##############################
### Table with differences ###
##############################
RES1 <- sim1_test %>%
  mutate(C_diff_FL_POIS=C_index_FL_pois-C_index_TRUE,
         C_diff_Pool_pois=C_index_Pool_pois- C_index_TRUE,
         C_diff_FL_Cox=C_index_FL_COX - C_index_TRUE,
         C_diff_Pool_cox= C_index_Pool_COX - C_index_TRUE,
         
         IBS_diff_FL_POIS=IBS_FL_pois -IBS_true_model,
         IBS_diff_Pool_Pois=IBS_pool_pois -IBS_true_model,
         IBS_diff_FL_Cox=IBS_FL_Cox-IBS_true_model,
         IBS_diff_Pool_Cox=IBS_Pool_COX-IBS_true_model,
         Simulation="Simulation 1"
  )
RES2 <- sim2_test %>%
  mutate(C_diff_FL_POIS=C_index_FL_pois-C_index_TRUE,
         C_diff_Pool_pois=C_index_Pool_pois- C_index_TRUE,
         C_diff_FL_Cox=C_index_FL_COX - C_index_TRUE,
         C_diff_Pool_cox= C_index_Pool_COX - C_index_TRUE,
         
         IBS_diff_FL_POIS=IBS_FL_pois -IBS_true_model,
         IBS_diff_Pool_Pois=IBS_pool_pois -IBS_true_model,
         IBS_diff_FL_Cox=IBS_FL_Cox-IBS_true_model,
         IBS_diff_Pool_Cox=IBS_Pool_COX-IBS_true_model,
         Simulation="Simulation 2"
  )
RES3 <- sim3_test %>%
  mutate(C_diff_FL_POIS=C_index_FL_pois-C_index_TRUE,
         C_diff_Pool_pois=C_index_Pool_pois- C_index_TRUE,
         C_diff_FL_Cox=C_index_FL_COX - C_index_TRUE,
         C_diff_Pool_cox= C_index_Pool_COX - C_index_TRUE,
         
         IBS_diff_FL_POIS=IBS_FL_pois -IBS_true_model,
         IBS_diff_Pool_Pois=IBS_pool_pois -IBS_true_model,
         IBS_diff_FL_Cox=IBS_FL_Cox-IBS_true_model,
         IBS_diff_Pool_Cox=IBS_Pool_COX-IBS_true_model,
         Simulation="Simulation 3"
  )
RES4 <- sim4_test %>%
  mutate(C_diff_FL_POIS=C_index_FL_pois-C_index_TRUE,
         C_diff_Pool_pois=C_index_Pool_pois- C_index_TRUE,
         C_diff_FL_Cox=C_index_FL_COX - C_index_TRUE,
         C_diff_Pool_cox= C_index_Pool_COX - C_index_TRUE,
         
         IBS_diff_FL_POIS=IBS_FL_pois -IBS_true_model,
         IBS_diff_Pool_Pois=IBS_pool_pois -IBS_true_model,
         IBS_diff_FL_Cox=IBS_FL_Cox-IBS_true_model,
         IBS_diff_Pool_Cox=IBS_Pool_COX-IBS_true_model,
         Simulation="Simulation 4"
  )

df_combined <- rbind(RES1,RES2,RES3,RES4) %>% as.data.frame()

vars <- c("C_index_TRUE","C_index_Pool_pois","C_index_FL_COX","C_index_Pool_COX","C_index_FL_pois",
          "C_diff_FL_POIS","C_diff_Pool_pois","C_diff_FL_Cox","C_diff_Pool_cox",
          "IBS_true_model","IBS_FL_pois","IBS_pool_pois","IBS_FL_Cox","IBS_Pool_COX",
          "IBS_diff_FL_POIS","IBS_diff_Pool_Pois","IBS_diff_FL_Cox","IBS_diff_Pool_Cox",
          "Avg_Surv_diff_FL_POIS","Avg_Surv_diff_Pool_Pois","Avg_Surv_diff_FL_Cox","Avg_Surv_diff_Pool_cox"
          
          
          )


tab1 <- CreateTableOne(vars = vars,data = df_combined,strata = "Simulation")
tab1 <- print(tab1,minMax = T,nonnormal = vars,contDigits = 4,test = F)%>% as.data.frame()

write.csv(x = tab1,file = "Q:/AalbUH-Haema-AnnonymisationOfPredictiveModels/Results/Tables/C_diff.csv",row.names = T)











