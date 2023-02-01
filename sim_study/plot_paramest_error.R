#' Just plotting the errors from the error assessment scheme 
library(tidyverse)
library(gridExtra)
suffix <- "" #this is an artifact of a previous version. We keep it in case
             #we need to use it for running legacy code (probably team-intern uses) 
source(paste0("param_set",suffix,".R"))
#' General parameters as for the simulations

#' To set by hand
s<- 100#100
coal_index <- 2 #beta = 1 psi = 2
param_to_write <- c("coalp","r","e")
#width_fp <- 40


#' Subsequent names/params
coal_name2 <- c("beta","psi")[coal_index]
load(paste0("results4plotstables/estimR2err_",coal_name2,"_s",s,
            suffix,".RData"))
params1 <- list(params_beta_exp,params_dirac_exp)[[coal_index]]
grid1 <- est_R2_list
coal_name <- c("alpha","psi")[coal_index]
name1 <- paste0("_",coal_name,"_",param_to_write)
coalp <- list(coal_beta,coal_dirac)[[coal_index]]  


#' Parameter to produce FALSE=plots, TRUE=L2 distances 
#stats_out <- FALSE#TRUE
#par(mfrow=c(2,2))
#if (stats_out){res1 <- list(paste0("L2err",name1),NULL)}

#pdf(paste0("err_r_coal_n",s,name1,".pdf"))

#if(!(stats_out)){pdf(paste0("err_n",s,"_",name1,".pdf"))}
#' Order values into tidy object
data1 <- NULL

params_temp <- as_tibble(params1)
params2 <- mutate(params_temp,e=misiden[1],
                  .keep = "all")
for (i in 2:length(misiden)){
params2 <- bind_rows(params2,
                     mutate(params_temp,e=misiden[i],
                     .keep = "all"))}
colnames(params2)[2] <- "coalp"

for (i in seq(along=est_R2_list)){
temp1 <- as_tibble(t(est_R2_list[[i]]))
colnames(temp1)[2:3] <- paste0("est_",c("r","coalp"))
temp1 <- temp1 %>% mutate(params2[i,]) 
data1 <- bind_rows(data1,temp1)
}

data1 <- data1 %>% mutate(coal_e = est_coalp - coalp,
                          r_e = est_r - r,
                          e_e = est_misiden - e)

#Plot overall errors for alpha, r, e
#Now in separate plots
#Make text sizes bigger
theme_set(theme_gray(base_size = 14))

#alpha
if (coal_index==1){
data_bp <- bind_cols(param=rep(coal_name,
                      nrow(data1)),
                     data1[,9])
colnames(data_bp)[2] <- "val" 

plot1 <- ggplot(data_bp, aes(param,val)) + 
  geom_boxplot(orientation = "x") + 
  labs( x = expression(alpha),
    y = expression(widehat(alpha)-alpha)) +
  theme(axis.ticks.x=element_blank(),
         axis.text.x=element_blank())

plot1a <- ggplot(data_bp, aes(param,val)) + 
  geom_count(aes(size = after_stat(prop))) + 
  labs( x = expression(alpha),
        y = expression(widehat(alpha)-alpha)) +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank())
#ggsave(paste0("pics/err_n",s,"_",coal_name2,"_bp_coalp.pdf"), 
#       width = 10, height = 10)
}

#PSI
if (coal_index==2){
  data_bp <- bind_cols(param=rep(coal_name,
                                 nrow(data1)),
                       data1[,9])
  colnames(data_bp)[2] <- "val" 
  
plot1 <- ggplot(data_bp, aes(param,val)) + 
        geom_boxplot(orientation = "x") + 
    labs( x = expression(Psi),
          y = expression(widehat(Psi)-Psi)) +
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank())

plot1a <- ggplot(data_bp, aes(param,val)) + 
  geom_count(aes(size = after_stat(prop))) + 
  labs( x = expression(Psi),
        y = expression(widehat(Psi)-Psi)) +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank())
#  ggsave(paste0("pics/err_n",s,"_",coal_name2,"_bp_coalp.pdf"), 
#         width = 10, height = 10)
}

#r
data_bp <- bind_cols(param=rep("rho",
                               nrow(data1)),
                     data1[,10])
colnames(data_bp)[2] <- "val" 

plot2 <- ggplot(data_bp, aes(param,val)) + 
  geom_boxplot(orientation = "x") + 
  labs( x = expression(g),
        y = expression(hat(g)-g)) +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank())

plot2a <- ggplot(data_bp, aes(param,val)) + 
  geom_count(aes(size = after_stat(prop))) + 
  labs( x = expression(g),
        y = expression(hat(g)-g)) +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank())
#ggsave(paste0("pics/err_n",s,"_",coal_name2,"_bp_r.pdf"), 
#       width = 10, height = 10)

#e
data_bp <- bind_cols(param=rep("e",
                               nrow(data1)),
                     data1[,11])
colnames(data_bp)[2] <- "val" 

plot3 <- ggplot(data_bp, aes(param,val)) + 
  geom_boxplot(orientation = "x") + 
  labs( x = expression(e),
        y = expression(widehat(e)-e)) +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank())

plot3a <- ggplot(data_bp, aes(param,val)) + 
  geom_count(aes(size = after_stat(prop))) + 
  labs( x = expression(e),
        y = expression(widehat(e)-e)) +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank())

g1 <- arrangeGrob(plot1,plot2,plot3, ncol=3, nrow = 1)
#ggsave(paste0("pics/err_n",s,"_",coal_name2,"_bp",
#              suffix,".pdf"),
#       plot = g1,  
#       width = 10, height = 10)

g2 <- arrangeGrob(plot1a,plot2a,plot3a, ncol=3, nrow = 1)
ggsave(paste0("pics/err_n",s,"_",coal_name2,"_cp",
              suffix,".pdf"),
       plot = g2,
       width = 10, height = 6)

#Reset size globally
theme_set(theme_gray(base_size = 11))

#Errors for each parameter combo
#Plot hat(alpha)-alpha
if (coal_index==1){
ggplot(data1, aes(factor(coalp),coal_e)) + 
  geom_count(aes(size = after_stat(prop))) +
  scale_size_area(max_size = 5) +
  facet_grid(rows = vars(r), cols = vars(e),
             labeller = "label_both") + 
  labs( x = expression(alpha), 
        y = expression(hat(alpha)-alpha)) +
  theme_minimal()
  } else {
    ggplot(data1, aes(factor(coalp),coal_e)) + 
      geom_count(aes(size = after_stat(prop))) +
      scale_size_area(max_size = 5) +
      facet_grid(rows = vars(r), cols = vars(e),
                 labeller = "label_both") + 
      labs( x = expression(Psi), 
            y = expression(hat(Psi)-Psi)) +
      theme_minimal()
    }
ggsave(paste0("pics/err_n",s,name1,suffix,".pdf")[1], 
       width = switch(coal_index,22,22),height = 10)

#Plot hat(r)-r
if (coal_index==1){
  ggplot(data1, aes(factor(coalp),r_e)) + 
    geom_count(aes(size = after_stat(prop))) +
    scale_size_area(max_size = 5) +
    facet_grid(rows = vars(r), cols = vars(e),
               labeller = "label_both") + 
    labs( x = expression(alpha), 
          y = expression(hat(r)-r)) +
    theme_minimal()
} else {
  ggplot(data1, aes(factor(coalp),r_e)) + 
    geom_count(aes(size = after_stat(prop))) +
    scale_size_area(max_size = 5) +
    facet_grid(rows = vars(r), cols = vars(e),
               labeller = "label_both") + 
    labs( x = expression(Psi), 
          y = expression(hat(r)-r)) +
    theme_minimal()
}
ggsave(paste0("pics/err_n",s,name1,suffix,".pdf")[2], 
       width = switch(coal_index,22,22), height = 10)

#Plot hat(e)-e

if (coal_index==1){
  ggplot(data1, aes(factor(coalp),e_e)) + 
    geom_count(aes(size = after_stat(prop))) +
    scale_size_area(max_size = 5) +
    facet_grid(rows = vars(r), cols = vars(e),
               labeller = "label_both") + 
    labs( x = expression(alpha), 
          y = expression(hat(e)-e)) +
    theme_minimal()
} else {
  ggplot(data1, aes(factor(coalp),e_e)) + 
    geom_count(aes(size = after_stat(prop))) +
    scale_size_area(max_size = 5) +
    facet_grid(rows = vars(r), cols = vars(e),
               labeller = "label_both") + 
    labs( x = expression(Psi), 
          y = expression(hat(e)-e)) +
    theme_minimal()
}
ggsave(paste0("pics/err_n",s,name1,suffix,".pdf")[3], 
       width = switch(coal_index,22,22),height = 10)