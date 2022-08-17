#Set Working Directory 
setwd("")


#Load Packages
library(Matrix)
library(gurobi)

# set global parameters -------------------------------
options(max.print = 99999)
set.seed(594)


# define restrictions to experiment with -------------------------------
# civilian cases in the paper (mission time limitations:50,24,8; variance limitations:1.5,0.5,0.5)
mission_time <- 8 # mission time --> 50, 24, 8 in the paper
variance_restriction <- 0.5 # mission time variance restriction --> 1.5 and 0.5 in the paper

  

# define problem parameters -------------------------------
n_target <- 28 # number of targets excluding origin and destination
n_node <- (n_target+2) # number of nodes including origin and destination
discount_percentage <- 30 # maximum discount in terms of %
epsilon_coef <- 0.1 # epsilon coefficient used in objective function for augmentation
n_revisit <- 3 # number of allowed revisits
coefficient_of_variation <- 0.2 
flight_speed <- 600 # flight speed (km per hour)
sensor_range <- 20 # sensor range
search_time_allocated <- 10/60 # search time allocation to a sub-region
collection_time_allocated <- 5/60 # collection time allocation to a sub-region
beta_f <- 0.9 # risk of being detected
mission_late_probability <- 0.3 #the allowed risk of not completing mission on time 

# define problem sets -------------------------------
target_set <- (2:(n_target+1))
node_set <- (1:n_node)
time_set <- (1:((n_target*n_revisit)+1)) # the number of movements the UAV is allowed is set to the maximum possible value (see the paper)
k_set <- (1:n_revisit) # set of allowed revisits


# read input files -------------------------------
mu_pIj <- read.table("information_probability_mean.txt",header = TRUE) # probability that a target has information -> mu_{pI_j}
reservoir_capacities <- read.table("reservoir_capacities.txt", header = TRUE)
coordinate_matrix <- as.matrix(read.table("coordinates.txt",header = TRUE))
coordinate_matrix <- cbind(node_set,coordinate_matrix); colnames(coordinate_matrix)<-NULL # adjust coordinate_matrix


# calculate the effectiveness of the search and collection sensors -------------------------------
target_edge <- t(reservoir_capacities) # edge lengths of targets
colnames(target_edge) <- NULL; rownames(target_edge) <- NULL
target_area <- (target_edge*2)^2 
number_of_sub_regions <- unlist(lapply(target_edge, function(x){if(x<=30){return(3)}else if(x<=40){return(4)} else {return(5)}}))
search_sensor_effectiveness <- (1-exp(-(flight_speed*sensor_range)/((target_area)/number_of_sub_regions))^search_time_allocated) #effectiveness of search sensor - by Xia's paper--28km is based on Global Hawk
collection_sensor_effectiveness <- 0.5 # fixed (see the manuscript for details)



# generate all possible connections -------------------------------
connection_list <- c()
for (d in 2:(n_node-1)) {
  connection_temp <- cbind(rep(d,( n_node-1)),   (2:(n_node)))
  connection_list <- rbind(connection_list,connection_temp)
}
connections_all <- rbind(cbind(rep(1,(n_node-2)), 2:(n_node-1)),connection_list)


# calculate travel time t_{ij} of connections -------------------------------
# calculate euclidean distances
euclidean_distance <- apply(connections_all,1,function (x) sqrt((coordinate_matrix[x[1],2] - coordinate_matrix[x[2],2])^2 + (coordinate_matrix[x[1],3] - coordinate_matrix[x[2],3])^2))
connections_all <- cbind(connections_all,euclidean_distance)

# update distances with discounted euclidean distances. discount rates linearly increase from 0 to discount_percentage (%30 in the paper)
discount_rates <- unlist(lapply(connections_all[,3],function (x) (discount_percentage/(max(connections_all[,3])-min(connections_all[,3])))*(x - min(connections_all[,3]))))
discount_amounts <- (connections_all[,3] * discount_rates)/100
connections_all[,3] <- (connections_all[,3] - discount_amounts)
# connections_all_DiscRate<-cbind(connections_all[,1:2],discount_rates)

connections_all[,3] <- connections_all[,3]/flight_speed # convert distance to time
tij <- connections_all[,3] # travel time t_{ij}
var_tij <- (connections_all[,3]*coefficient_of_variation)^2 # travel cost variances

connections_all <- cbind(connections_all,var_tij) # connections_all: # 3: travel time, # 4:travel time variance


#---------------------------------------------------------------------------------------------------------------
# calculate information search and collection time t_{j,k} for targets
# generate q_a -  the probability of detecting the information in the ath searched sub-region 
# generate m_{I_j} - the expected value of the information collection from target j
# probability that a target has information -> mu_{pI_j}
mu_pIj <- as.vector(t(mu_pIj))
stdev_pIj <- apply(cbind((1-mu_pIj)/3, mu_pIj/3),1, min)

# probability that information is detected in the ath searched sub-region q_a
# conditional prob of information exists in ath searched sub-region given that there is no information in previously searched (a-1) sub-regions.
e_a <- lapply(number_of_sub_regions, function (n){unlist(lapply((1:n), function(x) (1/(n-(x-1)))))})
e_a_2 <- lapply(e_a, function(n){c(1,(1-n))[1:length(n)]})
q_a <- sapply(1:length(number_of_sub_regions), function(n){
  sapply(1:number_of_sub_regions[n], function (x){
    row_x<-c(1)
    for(i in (1:x)){row_x=row_x * e_a_2[[n]][i]}
    row_x=row_x * e_a[[n]][x] * search_sensor_effectiveness[n]
    return(row_x)  
  })
})

# mean of the time spent at first visit
mu_tj1 <- sapply(1:length(number_of_sub_regions), function(n){
  sum(unlist(lapply(1:number_of_sub_regions[n], function(a){(a*search_time_allocated + collection_time_allocated)*mu_pIj[n]*q_a[[n]][a]}))) +
    (number_of_sub_regions[[n]]*search_time_allocated*(1-(mu_pIj[n]*sum(q_a[[n]]))))
}, simplify = TRUE)

# variance of the time spent at first visit
mu_tj1_2 <- sapply(1:length(number_of_sub_regions), function(n){
  sum(unlist(lapply(1:number_of_sub_regions[n], function(a){((a*search_time_allocated + collection_time_allocated)^2)*mu_pIj[n]*q_a[[n]][a]})))+ ((number_of_sub_regions[[n]]*search_time_allocated)^2*(1-(mu_pIj[n]*sum(q_a[[n]]))))
}, simplify = TRUE)
var_tj1 <-  mu_tj1_2 - mu_tj1^2

# mean of the time spent at revisits
mu_tjk <- mu_tj1

# variance of the time spent at revisits
var_tjk <- var_tj1


# Construct target attribute matrix---------------------

#organize mean of search/collection time of targets
visit_mu <- mu_tj1

for (i in 2:n_revisit){
  visit_mu <- cbind(visit_mu,mu_tjk)
}

#organize variance of search/collection time of targets
visit_var <- var_tj1

for (i in 2:n_revisit){
  visit_var <- cbind(visit_var,var_tjk)
}

#Target, Visit No, Mean S/C time, Var S/C time
target_attributes_all <- cbind(rep(target_set,each=n_revisit), rep(k_set,n_target),as.vector(t(visit_mu)),as.vector(t(visit_var)))



#---------------------------------------------------------------------------------------------------------------
# expected Information Collection From Targets I_{jk}-> information collected from target j at kth visit

# mean information collection at visits
i_gone <- 0
i_remain <- 1
i_revisists <- c()
for (i in 1:n_revisit){
  i_collect <- (1-i_gone)*collection_sensor_effectiveness
  i_gone <- i_gone+i_collect
  i_remain <- i_remain-i_collect
  i_revisists <- c(i_revisists,i_collect)
}

mu_Ijk <- c()
for (r in 1:n_revisit){
  mu_Ijk <- cbind(mu_Ijk, (mu_pIj * search_sensor_effectiveness * i_revisists[r]))
}

# Target, Visit No, Mean Information Collection, Var Information Collection
target_attributes_new_unord <- cbind(rep(target_set,n_revisit),rep(k_set, each=n_target),as.vector(t(mu_Ijk)))
target_attributes_new_ord <- target_attributes_new_unord[order(target_attributes_new_unord[,1], target_attributes_new_unord[,2]),]
target_attributes_new <- target_attributes_new_ord

#3:Mean S/C time, #4: Var S/C time #5:Mean Information Collection, #6Var Information Collection
target_attributes_all <- cbind(target_attributes_all,target_attributes_new[,3])



#---------------------------------------------------------------------------------------------------------------
#Risk of being detected
x_radar <- cbind(connections_all[,1:2],0)
y_radar <- cbind(rep(target_set,each=n_revisit),0)

#connections_all: #3: Travel Time, #4:Travel Time Variance, #5: Radar Detection Rate
connections_all <- cbind(connections_all,0)

#3:Mean S/C time, #4: Var S/C time #5:Mean Information Collection, #6 Var Information Collection #7Radar detection rate
target_attributes_all <- cbind(target_attributes_all,0)



# civilian cases in the paper (mission time limitations:50,24,8; variance limitations:1.5,0.5,0.5)
time_restriction_set <- c(50,24,8)
variance_restriction_set <- c(1.5,0.5,0.5)



#---------------------------------------------------------------------------------------------------------------
  # START MIP MODEL --------------------------------------
 
  # variable construction --------------------------------------
  # xpij Binary Variables: Variables for Information obtained at location j (excluding destination)

  xpij_matrix <- rbind(
    cbind(1,connections_all[1:(n_node-2),]),
    cbind(rep((2:tail(time_set,2)[1]),each=nrow(connections_all[(n_node-1):(nrow(connections_all)),])),
          do.call("rbind",replicate((length(time_set)-2),connections_all[(n_node-1):(nrow(connections_all)),],simplify=FALSE))),
    cbind(tail(time_set,1),connections_all[which(connections_all[,2]==n_node),])
  )
  
  rownames(xpij_matrix) <- NULL
  
  yjk_matrix <- target_attributes_all
  rownames(yjk_matrix) <- NULL
  
  xpij_ord <- as.matrix(cbind(paste0(xpij_matrix[,1],",",xpij_matrix[,2],",",xpij_matrix[,3] )))
  xpij_names <- apply(xpij_matrix[,1:3],1, function (x) paste0('x',",",x[1],",",x[2],",",x[3]))  
  xpij_obj <- rep(0,length(xpij_names))
  
  
  # yjk Binary Variables: Indicates  kth information collection from vertex j.
  yjk_names <- apply(target_attributes_all[,1:2],1, function (y) paste0('y',",",y[1],",",y[2]))
  yjk_obj <- target_attributes_all[,5]*1000 #00000
  
  
  # mu_Tf - total flight duration of the route f
  mu_Tf_max <- mission_time-(qnorm(1-mission_late_probability) * 1)
  mu_Tf_min <- 0
  mu_Tf_range <- (mu_Tf_max-mu_Tf_min)
  mu_Tf_names <- c('mu_Tf')
  mu_Tf_obj <- c(-epsilon_coef*(1/mu_Tf_range))
  
  
  # theta_f - average detection on the route f
  theta_f_max <- (-log(1-beta_f))
  theta_f_min <- 0
  theta_f_range <- (theta_f_max-theta_f_min)
  theta_f_names <- c('theta_f')
  theta_f_obj<- c(-epsilon_coef*(1/theta_f_range))
  
  
  # var_Tf - total Flight Duration of the route f
  var_Tf_max <- variance_restriction
  var_Tf_min <- 0
  var_Tf_range <- (var_Tf_max-var_Tf_min)
  var_Tf_names <- c('var_Tf')
  var_Tf_obj <- c(-epsilon_coef*(1/var_Tf_range))
  
  # sigma_Tf_dummy
  sigma_Tf_dummy_names <- c('sigma_Tf_dummy')
  sigma_Tf_dummy_obj <- c(0)
  
  # w_u
  u_vec <- seq(1,ceiling(mission_time/2),1)
  
  wu_names <- unlist(lapply(u_vec, function (w) paste0('w',",",w)))
  wu_obj <- rep(0,length(wu_names))
  
  
  # all variables - ordered
  all_names <- c(xpij_names,yjk_names,wu_names,sigma_Tf_dummy_names,theta_f_names,var_Tf_names,mu_Tf_names)
  all_obj <- c(xpij_obj,yjk_obj,wu_obj,sigma_Tf_dummy_obj,theta_f_obj,var_Tf_obj,mu_Tf_obj)
  
  
  # constraint Construction --------------------------------------
  
  # Constraint (12) in the paper - Probability of being late is less than or equal to a given probability
  C2_1_xpij <- spMatrix(1,(length(xpij_names)), i=rep(1, (length(xpij_names))), j=(1:(length(xpij_names))), x= c(xpij_matrix[,4]))
  C2_1_yjk <- spMatrix(1,(length(yjk_names)), i=rep(1, (length(yjk_names))), j=(1:(length(yjk_names))), x= c(yjk_matrix[,3]))
  C2_1_wu <- spMatrix(1, length(wu_names))
  C2_1_sigma_Tf_dummy <- spMatrix(1, length(sigma_Tf_dummy_names))
  C2_1_theta_f <- spMatrix(1, length(theta_f_names))
  # C2_1_var_If<-spMatrix(1, length(var_If_names))
  C2_1_var_Tf <- spMatrix(1, length(var_Tf_names))
  C2_1_mu_Tf <- spMatrix(1, length(mu_Tf_names),i=1,j=1,x=(-1))
  C2_1<- cbind(C2_1_xpij,C2_1_yjk,C2_1_wu,C2_1_sigma_Tf_dummy,C2_1_theta_f,C2_1_var_Tf,C2_1_mu_Tf)
  C2_1_rhs <- c(0)
  C2_1_dir <- c("=")
  
  C2_2_xpij <- spMatrix(1,(length(xpij_names)), i=rep(1, (length(xpij_names))), j=(1:(length(xpij_names))), x= c(xpij_matrix[,5]))
  C2_2_yjk <- spMatrix(1,(length(yjk_names)), i=rep(1, (length(yjk_names))), j=(1:(length(yjk_names))), x= c(yjk_matrix[,4]))
  C2_2_wu <-spMatrix(1, length(wu_names))
  C2_2_sigma_Tf_dummy <- spMatrix(1, length(sigma_Tf_dummy_names))
  C2_2_theta_f <- spMatrix(1, length(theta_f_names))
  C2_2_var_Tf <- spMatrix(1, length(var_Tf_names),i=1,j=1,x=(-1))
  C2_2_mu_Tf <- spMatrix(1, length(mu_Tf_names))
  C2_2 <- cbind(C2_2_xpij,C2_2_yjk,C2_2_wu,C2_2_sigma_Tf_dummy,C2_2_theta_f,C2_2_var_Tf,C2_2_mu_Tf)
  C2_2_rhs <- c(0)
  C2_2_dir <- c("=")
  
  
  C2_3_xpij <- spMatrix(1,(length(xpij_names)))
  C2_3_yjk <- spMatrix(1,(length(yjk_names)))
  C2_3_wu <- spMatrix(1, length(wu_names))
  C2_3_sigma_Tf_dummy <- spMatrix(1, length(sigma_Tf_dummy_names),i=1,j=1,x=qnorm(1-mission_late_probability))
  C2_3_theta_f <- spMatrix(1, length(theta_f_names))
  C2_3_var_Tf <- spMatrix(1, length(var_Tf_names))
  C2_3_mu_Tf <- spMatrix(1, length(mu_Tf_names),i=1,j=1,x=1)
  C2_3 <- cbind(C2_3_xpij,C2_3_yjk,C2_3_wu,C2_3_sigma_Tf_dummy,C2_3_theta_f,C2_3_var_Tf,C2_3_mu_Tf)
  C2_3_rhs <- c(mission_time)
  C2_3_dir <- c("<=")
  
  
  C2_4_xpij <- spMatrix(1,(length(xpij_names)))
  C2_4_yjk <- spMatrix(1,(length(yjk_names)))
  C2_4_wu <- spMatrix(1, length(wu_names),i=rep(1,length(wu_names)), j=1:length(wu_names), x=u_vec^2)
  C2_4_sigma_Tf_dummy <- spMatrix(1, length(sigma_Tf_dummy_names))
  C2_4_theta_f <- spMatrix(1, length(theta_f_names))
  C2_4_var_Tf <- spMatrix(1, length(var_Tf_names),i=1,j=1,x=-1)
  C2_4_mu_Tf <- spMatrix(1, length(mu_Tf_names))
  C2_4 <- cbind(C2_4_xpij,C2_4_yjk,C2_4_wu,C2_4_sigma_Tf_dummy,C2_4_theta_f,C2_4_var_Tf,C2_4_mu_Tf)
  C2_4_rhs <- c(0)
  C2_4_dir <- c(">=")
  
  
  C2_5_xpij <- spMatrix(1,(length(xpij_names)))
  C2_5_yjk <- spMatrix(1,(length(yjk_names)))
  C2_5_wu <- spMatrix(1, length(wu_names),i=rep(1,length(wu_names)), j=1:length(wu_names), x=u_vec)
  C2_5_sigma_Tf_dummy <- spMatrix(1, length(sigma_Tf_dummy_names),i=1,j=1,x=-1)
  C2_5_theta_f <- spMatrix(1, length(theta_f_names))
  C2_5_var_Tf <- spMatrix(1, length(var_Tf_names))
  C2_5_mu_Tf <- spMatrix(1, length(mu_Tf_names))
  C2_5 <- cbind(C2_5_xpij,C2_5_yjk,C2_5_wu,C2_5_sigma_Tf_dummy,C2_5_theta_f,C2_5_var_Tf,C2_5_mu_Tf)
  C2_5_rhs <- c(0)
  C2_5_dir <- c("=")
  
  
  C2_6_xpij <- spMatrix(1,(length(xpij_names)))
  C2_6_yjk <- spMatrix(1,(length(yjk_names)))
  C2_6_wu <- spMatrix(1, length(wu_names),i=rep(1,length(wu_names)), j=1:length(wu_names), x=rep(1,length(wu_names)))
  C2_6_sigma_Tf_dummy <- spMatrix(1, length(sigma_Tf_dummy_names))
  C2_6_theta_f <- spMatrix(1, length(theta_f_names))
  C2_6_var_Tf <- spMatrix(1, length(var_Tf_names))
  C2_6_mu_Tf <- spMatrix(1, length(mu_Tf_names))
  C2_6 <- cbind(C2_6_xpij,C2_6_yjk,C2_6_wu,C2_6_sigma_Tf_dummy,C2_6_theta_f,C2_6_var_Tf,C2_6_mu_Tf)
  C2_6_rhs <- c(1)
  C2_6_dir <- c("=")
  
  C2 <- rbind(C2_1,C2_2,C2_3,C2_4,C2_5,C2_6)
  C2_rhs <- rbind(C2_1_rhs,C2_2_rhs,C2_3_rhs,C2_4_rhs,C2_5_rhs,C2_6_rhs)
  C2_dir <- rbind(C2_1_dir,C2_2_dir,C2_3_dir,C2_4_dir,C2_5_dir,C2_6_dir)
  
  
  # Constraint (13) in the paper, probability of number of detection is larger than 1 must be smaller than or equal to a given probability
  C3_xpij <- spMatrix(1,(length(xpij_names)), i=rep(1, (length(xpij_names))), j=(1:(length(xpij_names))), x= c(xpij_matrix[,6]))
  C3_yjk <- spMatrix(1,(length(yjk_names)), i=rep(1, (length(yjk_names))), j=(1:(length(yjk_names))), x= c(yjk_matrix[,6]))
  C3_wu <- spMatrix(1, length(wu_names))
  C3_sigma_Tf_dummy <- spMatrix(1, length(sigma_Tf_dummy_names))
  C3_theta_f <- spMatrix(1, length(theta_f_names),i=1,j=1,x=(-1))
  C3_var_Tf <- spMatrix(1, length(var_Tf_names))
  C3_mu_Tf <- spMatrix(1, length(mu_Tf_names))
  C3 <- cbind(C3_xpij,C3_yjk,C3_wu,C3_sigma_Tf_dummy,C3_theta_f,C3_var_Tf,C3_mu_Tf)
  C3_rhs <- c(0)
  C3_dir <- c("=")
  
  # Constraint (14) in the paper: upper bound  is defined when defining the upper bounds of the variables of the model
  
  # Constraint (15) in the paper - Number of collections at a node must be smaller than or equal to the number of times the node appears on the route
  C6 <- do.call("rbind",sapply(target_set, function (a){
    C6_xpij <- spMatrix(1,(length(xpij_names)))
    temp_6_x <- which(xpij_matrix[,3]==a)
    temp_6_y <- which(yjk_matrix[,1]==a)
    
    C6_xpij <- spMatrix(1, length(xpij_names), i=rep(1,length(temp_6_x)), j=temp_6_x, x=rep(-1, length(temp_6_x)))
    C6_yjk <- spMatrix(1,(length(yjk_names)), i=rep(1, length(temp_6_y)), j=temp_6_y, x= rep(1, length(temp_6_y)))
    C6_wu <- spMatrix(1, length(wu_names))
    C6_sigma_Tf_dummy <- spMatrix(1, length(sigma_Tf_dummy_names))
    C6_theta_f <- spMatrix(1, length(theta_f_names))
    C6_var_Tf <- spMatrix(1, length(var_Tf_names))
    C6_mu_Tf <- spMatrix(1, length(mu_Tf_names))
    C6_temp <- cbind(C6_xpij,C6_yjk,C6_wu,C6_sigma_Tf_dummy,C6_theta_f,C6_var_Tf,C6_mu_Tf)
    return(C6_temp)
  }))
  C6_rhs <- rep(0,length(target_set))
  C6_dir <- rep("<=",length(target_set))
  
  
  # Constraint (16) in the paper - Defines the revisits. Target j must be visited before a revisit
  #  temproray matrix to generate order of yjk variables
  temp7 <- (length(target_set)*(length(k_set)-1)) #number of constraints
  temp7_2 <- cbind(1:(length(yjk_names)-1), 2:length(yjk_names))
  temp7_2 <- temp7_2[-seq(length(k_set),length(yjk_names)-length(k_set),by=length(k_set)),]
  
  C7_xpij <- spMatrix(temp7, length(xpij_names))
  C7_yjk <- spMatrix(temp7,length(yjk_names), i=rep(1:temp7,each=(length(k_set)-1)), j=as.vector(t(temp7_2)), x=rep(c(1,-1),temp7))
  C7_wu <- spMatrix(temp7, length(wu_names))
  C7_sigma_Tf_dummy <- spMatrix(temp7, length(sigma_Tf_dummy_names))
  C7_theta_f <- spMatrix(temp7, length(theta_f_names))
  C7_var_Tf <- spMatrix(temp7, length(var_Tf_names))
  C7_mu_Tf <- spMatrix(temp7, length(mu_Tf_names))
  C7 <- cbind(C7_xpij,C7_yjk,C7_wu,C7_sigma_Tf_dummy,C7_theta_f,C7_var_Tf,C7_mu_Tf)
  C7_rhs <- rep(0,temp7)
  C7_dir <- rep(">=",temp7)
  
  # Constraint (18) in the paper: flow balance equations --> if there is an incoming arc to a target, there has to be an outgoing arc as well
  C8_xpij <- do.call("rbind",lapply(head(tail(time_set,-1),-1), function (p){
    
    C8_xpij_temp <- do.call("rbind",lapply(target_set, function (j){
      temp8_1 <- which(xpij_matrix[,1]==p & xpij_matrix[,3]==j) #incoming arcs to j at p
      temp8_2 <- which(xpij_matrix[,1]==(p+1) & xpij_matrix[,2]==j) #outgoing arcs from j at p+1
      C8_xpij_1 <- spMatrix(1, length(xpij_names), i=rep(1,length(temp8_1)), j=temp8_1, x=rep(1,length(temp8_1))) #incoming arcs to j at p
      C8_xpij_2 <- spMatrix(1, length(xpij_names), i=rep(1,length(temp8_2)), j=temp8_2, x=rep(-1,length(temp8_2))) #outgoing arcs from j at p+1
      C8_xpij_j <- (C8_xpij_1 + C8_xpij_2)
      return(C8_xpij_j)
    }))
    return(C8_xpij_temp)
  }))
  
  C8_yjk <- spMatrix(nrow(C8_xpij),length(yjk_names))
  C8_wu <- spMatrix(nrow(C8_xpij), length(wu_names))
  C8_sigma_Tf_dummy <- spMatrix(nrow(C8_xpij), length(sigma_Tf_dummy_names))
  C8_theta_f <- spMatrix(nrow(C8_xpij), length(theta_f_names))
  C8_var_Tf <- spMatrix(nrow(C8_xpij), length(var_Tf_names))
  C8_mu_Tf <- spMatrix(nrow(C8_xpij), length(mu_Tf_names))
  C8 <- cbind(C8_xpij,C8_yjk,C8_wu,C8_sigma_Tf_dummy,C8_theta_f,C8_var_Tf,C8_mu_Tf)
  C8_rhs <- rep(0,nrow(C8_xpij))
  C8_dir <- rep("=",nrow(C8_xpij))
  
  # flow balance for time period 1 and 2
  C9_xpij_1 <- spMatrix(n_target, length(xpij_names), i=(1:n_target), j=which(xpij_matrix[,1]==1), x=rep(1,n_target))
  C9_xpij_2 <- do.call("rbind",lapply(target_set,function(s){
    spMatrix(1, length(xpij_names), i=rep(1,(n_target+1)), j=which(xpij_matrix[,1]==2 & xpij_matrix[,2]==s ), x=rep(-1,(n_target+1)))
  }))
  C9_xpij <- C9_xpij_1 + C9_xpij_2
  C9_yjk <- spMatrix(n_target,length(yjk_names))
  C9_wu <- spMatrix(n_target, length(wu_names))
  C9_sigma_Tf_dummy <- spMatrix(n_target, length(sigma_Tf_dummy_names))
  C9_theta_f <- spMatrix(n_target, length(theta_f_names))
  # C9_var_If<-spMatrix(n_target, length(var_If_names))
  C9_var_Tf <- spMatrix(n_target, length(var_Tf_names))
  C9_mu_Tf <- spMatrix(n_target, length(mu_Tf_names))
  C9 <- cbind(C9_xpij,C9_yjk,C9_wu,C9_sigma_Tf_dummy,C9_theta_f,C9_var_Tf,C9_mu_Tf)
  C9_rhs <- rep(0,n_target)
  C9_dir <- rep("=",n_target)
  
  # Constraint (19) in the paper: flow balance for time period 1-must leave
  C10_xpij <- spMatrix(1, length(xpij_names), i=rep(1,n_target), j=which(xpij_matrix[,1]==1), x=rep(1,n_target))
  C10_yjk <- spMatrix(1,length(yjk_names))
  C10_wu <- spMatrix(1, length(wu_names))
  C10_sigma_Tf_dummy <- spMatrix(1, length(sigma_Tf_dummy_names))
  C10_theta_f <- spMatrix(1, length(theta_f_names))
  C10_var_Tf <- spMatrix(1, length(var_Tf_names))
  C10_mu_Tf <- spMatrix(1, length(mu_Tf_names))
  C10 <- cbind(C10_xpij,C10_yjk,C10_wu,C10_sigma_Tf_dummy,C10_theta_f,C10_var_Tf,C10_mu_Tf)
  C10_rhs <- c(1)
  C10_dir <- c("=")
  
  # Constraint (20) in the paper: flow balance for time period 1-must arrive
  temp11 <- which(xpij_matrix[,3]==n_node)
  C11_xpij <- spMatrix(1, length(xpij_names), i=rep(1,length(temp11)), j=temp11, x=rep(1,length(temp11)))
  C11_yjk <- spMatrix(1,length(yjk_names))
  C11_wu <- spMatrix(1, length(wu_names))
  C11_sigma_Tf_dummy <- spMatrix(1, length(sigma_Tf_dummy_names))
  C11_theta_f <- spMatrix(1, length(theta_f_names))
  C11_var_Tf <- spMatrix(1, length(var_Tf_names))
  C11_mu_Tf <- spMatrix(1, length(mu_Tf_names))
  C11 <- cbind(C11_xpij,C11_yjk,C11_wu,C11_sigma_Tf_dummy,C11_theta_f,C11_var_Tf,C11_mu_Tf)
  C11_rhs <- c(1)
  C11_dir <- c("=")
  
  
  # model solve --------------------------
  model <- list()
  model$modelname <- 'UAV Prize Collection'
  model$modelsense <- 'max'
  
  # initialize data for variables
  model$lb       <- 0
  model$ub       <- c(rep(1, length(xpij_names)),rep(1, length(yjk_names)), rep(1, length(wu_names)),ceiling(mission_time/2),(-log(1-beta_f)),variance_restriction,mu_Tf_max)
  model$vtype    <- c(rep("B", length(xpij_names)),rep("B", length(yjk_names)), rep("B", length(wu_names)),"C","C","C","C")
  model$obj      <- all_obj
  model$varnames <- all_names
  
  # build constraint matrix
  model$A        <- rbind(C2,C3,C6,C7,C8,C9,C10,C11)
  model$rhs      <- c(C2_rhs,C3_rhs,C6_rhs,C7_rhs,C8_rhs,C9_rhs,C10_rhs,C11_rhs)
  model$sense    <- c(C2_dir,C3_dir,C6_dir,C7_dir,C8_dir,C9_dir,C10_dir,C11_dir)
  
  

  # set parameters
  params <- list()
  params$MipFocus = 1
  
  
  ptm<-proc.time()
  result <- gurobi(model,params)
  runTime<-(proc.time()-ptm)
  runTime<-sum(runTime[1:2])
  

# extract the results --------------------------
decvr <- result$x

# expected information collection
expected_info_collected <- (result$objval+(-1*(tail(decvr,3)%*% c(theta_f_obj, var_Tf_obj,mu_Tf_obj))))/1000

# xp,i,j results
xpij_results <- decvr[1:length(xpij_names)]
xpij_variables <- all_names[1:length(xpij_names)]
xpij_variables <- xpij_variables[which(xpij_results>0.0001)] # >0.0001 is due to avoiding gurobi rounding errors
xpij_results <- xpij_results[which(xpij_results>0.0001)]


#yjk results
yjk_results <- decvr[(length(xpij_names)+1):(length(xpij_names)+length(yjk_names))]
yjk_variables <- all_names[(length(xpij_names)+1):(length(xpij_names)+length(yjk_names))]
yjk_variables <- yjk_variables[which(yjk_results>0.0001)] # >0.0001 is due to avoiding gurobi rounding errors
yjk_results <- yjk_results[which(yjk_results>0.0001)]

# mu_Tf_names, theta_f_names, var_Tf_names
other_results <- tail(decvr,3)[c(3,1,2)]


cat(c("mission_time", "variance_restriction", "expected_info_collected", "mu_Tf_names", "theta_f_names", "var_Tf_names", "result$objval", "result$objbound","result$mipgap","result$runtime", "runTime"), file = "civilian_application_summary.txt", append = TRUE)
cat("\n", file ="civilian_application_summary.txt", append = TRUE)

cat(c(mission_time,variance_restriction, expected_info_collected, other_results, result$objval, result$objbound,result$mipgap,result$runtime, runTime), file = "civilian_application_summary.txt", append = TRUE)
cat("\n", file ="civilian_application_summary.txt", append = TRUE)


sink("civilian_application_dec_space.txt", append=TRUE)
cat("Mission Time:", mission_time)
cat("\n")
cat("----------------------------------------")
cat("\n")
cat("Status:",result$status )
cat("\n")
cat("Optimality gap:",(result$mipgap)*100)
cat("\n")
cat("Run Time given by Gurobi:", result$runtime)
cat("\n")
cat("Run Time given by R:", runTime)
cat("\n")
cat("Objective value:",result$objval)
cat("\n")
cat("Best available upper bound:",result$objbound)
cat("\n")
cat("Expected Inf Collection:",expected_info_collected)
cat("\n")
cat("Other Variables:",other_results)
cat("\n")
cat("xpij Result: ",  xpij_variables)
cat("\n")
cat("xpij Result: ",xpij_results)
cat("\n")
cat("yjk Result: ",  yjk_variables)
cat("\n")
cat("yjk Result: ",yjk_results)
cat("\n")
cat("\n")
cat("\n")
sink()

