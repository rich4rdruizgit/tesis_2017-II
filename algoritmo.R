# clear objects from the workspace ------
rm(list=ls(all=TRUE))
#options(OutDec = ".")
set.seed(3)
#==============================
# Library ----------

#install.packages("foreach")
#install.packages("doParallel")
#install.packages("parallel")
#install.packages("iterators")
#install.packages("modeest")


library(foreach)
library(doParallel)
library(parallel)
library(iterators)
library(modeest)
nameProject = "p3_x2_output4" #folder name
path = "~/clr/problem3/p3_x2_1207/"
pathCSV = paste("~/clr/data/data.csv",sep="")
pathTest = paste("~/clr/data/test_data.csv",sep="") #location of the 2 last years database 
pathProject = paste(path,nameProject, sep="")
print("Iniciando  Programa")
#Path modificados a mi pc
#D:\Universidad\Tesis\Ante- Proyecto\Paper_3


#directory ---------
dir.create(pathProject)

#Data ---------
#If there is a problem to read data, remove separator for both lines. For example : read.csv(pathCSV) instead of read.csv(pathCSV,sep=";")
mydata = read.csv(pathCSV) # the big database
mydatatest = read.csv(pathTest) # the 2 last years database

variable_names = c("age","aadt","trucks","elevation","precip", "min_temp", "max_temp", "wet_days", "freeze_thaw", "rut_depth", "factor(number_of_lanes)", "factor(sys_id)","factor(f_class)","factor(category)")
psi_column_position = which(colnames(mydata) == "psi")
numvariables = length(variable_names)
numsegments = max(mydata[,1])
numobservations = nrow(mydata)

#Simulated annealing parameters ---------
temp = 10
temp_min = 1E-5
alpha = 0.97
number_of_neighbors = 5
Boltz_const = 80

#Clustering parameters --------
min_clusters = 4
max_clusters = min_clusters
min_observations_cluster = 800
num_changes = 80

#Statistical parameters --------
level_of_significance = 0.05
max_vif = 10

#Script variables----------
psi=list()
variables <- list()
variables_promediadas = list()
#file
output_numcluster_file = paste(pathProject,"/Number of cluster.txt",sep="")
output_final_plot = paste(pathProject,"/allclusters.jpg",sep="")
all_variables_clusters <- list()

cl<-makeCluster(detectCores()-1)
registerDoParallel(cl)
print("inicializo las variables")

#---------------Data pre processing--------------
# never change these part of data pre processing
#Fill variables and psi with mydata
#Llenar variables y psi con mis datos , datos = data.csv
# ----------------------Normalizacion del dataset --------
normalizacion <- function(mdataset){
  for(i in 3:(NCOL(mdataset)-4)){
    mdataset[1:NROW(mdataset),i]=(mdataset[1:NROW(mdataset),i]-min(mdataset[,i]))/(max(mdataset[,i])-min(mdataset[,i]))
  }
  return(mdataset)
}
mydata <- normalizacion(mydata)
print("normalizo")

#-----------------------Promedio de cada Segmento con su respectivo identificador con modas -------
promedio_segmento <- function(variables, segment,id){
  promed<-colMeans(variables[[segment]][,1:10])
  modas<-apply(variables[[segment]][,11:14], 2, mlv,  method = "mfv")
  seg_pro_mod<-as.array(promed)
  for(i in 1:length(modas)){
    seg_pro_mod <- c(seg_pro_mod,modas[[i]]$M)
  }
  seg_pro_mod <- c(seg_pro_mod,id)
  return(seg_pro_mod)
}

print("cargo promedio")

#----------------Formula de distancia -------------------
formulirri <- function(seg_pro_mod,seg_pro_mod2){
  sum=0
  for(i in 1:(length(seg_pro_mod)-4)){
    sum= sum + ((seg_pro_mod[i]-seg_pro_mod2[i])^2)
  }
  sum2=0
  for(i in 11:length(seg_pro_mod)){
    if(seg_pro_mod[i]!=seg_pro_mod2[i])
    {
      sum2=sum2+1
    }
  }
  dist_euquidiana<-sqrt((sum+sum2))
  return (dist_euquidiana)
}

l=1
for (segment in 1:numsegments){
  #browser("Llenando psi")
  numyears = mydata[l,2]
  p=array(NA,dim=c(numyears))
  v=array(NA,dim=c(numyears,numvariables))
  for (year in 1:numyears){
    p[year]=mydata[l,psi_column_position]
    for (variable in 1:numvariables){
      v[year,variable] = mydata[l,variable+psi_column_position]
    }
    l=l+1
  }
  variables[[segment]] = v
  id_seg<-mydata[l-1,1]
  variables_promediadas[[segment]] <- promedio_segmento(variables,segment,id_seg)
  psi[[segment]] = p
}

#------------------Matriz de similitud --------------
matriz_distancia<-matrix(1:13498276,nrow=3674,ncol=3674,byrow = TRUE)
for(i in 1:NROW(matriz_distancia)){
  for(j in 1:NCOL(matriz_distancia)){
    matriz_distancia[i,j] = formulirri(variables_promediadas[[i]],variables_promediadas[[j]])
  } 
}

print("cargo funcion")
#------------------Generando orden en la asignacion de cluster-----------
generate_order_cluster <- function(){
  j=1
  clusters = rep(0,numsegments)
  for(i in 1:numsegments){
    clusters[i] <- j
    j <- j+1
    if(j==numclusters+1){
      j<-1
    }
  }
  return (clusters)
}

#------------------metodo para asignar clusters en baraja------
generate_random_cluster <- function(segmentos){
  clusters = rep(0,numsegments)
  positions <- sample(1:numsegments, numsegments, replace=F)
  for(i in 1:numsegments){
    clusters[positions[i]]<- segmentos[i]
  }
  return (clusters)
}

print("cargo baraja")
#------------------ Lo de ellos ------------
#-------------calculation for kmax---------------------
cluster_max <- function(matrix_data,n){
  browser("K max")
  if(n>nrow(matrix_data)){
    k_max <- 0
    result <- k_max
  }else{
    k_max <- 0
    max_observations <- max(table(matrix_data[,1])) #return the maximum number of observation for a unique segment
    m <- matrix(0,max_observations,2) #we create a matrix with max_observation rows and 2 columns.
    
    # we put into matrix the number of observation in the first column and the number of segment in the second column
    for(i in c(1:max_observations)){
      m[i,1] <- i
      m[i,2] <- sum(table(matrix_data[,1])==i)
    }
    #first of all we count how many segments have more than (n) observations
    if(n<=nrow(m)){
      for(i in c(n:nrow(m))){
        k_max <- k_max+m[i,2]
        m[i,2] <- 0
      }
    }
    
    result <- k_max #useful if n is 1
    browser("Hola k max")
    #if it still some segments, we create a matrix with the rest of the data (without 0)
    while(sum(m[,2])!=0){
      m <- subset(m, m[,2] != 0) #delete lines with 0
      alpha <- nrow(m)
      beta <- alpha
      m_save <- m
      if(sum(m[,1]*m[,2])<n){ #if there is not enough observations
        for(i in c(1:nrow(m))){
          m[i,2] <- 0
        }
        result <- k_max
      }else{
        s <- m[alpha,1]
        m[alpha,2] <- m[alpha,2]-1
        while(s != n){
          if(m[beta,2]==0){
            beta <- beta-1
            if(beta==0){
              m <- m_save
              beta <- nrow(m)
              n <- n+1
              s <- 0
            }
          }
          if(s>n){
            s <- s-m[beta,1]
            m[beta,2] <- m[beta,2]+1
            beta <- beta-1
            if(beta==0){
              alpha <- alpha-1
              m <- m_save
              beta <- nrow(m)
              s <- 0
              if(alpha==0){
                n <- n+1
                alpha <- nrow(m)
                beta <- alpha
              }else{
                s <- m[alpha,1]
                m[alpha,2] <- m[alpha,2]-1
              }
            }
          }
          s <- s+m[beta,1]
          m[beta,2] <- m[beta,2]-1
        }
        k_max <- k_max+1
        result <- k_max
      }
    }
  }
  return (result)
}

#------------Coeff function-----------------
my_coeff <- function(x){
  coeff <-  coef(summary(x))[,1]
}

# Description:  The function combination determines all combinations of a list of objects
combination <- function(objects){
  n <- length(objects) #length of the vector in parameter
  ind <- c(1:n) #vector of index
  number_combination <- (2^n)-1 #calculation of the number of combination
  comb <- list() # we create a list which will contain every combinations
  length(comb) <- 1 # length of the list comb
  v=-1
  #we create the list with number combination
  for(i in ind){
    matrix<-t(combn(n,i))
    v<-v+choose(n,(i-1))
    for(j in (1:nrow(matrix))){
      comb[[j+v]] <- matrix[j,]
    }
  }
  #we replace numbers by elements
  for(i in (1:length(comb))){
    vect_number <- comb[[i]]
    for(j in (1:length(vect_number))){
      comb[[i]][j]<-objects[as.numeric(vect_number[j])]
    }
  }
  comb
}


#==============================
# Vif functions
#==============================
# never change these 3 vif functions

"vif" <-
  function(xx, ...)
    UseMethod("vif")

"vif.default" <-
  function(xx, y.name, na.action=na.exclude, ...) {
    nnames <- names(xx)
    nn.x <- seq(along=nnames)
    if (missing(y.name))
      y.number <- 0
    else {
      y.number <- match(y.name, nnames, 0)
      nn.x <-  nn.x[-y.number]
    }
    r2 <- nn.x
    names(r2) <- nnames[nn.x]
    if (length(r2) < 2) stop("vif requires two or more X-variables.")
    for (i in nn.x) {
      tmp.lm <- lm(xx[,i] ~
                     data.matrix(xx[,-c(y.number, i)]),
                   na.action=na.action)
      r2[nnames[i]] <- summary(tmp.lm)$r.squared
    }
    1/(1-r2)
  }

"vif.lm" <-
  function(xx, na.action=na.exclude, ...) {
    xxx <- xx  ## deal with scope problem
    if(length(xxx$x)==0 ||
       !(class(xxx$x) == "model.matrix" || class(xxx$x) == "matrix")) {
      xxx <- try(update(xxx, x = TRUE), silent=TRUE)
      if (class(xxx) == "Error" || class(xx)=="try-error") ## S-Plus || R
        stop("Please recompute the 'lm' object with 'x=TRUE'.")
    }
    xx <- as.data.frame(unclass(xxx$x))[-1]
    vif(xx, na.action=na.action)
  }

#-------------create_my_data_cluster---------------
create_my_data_cluster <- function(c, member){
  #browser("My data cluster")
  m <- matrix(data=NA , nrow = 0, ncol = ncol(mydata))
  for(segment in 1:numsegments){
    if(member[segment]==c){ 
      segment_data <- mydata[mydata[, 1] == segment,] 
      m <- rbind(m,segment_data) 
    }
  }
  return(m)
}

#-------------create_my_data_cluster_log---------------
create_my_data_cluster_log <- function(m_data_cluster){
  m <- m_data_cluster
  for(i in 1:nrow(m)){
    for(j in 1:ncol(m)){
      if(j>=psi_column_position && j<psi_column_position+11){
        m[i,j]<-log(m[i,j])
      }
      if(m[i,j]==-Inf){
        m[i,j]<-0
      }
    }  
  }
  return(m)
}

#---------------remove_categorical_variable_one_value---------------
remove_categorical_variable_one_value <- function(data_c){
  #browser("Remover variable categorica")
  positions_factor <- grep("factor",variable_names) + psi_column_position
  thenames <- variable_names
  for(i in positions_factor){
    if(length(table(data_c[i])) < 2){
      print("Segun pacho nunca entra aca")
      thenames[grep(names(data_c[i]), thenames)] <- "NA"
    }
  }
  return(thenames)
}

#---------Remove high value of vif----------
remove_high_vif <- function(fit,newnames){
  #browser(" remueve valors mayor de vif")
  if(length(grep("NA",newnames))<(length(newnames)-2)){
    vif_fit = vif(fit)
    if(!is.nan(max(vif_fit))){ #si , son diferentes de null
      if(max(vif_fit)>max_vif){ # si el valor maximo de los vif_fit> max_vif
        max_position = which.max(vif_fit)
        if(grepl('^factor',names(max_position))){#if max is categorical variable
          newnames[grep(sub("factor\\(","",sub("\\).*","",names(max_position))),newnames)] = NA
        }else{
          newnames[grep(names(max_position),newnames)] = NA
        }
      }
    }
  }
  return (newnames)
}

#-----------Generate a neighbor clusters-----------
generate_neighbor_clusters <- function(clusters){
  # browser("Generar grupos vecinos")
  repeat{
    all_nb_neighbors <<- all_nb_neighbors + 1
    positions <- sample(1:length(clusters), num_changes, replace=F)
    for(i in 1:length(positions) ){
      position = positions[i]
      old_cluster = clusters[position]
      repeat{ # repeat if new cluster is still the same old cluster 
        new_cluster = sample(1:numclusters, 1)
        if (new_cluster != old_cluster) {break}
      }
      clusters[position] = new_cluster
    }
    if(valid_clusters(clusters)){break}
  }
  return(clusters)
}

#-------------Calculate BIC--------------------
calculate_BIC_overall <- function(fits){
  #browser("Calculo del BIC")
  error = 0
  for(i in 1:length(fits)){
    error = error + sum(resid(fits[[i]])^2)
  }
  var_T = 0
  for(i in 1:numclusters){
    var_T = var_T + length(my_coeff(fits[[i]]))
  }
  var_T = var_T + numclusters - 1
  BIC = numobservations + numobservations*log(2*pi) + numobservations*log(error/numobservations) + var_T*log(numobservations)
  W_BIC = BIC
  eval<-c(BIC, W_BIC, error)
  return(eval)
}

#-------------Check acceptance probability-----------
acceptance_probability <- function(actual_eval, neighbor_eval, temp){
  prob = 0 
  if(neighbor_eval < actual_eval){
    prob = 1
  }else{
    prob = exp(-(abs((neighbor_eval-actual_eval)/Boltz_const)/temp))
  }
  return (prob)
}

#-------------Estimate psi functions-----------
estimate_psi<- function(fits, clusters, bData){
  if(bData){
    data<-mydata
  }else{
    data<-mydatatest
  }
  estm_psi_plot <- rep(0,nrow(data))
  for (i in 1:nrow(data)){ #for each observation of data
    current_sample_id = data[i,1] #take the ID
    current_cluster = clusters[current_sample_id] #take the membership
    betas = cbind(names(fits[[current_cluster]]$coefficients),as.double(fits[[current_cluster]]$coefficients))
    betas[is.na(betas)] <- 0
    estm_psi_plot[i] <- exp(as.double(betas[1,2])) #intercept
    unique_var <- unique(substr(betas[1:nrow(betas),1],1,nchar(betas[1:nrow(betas),1])-1))
    nb_var <- length(unique_var)
    for(j in 2:nb_var){
      k=1
      var = gsub("\\(|\\)","",unique_var[j])
      for(k in 1:length(variable_names)){
        if(grepl(var, gsub("\\(|\\)","",variable_names[k]))){
          if(grepl('^factor',variable_names[k])){ # if categorical variable
            if(length(which(betas[,1]==paste(variable_names[k],data[i,k+psi_column_position],sep=""))) != 0){
              estm_psi_plot[i] <- estm_psi_plot[i] * exp(as.double(betas[which(betas[,1]==paste(variable_names[k],data[i,k+psi_column_position],sep="")),2]))
            }
          }else{
            if(data[i,k+psi_column_position]!=0){
              estm_psi_plot[i] <- estm_psi_plot[i] * (data[i,k+psi_column_position]**as.double(betas[which(betas[,1]==variable_names[k]),2]))
            }
          }
        }
      }
    }
  }
  return(estm_psi_plot)
}

#-------------Estimate log psi functions-----------
estimate_log_psi<- function(fits, clusters, bData){
  if(bData){
    data<-mydata
  }else{
    data<-mydatatest
  }
  estm_psi_plot <- rep(0,nrow(data))
  for (i in 1:nrow(data)){ #for each observation of data
    current_sample_id = data[i,1] #take the ID
    current_cluster = clusters[current_sample_id] #take the membership
    betas = cbind(names(fits[[current_cluster]]$coefficients),as.double(fits[[current_cluster]]$coefficients))
    betas[is.na(betas)] <- 0
    estm_psi_plot[i] <- as.double(betas[1,2]) #intercept
    unique_var <- unique(substr(betas[1:nrow(betas),1],1,nchar(betas[1:nrow(betas),1])-1))
    nb_var <- length(unique_var)
    for(j in 2:nb_var){
      k=1
      var = gsub("\\(|\\)","",unique_var[j])
      for(k in 1:length(variable_names)){
        if(grepl(var, gsub("\\(|\\)","",variable_names[k]))){
          if(grepl('^factor',variable_names[k])){ # if categorical variable
            if(length(which(betas[,1]==paste(variable_names[k],data[i,k+psi_column_position],sep=""))) != 0){
              estm_psi_plot[i] <- estm_psi_plot[i] + as.double(betas[which(betas[,1]==paste(variable_names[k],data[i,k+psi_column_position],sep="")),2])
            }
          }else{
            if(data[i,k+psi_column_position]!=0){
              estm_psi_plot[i] <- estm_psi_plot[i] + as.double(betas[which(betas[,1]==variable_names[k]),2])*log(data[i,k+psi_column_position])
            }
          }
        }
      }
    }
  }
  return(estm_psi_plot)
}

#-------------calcul_BIC_combinaison------------
calcul_BIC_combinaison <- function(comb, data_c){
  #   browser("Cmbinacion del calculo bic")
  result <- matrix(0, length(comb), 2)
  for(k in 1:length(comb)){
    var_names = paste(comb[[k]],collapse = "+")
    fit <- lm(paste("psi ~ ",var_names), data=data_c)
    result[k,1]<-k
    result[k,2]<-BIC(fit)
  }
  result = result[order(result[,2]),]#order tab with BIC
  return(result)
}

#==============================
# Print results functions
#==============================
# All of these functions dont have to change except print_overfitting if we create a new model (this one is for power model)

#-------------Print parameters--------------------
print_parameters <- function(){
  write(paste("########### Configuration ###########"),file=output_opt_file,append=TRUE)
  write(paste("psi_column_position = ",psi_column_position),file=output_opt_file,append=TRUE)
  write(paste("numvariables = ",numvariables),file=output_opt_file,append=TRUE)
  write(paste("numsegments = ",numsegments),file=output_opt_file,append=TRUE)
  write(paste("numobservations = ",numobservations),file=output_opt_file,append=TRUE)
  write(paste("numclusters = ",numclusters),file=output_opt_file,append=TRUE)
  write(paste("max_clusters = ",max_clusters),file=output_opt_file,append=TRUE)
  write(paste("min_observations_cluster = ",min_observations_cluster),file=output_opt_file,append=TRUE)
  write(paste("num_changes = ",num_changes),file=output_opt_file,append=TRUE)
  write(paste("Boltzmann constant = ", Boltz_const),file=output_opt_file,append=TRUE)
  write(paste("temp = ",temp),file=output_opt_file,append=TRUE)
  write(paste("temp_min = ",temp_min),file=output_opt_file,append=TRUE)
  write(paste("alpha = ",alpha),file=output_opt_file,append=TRUE)
  write(paste("number_of_neighbors = ",number_of_neighbors),file=output_opt_file,append=TRUE)
}

#-------------Print BIC values--------------------
print_optimization <- function(j, bic, wbic, sse){
  if(j == 0){
    write(paste("\n########### Optimization ############"),file=output_opt_file,append=TRUE)
    write(paste("Iteration              BIC                W_BIC                  SSE"),file=output_opt_file,append=TRUE)
  }
  print(paste("[",j,"]    ",bic,"     ",wbic,"      ",sse))
  write(paste("[",j,"]    ",bic,"     ",wbic,"      ",sse),file=output_opt_file,append=TRUE)
}

#-------------Print time values-------------------
print_time <- function(start_time, end_time){
  write(paste("\n################ End #################"),file=output_opt_file,append=TRUE)
  write(paste("Program start time = ", start_time),file=output_opt_file,append=TRUE)
  print(paste("Program start time = ", start_time))
  print(paste("Program end time = ", Sys.time()))
  write(paste("Program end time = ", Sys.time()),file=output_opt_file,append=TRUE)
  print(difftime(end_time,start_time))
  write(paste("Total = ", as.character.POSIXt(difftime(end_time,start_time))),file=output_opt_file,append=TRUE)
}

#-------------Print fit values--------------------
print_fits <- function(best_fits){
  write(paste("Number of clusters : ",numclusters),file=output_betas_file,sep="\n",append=TRUE)
  for(i in 1:length(best_fits)){
    cat("\n#####################################################################################################",file=output_betas_file,sep="\n",append=TRUE)
    cat(paste("Cluster [",i,"]"),file=output_betas_file,sep="\n",append=TRUE)
    out1<-capture.output(summary(best_fits[[i]]))
    cat(out1,file=output_betas_file,sep="\n",append=TRUE)
    cat(paste("BIC: ",BIC(best_fits[[i]])), file=output_betas_file, sep="\n", append=TRUE)
  } 
}

#-------------Overfitting test results -----------
print_overfitting <- function(fits, estm_psi, bData){
  TSSY = sum((log(mydata$psi) - mean(log(mydata$psi)))^2)
  RMSE = 0
  NRMSE = 0
  WCSSY = 0
  BCSSY = 0
  SSRY = 0
  SSEY = 0
  MAPE = 0
  MAE = 0
  
  for(k in 1:numclusters){ # for each cluster
    WCSSY = WCSSY + sum((data_clusters_log[[k]]$psi - mean(data_clusters_log[[k]]$psi))^2)
    BCSSY = BCSSY + nrow(data_clusters_log[[k]]) * (mean(data_clusters_log[[k]]$psi)-mean(log(mydata$psi)))^2
    SSRY = SSRY + sum((predict(fits[[k]]) - mean(data_clusters_log[[k]]$psi))^2)
    SSEY = SSEY + sum(resid(fits[[k]])^2)
    p = length(my_coeff(fits[[k]]))-1
    n = nrow(data_clusters[[k]])
  }
  
  overall_R2 = 1 - SSEY/TSSY
  overall_BIC = calculate_BIC_overall(fits)
  
  if(bData){
    write(paste("\n########################## Number of Clusters",numclusters,"###############################"),output_overfitting_file,append=TRUE)
    write(paste("\nTSSY =",TSSY),output_overfitting_file,append=TRUE)
    write(paste("WCSSY =",WCSSY),output_overfitting_file,append=TRUE)
    write(paste("BCSSY =",BCSSY),output_overfitting_file,append=TRUE)
    write(paste("SSRY =",SSRY),output_overfitting_file,append=TRUE)
    write(paste("SSEY =",SSEY),output_overfitting_file,append=TRUE)
    write(paste("\nOverall R-squared =",overall_R2),output_overfitting_file,append=TRUE)
    write(paste("Overall W_BIC =",overall_BIC[2]),output_overfitting_file,append=TRUE)
    RMSE = sqrt(sum((log(mydata$psi) - estm_psi)^2, na.rm=TRUE)/nrow(mydata))
    NRMSE = RMSE/(max(log(mydata$psi))-min(log(mydata$psi)))
    MAPE = (100/nrow(mydata))*sum(abs((log(mydata$psi) - estm_psi)/log(mydata$psi)), na.rm=TRUE)
    MAE = (1/nrow(mydata))*sum(abs(log(mydata$psi) - estm_psi), na.rm=TRUE)
    write(paste("\nRMSE_training =",RMSE),output_overfitting_file,append=TRUE)
    write(paste("NRMSE_training =",NRMSE),output_overfitting_file,append=TRUE)
    write(paste("MAPE_training =",MAPE),output_overfitting_file,append=TRUE)
    write(paste("MAE_training =",MAE),output_overfitting_file,append=TRUE)
    cat("\n\n#####################################################################################################",file=output_overfitting_file,sep="\n",append=TRUE)
    cat("Statistics using test dataset",file=output_overfitting_file,sep="\n",append=TRUE)
    cat("Overall",file=output_overfitting_file,sep="\n",append=TRUE)
    jpeg(output_overfitting_plot_file_predict)
    plot(log(mydata$psi),estm_psi,xlab="Observed ln(PSI)", ylab="Predicted ln(PSI)",xlim=c(1,log(5)), ylim=c(1,log(5)),pch=2,main=paste("Number of clusters = ",numclusters))
    abline(a=0,b=1, col="red")
    legend(3.5,2,c(paste("\nRMSE =", round(RMSE, digits = 3)),paste("n =",format(nrow(mydata), big.mark=",", scientific = FALSE))))
    dev.off()
  }else{
    RMSE = sqrt(sum((log(mydatatest$psi) - estm_psi)^2, na.rm=TRUE)/nrow(mydatatest))
    NRMSE = RMSE/(max(log(mydatatest$psi))-min(log(mydatatest$psi)))
    MAPE = (100/nrow(mydatatest))*sum(abs((log(mydatatest$psi) - estm_psi)/log(mydatatest$psi)), na.rm=TRUE)
    MAE = (1/nrow(mydatatest))*sum(abs(log(mydatatest$psi) - estm_psi), na.rm=TRUE)
    write(paste("\nRMSE =",RMSE),output_overfitting_file,append=TRUE)
    write(paste("NRMSE =",NRMSE),output_overfitting_file,append=TRUE)
    write(paste("MAPE =",MAPE),output_overfitting_file,append=TRUE)
    write(paste("MAE =",MAE),output_overfitting_file,append=TRUE)
    jpeg(output_overfitting_plot_file)
    plot(log(mydatatest$psi),estm_psi,xlab="Observed ln(PSI)", ylab="Predicted ln(PSI)",xlim=c(1,log(5)), ylim=c(1,log(5)),pch=2,main=paste("Number of clusters = ",numclusters))
    abline(a=0,b=1, col="red")
    legend(3.5,2,c(paste("\nRMSE =", round(RMSE, digits = 3)),paste("n =",format(nrow(mydatatest), big.mark=",", scientific = FALSE))))
    dev.off()
  }
}

#------------- print rmse for each cluster--------------
print_rmse <- function(best_fits, best_clusters, estm_psi){
  for(i in 1:length(best_fits)){
    psicluster<-NULL
    psipredicted<-NULL
    sample_id = which(best_clusters==i)
    for(j in 1:length(sample_id)){
      psicluster <- c(psicluster, mydatatest$psi[which(mydatatest$sample_id==sample_id[j])])
      psipredicted <- c(psipredicted, estm_psi[which(mydatatest$sample_id==sample_id[j])])
    }
    RMSE = sqrt(sum((log(psicluster) - psipredicted)^2, na.rm=TRUE)/length(psicluster))
    NRMSE = RMSE/(max(log(psicluster))-min(log(psicluster)))
    MAPE = (100/length(psicluster))*sum(abs((log(psicluster) - psipredicted)/log(psicluster)), na.rm=TRUE)
    MAE = (1/length(psicluster))*sum(abs(log(psicluster) - psipredicted), na.rm=TRUE)
    cat("\n#####################################################################################################",file=output_overfitting_file,sep="\n",append=TRUE)
    cat(paste("Cluster [",i,"]"), file=output_overfitting_file,sep="\n",append=TRUE)
    cat(paste("\nRMSE = ",RMSE), file=output_overfitting_file, sep="\n", append=TRUE)
    cat(paste("NRMSE = ",NRMSE), file=output_overfitting_file, sep="\n", append=TRUE)
    cat(paste("MAPE = ",MAPE), file=output_overfitting_file, sep="\n", append=TRUE)
    cat(paste("MAE = ",MAE), file=output_overfitting_file, sep="\n", append=TRUE)
  }
}

#-------------Create comprehensive tab------------
print_comprehensive_tab <- function(fits, all_variables_clusters){
  num_coef = 0
  vect <- NULL
  #we group in a vector every names of coefficients of every clusters
  for(j in 1:numclusters){
    vect <- c(vect,all_variables_clusters[[j]])
  }
  #then we look how many different element there are
  vect <- unique(vect)
  num_coef <- length(vect)
  tab = matrix(0,num_coef+3,numclusters)
  dimnames(tab) = list(c(vect,"BIC","R-squared","Adjusted R-squared"),paste("K=",1:numclusters,sep=""))
  for(i in 1:numclusters){
    p_values <- coef(summary(fits[[i]]))[,4] #check pvalue
    j <- 1
    for(nb in 1:num_coef){
      if(j<=length(fits[[i]]$coefficients)){
        if(names(tab[,1][nb])==(names(fits[[i]]$coefficients[j]))){
          if(!is.na(p_values[j])){
            if(p_values[j] > level_of_significance){
              tab[nb,i]=fits[[i]]$coefficients[j]
            }else{
              tab[nb,i]=paste(fits[[i]]$coefficients[j],"*")
            }
          }else{
            tab[nb,i]=fits[[i]]$coefficients[j]
          }
          j <- j+1
        }
      }
    }
    tab[num_coef+1,i] = BIC(fits[[i]])
    tab[num_coef+2,i] = summary(fits[[i]])$r.squared
    tab[num_coef+3,i] = summary(fits[[i]])$adj.r.squared
  }
  write.csv(tab,file=output_comprehensive_betas_file)
}

#-------------print predicted psi------------
print_predicted_psi<-function(predicted_psi, best_clusters, bData){
  if(bData){
    data<-mydata
  }else{
    data<-mydatatest
  }
  data<-cbind(data,predicted_psi)
  membership<-NULL
  j=1
  for(i in 1:nrow(data)){
    repeat{
      if(data$sample_id[i]==j){
        membership<-c(membership,best_clusters[j])
        break
      }
      j <- j+1
    }
  }
  data<-cbind(data,membership)
  if(bData){
    write.csv(data,output_data_file_predicted_psi)
  }else{
    write.csv(data,output_datatest_file_predicted_psi)
  }
}

#-------------Create plot best number of cluster----------------
create_plot_best_nb_clusters <- function(){
  jpeg(output_final_plot)
  tab = read.table(file = output_numcluster_file,sep="=")
  plot(tab,type ="l",xlab="Number of clusters, K", ylab="BIC")
  for(i in 1:nrow(tab)){
    if(tab[i,2] == min(tab[2])){
      points(min(tab[i,1]),min(tab[2]),col=2,pch=1)
    }
  }
  legend(2, min(tab[2])+10, "Optimum number of clusters", col=2, pch=1)
  dev.off()
}

#-------------Rename best folder---------------
rename_best_folder <- function(){
  tab = read.table(file = output_numcluster_file,sep="=")
  for(i in 1:nrow(tab)){
    if(tab[i,2] == min(tab[2])){
      src.dir = paste(pathProject,"/K=",i+1,"/",sep="")
      dest.dir = paste(pathProject,"/best_K=",i+1,"/",sep="") 
      file.rename(from=src.dir, to=dest.dir) 
    }
  }
}

#==============================
# RUN CALIBRATION-----------
#==============================

#max_clusters = cluster_max(mydata,min_observations_cluster)

print("inicia algoritmo")
for (numclusters in min_clusters:max_clusters){
  all_nb_neighbors <<- 0
  dir.create(paste(pathProject,"/K=",numclusters,sep=""))
  output_opt_file = paste(pathProject, "/K=",numclusters,"/Optimization Info, numclusters ",numclusters,".txt",sep="")
  output_members_file = paste(pathProject, "/K=",numclusters,"/Membership, numclusters ",numclusters,".txt",sep="")
  output_betas_file = paste(pathProject, "/K=",numclusters,"/Betas, numclusters ",numclusters,".txt",sep="")
  output_comprehensive_betas_file = paste(pathProject, "/K=",numclusters,"/Comprehensive file for betas, numclusters ",numclusters,".csv",sep="")
  output_plot_file = paste(pathProject, "/K=",numclusters,"/BIC plot, numclusters ",numclusters,".jpg",sep="")
  output_overfitting_file = paste(pathProject,"/K=",numclusters,"/Model assessment statistics, numclusters ",numclusters,".txt",sep="")
  output_overfitting_plot_file = paste(pathProject,"/K=",numclusters,"/Predicted vs observed with test dataset, numclusters ",numclusters,".jpg",sep="")
  output_overfitting_plot_file_predict = paste(pathProject,"/K=",numclusters,"/Predicted vs observed with training dataset ",numclusters,".jpg",sep="")
  output_datatest_file_predicted_psi <- paste(pathProject, "/K=",numclusters,"/Test dataset with predicted psi and memberships, numclusters ",numclusters,".csv",sep="")
  output_data_file_predicted_psi <- paste(pathProject, "/K=",numclusters,"/Training dataset with predicted psi and memberships, numclusters ",numclusters,".csv",sep="")
  
  start_time = Sys.time()   
  print_parameters()
  profile <- data.frame()
  
  repeat{ # Repeat while significance is not OK
    best_clusters <- generate_order_cluster()
    
    best_clusters <- generate_random_cluster(best_clusters)
    readline("pausa alg")
    
    #best_clusters <- generate_random_cluster() #generate clusters initial randomly
    #this loop is to get the fit with remove vif function
    #======================================
    fits<-list()
    for(cluster in 1:numclusters){
      #browser("Primer for")
      #funciones de poweeeeeeer
      mydata_cluster <- create_my_data_cluster(cluster, best_clusters)
      out_clusters_file<-capture.output(mydata_cluster)
      cat(out_clusters_file,file="D:/Universidad/Tesis/Ante- Proyecto/Paper_3/saluda_cluster.txt",sep="", append = FALSE,fill = TRUE)
      #print(best_clusters)
      mydata_cluster_log <- create_my_data_cluster_log(mydata_cluster)
      out_clusters_file_log<-capture.output(mydata_cluster_log)
      cat(out_clusters_file_log,file="D:/Universidad/Tesis/Ante- Proyecto/Paper_3/saluda_cluster_log.txt",sep="", append = FALSE,fill = TRUE)
      
      
      newnames <- remove_categorical_variable_one_value(mydata_cluster_log)
      #print(newnames)
      repeat{
        #browser("primer repeat")
        #out<-capture.output (newnames)
        #cat(out,"richasr.txt",append = FALSE,fill = TRUE)
        names_without_na <- gsub("NA + ","",gsub(" + NA","",paste(newnames, collapse=" + "),fixed = TRUE),fixed = TRUE)
        #print(names_without_na)
        tryCatch({
          fit <- lm(paste("psi ~ ",names_without_na), data=mydata_cluster_log)
          #print(fit)
        }, error = function(ex) {return(NULL)})
        if(!exists("fit")){
          return(NULL)
        }
        old_newnames <- newnames
        newnames <- remove_high_vif(fit, newnames)
        if(length(grep("NA",paste(newnames))) == length(grep("NA",paste(old_newnames)))){break}
      } #termina segundo repeat
      comb_name <- combination(newnames[grep("NA",paste(newnames), invert=TRUE)])
      tab_best_BIC <- calcul_BIC_combinaison(comb_name, mydata_cluster_log)
      
      #check level of signficance
      valid = FALSE
      p=1
      repeat{
        #browser("segundo repeat")
        var_names <- paste(comb_name[[tab_best_BIC[p,1]]],collapse = "+")
        fit <- lm(paste("psi ~ ",var_names), data=mydata_cluster_log)
        p_values <- coef(summary(fit))[,4] #check pvalue
        if(max(p_values) > level_of_significance){
          valid = FALSE
        }else{
          valid = TRUE
        }
        if(valid){break}
        p = p + 1
        if(p > nrow(tab_best_BIC)){ return (NULL) }
      }
      fits[[cluster]] <- fit
    }
    #======================================
    if(!is.null(fits)){break}
  }
  #Aqui terminar el loop de arriba :;
  
  best_eval <- calculate_BIC_overall(fits)
  print_optimization(0, best_eval[1], best_eval[2], best_eval[3]) #aqui voy
  j = 1
  while (temp > temp_min){
    #browser("primer while")
    i = 1
    while (i <= number_of_neighbors){ #number of neighbors
      repeat{ # Repeat while significance is not OK
        neighbor_clusters <- generate_neighbor_clusters(best_clusters)
        #this loop is to get the fit with remove vif function
        #======================================
        fits<-list()
        for(cluster in 1:numclusters){
          mydata_cluster <- create_my_data_cluster(cluster, neighbor_clusters)
          mydata_cluster_log <- create_my_data_cluster_log(mydata_cluster)
          newnames <- remove_categorical_variable_one_value(mydata_cluster_log)
          repeat{
            names_without_na <- gsub("NA + ","",gsub(" + NA","",paste(newnames, collapse=" + "),fixed = TRUE),fixed = TRUE)
            tryCatch({
              fit <- lm(paste("psi ~ ",names_without_na), data=mydata_cluster_log)
            }, error = function(ex) {return(NULL)})
            if(!exists("fit")){
              return(NULL)
            }
            old_newnames <- newnames
            newnames <- remove_high_vif(fit, newnames)
            if(length(grep("NA",paste(newnames))) == length(grep("NA",paste(old_newnames)))){break}
          }
          comb_name <- combination(newnames[grep("NA",paste(newnames), invert=TRUE)])
          tab_best_BIC <- calcul_BIC_combinaison(comb_name, mydata_cluster_log)
          
          #check level of significance
          valid = FALSE
          p=1
          repeat{
            var_names <- paste(comb_name[[tab_best_BIC[p,1]]],collapse = "+")
            fit <- lm(paste("psi ~ ",var_names), data=mydata_cluster_log)
            p_values <- coef(summary(fit))[,4] #check pvalue
            if(max(p_values) > level_of_significance){
              valid = FALSE
            }else{
              valid = TRUE
            }
            if(valid){break}
            p = p + 1
            if(p > nrow(tab_best_BIC)){ return (NULL) }
          }
          fits[[cluster]] <- fit
        }
        #======================================
        if(!is.null(fits)){break}
      }
      neighbor_eval = calculate_BIC_overall(fits)
      ap = acceptance_probability(best_eval[2],neighbor_eval[2],temp)
      if(ap > runif(1, 0, 1)){
        best_clusters = neighbor_clusters
        best_fits = fits
        #this loop is to get the fit with remove vif function
        #======================================
        fits<-list()
        data_clusters<-list()
        data_clusters_log<-list()
        for(cluster in 1:numclusters){
          mydata_cluster <- create_my_data_cluster(cluster, neighbor_clusters)
          mydata_cluster_log <- create_my_data_cluster_log(mydata_cluster)
          data_clusters[[cluster]] <- mydata_cluster
          data_clusters_log[[cluster]] <- mydata_cluster_log
          newnames <- remove_categorical_variable_one_value(mydata_cluster_log)
          repeat{
            names_without_na <- gsub("NA + ","",gsub(" + NA","",paste(newnames, collapse=" + "),fixed = TRUE),fixed = TRUE)
            tryCatch({
              fit <- lm(paste("psi ~ ",names_without_na), data=mydata_cluster_log)
            }, error = function(ex) {return(NULL)})
            if(!exists("fit")){
              return(NULL)
            }
            old_newnames <- newnames
            newnames <- remove_high_vif(fit, newnames)
            if(length(grep("NA",paste(newnames))) == length(grep("NA",paste(old_newnames)))){break}
          }
          comb_name <- combination(newnames[grep("NA",paste(newnames), invert=TRUE)])
          tab_best_BIC <- calcul_BIC_combinaison(comb_name, mydata_cluster_log)
          
          #check level of significance
          valid = FALSE
          p=1
          repeat{
            var_names <- paste(comb_name[[tab_best_BIC[p,1]]],collapse = "+")
            fit <- lm(paste("psi ~ ",var_names), data=mydata_cluster_log)
            p_values <- coef(summary(fit))[,4] #check pvalue
            if(max(p_values) > level_of_significance){
              valid = FALSE
            }else{
              valid = TRUE
            }
            if(valid){break}
            p = p + 1
            if(p > nrow(tab_best_BIC)){ return (NULL) }
          }
          fits[[cluster]] <- fit
        }
        #======================================
        best_eval = neighbor_eval
      }
      i = i + 1
    }
    print_optimization(j, best_eval[1], best_eval[2], best_eval[3])
    profile <- rbind(profile, c(j, best_eval[2]))
    plot(profile[,2], type ="l", xlab="Iteration", ylab="BIC")
    jpeg(output_plot_file)
    plot(profile[,2], type ="l", col= "red", xlab="Iteration", ylab="BIC", main=paste("Number of clusters = ",numclusters))
    dev.off()
    temp = temp * alpha
    j = j + 1
    
  }
  end_time = Sys.time()
  write(paste(numclusters,"=",best_eval[2]),file=output_numcluster_file,append=TRUE)
  print_time(start_time, end_time)
  print_fits(best_fits)
  write.csv(matrix(cbind(1:length(best_clusters),best_clusters), ncol =2, dimnames = list(NULL,c("Segment ID", "Cluster Membership"))),file=output_members_file)
  write(paste("\nTotal number of valid neighbors tested =",(i-1)*(j-1)),output_opt_file,append=TRUE)
  write(paste("\nTotal number of neighbors tested =",all_nb_neighbors),output_opt_file,append=TRUE)
  
  best<-best_fits
  #this loop is to get the fit with remove vif function
  #======================================
  fits<-list()
  data_clusters <- list()
  data_clusters_log<-list()
  for(cluster in 1:numclusters){
    mydata_cluster <- create_my_data_cluster(cluster, best_clusters)
    mydata_cluster_log <- create_my_data_cluster_log(mydata_cluster)
    data_clusters[[cluster]] <- mydata_cluster
    data_clusters_log[[cluster]] <- mydata_cluster_log
    newnames <- remove_categorical_variable_one_value(mydata_cluster_log)
    temp=1
    repeat{
      names_without_na <- gsub("NA + ","",gsub(" + NA","",paste(newnames, collapse=" + "),fixed = TRUE),fixed = TRUE)
      tryCatch({
        fit <- lm(paste("psi ~ ",names_without_na), data=mydata_cluster_log)
      }, error = function(ex) {return(NULL)})
      if(!exists("fit")){
        return(NULL)
      }
      if(temp==1){
        all_variables_clusters[[cluster]]<-names(fit$coefficients)
      }
      temp<-temp+1
      old_newnames <- newnames
      newnames <- remove_high_vif(fit, newnames)
      if(length(grep("NA",paste(newnames))) == length(grep("NA",paste(old_newnames)))){break}
    }
    comb_name <- combination(newnames[grep("NA",paste(newnames), invert=TRUE)])
    tab_best_BIC <- calcul_BIC_combinaison(comb_name, mydata_cluster_log)
    
    #check level of significance
    valid = FALSE
    p=1
    repeat{
      var_names <<- paste(comb_name[[tab_best_BIC[p,1]]],collapse = "+")
      fit <- lm(paste("psi ~ ",var_names), data=mydata_cluster_log)
      p_values <- coef(summary(fit))[,4] #check pvalue
      if(max(p_values) > level_of_significance){
        valid = FALSE
      }else{
        valid = TRUE
      }
      if(valid){break}
      p = p + 1
      if(p > nrow(tab_best_BIC)){ return (NULL) }
    }
    fits[[cluster]] <- fit
  }
  #======================================
  best_fits <- fits
  print("###############################")
  print(best_fits)
  print_comprehensive_tab(best_fits, all_variables_clusters)
  estm_log_psi_datatest <- estimate_log_psi(best_fits,best_clusters,FALSE)
  estm_log_psi_data <- estimate_log_psi(best_fits,best_clusters,TRUE)
  print_overfitting(best_fits, estm_log_psi_data, TRUE)
  print_overfitting(best_fits, estm_log_psi_datatest, FALSE)
  print_rmse(best_fits, best_clusters, estm_log_psi_datatest)
  print_predicted_psi(estm_log_psi_datatest, best_clusters, FALSE)
  print_predicted_psi(estm_log_psi_data, best_clusters, TRUE)
}
create_plot_best_nb_clusters()
rename_best_folder()
pdf(NULL)
stopCluster(cl)
