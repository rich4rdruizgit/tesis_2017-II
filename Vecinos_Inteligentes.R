#NOTA IMPORTANTE este algoritmo esta MODIFICADO para un optimizar el rendimiento 

#limpia objectos del espacio de trabajo
rm(list=ls(all=TRUE))

#asignacio de semilla 
set.seed(3)


#librerias 
library(foreach)
library(doParallel)
library(parallel)
library(iterators)
library(modeest)

#directorios
nameProject = "p3_x2_output4" #folder name
path = "~/clr/problem3/p3_x2_1207/"
pathProject = paste(path,nameProject, sep="")
dir.create(pathProject)

# Cargando los datos Data ---------

load("datos.rda")
load("datos_t.rda")
load("datos_nor.rda")




# inicializando variables ---------
variable_names = c("age","aadt","trucks","elevation","precip", "min_temp", "max_temp", "wet_days", "freeze_thaw", "rut_depth", "factor(number_of_lanes)", "factor(sys_id)","factor(f_class)","factor(category)")
psi_column_position = which(colnames(mydata) == "psi")
numvariables = length(variable_names)
numsegments = max(mydata[,1])
numobservations = nrow(mydata)

#Simulated annealing parameters ---------
temp = 10
temp_min = 1E-25
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

matriz_distancia<-matrix(1:13498276,nrow=3674,ncol=3674,byrow = TRUE)
for(i in 1:NROW(matriz_distancia)){
  for(j in 1:NCOL(matriz_distancia)){
    matriz_distancia[i,j] = formulirri(variables_promediadas[[i]],variables_promediadas[[j]])
  } 
}

print("fin")



