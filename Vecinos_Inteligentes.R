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

load("datos_Ent.rdata")
load("datos_test.rda")



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




normalizacion <- function(mdataset){
  print(NCOL(mdataset))
  print("numero de cool")
  for(i in 3:(NCOL(mdataset)-4)){
    mdataset[1:NROW(mdataset),i]=(mdataset[1:NROW(mdataset),i]-min(mdataset[,i]))/(max(mdataset[,i])-min(mdataset[,i]))
  }
  return(mdataset)
}
mydatan <- normalizacion(mydata)

