library(foreach)
library(doParallel)
library(parallel)
library(iterators)
library(modeest)


#Semilla
set.seed(3)


# install.packages("modeest")

pathCSV = paste("~/clr/data/data.csv",sep="")
mydata = read.csv(pathCSV) #lee los archivos del dataset
#inicializacion de variables

variable_names = c("age","aadt","trucks","elevation","precip", "min_temp", "max_temp", "wet_days", "freeze_thaw", "rut_depth", "factor(number_of_lanes)", "factor(sys_id)","factor(f_class)","factor(category)")
psi_column_position = which(colnames(mydata) == "psi")
numvariables = length(variable_names)
numsegments = max(mydata[,1])
numobservations = nrow(mydata)
variables=list()
#nuestras variables

variables_promediadas=list()
psi=list()
ord<-c(11,12) #ejemplos de las columas ordinales para el metodo de promedio
var_exp_ord =c(1:10, ord)
var_cat=c(11,14)#ejemplos de las columas ordinales para el metodo de promedio


#1
#normalizacioncirri
tabla_normalizada<-normalizacion(mydata)

normalizacion <- function(mdataset){
  
  for(i in 3:(NCOL(mdataset)-4)){
    
    mdataset[1:NROW(mdataset),i]=(mdataset[1:NROW(mdataset),i]-min(mdataset[,i]))/(max(mdataset[,i])-min(mdataset[,i]))
  }
  
  return(mdataset)
}

mydata <-normalizacion(mydata)



# max(mydata[,10])
# min(mydata[,10])
#promedio de cada segmento
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



#formulirri para la distancirri
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
print("inicializando los segmentirris")
for (segment in 1:numsegments){
  
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



#richard
set.seed(1)
positions <- sample(1:numsegments, numsegments, replace=F)

