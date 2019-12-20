#ParametersPertinents
ParametersPertinents=function(graph,pertinent=TRUE,Round=6)
{
#Size
Size=vcount(graph)

#Inverse centrality
AShortest=NULL
sh=shortest.paths(graph)
Inverse=NULL
path=1/sh
Inverse=mean(path[!is.infinite(path)])

#Average shortest path length
# AShortest=mean(sh[!is.infinite(sh)])
sh[sh==0] <- Inf
pp=rowMeans(sh*is.finite(sh),na.rm=TRUE)
AShortest=sum(pp*is.finite(pp),na.rm=TRUE)/Size

#Betweenness
bet=betweenness(graph,normalized = TRUE) 
BTW=mean(bet) 

#Degree
Degreee=mean(degree(graph))

#Clustering coefficient
clus=NULL
clus=transitivity(graph,type = c("average"),isolates = c("zero"))

resultsPara=round(c(AShortest,BTW,Degreee,Inverse,clus,Size),Round)
names(resultsPara)=c("ShortestPath","Betweenness","Degree","Inverse","ClusteringCoef","Size")
if(!pertinent)
{
#Neighborhood connectivity 
voisins=NULL
Neib=NULL
NBC=NULL
for(i in 1:Size)
{
voisins=neighbors(graph, i)
Neib[i]=mean(degree(graph)[voisins])
}
NBC=mean(Neib,na.rm=TRUE)

#Radiality
rad=NULL	
sh[is.infinite(sh)] <- NA
diam=diameter(graph) 	 
for (v in V(graph))
rad[v] <- sum(diam+1-sh[v,],na.rm=TRUE)/(Size-1)/diam
RAD=mean(rad)

#Stress
library(sna)
adjan=get.adjacency(graph)
adjan=as(adjan,"matrix") 
stress=stresscent(adjan)
STR=mean(stress) 
detach("package:sna", unload=TRUE)

#Topoligical coefficient
####################
TOPO=function(gra)
{
fe=get.adjacency(gra) # Adjacency matrix
fe=as.matrix(fe)

if (nrow(fe) != Size)
{
tra=matrix(0,Size,Size)
nplus=0
nplus=Size-nrow(fe)
tra= bdiag(fe, diag(nplus))
pmal=Size-seq(1:nplus)+1
tra[pmal,pmal]=0
fe=as.matrix(tra)
}

J=matrix(0,Size,Size)
T=NULL

for ( i in 1:(Size-1))
for ( j in (i+1):Size)
{
wh1=which(fe[i,]==fe[j,])
if (length(which(fe[i,wh1]==1)) !=0 )
{
J[i,j]=length(which(fe[i,wh1]==1)) + fe[i,j]
J[j,i]=J[i,j]
}
}
for ( i in 1:Size)
{ 
if (sum(fe[i,]) > 2 )
T[i]=(sum(J[i,])/length(which(J[i,]!=0)))/sum(fe[i,])
if (sum(fe[i,]) <= 2 )
T[i]=0
T[is.na(T)] <- 0
}
return(mean(T))
}

####################
TOPOL=TOPO(graph)

#Eccentricity
ecc=eccentricity(graph) #eccentricity pour chaque gene
ECC=mean(ecc)

#Closeness
Closeness=NULL
Ipclo=NULL
Ipclo=1/sh
Ipclo[is.infinite(Ipclo)]=0 
Closeness=mean(Ipclo,na.rm=TRUE)

resultsPara=round(c(AShortest,Closeness,BTW,Degreee,ECC,NBC,STR,RAD,TOPOL,Inverse,clus,Size),Round)
names(resultsPara)=c("ShortestPath","Closeness","Betweenness","Degree","Eccentricity","Neighborhood","Stress","Radiality","TopoligicalCoef","Inverse","ClusteringCoef","Size")
}
return(resultsPara)
}

