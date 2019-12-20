summaryTopoFunc=function(FMinitial,FMfinal)
{
summaryy=matrix(0,2,4)

#FM
FM=genesGO
grafitness=induced.subgraph(graphALL, FMinitial)
descriptorsfitness=FunctionParam(grafitness) #descriptors
p1=predict(FunctionLDA,newdata=descriptorsfitness)
LDAy=LDAnorma(p1$x)
LDAy=ifelse(LDAy<0,0,LDAy)

nsize=length(FMinitial)
d=distGO[FMinitial,FMinitial]
d[is.na(d)]=0
betta=2*sum(d[lower.tri(d, diag = FALSE)])/(nsize*(nsize-1))
PercentageSize= length(intersect(genesGO,FMinitial))/length(genesGO)

colnames(summaryy)=c("Size","ScoreTopo","ScoreFunc","Criteria")
rownames(summaryy)=c("Initial_Module","Final_Module")

summaryy[1,1]=nsize
summaryy[1,2]=LDAy
summaryy[1,3]=betta
summaryy[1,4]=Criteria(LDAy,betta,PercentageSize)

#FMprime
grafitness=induced.subgraph(graphALL, FMfinal)
descriptorsfitness=FunctionParam(grafitness) #descriptors
p1=predict(FunctionLDA,newdata=descriptorsfitness)
LDAy=LDAnorma(p1$x)
LDAy=ifelse(LDAy<0,0,LDAy)

nsize=length(FMfinal)
d=distGO[FMfinal,FMfinal]
d[is.na(d)]=0
betta=2*sum(d[lower.tri(d, diag = FALSE)])/(nsize*(nsize-1))
PercentageSize= length(intersect(genesGO,FMfinal))/length(genesGO)

summaryy[2,1]=nsize
summaryy[2,2]=LDAy
summaryy[2,3]=betta
summaryy[2,4]=Criteria(LDAy,betta,PercentageSize)
return(summaryy)
}


