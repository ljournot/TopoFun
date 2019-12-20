
#Replace "path to TopoFunc functions" by the path to the folder that stores the TopoFunc functions.
setwd("path to TopoFunc functions")

#TopoFunc

TopoFunc<-function(	TrueGenesGO,
					TrueGenesALL,
					graphALL,
					ReduceSpace=FALSE,
					OptimiseInitialisation=FALSE,
					Ngene=100,
					aleatoire=FALSE,
					Tpop_A=100,
					Tpop_B=100,
					Tpop_C=300,
					Pmax=500,
					FunctionParam=ParametersPertinents,
					FunctionLDA=ResultsLDA,
					distGO=MatriceGeneSim,
					LDAnorma=LDAnormalization,
					Criteria=fitness_function,
					selection="Tournament",
					tourn=2,
					pm=0.8,
					pc=0.5,
					AgeMax=10,
					plocImmig=1,
					conv=100
					)
{
#reswaves=res02
## random choice of genes
choicevar<-function(vec,Pmaxloc,Posis,nreste)
{
	vec=rep(NA,Pmaxloc+6)
	nbvar=sample(1:nreste,1)
	varsel=sort(sample(Posis,nbvar))
	vec[7:(6+nbvar)]=varsel
	vec[4]=nbvar
	vec[3]=0
	return(vec)
}

## deleting a gene
deletevar<-function(vec)
		{
			choixsup=sample(1:vec[4],1)				
			vec[6+choixsup]=NA
			vec[7:length(vec)]=sort(vec[7:length(vec)],na.last=TRUE)
			vec[4]=vec[4]-1
			vec[3]=0
			return(vec)
		}

## adding a gene
addvar<-function(vec,Posis)
		{
			poss=Posis
			#poss=poss[-vec[7:(6+vec[4])]]
			pola=vec[7:(6+vec[4])]
			poss=poss[!is.element(poss,pola)]
			if(length(poss)!=0)
			{
			choixadd=sample(poss,1)
			vec[7+vec[4]]=choixadd
			vec[4]=vec[4]+1
			vec[7:length(vec)]=sort(vec[7:length(vec)],na.last=TRUE)
			vec[3]=0
			}
			return(vec)
		}

## change of a gene
changevar<-function(vec,Posis)
		{
			choixsup=sample(1:vec[4],1)
			##choice of the variable to add among those which are not in the selection
			poss=Posis
			#poss=poss[-vec[7:(6+vec[4])]]
			pola=vec[7:(6+vec[4])]
			poss=poss[!is.element(poss,pola)]
			if(length(poss)!=0)
			{
			choixadd=sample(poss,1)
			vec[6+choixsup]=choixadd	
			vec[7:length(vec)]=sort(vec[7:length(vec)],na.last=TRUE)
			vec[3]=0
			}
			return(vec)
		}

## mutation
	mutation<-function(vec,Posis=posis,Pmaxloc=Pmax,Fbarre=Fbarre,Fmax=Fmax,pm=pm)
		{
			enrsol=vec
			
			pmloc=1 #if Fmax = Fbarre then all the solutions have the same fitness; we must make mutations
			if (enrsol[1]< Fbarre) pmloc=pm  
			
			if (enrsol[1]>= Fbarre & Fmax!=Fbarre)
			pmloc= pm * ((Fmax - enrsol[1])/(Fmax-Fbarre))
			
			mut<-runif(1,0,1)
			if (mut<pmloc)
			   {
				typemut<-runif(1,0,1)
				if (vec[4]>1 & vec[4]<Pmaxloc)
				    {
					if (typemut<0.5) return(changevar(vec,Posis)) else
					    {
						if (typemut<0.75) return(deletevar(vec)) else return(addvar(vec,Posis))
					    }
				    }
				else
				    {
					if (vec[4]==1)
					    {
						if (typemut<0.5) return(changevar(vec,Posis)) else return(addvar(vec,Posis))
					    }
					else
					    {
						if (typemut<0.5) return(changevar(vec,Posis)) else return(deletevar(vec))
					    }
				    }
			    }
			else return(vec)
		}

fitness<-function(vec)
{
genesfitness=unique(c(genesCliqGO,vec[7:(6+vec[4])]))
grafitness=induced.subgraph(graphALL, genesfitness)
descriptorsfitness=FunctionParam(grafitness) #descripteurs
p1=predict(FunctionLDA,newdata=descriptorsfitness)
LDAy=LDAnorma(p1$x)
LDAy=ifelse(LDAy<0,0,LDAy)
###################################################
 truegenesfitness=as.character(TrueGenesALL[genesfitness])
 nsize=length(genesfitness)
 d=distGO[truegenesfitness,truegenesfitness]
 d[is.na(d)]=0
 betta=2*sum(d[lower.tri(d, diag = FALSE)])/(nsize*(nsize-1))

PercentageSize= length(intersect(genesGO,genesfitness))/length(genesGO)
vec[1]=Criteria(LDAy,betta,PercentageSize)
vec[2]=LDAy
vec[5]=betta
vec[6]=PercentageSize
#vec[2]=paste(LDAy,betta)
return(vec)
}

	
####Start
if(ReduceSpace) print(paste("Warning! Reducing SPACE prevent to identify some genes in the data."))
if(ReduceSpace) OptimiseInitialisation=FALSE
if(OptimiseInitialisation) ReduceSpace=TRUE
allthefiless99=TrueGenesALL

if(ReduceSpace) #To reduce the genes space
{

allthefiless00=TrueGenesALL
graphALL00=graphALL
genesGO=which(allthefiless00 %in% TrueGenesGO)

ng=length(genesGO)
gg=list(0,ng)
for(i in 1:ng)
gg[[i]]=neighbors(graphALL00,genesGO[i])
tbt=table(unlist(gg))
aa=which(tbt>=mean(tbt))
testgenes1200=noquote(names(aa))
testgenes1200=unique(c(testgenes1200,genesGO))
testgenes1200=as.numeric(testgenes1200)

testgenes=testgenes1200
OtherGenesGo=allthefiless00[testgenes]

TrueGenesALL=unique(c(TrueGenesGO,OtherGenesGo))
genesallthefiless=which(allthefiless00 %in% TrueGenesALL)
graphALLtest=induced.subgraph(graphALL00, genesallthefiless)
graphALL=graphALLtest
allthefiless99=allthefiless00[genesallthefiless]
TrueGenesALL=allthefiless99


}

genesALL=1:length(TrueGenesALL) #position de tous les 20959 genes
genesGO=which(TrueGenesALL %in% TrueGenesGO) 
graphGO=induced.subgraph(graphALL, genesGO)
genesCliqGO=genesGO[largest.cliques(graphGO)[[1]]]
genesAntiCliqGO=genesGO[!is.element(genesGO,genesCliqGO)]
genesAntiGO=genesALL[!is.element(genesALL,genesGO)]

genesAntiCliq=genesALL[!is.element(genesALL,genesCliqGO)]#positions de tous les genes sans cliq
posis=genesAntiCliq
nreste=length(genesAntiCliq) ##number of remaining genes (20959-28)

Tpop = Tpop_A + Tpop_B + Tpop_C

if(aleatoire) #we generate solutions by chance

{	
### choix des variables
pop=matrix(NA,Tpop,6+Pmax)
pop=t(apply(pop,1,choicevar,Pmaxloc=Pmax,Posis=posis,nreste=Pmax)) #Pmaxloc >= nreste
###calcul du critère (fitness)
pop=t(apply(pop,1,fitness))
}

if(!aleatoire)
{
### choix des variables
#Tpop_A=round(Tpop_A*Tpop)
 vecTN=rep(NA,6+Pmax) 
 vecTN[4]=length(genesAntiCliqGO)
 vecTN[7:(6+vecTN[4])]=genesAntiCliqGO
 pop_A=matrix(rep(vecTN,each=Tpop_A),nrow=Tpop_A)
 pop_A=matrix(rep(fitness(pop_A[1,]),each=Tpop_A),nrow=Tpop_A)
###
#Tpop_B=round(Tpop_B*Tpop)
 pop_B<-NULL
 for(PB in 1:Tpop_B)
 {
 vecTN=rep(NA,6+Pmax) 
 nbrgenerestant=length(genesAntiCliqGO) 
 vecTN[4]=round(nbrgenerestant*0.8)
 vecTN[7:(6+vecTN[4])]=sample(genesAntiCliqGO,round(nbrgenerestant*0.8)) #80% of GO genes
 pop_B=rbind(pop_B,vecTN)
 }
 pop_B=t(apply(pop_B,1,fitness))
###
#Tpop_C=round(Tpop_C*Tpop) 
 pop_C<-NULL
 for(PC in 1:Tpop_C)
 {
vecTN=rep(NA,6+Pmax)
nbrgenerestant=Pmax-length(genesGO)
addgenes=min(nbrgenerestant,10)
vecTN[4]=length(genesAntiCliqGO)+addgenes

vecTN[7:(6+vecTN[4])]=c(genesAntiCliqGO,sample(genesAntiGO,addgenes))
pop_C=rbind(pop_C,vecTN)
}
pop_C=t(apply(pop_C,1,fitness))

pop=rbind(pop_A,pop_B,pop_C)
Sabc=Tpop_A+Tpop_B+Tpop_C
pop=pop[1:Sabc,]
popINITIAL=pop
}

if(OptimiseInitialisation)
{
pp=which(allthefiless00 %in% TrueGenesALL)
for(i in 1:Tpop)
{
ccc=7:(6+pop[i,4])
vvv=pp[pop[i,ccc]]
pop[i,ccc]=vvv
}

TrueGenesALL=allthefiless00
graphALL=graphALL00
genesALL=1:length(TrueGenesALL)
genesGO=which(TrueGenesALL %in% TrueGenesGO) 
graphGO=induced.subgraph(graphALL00, genesGO)
genesCliqGO=genesGO[largest.cliques(graphGO)[[1]]]
genesAntiCliqGO=genesGO[!is.element(genesGO,genesCliqGO)]
genesAntiGO=genesALL[!is.element(genesALL,genesGO)]

genesAntiCliq=genesALL[!is.element(genesALL,genesCliqGO)]
posis=genesAntiCliq
nreste=length(genesAntiCliq)
allthefiless99=allthefiless00
}

#####graphique

			# plot(c(0,Ngene),c(-0.05,0.6),type="n")
			# points(0,mean(pop[,1])/(diff(LDAb)*Pmax),pch=16,col="red") #mean of fitness
			# points(0,max(pop[,1])/(diff(LDAb)*Pmax),pch="+",col="red") #max of fitness
			# points(0,mean(pop[,2])/diff(LDAb),pch=16,col="green") #mean of LDA
			# points(0,max(pop[,2])/diff(LDAb),pch="+",col="green") #max of LDA
			# points(0,mean(pop[,4])/Pmax,pch=16,col="blue") #mean of size
			# points(0,mean(pop[,4])/Pmax,pch="+",col="blue") #max of size
			
			ab1=min(pop[1,1:2])-1
			ab2=max(pop[1,1:2])+2
			col = c("red", "green", "blue", "black", "pink") #adjustcolor("green3", alpha.f = 0.1)
			plot(c(0,Ngene),c(ab1,ab2),type="n",ylim=c(-0.5,1),xlab="Number of generations",ylab="Values",cex.lab=1.5)
			points(0,mean(pop[,1]),pch=1,col="red",type = "o") #mean of fitness
			points(0,max(pop[,1]),pch=16,col="red") #max of fitness
			points(0,mean(pop[,2]),pch=1,col="green") #mean of LDA
			points(0,max(pop[,2]),pch=16,col="green") #max of LDA
			points(0,mean(pop[,4]+length(genesCliqGO))/Pmax,pch=1,col="blue") #mean of size
			points(0,max(pop[,4]+length(genesCliqGO))/Pmax,pch=16,col="blue") #max of size

			points(0,mean(pop[,5]),pch=1,col="black") #mean of betta
			points(0,max(pop[,5]),pch=16,col="black") #max of betta
			
			points(0,mean(pop[,6]),pch=1,col="pink") #mean of PercentageSize
			points(0,max(pop[,6]),pch=16,col="pink") #max of PercentageSize


      legend("bottomright", 
             legend = c("BestFitness_max", "ScoreTopo_max", "Size_max","ScoreFunc_max","PercentageSize_max"), 
             col = col, pch = 16, 
             lty = 0,
             pt.cex = c(rep(1,5)), 
             inset = 0.02)
			 
		legend("bottomleft", 
             legend = c("BestFitness_mean", "ScoreTopo_mean", "Size_mean","ScoreFunc_mean","PercentageSize_mean"), 
             col = col, pch = 1, 
             lty = 0,
             pt.cex = c(rep(1,5)), 
             inset = 0.02)
			 
 cat(paste("Number of genes in space =",length(TrueGenesALL))) 
 cat("\n")
 cat(paste("GA | iteration =", 0, "   |    MeanFitness =", noquote(format(round(mean(pop[,1]),3), nsmall = 3)), "   |    BestFitness =", noquote(format(round(max(pop[,1]),3), nsmall = 3))))
 cat("\n")
 
	# else	
		# {
			# plot(c(0,p),c(0,Ngene),type="n")
			# points(pop[which.max(pop[,1]),4:ncol(pop)],rep(0,Pmax),pch="|")
		# }
			
	# cat(c(mean(pop[,1]),max(pop[,1]),mean(pop[,2]),max(pop[,2]),mean(pop[,4]),max(pop[,4]),mean(pop[,5]),max(pop[,5]),mean(pop[,6]),max(pop[,6])),"\n")
    list.pop1=list(0,Ngene)
	list.pop2=list(0,Ngene)
	list.pop4=list(0,Ngene)
	list.pop5=list(0,Ngene)
	list.pop6=list(0,Ngene)
	list.elitisme=list(0,Ngene)
	pop[,3]=0
	indivpos0=1
	elitisme0=pop[indivpos0,]
	gene=1
	popFINAL=NULL
	bestfitness=NULL
	elitismeM1=0
	conver=0


# pop=test$popFINAL
# elitisme0=test$elitisme[[1000]]
	while (gene <= Ngene)
	   {
		t1=Sys.time()
		
############################################################################################################################################	
		##check stopping criteria
			indivpos<-which.max(pop[,1])
			elitisme<-pop[indivpos,]

			ze0=elitisme0[7:(6+elitisme0[4])]
			ze1=elitisme[7:(6+elitisme[4])]
			su=(sum(ze0 %in% ze1)+sum(ze1 %in% ze0))/(length(ze0)+length(ze1))
			if(su==1) { 
			conver=conver+1
			elitisme[3]=elitisme[3]+1
			}
			if(su<1) conver=0
		
			if(conver==conv) #conv=100 par defaut BREAK
			gene=Ngene+1
		
		elitisme0=elitisme
		pop[,3]=pop[,3] + 1 
############################################################################################################################################
		
		if(selection=="Tournament")
		{
		##selection Tournament
		newpop=matrix(NA,0,Pmax+6)
		for ( i in 1:Tpop){
		particip=sample(1:Tpop,tourn,replace=F)
		newpop=rbind(newpop,pop[particip[which.max(pop[particip,1])],])}
		pop<-newpop
		}
		
		if(selection!="Tournament")
		{
		##selection Average
		M<-median(1:Tpop)
		a1<-2/(Tpop*(3*Tpop-4*M+1))
		a2<-2*(Tpop-2*M)/(Tpop*(3*Tpop-4*M+1))
		rang<-rank(pop[,1],ties.method="average")
		proba<-a2+a1*rang
		selec=sample(1:Tpop,Tpop,replace=T,prob=proba)
		pop<-pop[selec,]
		}
############################################################################################################################################	
		
		###crossover
	      	nvpop=pop
			Fbarre=mean(nvpop[,1])
			Fmax=max(nvpop[,1])    
        	seqce=sample(1:Tpop,Tpop)
      		for (i in 1:(length(seqce)/2))
		    {
				soloc=nvpop[seqce[(2*i-1):(2*i)],]
				floc=soloc[,1]
				Fprim=max(floc)
				
				pcloc=1
				if (Fprim < Fbarre) pcloc=pc 
				
				if (Fprim >= Fbarre & Fmax!=Fbarre)
				pcloc= pc * ((Fmax - Fprim)/(Fmax-Fbarre))
				
			cross=runif(1,0,1)
		        if (cross<pcloc)
			    {
					#soloc=nvpop[seqce[(2*i-1):(2*i)],]
					###croisement à un point
					ptscrois<-sample(7:(7+(min(soloc[,4])-1)),1)#!?
					fils1<-soloc[1,]
					fils2<-soloc[2,]
					fils1[7:ptscrois]<-soloc[2,7:ptscrois]
					repdup=duplicated(fils1[7:length(fils1)],incomparables=NA)#positions des duplic. sans les NA.
					fils1[7:length(fils1)]=sort(c(unique(fils1[7:length(fils1)],incomparables=NA),rep(NA,sum(repdup))),na.last=TRUE)
					fils2[7:ptscrois]<-soloc[1,7:ptscrois]
					repdup=duplicated(fils2[7:length(fils2)],incomparables=NA)
					fils2[7:length(fils2)]=sort(c(unique(fils2[7:length(fils2)],incomparables=NA),rep(NA,sum(repdup))),na.last=TRUE)
					###calcul des nouveaux paramètres
					fils1[4]=sum(!is.na(fils1[7:length(fils1)]))
					fils2[4]=sum(!is.na(fils2[7:length(fils2)]))
					nvpop[seqce[(2*i-1):(2*i)],]<-rbind(fils1,fils2)
					nvpop[seqce[(2*i-1):(2*i)],3]=rep(0,2) #migration
			     }
		    }

############################################################################################################################################

		###mutation
		nvpop1=t(apply(nvpop,1,mutation,Posis=posis,Pmaxloc=Pmax,Fbarre=Fbarre,Fmax=Fmax,pm=pm))
		pop<-nvpop1
		pop=t(apply(pop,1,fitness))
		
############################################################################################################################################	
		
		###élitisme
		pop=rbind(pop,elitisme)
		indivpos<-which.max(pop[,1])
		elitisme<-pop[indivpos,]
		pop=pop[1:Tpop,]
		
############################################################################################################################################	    
	
		###Immigration 
		pop=rbind(pop,elitisme) #Ajouté
		agemax=which(pop[,3]>=AgeMax)	
		#print(sprintf("Il y a %d solutions superieures à Agemax",length(agemax)))
		genesElitisme=elitisme[7:(6+elitisme[4])]
		Pimmig=Pmax-length(genesCliqGO)-elitisme[4]
		genesAntiElitisme=genesAntiCliq[!is.element(genesAntiCliq,genesElitisme)]
		if(length(agemax)>0)
		{
		for (i in agemax){	
		#ploc=sample(Pimmig,1)
		ploc=plocImmig
		newgenes=sample(genesAntiElitisme,ploc)
			pop[i,7:(elitisme[4]+ploc+6)]=c(genesElitisme,newgenes)
			pop[i,4]=length(c(genesElitisme,newgenes))
			pop[i,3]=0}			
		
		newpop=pop[agemax,]
		if(length(agemax)==1)
		newpop=as.matrix(t(newpop))
		newpop=t(apply(newpop,1,fitness))
		pop[agemax,]=newpop
		}
############################################################################################################################################
		
		sup= which(pop[,1]==min(pop[,1]))
		#sup=which.min(pop[,1])
		#sup= which(pop[,1]<=mean(pop[,1]))
		for(elit in 1:length(sup))
		pop[sup[elit],]=elitisme
		list.elitisme[[gene]]=elitisme	
		
		
	

		###graphiques
		#list.pop[[gene]]=pop
		list.pop1[[gene]]=mean(pop[,1])
		list.pop2[[gene]]=mean(pop[,2])
		list.pop4[[gene]]=mean(pop[,4])
		list.pop5[[gene]]=mean(pop[,5])
		list.pop6[[gene]]=mean(pop[,6])

			
			points(gene,mean(pop[,1]),pch=1,col="red") #mean of fitness
			points(gene,max(pop[,1]),pch=16,col="red") #max of fitness
			points(gene,mean(pop[,2]),pch=1,col="green") #mean of LDA
			points(gene,max(pop[,2]),pch=16,col="green") #max of LDA
			points(gene,mean(pop[,4]+length(genesCliqGO))/Pmax,pch=1,col="blue") #mean of size
			points(gene,mean(pop[,4]+length(genesCliqGO))/Pmax,pch=16,col="blue") #max of size
		
			points(gene,mean(pop[,5]),pch=1,col="black") #mean of betta
			points(gene,max(pop[,5]),pch=16,col="black") #max of betta
			
			points(gene,mean(pop[,6]),pch=1,col="pink") #mean of PercentageSize
			points(gene,max(pop[,6]),pch=16,col="pink") #max of PercentageSize

		#else points(pop[which.max(pop[,1]),4:ncol(pop)],rep(gene,Pmax),pch="+")
	    #gene=gene+1
				
		# cat("abc...")
# cat(c(gene,mean(pop[,1]),max(pop[,1]),mean(pop[,2]),max(pop[,2]),mean(pop[,4]),max(pop[,4]),mean(pop[,5]),max(pop[,5]),mean(pop[,6]),max(pop[,6])),"\n")
		# cat("xyz...")
		
t2=Sys.time()	

 cat(paste("GA | iteration =", gene, "   |    MeanFitness =", noquote(format(round(mean(pop[,1]),3), nsmall = 3)), "   |    BestFitness =", 
 noquote(format(round(max(pop[,1]),3), nsmall = 3)),"   |    Size =",elitisme[4]+length(genesCliqGO)))
 

# , "| Time =",print(sprintf("%f",t2-t1)))
 cat("\n")
 elitismeM1=elitisme
		# indivpos0=which(pop[,1]==max(pop[,1]))
		# indivpos0=indivpos0[length(indivpos0)]
		#cat(elitisme[1:6],"\n")
		
		#print(sprintf("Time=%f",t2-t1))
		bestfitness=c(bestfitness,elitisme[1])
		popFINAL=pop
		# if(gene %% 10==0)
		# {
		

		#save(outputp,file=paste("AG",ip,".abcd.",repetition,'.RData',sep="")) #essainum=30
		# }
		iteration=gene
		gene=gene+1
		if(difftime(t2,t1, unit="sec")>120) #21/11/2018
	    gene=Ngene+1
}

##Results
#FMprime		
FMprime=elitisme
FMprime=FMprime[-c(1:6)]
FMprime <- FMprime[!is.na(FMprime)]
FMprime=c(genesCliqGO,FMprime)
FMprime=which(allthefiless00 %in% allthefiless99[FMprime])

PositionCorrection=function(value) #Position Correction
{
solution1=value
solution1=solution1[-c(1:6)]
solution1 <- solution1[!is.na(FMprime)]
solution1=which(allthefiless00 %in% allthefiless99[solution1])
return(solution1)
}

#Summary
summaryTopoFunc=function(FMinitial,FMfinal)
{
summaryy=matrix(0,2,4)

#FM
FM=genesGO
grafitness=induced.subgraph(graphALL00, FMinitial)
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
grafitness=induced.subgraph(graphALL00, FMfinal)
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

genesGO=which(allthefiless00 %in% TrueGenesGO) 
fitnessSummary=summaryTopoFunc(FMinitial=genesGO,FMfinal=FMprime)


	    outputp=NULL
		#outputp$popINITIAL=popINITIAL
		#outputp$popFINAL=popFINAL
		outputp$pop1=list.pop1
		outputp$pop2=list.pop2
		outputp$pop4=list.pop4
		outputp$pop5=list.pop5
		outputp$pop6=list.pop6
		outputp$genesCliqGO=genesCliqGO
		outputp$Pmax=Pmax
		outputp$iteration=iteration
		#outputp$elitisme=list.elitisme
		outputp$bestfitness=bestfitness
		outputp$FMprime=FMprime
		outputp$fitnessSummary=fitnessSummary

#print(FMprime)
cat("\n")
cat("\n")
  cat("+-----------------------------------+\n")
  cat("|             TopoFunc     	       |\n")
  cat("+-----------------------------------+\n\n")
  
  cat("+-----------------------------------+\n") 
  cat("\nTopoFunc results: \n")
   cat("\n")
  cat(paste("Population size       = ", Tpop, "\n"))
  cat(paste("Number of generations = ", Ngene, "\n"))
  cat(paste("Iterations=", iteration, "\n"))
  cat(paste("Best Fitness =", round(fitnessSummary[2,4],7), "\n"))

  cat("+-----------------------------------+\n")  
cat("\nTopoFunc summary: \n")
 cat("\n")
print(fitnessSummary)		
return(outputp)
}
