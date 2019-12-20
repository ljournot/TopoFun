#ConstructionModule
ConstructionModule=function(gene_ids, threshold, path)
 {

 gene_ids=unique(gene_ids)
 ng=length(gene_ids) #Number of genes

 #Read the genes files
 setwd(path)
 readata <- lapply(as.character(gene_ids), read.table)
 readata=data.frame(readata)
 v3=c("ID","Mutual_Rank","COR")
 colnames(readata)=rep(v3,ng)

 #Matrix of Relation between genes set.
 matric=NULL
 AMR=NULL
 for(i in 1:ng)
 {
 v12=readata[,c(3*i-2,3*i-1)]
 mat=cbind(gene_ids[i],v12[,c(2,1)])
 mat=mat[which(mat[,2]<=threshold),]

 mat=mat[-1,] #length=804-1=803
 maty=mat[which(mat[,3] %in% gene_ids),] #dim=5,3
 AMR[i]=NA
 if(length(maty)>=3)
 AMR[i]=mean(maty[,2])

 matric=rbind(matric,maty)

 }

 matric=unique(t(apply(matric,1,function(x) c(min(x[1],x[3]),x[2],max(x[1],x[3])))))
 row.names(matric) <- NULL
 colnames(matric)=c("ID","Mutual_Rank","ID")

 #Transform genes ID to genes in the Matrix of Relation
 GN=matric
 GN[,1]=sapply(matric[,1],function(z) which(gene_ids %in% z))
 GN[,3]=sapply(matric[,3],function(z) which(gene_ids %in% z))
GN=unique(t(apply(GN,1,function(x) c(min(x[1],x[3]),x[2],max(x[1],x[3])))))
colnames(GN)=colnames(matric)
row.names(GN) <- NULL

#Adjacency matrix
 adj=GN[,c(1,3)]
 adj=rbind(adj,cbind(1:ng,1:ng))

 #Construction of Graph
 tre=as.vector(t(adj))
 gra=graph(tre,directed=FALSE)
 gra=simplify(gra, remove.loops = TRUE)

 ResultsConstructionModule=NULL
 ResultsConstructionModule$RelationMatrix=matric #Relation Matrix
 ResultsConstructionModule$AverageMutualRank=AMR #Average Mutual Rank
 ResultsConstructionModule$Graph=gra #Graph

 return(ResultsConstructionModule)

 }
 