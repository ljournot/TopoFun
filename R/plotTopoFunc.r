#plot TopoFunc
plotTopoFunc=function(resAG)
{
col=c("red","green","blue","black","pink")
plot(x=1:resAG$iteration,y=resAG$bestfitness,col=col[1],type ="o",pch=16,ylim=c(-0.5,1),xlab="Number of generations",ylab="Values",cex.lab=1.5)
points(x=1:resAG$iteration,y=resAG$pop2,col=col[2],type ="o",pch=16)
points(x=1:resAG$iteration,y=lapply(resAG$pop4,function(x) (x+length(resAG$genesCliqGO))/resAG$Pmax),col=col[3],type ="o",pch=16)
points(x=1:resAG$iteration,y=resAG$pop5,col=col[4],type ="o",pch=16)
points(x=1:resAG$iteration,y=resAG$pop6,col=col[5],type ="o",pch=16)

      legend("bottomright", 
             legend = c("BestFitness", "ScoreTopo", "Size","ScoreFunc","PercentageSize"), 
             col = col, pch = 16, 
             lty = 1,
             pt.cex = c(rep(1,5)), 
             inset = 0.02)
			 
}

#plot.TopoFunc(resAG=test)
