#This script exemplifies the use of TopoFunc on GO:0006413 ('translational initiation').
#GO:0006413 contains 52 genes whose GeneIDs are in column 27 of the matrix 'matrixGenesGOavant'.

#Loading packages
library(igraph)
library(MASS)

#Loading data
#Replace "path to TopoFunc data" by the path to the folder that stores the TopoFunc data.
setwd("path to TopoFunc data")
load("TrueGenesALLdata.RData")
load("ResultsLDA.RData")
load("graphALLdata.RData")
load("MatriceGeneSim.RData")
load("matrixGenesGOavant.RData")

#Loading functions
#Replace "path to TopoFunc functions" by the path to the folder that stores the TopoFunc functions.
setwd("path to TopoFunc functions")
source("ParametersPertinents.r")
source("LDAnormalization.r")
source("summaryTopoFunc.r")
source("plotTopoFunc.r")
source("fitness_function.r")
source("TopoFunc.r")

#################################################################

#Get the M0 module (=GO:0006413)
ip=27
M0=TrueGenesALLdata[which(matrixGenesGOavant[,ip]==1)]

#Run TopoFunc on the M0 module
x11()
results=TopoFunc(   TrueGenesGO=M0, 
                    graphALL=graphALLdata,
                    TrueGenesALL=TrueGenesALLdata,
                    ReduceSpace=FALSE,
                    OptimiseInitialisation=TRUE,
                    AgeMax=5,
                    Ngene=10,
                    Tpop_A=100,
                    Tpop_B=100,
                    Tpop_C=300,
                    aleatoire=FALSE,
                    Pmax=500,
                    tourn=10,
                    pm=0.5,
                    pc=1,
                    FunctionParam=ParametersPertinents,
                    FunctionLDA=ResultsLDA,
                    distGO=MatriceGeneSim,
                    LDAnorma=LDAnormalization,
                    Criteria=fitness_function,
                    plocImmig=1
                )

#Plot results
x11()
plotTopoFunc(results)
print(results$fitnessSummary) 

#Get the new M1 module
TrueGenesGOprime=TrueGenesALLdata[results$FMprime]