## Plots and Descriptive Statistics of Social Networks ##
## James O'Malley 12/23/2022 ##

#Link relevant libraries
library(igraph) #Caution: Masks functions in SNA
library(sna)

#Load network ("relational") data
reldir <- scan("ICHPS_PhysPractBin.txt")
nr <- sqrt(length(reldir)) #Number of physicians
reldir <- matrix(reldir,ncol=nr,nrow=nr,byrow=T)

#Form binary directed network
pnet <- network(reldir,directed=TRUE,matrixtype="adjacency")
par(mfrow=c(1,1))
gplot(pnet,gmode="digraph",mode="fruchtermanreingold",displaylabels=T)

#Node 2 is an outward tie isolate, but is quite popular (in-degree of 9)
reldir[,2]
reldir[2,]

#Binary undirected network
relmut <- ifelse(reldir+t(reldir)>0,1,0) #Can use either the "or" or the "and" rules
mutnet <- network(relmut,directed=FALSE,matrixtype="adjacency")
gplot(mutnet,gmode="graph",mode="fruchtermanreingold",displaylabels=T) #Use graph instead of digraph for undirected

#Aside: Form edgelist from adjacency matrix of directed network: useful for using igraph
idto <- seq(1,nr)
edgelist=c()
for (i in 1:nr) {
  edge <- idto[(reldir[i,]==1)]
  nc=length(edge)
  edgelist=rbind(edgelist,cbind(rep(i,nc),edge))
}
edgelist <- data.frame(edgelist)
names(edgelist) <- c("source","target")
write.csv(edgelist,file="ICHPS_ClinEdgeList.csv",row.names=FALSE)

## If network already stored in edgelist form ##
#reldir <- read.csv("ICHPS_ClinEdgeList.csv")
#pnet <- network(reldir,directed=TRUE,matrixtype="edgelist")

#Load node covariate ("attribute") data
covdata <- scan("ICHPS_nodecov.dat")
covdata <- matrix(covdata,nrow=nr,byrow=T)
colnames(covdata) <- c("id","male","whexpert","pctwom",
                       "numsess","practice","sumhrt")
covdata <- data.frame(covdata)
nodecov <- list(male=covdata[,2], whexpert=covdata[,3], pctwom=covdata[,4],
                numsess=covdata[,5], practice=covdata[,6])

#Make categorical variables of the attribute variables
nodecov$bcma <- ifelse(nodecov$practice==1,1,0)
nodecov$bima <- ifelse(nodecov$practice==2,1,0)
nodecov$bpp <- ifelse(nodecov$practice==3,1,0)
nodecov$wnhlth <- ifelse(nodecov$practice==4,1,0)
nodecov$numcat <- 1+ifelse(nodecov$numsess>2,1,0)+ifelse(nodecov$numsess>6,1,0)
nodecov$pctcat <- 1+ifelse(nodecov$pctwom>50,1,0)+ifelse(nodecov$pctwom>75,1,0)
nodecov<-data.frame(nodecov)

#igraph plots of network
nodes <- unique(c(edgelist[,1],edgelist[,2]))
gnet <- graph_from_data_frame(d=edgelist, vertices=nodes, directed=TRUE) #Make graph object
print(gnet, e=TRUE, v=TRUE)

V(gnet)$color <- covdata$whexpert+1
par(mar=c(0,0,0,0))
plot(gnet,
     vertex.color = V(gnet)$color, # change color of nodes
#     vertex.size = V(gnet)$size,
     vertex.label.color = "black", # change color of labels
     vertex.label.cex = .75, # change size of labels to 75% of original size
     edge.curved=0.25, # add a 25% curve to the edges
     edge.color="grey20", # change edge color to grey
     edge.arrow.size=0.3)
dev.copy2pdf(file="PhysNetPractPlot.pdf", width=6, height=6) #to file

V(gnet)$color <- covdata$practice
V(gnet)$size <- 10*(covdata$whexpert+1)
par(mar=c(0,0,0,0))
plot(gnet,
     vertex.color = V(gnet)$color, # change color of nodes
     vertex.size = V(gnet)$size,
     vertex.label.color = "black", # change color of labels
     vertex.label.cex = .75, # change size of labels to 75% of original size
     edge.curved=.25, # add a 25% curve to the edges
     edge.color="grey20", # change edge color to grey
     edge.arrow.size=0.3)
dev.copy2pdf(file="PhysNetNicePlot.pdf", width=6, height=6) #to file

reldir[,9]
reldir[9,]

## Descriptive analyses using SNA: Directed network ##
detach("package:igraph", unload = TRUE) #Unload igraph as it masks functions in SNA
#Or when calling functions, precede function name with package name

#Dyad and Triad census
dyadc=dyad.census(reldir)
recip=grecip(reldir,measure="edgewise") #Proportion of reciprocated ties from dyad census (=2*26/(2*26+111))
recipd=grecip(reldir,measure="dyadic") #Proportion of dyads that are symmetric from dyad census (=(26+391)/(26+111+391))
triadc_mut=triad.census(relmut,g="graph")
triadc_dir=triad.census(reldir,g="digraph")
trans=gtrans(reldir,mode="digraph")

#Degree distributions of directed network
idegree=degree(reldir,cmode="indegree")
odegree=degree(reldir,cmode="outdegree")
central=centralization(reldir,degree)

#Centrality measures
reldir <- t(reldir) #Take transpose if want centrality for directed network to be based on inbound ties 
closecent=closeness(reldir,gmode="digraph")
bcent=betweenness(reldir,gmode="digraph")
eigcent=evcent(reldir,gmode="digraph",use.eigen=FALSE)
powcent=bonpow(reldir,gmode="digraph",exponent=1,rescale=FALSE) #exponent is our alpha parameter
powcent=bonpow(reldir,gmode="digraph",exponent=0,rescale=FALSE) #Get measure proportional to degree

#Cycle and path census' (statistics beyond those described in the slides)
#reldir <- t(reldir) #Transpose back
kcycle=kcycle.census(reldir)
kpath=kpath.census(reldir)

#Connectedness
connectstr=is.connected(reldir,connected="strong")
connectwk=is.connected(reldir,connected="weak")


