## Plots and Descriptive Statistics of Social Networks ##
## Using physician network in the grandpa package (Bobak at al, 2022) ##
## James O'Malley 12/23/2022 ##

#Install relevant packages: For the first time you want to use them
firsttime <- 0
if (firsttime==1) {
 options(timeout = max(1000, getOption("timeout")))
 install.packages("devtools")
 install.packages("stringi")
 devtools::install_github("https://github.com/CarlyBobak/grandpa")
}

#Link relevant libraries
library(devtools)
library(stringi)
library(grandpa)
library(igraph) #Caution: Masks functions in SNA
library(latentnet)
library(statnet)
library(amen)
library(eigenmodel)
library(coda)

#Learn about grandpa package
?physNet #An undirected simulated physician network in grandpa

#Attribute data
covdata <- vertex_attr(physNet)
#Label1 = Speciality of physician
#CommunityLabel = Hospital-like latent cluster of nodes determined from the data
#linchLabel = Measure of node-uniqueness is linking communities
#...

## Visualization and Descriptive Analyses of physNet network ##

#Network object comes with in-built color and size attributes for plotting
plot(physNet)

#However, we can reset these and make a new plot
V(physNet)$color2 <- as.numeric(covdata$CommunityLabel) #Color nodes by community membership
V(physNet)$size2 <- 5*as.numeric(covdata$linchLabel) #Size nodes by linchLabel status
par(mar=c(0,0,0,0))
plot(physNet,
     vertex.color = V(physNet)$color, #Change color of nodes by exchanging color with color2
     vertex.size = V(physNet)$size, #Change size of nodes by exchanging size with size2
     vertex.label.color = "black", # change color of labels
     vertex.label.cex = .25, # change size of labels to 25% of original size
     edge.curved=.25, # add a 25% curve to the edges
     edge.color="grey20", # change edge color to grey
     edge.arrow.size=0.3)
dev.copy2pdf(file="PhysNet_grandpa.pdf", width=6, height=6) #Output plot of network

par(mar=c(0,0,0,0))
plot(physNet,
     vertex.color = V(physNet)$color2, # change color of nodes
     vertex.size = V(physNet)$size2,
     vertex.label.color = "black", # change color of labels
     vertex.label.cex = .25, # change size of labels to 25% of original size
     edge.curved=.25, # add a 25% curve to the edges
     edge.color="grey20", # change edge color to grey
     edge.arrow.size=0.3)
dev.copy2pdf(file="PhysNet_grandpa2.pdf", width=6, height=6) #to file

#summarizing the network
size <- length(physNet$name) #Number of physicians in network (608)
degree <- degree(physNet) #There is a long degree as network is undirected
centralization <- centralize(degree,theoretical.max=608,normalized=TRUE)
trans <- transitivity(physNet) #Transitivity for network as a whole
loc_trans <- transitivity(physNet,type="local") #Get node specific transitivity
cent_bet <- betweenness(physNet) #Betweenness centrality
cent_eig <- eigen_centrality(physNet)$value #Eigenvector centrality; multiple components might cause problems?
cent_close <- closeness(physNet) #Closeness centrality
cent_bonpow <- power_centrality(physNet) #Bonachich power centrality
triadc <- motifs(physNet,size=3)

## Relational Models estimated on physNet network ##

#Deconstruct network object into key components to define object for estimating models
netedges <- as_edgelist(physNet)
covdata <- vertex_attr(physNet) #Labels as above
covdatared <- list(Community=as.numeric(covdata$CommunityLabel),
                   linchphys=as.numeric(covdata$linchLabel)) #Reduced set of attributes for illustration
covnames <- names(covdatared)
covdatamat <- matrix(cbind(covdatared$Community,covdatared$linchphys),ncol=2,byrow=TRUE)
covdatamat <- data.frame(covdatamat)
names(covdatamat) <- c("Community","Linchphys")

#Make network object using network function (tailored for ERGM and latent-space models)
pnet <- network(netedges,directed=FALSE,matrix.type="edgelist",
                vertex.attr=covdatared,vertex.attrnames=covnames)
par(mfrow=c(1,1))
plot(pnet,mode="fruchtermanreingold",displaylabels=T) #Rough plot (not as nice as igraph plots)

#Comparison of simple ERGM and latent-space models (emulating illustration in slides but for a larger network)
model0 <- ergm(pnet ~ edges + nodematch("Community",diff=FALSE))
lmodel0 <- ergmm(pnet ~ nodematch("Community",diff=FALSE)) #latent-space model equivalent to the ERGM on above line

model1 <- ergm(pnet~edges + nodecov("linchphys") +
                 nodecov("Community") + nodematch("linchphys",diff=F) + 
                 nodematch("Community",diff=F))

model2 <- ergm(pnet~edges + nodecov("linchphys") +
                 nodecov("Community") + nodematch("linchphys",diff=T) + 
                 nodematch("Community",diff=F))

#model3 <- ergm(pnet~edges + triangles + nodecov("linchphys") +
#                 nodecov("Community") + nodematch("linchphys",diff=F) + 
#                 nodematch("Community",diff=F))
#Takes a very long time to run. Solution is predicted to be degenerate!

lmodel2 <- ergmm(pnet ~ nodecov("linchphys") +
                   nodecov("Community") + nodematch("linchphys",diff=T) + 
                   nodematch("Community",diff=F))
#Takes a very long time to run!

#Checking MCMC diagnostics
mcmc.diagnostics(lmodel2)

#Assessment of Goodness of fit
lmodel2.gof <- gof(lmodel2,GOF=~idegree,control=ergmm.control(nsim=100),verbose=T)
plot(lmodel2.gof)

## Peer Association Models: Cross-sectional analyses using functions in SNA ##

#For illustration evaluate peer-association of community membership (this is tautological, but useful for illustration)

#First, form adjacency matrix
nr <- size
relmut <- matrix(0,nrow=nr,ncol=nr)
for (i in 1:ncol(edges)) {
 relmut[edges[i,1],edges[i,2]] <- 1
 relmut[edges[i,2],edges[i,1]] <- 1
}

#Scale rows of adjacency matrix so that row sums = 1
numalters <- apply(relmut,1,sum) #Compute row sums (the physicians' outdegrees)
scale=as.vector(numalters^(-1))
noalters <- ifelse(is.infinite(scale)==1,1,0) #Indicate of physician having degree of 0 (an isolate)
scale[noalters*seq(1,nr)] <- 1/(nr-1) #Get's rid of infinity for isolates; has not impact if isoequal=0
on=matrix(1,ncol=nr,nrow=1)
wtrelmut <- relmut * (scale %*% on)

#Compute mean value of dependent variable (community membership) for directly connected
# actors - i.e., based on the adjacency matrix
Communityalt <- wtrelmut %*% as.vector(covdatamat$Community)
regdata <- data.frame(covdatamat,Communityalt=Communityalt,noalters=noalters,numalters=numalters)

#Setup contrasts for predictor variables
on <- as.vector(rep(1,nr))
x <- as.matrix(cbind(on,regdata[,c("numalters")]))

#Simplistic but incorrect estimation of outcome and network autocorrelation models
reg.adj <- lm(Community~Communityalt+numalters-1, data=regdata)

#Using lnam package
lnam1.adj <- lnam(regdata$Community,x,wtrelmut)
lnam2.adj <- lnam(regdata$Community,x,NULL,wtrelmut)

#Print out all results
summary(reg.adj) #Naive regression
summary(lnam1.adj) #Network autocorrelated outcome model
summary(lnam2.adj) #Network autocorrelation model
