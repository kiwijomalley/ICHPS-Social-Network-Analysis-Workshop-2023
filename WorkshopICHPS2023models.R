## Relational and Peer Association Models of Social Networks ##
## James O'Malley 12/23/2022 ##

library(latentnet)
library(statnet)
library(amen)
library(eigenmodel)
library(sna)
library(coda)
library(numDeriv)

#Get information about in-built data sets
data(package="latentnet")

#Load network ("relational") data
reldir <- scan('ICHPS_PhysPractBin.txt')
nr=sqrt(length(reldir))
reldir=matrix(reldir,ncol=nr,nrow=nr,byrow=T)

#Load node covariate ("attribute") data
covdata <- scan("ICHPS_nodecov.dat",sep="")
covdata <- matrix(covdata,nrow=nr,byrow=T)
colnames(covdata) <- c("id","male","whexpert","pctwom",
                       "numsess","practice","sumhrt")
covdata <- data.frame(covdata)

#Form list of attributes for help with forming network object later
nodecov <- list(male=covdata[,2], whexpert=covdata[,3], pctwom=covdata[,4],
 numsess=covdata[,5], practice=covdata[,6])

#Make categorical variables of attribute variables
nodecov$bcma <- ifelse(nodecov$practice==1,1,0)
nodecov$bima <- ifelse(nodecov$practice==2,1,0)
nodecov$bpp <- ifelse(nodecov$practice==3,1,0)
nodecov$wnhlth <- ifelse(nodecov$practice==4,1,0)
nodecov$numcat <- 1+ifelse(nodecov$numsess>2,1,0)+ifelse(nodecov$numsess>6,1,0)
nodecov$pctcat <- 1+ifelse(nodecov$pctwom>50,1,0)+ifelse(nodecov$pctwom>75,1,0)

#Make network object
pnet <- network(reldir,directed=TRUE,matrixtype="adjacency",
 vertex.attr=nodecov,
 vertex.attrnames=c("male","whexpert","pctwom","numsess","practice",
  "bcma","bima","bpp","wnhlth","numcat","pctcat"))
par(mfrow=c(1,1))
plot(pnet,mode="fruchtermanreingold",displaylabels=T)

## Models of relational data ##

#Comparison of simple ERGM and latent-space models to show equivalence between them
model0 <- ergm(pnet ~ edges + nodematch("male",diff=FALSE))
lmodel0 <- ergmm(pnet ~ nodematch("male",diff=FALSE))

#More advanced ERGM models

#Exponential random graph models
model1 <- ergm(pnet~edges + mutual + nodeocov("whexpert") +
                  nodeocov("pctwom") + nodeocov("numsess") + nodematch("male",diff=F) +
                  nodematch("bcma",diff=F) + nodematch("bima",diff=F) +
                  nodematch("bpp",diff=F) + nodematch("wnhlth",diff=F))
model2 <- ergm(pnet~edges + mutual + nodeocov("whexpert") +
                  nodeocov("pctwom") + nodeocov("numsess") + nodematch("male",diff=T) +
                  nodematch("practice",diff=T))
#Note: Models 1 and 2 differ because the homophily statistics for Model 1 add up the number of nodes
#for which the value of the attribute in both nodes = 0 or both nodes = 1, whereas the separate network 
#statistics for Model 2 add up the number of nodes for which both nodes = k.

#Assess Goodness of Fit of the fitted model for Model 1
model1.gof <- gof(model1~distance+espartners+dspartners+idegree+odegree,
                   control=control.gof.formula(nsim=100),verbose=T)
model1.gof <- gof(model1~idegree,control=control.gof.formula(nsim=100),verbose=T)
par(mfrow=c(2,1))
#plot(model1.gof) #Something seems wrong with this function as GOF is terrible

#More advanced latent space models

#Firstly, consider latent-space model with euclidean distance between latent node positions
lmodel1 <- ergmm(pnet ~ rreceiver+rsender+euclidean(d=3,G=2)+nodematch("male",diff=T))
plot(lmodel1)

#Checking MCMC diagnostics
mcmc.diagnostics(lmodel1)

#Assess Goodness of fit
#For in-degree
lmodel1i.gof <- gof(lmodel1,GOF=~idegree,control=ergmm.control(nsim=100),verbose=T)
plot(lmodel1i.gof)
#For out-degree
lmodel1o.gof <- gof(lmodel1,GOF=~odegree,control=ergmm.control(nsim=100),verbose=T)
plot(lmodel1o.gof)

#Now use bilinear model: adds products of latent positions as a latent predictor
bmodel1 <- ergmm(pnet ~ rreceiver+rsender+bilinear(d=3,G=2)+nodematch("male",diff=T))
plot(bmodel1)

#Assess Goodness of fit
bmodel1i.gof <- gof(bmodel1,GOF=~idegree,control=ergmm.control(nsim=100),verbose=T)
plot(bmodel1i.gof)
bmodel1o.gof <- gof(bmodel1,GOF=~odegree,control=ergmm.control(nsim=100),verbose=T)
plot(bmodel1o.gof)

## Cross-sectional social influence analyses using functions in SNA package ##

#First, standardize adjacency matrix so that row sums = 1
numalters <- apply(reldir,1,sum) #Compute row-sums (the physicians' outdegrees)
scale=as.vector(numalters^(-1))
noalters <- ifelse(is.infinite(scale)==1,1,0) #Indicator of physician having an out-degree of 0 (an outbound isolate)
scale[noalters*seq(1,nr)] <- 1/(nr-1) #Get rid of infinity for isolates; has no impact if isoequal=0
on=matrix(1,ncol=nr,nrow=1)
wtreldir <- reldir * (scale %*% on)

#Compute mean value of dependent variable (hormone replacement therapy, hrt) for directly connected
# actors - i.e., based on the adjacency matrix
hrtalt <- wtreldir %*% as.vector(covdata$sumhrt)
regdata <- data.frame(covdata,hrtalt=hrtalt,noalters=noalters,numalters=numalters)

#Setup contrasts for predictor variables
on <- as.vector(rep(1,nr))
x <- as.matrix(cbind(on,regdata[,c("male","pctwom","numalters")]))

#Simplistic but incorrect estimation of outcome and network autocorrelation models
reg.adj <- lm(sumhrt~x+hrtalt-1, data=regdata)

#Using lnam package
lnam1.adj <- lnam(regdata$sumhrt,x,wtreldir)
lnam2.adj <- lnam(regdata$sumhrt,x,NULL,wtreldir)

#Print out all results
summary(reg.adj) #Naive regression
summary(lnam1.adj) #Network autocorrelated outcome model
summary(lnam2.adj) #Network autocorrelation model
