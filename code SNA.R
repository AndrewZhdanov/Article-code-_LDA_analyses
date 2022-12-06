library('rmarkdown')
library('RColorBrewer')
library('network')
library('sna')
library('igraph')
library('readxl')
library(tidyr)
library(intergraph)

library(ergm)
library(coda)
library(RColorBrewer)
library(rgl)
library(sna)

library(leaps)
library(stargazer)


SWFile <- "/Candidates/bidentrumpsanders.tsv" #this is the name of the file
wd<-getwd()
SWFilePath <- paste(wd, SWFile, sep = "")
SWrawdata<-read.table(SWFilePath, sep = "\t", header = TRUE, row.names = 1) #цветные не распознаёт
#SWrawdata

data <- SWrawdata
data[is.na(data)] <- 0
#data
data <- t(data)

SWnet<-graph_from_incidence_matrix(data)
#par(mar=c(0,0,0,0))
#plot(SWnet)


colors<-brewer.pal(8, 'Accent') #select and load a palette
# Next, change attributes we want.
# Doing it one at a time, as we are sure you've figured out,
# is much easier than inside the plot function
V(SWnet)$color <- c(colors[1],colors[6])[V(SWnet)$type+1] #set two different colors for two node types
V(SWnet)$shape <- c("square", "circle")[V(SWnet)$type+1] #change shapes based on node types
V(SWnet)$label.color<-c("black", "transparent")[V(SWnet)$type+1]
V(SWnet)$label.cex<-c(0.5, 0.7)[V(SWnet)$type+1]
V(SWnet)$label.font=2
# Calculate the indegree of events:
V(SWnet)$indegree <- degree(SWnet, mode = "in") #calculate indegree

V(SWnet)$size<-ifelse(V(SWnet)$type==TRUE,V(SWnet)$indegree*0,000015) #size of circle

#square label
for (i in 0:98){
  V(SWnet)$label.cex[i]<-V(SWnet)$indegree[i]*0.005
}

#square size
for (i in 0:98){
  V(SWnet)$size[i]<-V(SWnet)$indegree[i]*0.05
}


V(SWnet)$label.cex[1] <-  0.9
V(SWnet)$label.cex[22] <-  0.9
V(SWnet)$label.cex[52] <-  0.9

##################
####Appendix 4 (1)
##################

par(mar=c(0,0,0,0))
plot(SWnet)

df <- cbind(newColName = rownames(data), data)
rownames(df) <- 1:nrow(df)


char_mes <- read_excel('characteristic of all messenges.xlsx')

SWnet.sep<-bipartite.projection(SWnet)
#par(mar=c(0,0,0,0), mfrow=c(1,2))
#plot(SWnet.sep$proj1)
#plot(SWnet.sep$proj2) 


q <- SWnet.sep$proj2
edge.attributes(q)$weight
q <- delete_edge_attr(q, "weight")
edge.attributes(q)$weight
messenges<-asNetwork(q)

#par(mar=c(0,0,0,0))
#plot(messenges)


messenges %v% "char_mes" <- c(char_mes$type)
messenges

###########
###Figure 1
###########


## All Candidates
flomodel.01 <- ergm(messenges~edges+nodecov('char_mes')) #чем выше враждебность, тем большей связей
summary(flomodel.01)

messenges %v% "low" <- c(char_mes$low)

flomodel.02 <- ergm(messenges~edges+nodecov('low')) #нейтральные сообщения уменьшают связи
summary(flomodel.02)

messenges %v% "medium" <- c(char_mes$medium)

flomodel.03 <- ergm(messenges~edges+nodecov('medium')) #средней враждебности сообщения уменьшают связи
summary(flomodel.03)

messenges %v% "hard" <- c(char_mes$hard)

flomodel.04 <- ergm(messenges~edges+nodecov('hard')) #средней враждебности сообщения уменьшают связи
summary(flomodel.04)

# Biden 

SWFile <- "/Candidates/biden.tsv" #this is the name of the file
wd<-getwd()
SWFilePath <- paste(wd, SWFile, sep = "")
SWrawdata<-read.table(SWFilePath, sep = "\t", header = TRUE, row.names = 1) #цветные не распознаёт
SWrawdata

data <- SWrawdata
data[is.na(data)] <- 0
data
data <- t(data)

SWnet<-graph_from_incidence_matrix(data)

SWnet.sep<-bipartite.projection(SWnet)
q <- SWnet.sep$proj2
edge.attributes(q)$weight
q <- delete_edge_attr(q, "weight")
edge.attributes(q)$weight
messenges_b<-asNetwork(q)

par(mar=c(0,0,0,0))
#plot(messenges_b)
messenges_b

length(messenges_b)

char_mes_b <- read_excel('characteristic of biden messenges.xlsx')
char_mes_b$type

#flomodel.01 <- ergm(messenges_b ~ edges)
#summary(flomodel.01)

messenges_b %v% "char_mes_b" <- c(char_mes_b$type)

flomodel.05 <- ergm(messenges_b~edges+nodecov('char_mes_b')) 
summary(flomodel.05)

messenges_b %v% "low" <- c(char_mes_b$low)

flomodel.06 <- ergm(messenges_b~edges+nodecov('low')) 
summary(flomodel.06)

messenges_b %v% "medium" <- c(char_mes_b$medium)

flomodel.07 <- ergm(messenges_b~edges+nodecov('medium'))
summary(flomodel.07)

messenges_b %v% "hard" <- c(char_mes_b$hard)

flomodel.08 <- ergm(messenges_b~edges+nodecov('hard'))
summary(flomodel.08)

stargazer(flomodel.01, flomodel.02, flomodel.03, flomodel.04, 
          flomodel.05, flomodel.06, flomodel.07, flomodel.08, type = 'text')


###########
###Figure 2
###########

#Trump 
SWFile <- "/Candidates/trump1.tsv" #this is the name of the file
wd<-getwd()
SWFilePath <- paste(wd, SWFile, sep = "")
SWrawdata<-read.table(SWFilePath, sep = "\t", header = TRUE, row.names = 1) #цветные не распознаёт + нужны разные строки
#SWrawdata

data <- SWrawdata
data[is.na(data)] <- 0
#data
data <- t(data)

SWnet<-graph_from_incidence_matrix(data)
#par(mar=c(0,0,0,0))
#plot(SWnet)


SWnet.sep<-bipartite.projection(SWnet)
q <- SWnet.sep$proj2
edge.attributes(q)$weight
q <- delete_edge_attr(q, "weight")
edge.attributes(q)$weight
messenges_t<-asNetwork(q)

#par(mar=c(0,0,0,0))
#plot(messenges_t)
messenges_t

length(messenges_t)

char_mes_t <- read_excel('characteristic of trump messenges.xlsx')
#char_mes_t$type

flomodel.01 <- ergm(messenges_t ~ edges)
summary(flomodel.01)

messenges_t %v% "char_mes_t" <- c(char_mes_t$type)

flomodel.01 <- ergm(messenges_t~edges+nodecov('char_mes_t')) #чем выше враждебность, тем большей связей
summary(flomodel.01)

messenges_t %v% "low" <- c(char_mes_t$low)

flomodel.02 <- ergm(messenges_t~edges+nodecov('low')) #нейтральные сообщения уменьшают связи
summary(flomodel.02)

messenges_t %v% "medium" <- c(char_mes_t$medium)

flomodel.03 <- ergm(messenges_t~edges+nodecov('medium')) #средней враждебности сообщения уменьшают связи
summary(flomodel.03)

messenges_t %v% "hard" <- c(char_mes_t$hard)

flomodel.04 <- ergm(messenges_t~edges+nodecov('hard')) #средней враждебности сообщения уменьшают связи
summary(flomodel.04)

#Sanders

SWFile <- "/Candidates/sanders1.tsv" #this is the name of the file
wd<-getwd()
SWFilePath <- paste(wd, SWFile, sep = "")
SWrawdata<-read.table(SWFilePath, sep = "\t", header = TRUE, row.names = 1) #цветные не распознаёт + нужны разные строки
SWrawdata

data <- SWrawdata
data[is.na(data)] <- 0
data
data <- t(data)

SWnet<-graph_from_incidence_matrix(data)
#par(mar=c(0,0,0,0))
#plot(SWnet)



q <- SWnet.sep$proj2
edge.attributes(q)$weight
q <- delete_edge_attr(q, "weight")
edge.attributes(q)$weight
messenges_s<-asNetwork(q)

#par(mar=c(0,0,0,0))
#plot(messenges_s)
messenges_s

char_mes_s <- read_excel('characteristic of sanders messenges.xlsx')
char_mes_s$type

#flomodel.05 <- ergm(messenges_s ~ edges)
#summary(flomodel.05)

messenges_s %v% "char_mes_s" <- c(char_mes_s$type)

flomodel.05 <- ergm(messenges_s~edges+nodecov('char_mes_s'))
summary(flomodel.05)

messenges_s %v% "low" <- c(char_mes_s$low)

flomodel.06 <- ergm(messenges_s~edges+nodecov('low')) 
summary(flomodel.06)

messenges_s %v% "medium" <- c(char_mes_s$medium)

flomodel.07 <- ergm(messenges_s~edges+nodecov('medium')) 
summary(flomodel.07)

#Bernie Sanders didn't have a single message with the greatest affect on the campaign trail 
messenges_s %v% "hard" <- c(char_mes_s$hard)

flomodel.08 <- ergm(messenges_s~edges+nodecov('hard')) 
summary(flomodel.08)


stargazer(flomodel.01, flomodel.02, flomodel.03, flomodel.04, 
          flomodel.05, flomodel.06, flomodel.07, type = 'text')
