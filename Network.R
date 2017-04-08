otu<-read.table("otu.txt",,header=TRUE,sep="\t",row.names=1)
otu<-t(otu)
a<-otu
ra<-a[,colSums(a)/nrow(a)<0.00001]
abd<-a[,colSums(a)/nrow(a)>0.0005]
b<-a[,colSums(a)/nrow(a)>=0.00001]
int<-b[,colSums(b)/nrow(b)<=0.0005]

######generate Core matrix#####
h<-c(1:ncol(otu));f<-0
for (i in 1:ncol(otu))
{b<-otu[,i];c<-sum(b)
f[i]<-length(b[b>0])
if (f[i]==0){h[i]=0} 
else{h[i]=c/f[i]}}
tt<-t(rbind(h,f))
rownames(tt)<-colnames(otu)
b<-rownames(subset(tt,f>nrow(otu)*0.8))
core<-otu[,b]

library(Hmisc)
library(igraph)
library(agricolae)
library(fdrtool)
########network######
b<-rcorr(otu,type="spearman")
r<-b$r;p<-b$P
diag(p)<-1
pp<-as.vector(p)
qvalue<-fdrtool(pp,statistic="pvalue")
q<-qvalue$qval
q<-matrix(q,ncol(p),ncol(p))
q[q>0.01]<-0
q[q<=0.01&q>0]<-1
r[r<0.6]<-0
g<-r*q
g <- graph.adjacency(g, weighted=TRUE, mode="undirected")
g <- simplify(g)
g<-delete.vertices(g,names(degree(g)[degree(g)==0]))
d<-degree(g)
write.graph(g,"whole_network.gml", format="gml")#####generate network by Gephi#######

###### Topological features at node level#######
dd<-matrix(NA,ncol=5,nrow=length(d))
dd<-as.data.frame(dd)
rownames(dd)<-names(d)
colnames(dd)<-c("degree","betweenness","closeness","eigenvector","category")
dd[,1]<-d
dd[,5]<-"other"
dd[intersect(names(d),colnames(ra)),5]<-"rare"
dd[intersect(names(d),colnames(abd)),5]<-"abundant"
dd[intersect(names(d),colnames(int)),5]<-"intermediate"
btw<-betweenness(g);cls<-closeness(g);egv<-evcent(g)
dd[,2]<-btw;dd[,3]<-cls;dd[,4]<-egv$vector

dd_core<-matrix(NA,ncol=5,nrow=length(d))
dd_core<-as.data.frame(dd_core)
rownames(dd_core)<-names(d)
colnames(dd_core)<-c("degree","betweenness","closeness","eigenvector","category")
dd_core[,1]<-d
dd_core[,5]<-"other"
dd_core[intersect(names(d),colnames(core)),5]<-"core"
btw<-betweenness(g);cls<-closeness(g);egv<-evcent(g)
dd_core[,2]<-btw;dd_core[,3]<-cls;dd_core[,4]<-egv$vector
######Non-parametric test######
###rare vs.abundant####
dd<-subset(dd,category!="intermediate")
a<-with(dd,wilcox.test(degree~category))
fdrtool(a$p.value,statistic="pvalue")
a<-with(dd,kruskal(betweenness,category))
fdrtool(a$p.value,statistic="pvalue")
a<-with(dd,kruskal(closeness,category))
fdrtool(a$p.value,statistic="pvalue")
a<-with(dd,kruskal(eigenvector,category))
fdrtool(a$p.value,statistic="pvalue")
pdf("node_level boxplot.pdf")
par(mfrow=c(2,2))
with(dd,boxplot(degree~category,ylab="Degree",las=1,
col=c("darkorange","cornflowerblue")))
with(dd,boxplot(betweenness~category,ylab="Betweenness",las=1,
col=c("darkorange","cornflowerblue"),ylim=c(0,300000)))
with(dd,boxplot(closeness~category,ylab="Closeness",las=1,
col=c("darkorange","cornflowerblue"),ylim=c(0.0000024,0.0000027)))
with(dd,boxplot(eigenvector~category,ylab="Eigenvector",las=1,
col=c("darkorange","cornflowerblue")))
dev.off()
###core vs.other####
a<-with(dd_core,wilcox.test(degree~category))
fdrtool(a$p.value,statistic="pvalue")
a<-with(dd_core,kruskal(betweenness,category))
fdrtool(a$p.value,statistic="pvalue")
a<-with(dd_core,kruskal(closeness,category))
fdrtool(a$p.value,statistic="pvalue")
a<-with(dd_core,kruskal(eigenvector,category))
fdrtool(a$p.value,statistic="pvalue")
pdf("node_level boxplot-1.pdf")
par(mfrow=c(2,2))
with(dd_core,boxplot(degree~category,ylab="Degree",las=1,
col=c("darkorange","cornflowerblue")))
with(dd_core,boxplot(betweenness~category,ylab="Betweenness",las=1,
col=c("darkorange","cornflowerblue"),ylim=c(0,300000)))
with(dd_core,boxplot(closeness~category,ylab="Closeness",las=1,
col=c("darkorange","cornflowerblue"),ylim=c(0.0000024,0.0000027)))
with(dd_core,boxplot(eigenvector~category,ylab="Eigenvector",las=1,
col=c("darkorange","cornflowerblue")))
dev.off()
###### Topological features at network level#######
a<-subset(dd,category=="rare")
g_ra<-subgraph(g,rownames(a))
a<-subset(dd,category=="abundant")
g_abd<-subgraph(g,rownames(a))
a<-subset(dd_core,category=="core")
g_core<-subgraph(g,rownames(a))
a<-subset(dd_core,category=="other")
g_other<-subgraph(g,rownames(a))

net<-matrix(NA,ncol=6,nrow=4)
net<-as.data.frame(net)
net[1,]<-c(mean(degree(g_abd)),transitivity(g_abd),average.path.length(g_abd),
graph.density(g_abd),diameter(g_abd),modularity(walktrap.community(g_abd)))
net[2,]<-c(mean(degree(g_ra)),transitivity(g_ra),average.path.length(g_ra),
graph.density(g_ra),diameter(g_ra),modularity(walktrap.community(g_ra)))
net[3,]<-c(mean(degree(g_core)),transitivity(g_core),average.path.length(g_core),
graph.density(g_core),diameter(g_core),modularity(walktrap.community(g_core)))
net[4,]<-c(mean(degree(g_other)),transitivity(g_other),average.path.length(g_other),
graph.density(g_other),diameter(g_other),modularity(walktrap.community(g_other)))
colnames(net)<-c("AD","Cc","Apl","density","diameter","modula")
rownames(net)<-c("Abundant","Rare","Core","Other")
