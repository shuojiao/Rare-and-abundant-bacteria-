otu<-read.table("otu.txt",,header=TRUE,sep="\t",row.names=1)
env<-read.table("env.txt",,header=TRUE,sep="\t",row.names=1)####environmental factors######
a<-otu
ra<-a[,colSums(a)/nrow(a)<0.00001]
abd<-a[,colSums(a)/nrow(a)>0.0005]
b<-a[,colSums(a)/nrow(a)>=0.00001]
int<-b[,colSums(b)/nrow(b)<=0.0005]

library(ggplot2)
library(geosphere)
library(fdrtool)
library(vegan)
library(TSA)
library(gplots)
library(VennDiagram)

##############Distance decay relationship###########
col<-c("red2","darkorange","grey92","cornflowerblue")
aa<-distm(env[,15:14])#####location####
df<-data.frame(as.vector(1-vegdist(core)),
as.vector(1-vegdist(abd)),
as.vector(1-vegdist(int)),
as.vector(1-vegdist(ra)),
as.vector(as.dist(aa)/100000))
names(df)<-c("core","abundant","intermediate","rare","geodis")
pdf("DDR.pdf",height=6,width=8)
par(mar=c(5,5,2,2))
plot(df$abundant~df$geodis,type="n",cex.lab=1.5,cex.axis=1.3,ylim=c(0,1.1),
xlab="Geographic distance",ylab="Community similarity")
for(i in 1:3)
{reg<-lm(df[,i]~df$geodis)
abline(reg,col=col[i],lwd=3)}
for(i in 1:3)
{points(df[,i]~df$geodis,pch=20,cex=1.2,col=col)}
box(lwd=2)
dev.off()

#########################################
#########partioning beta########
library(betapart)
a1<-bray.part(ra)
a2<-bray.part(abd)
a3<-bray.part(int)
a4<-bray.part(core)
a<-bray.part(otu)
tu1<-a1$bray.bal/a1$bray
rh1<-a1$bray.gra/a1$bray
br1<-a1$bray
tu2<-a2$bray.bal/a2$bray
rh2<-a2$bray.gra/a2$bray
br2<-a2$bray
tu3<-a3$bray.bal/a3$bray
rh3<-a3$bray.gra/a3$bray
br3<-a3$bray
tu4<-a4$bray.bal/a4$bray
rh4<-a4$bray.gra/a4$bray
br4<-a4$bray
tu<-a$bray.bal/a$bray
rh<-a$bray.gra/a$bray
br<-a$bray
tu1<-a1$bray.bal
rh1<-a1$bray.gra
br1<-a1$bray
tu2<-a2$bray.bal
rh2<-a2$bray.gra
br2<-a2$bray
tu3<-a3$bray.bal
rh3<-a3$bray.gra
br3<-a3$bray
tu4<-a4$bray.bal
rh4<-a4$bray.gra
br4<-a4$bray
tu<-a$bray.bal
rh<-a$bray.gra
br<-a$bray
b<-c(as.vector(tu),as.vector(tu4),as.vector(tu1),as.vector(tu3),as.vector(tu2))
bb<-rep(c("whole","core","ra","int","abd"),each=703)
kruskal(b,bb,p.adj="fdr",group=T,console=T)
df<-matrix(NA,ncol=2,nrow=5)
rownames(df)<-c("Whole","Core","Rare","Intermediate","Abundant")
df[,1]<-c(mean(tu),mean(tu4),mean(tu1),mean(tu3),mean(tu2))
df[,2]<-c(mean(rh),mean(rh4),mean(rh1),mean(rh3),mean(rh2))
df<-df*100
pdf("beta-partion.pdf",height=5,width=5.5)
par(mar=c(5,7,2,2))
barplot(t(df),horiz=T,col=c("grey36","grey60"),las=1,border=NA
,cex.axis=1.2,cex.names=1.2,cex.lab=1.4,xlab="Proportion(%)")
dev.off()
cc<-c(mean(br),mean(br4),mean(br1),mean(br3),mean(br2))
names(cc)<-c("Whole","Core","Rare","Intermediate","Abundant")
col<-c("darkgreen","red2","darkorange","grey92","cornflowerblue")
pdf("Beta-diversity.pdf",height=5,width=5.5)
par(mar=c(5,7,2,2))
barplot(cc,horiz=T,col=col,las=1,border=NA
,cex.axis=1.2,cex.names=1.2,cex.lab=1.4,xlab="Beta-diversity")
dev.off()
b<-c(as.vector(br),as.vector(br4),as.vector(br1),as.vector(br3),as.vector(br2))
bb<-rep(c("whole","core","ra","int","abd"),each=703)
kruskal(b,bb,p.adj="fdr",group=T,console=T)

#####################################################
#########driving factors of ¦Â-diversity#########
aa<-env[,15:14]
pc<-pcnm(aa)
pc<-pc$vectors
env1<-cbind(pc,env[,1:13])
env1<-as.matrix(env1)
#####forward model selection#####
mod1<-cca(abd)~.,env1)
mod0<-cca(abd~1.,env1)
mod<-ordiR2step(mod0,scope=formula(mod1),perm.max=999)
mod$anova
mod_abd<-mod
mod1<-cca(ra)~.,env1)
mod0<-cca(ra~1.,env1)
mod<-ordiR2step(mod0,scope=formula(mod1),perm.max=999)
mod$anova
mod_ra<-mod

col<-brewer.pal(5,"Set1")
mycol<-rep(col,c(9,9,9,4,7))
si<-scores(mod_abd)$sites
sp<-mod_abd$CCA$biplot,ylim=c(-1.5,1)
pdf("cca_abd.pdf")
plot(si,pch=16,col=mycol,font=2,font.lab=2,
cex.lab=1.2,cex=1.2,xlab="",ylab="")
legend("bottomleft",c("XY","JB","YP","YC","JKH"),
pch=16,col=col,bty="n",
text.col= col,text.font=2,cex=1.2, pt.cex=1,inset=.02)
for(i in 1:nrow(sp))
{Arrows(0,0,sp[i,1],sp[i,2],lwd=2,lty=1,col=colours()[294],arr.length=0.2,
arr.type="triangle")
}
pointLabel(x=sp[,1],y=sp[,2],labels=rownames(sp),
col=colours()[450],cex=1.3,font=4)
abline(h=0,v=0,lty = 3)
dev.off()

si<-scores(mod_ra)$sites
sp<-mod_ra$CCA$biplot,ylim=c(-1.5,1)
pdf("cca_ra.pdf")
plot(si,pch=16,col=mycol,font=2,font.lab=2,
cex.lab=1.2,cex=1.2,xlab="",ylab="")
legend("bottomleft",c("XY","JB","YP","YC","JKH"),
pch=16,col=col,bty="n",
text.col= col,text.font=2,cex=1.2, pt.cex=1,inset=.02)
for(i in 1:nrow(sp))
{Arrows(0,0,sp[i,1],sp[i,2],lwd=2,lty=1,col=colours()[294],arr.length=0.2,
arr.type="triangle")
}
pointLabel(x=sp[,1],y=sp[,2],labels=rownames(sp),
col=colours()[450],cex=1.3,font=4)
abline(h=0,v=0,lty = 3)
dev.off()
######mantel test######
####abd__cca######
local<-c("Pb","TPH","TN","PH")
local<-env1[,local]
local.dis<-vegdist(local,"euclidean")
reg<-c("PCNM3","PCNM9","PCNM20")
reg<-env1[,reg]
reg.dis<-vegdist(reg,"euclidean")
mantel.partial(vegdist(abd),local.dis,reg.dis,"spear")
mantel.partial(vegdist(abd),reg.dis,local.dis,"spear")
mantel(vegdist(abd),local.dis,"spear")
mantel(vegdist(abd),reg.dis,"spear")
pdf("varp_abd.pdf")
otu.varpart= varpart(abd,local,reg)
plot(otu.varpart,digits=4,cex=1.2)
dev.off()
anova.cca(rda(abd,local,reg))
anova.cca(rda(abd,reg,local))
####ra__cca######
local<-c("Pb","Om")
local<-env1[,local]
local.dis<-vegdist(local,"euclidean")
reg<-c("PCNM2","PCNM3","PCNM6")
reg<-env11[,reg]
reg.dis<-vegdist(reg,"euclidean")
mantel.partial(vegdist(ra),local.dis,reg.dis,"spear")
mantel.partial(vegdist(ra),reg.dis,local.dis,"spear")
mantel(vegdist(ra),local.dis,"spear")
mantel(vegdist(ra),reg.dis,"spear")
pdf("varp_ra.pdf")
otu.varpart= varpart(ra,local,reg)
plot(otu.varpart,digits=4,cex=1.2)
dev.off()
anova.cca(rda(ra,local,reg))
anova.cca(rda(ra,reg,local))

###############################
########Map###########
library(mapdata)
library(maptools)
xy<-unique(env1[,15:14])
pdf("map.pdf")
map("china",panel.first=grid(),xlim=c(107,112),ylim=c(33,39),lwd=1.2,col="gray35")
axis(1,lwd=2)
axis(2,lwd=2,las=1)
axis(3,lwd=2)
axis(4,lwd=2,las=1)
box(lwd=2)
points(xy,pch=20,cex=1.2)
pointLabel(x=xy[,1],y=xy[,2],labels=c("Xianyang","Yulin","Yanchuan","Yanchang","Luochuan")
,cex=1,font=3)
legend("bottomleft",legend="Sampling site",pch=20,bty="n",inset=0.05,cex=1)
dev.off()

#####################################
#############NMDS###########
mantel(vegdist(abd),vegdist(otu1),method="spear")
mantel(vegdist(abd),vegdist(ra),method="spear")
mantel(vegdist(ra),vegdist(otu1),method="spear")
a<-rbind(otu,abd,ra)
b<-as.factor(rep(c("Whole","Abundant","Rare"),each=38))
col<-c("darkorange","cornflowerblue","darkgrey")
NMDS<-metaMDS(a)
si<-NMDS$points
pdf("NMDS.pdf")
dataEllipse(si,groups=b,,cex=1.1,cex.lab=1.2,font.lab=2,font=2,las=1,
add=F,levels=0.8,center.pch="+",lwd=2,center.cex=1,
sub="80% confidence ellipses are shown",col=col,pch=c(15,18,20,17))
text(0.7,0.9,"Stress=0.143",font=2)
dev.off()


b2<-rcorr(crt,env[,1:13],type="spearman")
r<-b2$r
p<-b2$P
r[abs(r)<0.5]<-0
diag(p)<-1
pp<-as.vector(p)
qvalue<-fdrtool(pp,statistic="pvalue")
q<-qvalue$qval
p<-matrix(q,ncol(p),ncol(p))
p[p>0.05]<-0
p[p<=0.05&p>0]<-1
c<-p*r
cc<-c[1:ncol(crt1),(ncol(crt1)+1):(ncol(crt1)+11)]
cc<-cc[rowSums(cc)!=0,]
ncc<-otuname1[rownames(cc),2:6]
library(corrplot)
pdf("crt-env.pdf",height=10,width=3.5)
corrplot(cc,"color",cl.pos="b",cl.length=5,mar=c(0,1,0,1),
tl.col = "black",tl.srt=30,tl.cex=1)
dev.off()
write.table(ncc,"crt-env.xls",sep="\t")
