#Alaska plaice: habitat associations across different size groups. Only juvenile and adults. Eggs and larval data are not included
library(mgcv)
library(maps)
library(mapdata)
library(spacetime)
library(fields)
library(date)
library(itsadug) 

setwd('/Users/lciannel/Documents/MyDocuments/...') #Change this as appropirate
bathy.dat<-read.table('BeringDepth.txt',sep='')#Depth dara are from NOAA repositories (1 min resolution)
names(bathy.dat)<-c('lon','lat','depth')
bathy.dat$depth[bathy.dat$depth>0]<-NA#Avoid points above water
head(bathy.dat)
bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]

#Adult and Juvenile data: summer bottom trawl survey. Data goes from 1982 to 2006
akp_adults_catch<-read.csv(file = 'akp_hauls.csv', header=TRUE, check.names=TRUE)
akp_adults_length<-read.csv(file = 'akp_flatfile.csv', header=TRUE, check.names=TRUE)
akp_adults_length$cpuelt<-log(akp_adults_length$CPUE_LENGTH)
akp_adults_length$join<-paste(akp_adults_length$VESSEL, akp_adults_length$YEAR, akp_adults_length$HAUL,sep='.')

#Cumulative distribution of catches as a function of size. This is needed to calculate size bins with approximately equal CPUE
length.seq<-seq(min(akp_adults_length$LENGTH),max(akp_adults_length$LENGTH),length=100)
length.cum<-length.seq*NA
for(i in 1:length(length.seq)){
length.cum[i]<-sum(akp_adults_length$CPUE_LENGTH[akp_adults_length$LENGTH<=length.seq[i]])	
}

#Size bins = each size bin cintains XX% of total numerical cacth
XX=0.125
perc.bin=seq(0,1,by=XX)
size.bin<-1:length(perc.bin)*NA
size.bin[1]<-min(akp_adults_length$LENGTH)
size.bin[length(size.bin)]<-max(akp_adults_length$LENGTH)
for(i in 2:(length(perc.bin)-1)){
size.bin[i]<-length.seq[(length.cum/sum(akp_adults_length$CPUE_LENGTH))>perc.bin[i]][1]
}

plot(length.seq,length.cum/sum(akp_adults_length$CPUE_LENGTH)*100,ylab="Cumulative biomass (%)")
abline(v=size.bin[2:length(size.bin)])
abline(h=(perc.bin*100)[2:length(perc.bin)])

#Find out the sample size for each size bin
sample.size<-1:(length(size.bin)-1)
for(i in 1:(length(size.bin)-1)){
tmp<-akp_adults_length[akp_adults_length$LENGTH<=size.bin[i+1]&akp_adults_length$LENGTH>size.bin[i],]
pk.l<-aggregate(tmp[,c("YEAR","LATITUDE","LONGITUDE")],list(tmp$join),mean)
names(pk.l)<-c('join','year','lat','lon')
pk.l$cpuelt<-aggregate(tmp$cpuelt,list(tmp$join),sum)$x
sample.size[i]<-nrow(pk.l)}

#Estimate habitat constraint (Rsq) and plot results on a map
dev.new(height=6,width=12)
par(mfrow=c(2,4),mai=c(0.3,0.3,0.2,0.2))

#Loop for Juveniles and Adults
rsq.ad<-1:(length(size.bin)-1)
for(i in 1:(length(size.bin)-1)){
tmp<-akp_adults_length[akp_adults_length$LENGTH<=size.bin[i+1]&akp_adults_length$LENGTH>size.bin[i],]
pk.l<-aggregate(tmp[,c("YEAR","LATITUDE","LONGITUDE")],list(tmp$join),mean)
names(pk.l)<-c('join','year','lat','lon')
pk.l$cpuelt<-aggregate(tmp$CPUE_LENGTH,list(tmp$join),sum)$x

#GAM of whole area
pk.gaml1<-gam(log(cpuelt)~factor(year)+s(lon,lat),data=pk.l)
rsq<-summary(pk.gaml1)$r.sq
rsq.ad[i]<-rsq

#Map
vis.gam(pk.gaml1,view=c('lon','lat'),plot.type="contour",color="topo",too.far=.03,xlim=range(akp_adults_length$LONGITUDE),ylim=c(53,max(akp_adults_length$LATITUDE)),xlab="",ylab="",main="",axes=F)
axis(1,at=-seq(175,160,by=-5),labels=as.character(seq(175,160,by=-5)))
axis(2,labels=T)
box()
text(-172.8,54,labels=paste(round(size.bin[i],0),"-",round(size.bin[i+1],0),'mm'),cex=1.3)
text(-173.25,54.6,labels=expression(paste('R'^2,' =   ','   %')),cex=1.3)
text(-172.8,54.53,labels=eval(round(rsq*100,0)),cex=1.3)
symbols(pk.l$lon,pk.l$lat,circles=log(pk.l$cpuelt+1),inches=0.08,bg=alpha('grey',f=0.1),fg=alpha('grey',f=0.3),add=T)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(200),labcex=0.4,col='black',add=T,)
map("worldHires",fill=T,col="lightblue4",add=T)}
#dev.copy(jpeg,'AP.jpg',height=6,width=14,res=100,units='in')
#dev.off()

