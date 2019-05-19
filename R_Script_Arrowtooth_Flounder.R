#ATF
library(mgcv)
library(maps)
library(mapdata)
library(spacetime)
library(fields)
library(date)
library(colorRamps)  ## to create color palletes
library(itsadug)    ##  to create color.legends ---- WATCH OUT: masks alpha function from scales (used for color transparency)
library(RColorBrewer)


setwd('/Users/lciannel/Documents/MyDocuments/...') #Change this as appropirate

#Bathymetry data from NOAA repositories
bathy.dat<-read.table('BeringDepth.txt',sep='')
names(bathy.dat)<-c('lon','lat','depth')
bathy.dat$depth[bathy.dat$depth>0]<-NA#Avoid points above water
head(bathy.dat)
bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]

#Import juvenile and adult survey data from summer groundfish survey (1982-2012)
atf_adults_catch<-read.csv(file = 'atf_cpue_upto2012_withPhi.csv', header=TRUE, check.names=TRUE)
atf_adults_length<-read.csv(file = 'atf_lengthcpue_sex_1992_2012.csv', header=TRUE, check.names=TRUE)
atf_adults_catch$join<-paste(atf_adults_catch$VESSEL,atf_adults_catch$YEAR,atf_adults_catch$HAUL,sep='.')
atf_adults_length$join<-paste(atf_adults_length$VESSEL, atf_adults_length$YEAR, atf_adults_length$HAUL,sep='.')
atf_adults_length$id_catch<-match(atf_adults_length$join,atf_adults_catch$join)
atf_adults_length$LATITUDE <-atf_adults_catch$LATITUDE[atf_adults_length$id_catch]
atf_adults_length$LONGITUDE <-atf_adults_catch$LONGITUDE[atf_adults_length$id_catch]
atf_adults_length$bt<-atf_adults_catch$BOTTOM_TEMP_CELCIUS[atf_adults_length$id_catch]

#Cumulative distribution of catches as a function of size
length.seq<-seq(min(atf_adults_length$LENGTH),max(atf_adults_length$LENGTH),length=100)
length.cum<-length.seq*NA
for(i in 1:length(length.seq)){
length.cum[i]<-sum(atf_adults_length$NO_HA[atf_adults_length$LENGTH<=length.seq[i]])	
}

#Size bins = XX% of total numerical cacth
XX=0.125
perc.bin=seq(0,1,by=XX)
size.bin<-1:length(perc.bin)*NA
size.bin[1]<-0
size.bin[length(size.bin)]<-max(atf_adults_length$LENGTH)
for(i in 2:(length(perc.bin)-1)){
size.bin[i]<-length.seq[(length.cum/sum(atf_adults_length$NO_HA))>perc.bin[i]][1]
}

plot(length.seq,length.cum/sum(atf_adults_length$NO_HA)*100,ylab="Cumulative biomass (%)")
abline(v=size.bin[2:length(size.bin)])
abline(h=(perc.bin*100)[2:length(perc.bin)])

#Estimate habitat constraint (Rsq) and plot results on a map
dev.new(height=6,width=12)
par(mfrow=c(2,4),mai=c(0.3,0.3,0.2,0.2))

#Juveniles and Adults
all.rsq<-1:(length(size.bin)-1)
for(i in 1:(length(size.bin)-1)){
tmp<-atf_adults_length[atf_adults_length$LENGTH<=size.bin[i+1]&atf_adults_length$LENGTH>size.bin[i],]
pk.l<-aggregate(tmp[,c("YEAR","LATITUDE","LONGITUDE")],list(tmp$join),mean)
names(pk.l)<-c('join','year','lat','lon')
pk.l$cpuelt<-aggregate(tmp$NO_HA,list(tmp$join),sum)$x

#GAM of whole area
pk.gaml1<-gam(log(cpuelt)~factor(year)+s(lon,lat),data=pk.l)
rsq<-summary(pk.gaml1)$r.sq
all.rsq[i]<-rsq

#Map
vis.gam(pk.gaml1,view=c('lon','lat'),plot.type="contour",color="topo",too.far=.03,xlim=range(atf_adults_length$LONGITUDE),ylim=c(53,max(atf_adults_length$LATITUDE)),xlab="",ylab="",main="",axes=F)
axis(1,at=-seq(175,160,by=-5),labels=as.character(seq(175,160,by=-5)))
axis(2,labels=T)
box()
text(-172.8,54,labels=paste(round(size.bin[i],0),"-",round(size.bin[i+1],0),'mm'),cex=1.3)
text(-173.25,54.6,labels=expression(paste('R'^2,' =   ','   %')),cex=1.3)
text(-172.8,54.53,labels=eval(round(rsq*100,0)),cex=1.3)
symbols(pk.l$lon,pk.l$lat,circles=log(pk.l$cpuelt+1),inches=0.08,bg=alpha('grey',f=0.1),fg=alpha('grey',f=0.3),add=T)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(200),labcex=0.4,col='black',add=T,)
map("worldHires",fill=T,col="lightblue4",add=T)}

dev.new(height=4,width=8)
allages<-paste(round(size.bin[1:(length(size.bin)-1)],0),round(size.bin[2:length(size.bin)],0),sep="-")
barplot(all.rsq,names.arg=allages,ylab=expression('R'^2),xlab='Size (mm)',main="Explained spatio-temporal variance by size group")
box()

