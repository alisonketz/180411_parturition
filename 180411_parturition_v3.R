###
### 1/31/2018 
### Alison Ketz
### Movement data preprocessing using adehabitatLT package
###


###
### Preliminaries
###

rm(list=ls())

library(geosphere)
library(lubridate)
library(Hmisc)
library(ggplot2)
library(adehabitatLT)
library(mvtnorm)
library(beepr)
library(stringr)

setwd("~/Documents/Parturition/180411_parturition")

###
### Load previous runs
###

load("anomaly_v3.Rdata")

###
### Load VIT data
###

d.vit=mdb.get('~/Documents/Data/SWDPPdeerDB.mdb',tables= "VIT")
names(d.vit)=tolower(gsub('[[:punct:]]',"",names(d.vit)))
names(d.vit)[10]="lowtag"
d.vit$datedropped=as.character(d.vit$datedropped)
d.vit$juliandropped=yday(mdy_hms(d.vit$datedropped))
d.vit$datefawnfound=as.character(d.vit$datefawnfound)
d.vit$julianfawnfound=yday(mdy_hms(d.vit$datefawnfound))
d.vit=d.vit[d.vit$datedropped!="",]
n.vit=length(d.vit$lowtag)

#reorder vit data by lowtag/lowtag
d.vit=d.vit[order(d.vit$lowtag),]

#extract the individual ids
individs=d.vit$lowtag


#manipulate the julian day for the doe 6865 with best guess based on movement rates

d.vit$juliandropped[9] = 147
d.vit$lowtag[9]

# d.vit=d.vit[-9,]
# individs=d.vit$lowtag
# n.vit=length(d.vit$lowtag)


#par(mfrow=c(1,1))
#plot(d$julian[d$lowtag==6865],d$distance[d$lowtag==6865])

###
### Load data GPS location data
###

d = matrix(NA,nr=1,nc=13)
#for loop for reading in data, using vector of lowtag's from the vit dataset
for(i in 1:n.vit){
    d.temp = read.table(paste("/home/aketz/Documents/Data/GPS_locations_ind/",d.vit$lowtag[i],".csv",sep=""),sep=",",header=TRUE,stringsAsFactors = FALSE)
    d.temp$lowtag = d.vit$lowtag[i]
    names(d.temp)=tolower(names(d.temp))
    d=rbind(d,as.matrix(d.temp))
}
d=d[-1,]
d=data.frame(d,stringsAsFactors = FALSE)

for(j in 1:dim(d)[2]){
    d[,j]=str_trim(d[,j],side="both")
}

class.indx=c(5:7,9:12)
for(j in class.indx){
    d[,j]=as.numeric(d[,j])
}

d$lowtag=as.factor(d$lowtag)

head(d)
###
### double checking for duplicate locations
###

summary(duplicated(d))


#calculating julian day and omitting outside of parturition window
d$julian=yday(mdy_hms(d$date_time_local))

###
### increased fixes parturition window
###

start=yday(mdy("03/20/2017")) # May 6 beginning of parturition period
end=yday(mdy("07/07/2017")) #end of parturition period

###
### subset entire dataset to parturition window
###

d=d[d$julian>start & d$julian <= end,]

#removing last day of parturition of 5004, outlier movements
rm5004 = (dim(d[d$lowtag==5004,])[1]-15):dim(d[d$lowtag==5004,])[1]
d=d[-rm5004,]
tail(d[d$lowtag==5004,])


#or keep all data, but exclude first couple hundred observations
#for plotting all data, not just window of parturition

# d=data.frame(d %>% group_by(lowtag) %>% slice(160:(n()-300)))
# n.obs=table(d$lowtag)
# n.obs
# 
# indxing=c(1,rep(NA,n.vit-1),dim(d)[1])
# for(j in 2:n.vit){
#     indxing[j]=indxing[j-1]+n.obs[j-1]
# }
# d[indxing,]
# 
# d[d$lowtag==individs[1],]


###
### double checking for duplicate locations
###

summary(duplicated(d))

###
### Adding Vit data to main GPS dataframe
###

d$dropped = 0

records=dim(d)[1]
records
for(i in 1:records){
    for(j in 1:n.vit){
        if(d$lowtag[i]==d.vit$lowtag[j]){
            if(d$julian[i]==d.vit$juliandropped[j])d$dropped[i]=1
        }
    }
}
sum(d$dropped)


###
### Converting date to POSIXct format
###

d$date_time_local=as.POSIXct(strptime(d$date_time_local,format="%m-%d-%Y %H:%M:%S"),tz="CST6CDT")
d$date_time_gmt=as.POSIXct(strptime(d$date_time_gmt,format="%m-%d-%Y %H:%M:%S"),tz="GMT")



#Create time lag between successive locations to censor data if needed.
time.diff <- diff(d$date_time_local)
d=d[-1,]
d$timediff <-round(as.numeric(abs(time.diff)))
rm=which(d$timediff>10)
d=d[-rm,]
names(d)[1]="devname"

###
### Dealing with missing data
###

#option 1, just remove
# d=d[!is.na(d$longitude),]

#option 2 impute
d$missing=0
for( i in 1:dim(d)[1]){
    if(is.na(d$longitude[i]))d$missing[i]=1
}

miss.per.ind=c()
for(j in 1:n.vit){
    miss.per.ind=c(miss.per.ind,sum(d$missing[d$lowtag==individs[j]]))
}
miss.per.ind
d=d[-c(1:3),]
for(i in 2:(dim(d)[1]-1)){
    if(is.na(d$longitude[i])){
        a=i-1
        while(is.na(d$longitude[a])){a=a-1}
        b=i+1
        while(is.na(d$longitude[b])){b=b+1}
        d[i,6:5] = midPoint(d[a,6:5],d[b,6:5])
    }
}

### 
### Without projection of datum into R, can use geospatial package to calculate distance and bearings
###

bearing.out=bearing(cbind(d$longitude,d$latitude))
d$bearing=bearing.out
dist.out = distHaversine(d[,6:5])
d=d[-1,]
d$distance = dist.out
d=d[-c(dim(d)[1]-1,dim(d)[1]),]#remove last 2 entries which are NA and NaN

###
### Projections! 
###

# setup coordinates
coords = cbind(d$longitude, d$latitude)
sp = SpatialPoints(coords)

# make spatial data frame
# spdf = SpatialPointsDataFrame(coords, d)
spdf = SpatialPointsDataFrame(sp, d)

# EPSG strings
latlong = "+init=epsg:4326"
proj4string(spdf) = CRS(latlong)

d.sp.proj = spTransform(spdf, CRS("+proj=tmerc +lat_0=0 +lon_0=-90 +k=0.9996 +x_0=520000
                                  +y_0=-4480000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
d=data.frame(d.sp.proj)

#20 and 21 are the coordinates in UTMs x y

###
### Animal paths
### Using the adelehabitatLT

d.traj <- as.ltraj(d[,20:21], date=d$date_time_local,id=d$lowtag)
plot(d.traj)

#Plot the trajectory for each individual - the numbers correspond to the ID in the d.traj object above
#blue is start
#red is end point

par(mfrow=c(2,2))
for(i in 1:n.vit){
    plot(d.traj[i])
}
###
### The last movements of individ1 is crazy
### removed those points up above
# plot(d.traj[1])

#converts traj object to data frame
dfdeer <- ld(d.traj)
dfdeer$id=as.character(dfdeer$id)
dfdeer=rbind(rep(NA,dim(dfdeer)[2]),dfdeer)
dfdeer=dfdeer[-dim(dfdeer)[1],]
d$rel.angle=dfdeer$rel.angle
d$dist.traj=dfdeer$dist
d$R2n=dfdeer$R2n
d$dx=dfdeer$dx
d$dy=dfdeer$dy
d$dt=dfdeer$dt

###
### Description of parameters returned by adehabitatLT object model. 
###

# • dx, dy, dt: these parameters measured at relocation i describe the increments
# of the x and y directions and time between the relocations i and
# i + 1. Such parameters are often used in the framework of stochastic differential
# equation modelling (e.g. Brillinger et al. 2004, Wiktorsson et al.
#                     2004);
# • dist: the distance between successive relocations is often used in animal
# movement analysis (e.g. Root and Kareiva 1984, Marsh and Jones 1988);

# • abs.angle: the absolute angle αi between the x direction and the step
# built by relocations i and i + 1 is sometimes used together with the parameter
# dist to fit movement models (e.g. Marsh and Jones 1988);

# • rel.angle: the relative angle βi measures the change of direction between
# the step built by relocations i − 1 and i and the step built by relocations
# i and i + 1 (often called “turning angle”). It is often used together with
# the parameter dist to fit movement models (e.g. Root and Kareiva 1984,
#                                            Marsh and Jones 1988);

# • R2n: the squared distance between the first relocation of the trajectory
# and the current relocation is often used to test some movements models
# (e.g. the correlated random walk, see the seminal paper of Kareiva and
# Shigesada, 1983).


########################################################################################################
###
### remove the initial index of each individual 
### because it makes no sense to calculate distances 
### between the locations of the individuals
###
#######################################################################################################

remove.indx=which(is.na(d$dist.traj))
d=d[-remove.indx,]



d$step = d$distance/d$timediff
# d$logstep = log(d$step)



##############################################################################################
###
### creating a three day (prior) based time step data frame, and running anomaly detection 
###
##############################################################################################
# features to calculate within day means and standard deviations, and assign to daily
#lowtag, julian, dropped,step,altittude, temp, bearing, rel.angle, R2n, dx,dy
#lowtag, julian, dropped,
# mean and sd of : step,altitude, temp, bearing, rel.angle, R2n, dx,dy
#matrix dim should  = 8*2+3 = 19


d$lowtag = as.character(d$lowtag)
d.dos.temp=matrix(NA,nr = 1,nc=19)#this will be too big, will trim after
for(j in 1:n.vit){
    d.temp=d[d$lowtag == individs[j],]
    julian.temp = unique(d.temp$julian)
    temp.mat = matrix(NA,nr = length(julian.temp),nc = 19)
    for(i in 1:length(julian.temp)){
        temp.mat[i,1] = d.temp$lowtag[1]
        temp.mat[i,2] = julian.temp[i]
        temp.mat[i,3] = ifelse(sum(d.temp$dropped[d.temp$julian == julian.temp[i]])>1,1,0)
        temp.mat[i,4] = mean(d.temp$step[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.step
        temp.mat[i,5] = sd(d.temp$step[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.step
        temp.mat[i,6] = mean(d.temp$altitude[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.altitude
        temp.mat[i,7] = sd(d.temp$altitude[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.altitude
        temp.mat[i,8] = mean(d.temp$temp[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.temp
        temp.mat[i,9] = sd(d.temp$temp[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.temp
        temp.mat[i,10] = mean(d.temp$bearing[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.bearing
        temp.mat[i,11] = sd(d.temp$bearing[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.bearing
        temp.mat[i,12] = mean(d.temp$rel.angle[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.rel.angle
        temp.mat[i,13] = sd(d.temp$rel.angle[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.rel.angle
        temp.mat[i,14] = mean(d.temp$R2n[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.R2n
        temp.mat[i,15] = sd(d.temp$R2n[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.R2n
        temp.mat[i,16] = mean(d.temp$dx[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.dx
        temp.mat[i,17] = sd(d.temp$dx[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.dx
        temp.mat[i,18] = mean(d.temp$dy[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#mu.dy
        temp.mat[i,19] = sd(d.temp$dy[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1)],na.rm=TRUE)#sig.dy
    }
    d.dos.temp = rbind(d.dos.temp,temp.mat)
}
d.dos.temp=d.dos.temp[-1,]
d.dos = data.frame(d.dos.temp,stringsAsFactors = FALSE)
names(d.dos)=c("lowtag","julian","dropped","mu.step","sig.step","mu.altitude","sig.altitude","mu.temp","sig.temp",
                "mu.bearing","sig.bearing","mu.rel.angle","sig.rel.angle","mu.R2n","sig.R2n","mu.dx","sig.dx","mu.dy","sig.dy")
for(i in 3:19){
    d.dos[,i]=as.numeric(d.dos[,i])
}
head(d.dos)
###
### center and scale variables for prediction
###

for(j in 1:n.vit){
    d.dos[d.dos$lowtag == individs[j],4:19]=apply(d.dos[d.dos$lowtag == individs[j],4:19],2,scale)
}

d.dos$julian = as.integer(d.dos$julian)

head(d.dos)
###
### Paturition window
###

part.window=yday("2017-05-06")

###
### Total number of possible anomalies to detect
###

possible.hits=rep(NA,n.vit)
for(j in 1:n.vit){
    possible.hits[j] = sum(d.dos$dropped[d.dos$lowtag==individs[j]])
}
possible.hits



##############################################################################################
###
### Expanding Anomaly Detection Training 
###
##############################################################################################
d.train=d.dos
ph=possible.hits
nInd.train=length(unique(d.train[,1]))

vitdrop = d.vit$juliandropped
# source("anomalyDetectUniv.R")
source("anomalyDetect.R")
source("training_detector.R")
source("evaluate.R")

#test 1 individual , 1 covariate
d.train=d.train[d.train[,1]==individs[1],1:4]
epsilon=.15
fit.train=training(d.train,eps=epsilon,pw=part.window,vd=vitdrop[1])
fit.train$alarm[[1]]
evaluate(alarm=fit.train$alarm,possible.hits=ph[1],nInd=1,vitdropday=d.vit$juliandropped[1])
# evaluate= function(alarm,possible.hits,nInd,vitdropday){
nCovs=1
fit.train2 = anomalyDetect(n.vit=length(unique(d.train$lowtag)),id=individs[1],d=d.train,eps=epsilon,covs.indx = 4:(3+nCovs))# fit model
fit.train2$alarm
evaluate(alarm=fit.train2$alarm,possible.hits=ph[1],nInd=1,vitdropday=d.vit$juliandropped[1])


fit.train$ind.mean
fit.train$ind.sd
fit.train2$ind.mean
fit.train2$ind.sd

fit.train$detect.quant
fit.train2$detect.quant

fit.train$threshold.density
fit.train2$threshold.density

fit.train$detect.density

fit.train2$detect.density

fit.train$alarm[[1]][,2]==fit.train2$alarm[[1]][,2]



#test 1 individual, multiple covariates
source("anomalyDetect.R")
source("training_detector.R")
source("evaluate.R")

d.train = d.dos
d.train=d.train[d.train[,1]==individs[1],]
epsilon=rep(.12,dim(d.train)[2]-3)
fit.train=training(d.train,eps=epsilon,pw=part.window,vd=vitdrop[1])
fit.train

nCovs=16

fit.train2 = anomalyDetect(n.vit=length(unique(d.train$lowtag)),id=individs[1],d=d.train,eps=epsilon,covs.indx = 4:(3+nCovs))# fit model
fit.train2
evaluate(alarm=fit.train$alarm,possible.hits=ph[1],nInd=1,vitdropday=d.vit$juliandropped[1])
evaluate(alarm=fit.train2$alarm,possible.hits=ph[1],nInd=1,vitdropday=d.vit$juliandropped[1])

fit.train$ind.mean==fit.train2$ind.mean

fit.train$ind.sd==fit.train2$ind.sd

fit.train$detect.quant
fit.train2$detect.quant

fit.train$threshold.density
fit.train2$threshold.density

fit.train$detect.density
fit.train2$detect.density

fit.train$alarm[[1]][,2]==fit.train2$alarm[[1]][,2]
fit.train$alarm[[1]][,2]
  
fit.train2$alarm[[1]][,2]



#test multiple individuals, 1 covariate
source("anomalyDetect.R")
source("training_detector.R")
source("evaluate.R")
d.train = d.dos
d.train=d.train[,1:4]
nInd.train=length(unique(d.train[,1]))
epsilon=rep(.12,dim(d.train)[2]-3)
fit.train=training(d.train,eps=epsilon,pw=part.window,vd=vitdrop)
fit.train$alarm
nCovs=1
fit.train2 = anomalyDetect(n.vit=length(unique(d.train$lowtag)),id=individs,d=d.train,eps=epsilon,covs.indx = 4:(3+nCovs))# fit model

evaluate(alarm=fit.train$alarm,possible.hits=ph,nInd=12,vitdropday=d.vit$juliandropped)
evaluate(alarm=fit.train2$alarm,possible.hits=ph,nInd=12,vitdropday=d.vit$juliandropped)

fit.train$ind.mean==fit.train2$ind.mean
fit.train$ind.mean
fit.train2$ind.mean

fit.train$ind.sd==fit.train2$ind.sd

fit.train$detect.quant
fit.train2$detect.quant

fit.train$threshold.density
fit.train2$threshold.density

fit.train$detect.density
fit.train2$detect.density

fit.train$alarm[[1]][,2]==fit.train2$alarm[[1]][,2]
fit.train$alarm[[1]][,2]

fit.train2$alarm[[1]][,2]




evaluate(alarm=fit.train$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=d.vit$juliandropped)


#test multiple individuals, multiple covariates
source("anomalyDetect.R")
source("training_detector.R")
source("evaluate.R")
d.train = d.dos
epsilon=rep(.15,dim(d.train)[2]-3)
aa=Sys.time()
fit.train=training(d.train,eps=epsilon,pw=part.window,vd=vitdrop)
Sys.time()-aa
fit.train$alarm[[1]]
nCovs=16
fit.train2 = anomalyDetect(n.vit=length(unique(d.train$lowtag)),id=individs,d=d.train,eps=epsilon,covs.indx = 4:(3+nCovs))# fit model

evaluate(alarm=fit.train$alarm,possible.hits=ph,nInd=12,vitdropday=d.vit$juliandropped)
evaluate(alarm=fit.train2$alarm,possible.hits=ph,nInd=12,vitdropday=d.vit$juliandropped)

fit.train$ind.mean==fit.train2$ind.mean
fit.train$ind.mean
fit.train2$ind.mean

fit.train$ind.sd==fit.train2$ind.sd

fit.train$detect.quant
fit.train2$detect.quant

fit.train$threshold.density
fit.train2$threshold.density

fit.train$detect.density
fit.train2$detect.density

fit.train$alarm[[1]][,2]==fit.train2$alarm[[1]][,2]
fit.train$alarm[[1]][,2]

fit.train2$alarm[[1]][,2]

evaluate(alarm=fit.train$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=d.vit$juliandropped)


###
### Running the leave one out cross validation training evaluation function
###
source("training_detector.R")
source("evaluate.R")
source("training_function.R")
source("anomalyDetect.R")

# starting=Sys.time()
fit.loo=loo.train(d.train=d.dos,part.window=126,ph=possible.hits,vdrop=d.vit$juliandropped)

# loo_train = function(d.train,part.window=126,ph,vdrop){
# ending=Sys.time()
# total.time=ending-starting
# total.time

fit.loo$decide
fit.loo$train.cov[,,2]
fit.loo$compare
fit.loo$loo.eval

fit.loo$train.cov[,,2]



for(j in 1:10){
  beep(sound=5)
}

fit.loo$loo.fittest

#fitting with different data for determining the quantiles for the density threshold calculation
# fit.loo.v1 = fit.loo#<(vd-2)
# fit.loo.v2 = fit.loo#=1:nT
# fit.loo.v3 = fit.loo#<pw
# fit.loo.v4 = fit.loo#>=pw:vd

detect.quant.test.v1=matrix(NA,nr=36,nc=16)
detect.quant.test.v2=matrix(NA,nr=36,nc=16)
detect.quant.test.v3=matrix(NA,nr=36,nc=16)
detect.quant.test.v4=matrix(NA,nr=36,nc=16)

for(k in 1:36){
  detect.quant.test.v1[k,]=fit.loo.v1$loo.fittest[[k]]$detect.quant
  detect.quant.test.v2[k,]=fit.loo.v2$loo.fittest[[k]]$detect.quant
  detect.quant.test.v3[k,]=fit.loo.v3$loo.fittest[[k]]$detect.quant
  detect.quant.test.v4[k,]=fit.loo.v4$loo.fittest[[k]]$detect.quant
}
detect.quant.test.v1
detect.quant.test.v2
detect.quant.test.v3
detect.quant.test.v4

fit.loo.v1$decide
nCovs=16
decide.versions = array(NA,c(nCovs,3,4))
decide.versions[,,1]=fit.loo.v1$decide
decide.versions[,,2]=fit.loo.v2$decide
decide.versions[,,3]=fit.loo.v3$decide
decide.versions[,,4]=fit.loo.v4$decide

decide.versions[,2,]
# fit.loo.v1 = fit.loo#<(vd-2)
# fit.loo.v2 = fit.loo#=1:nT
# fit.loo.v3 = fit.loo#<pw
# fit.loo.v4 = fit.loo#>=pw:vd

fit.recall.v1 = training(d.dos,eps=fit.loo.v1$decide[,2],pw=part.window,vd=vitdrop)# fit model v1 <(vd-2)
eval.recall.v1=evaluate(alarm=fit.recall.v1$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.v2 = training(d.dos,eps=fit.loo.v2$decide[,2],pw=part.window,vd=vitdrop)# fit model #=1:nT
eval.recall.v2=evaluate(alarm=fit.recall.v2$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.v3 = training(d.dos,eps=fit.loo.v3$decide[,2],pw=part.window,vd=vitdrop)# fit model <pw
eval.recall.v3=evaluate(alarm=fit.recall.v3$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)
fit.recall.v4 = training(d.dos,eps=fit.loo.v4$decide[,2],pw=part.window,vd=vitdrop)# fit model >=pw:vd
eval.recall.v4=evaluate(alarm=fit.recall.v4$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=vitdrop)

fit.recall = list(fit.recall.v1,fit.recall.v2,fit.recall.v3,fit.recall.v4)
hit.recall=rep(list(),4)
output=function(x){
  hit.recall=array(NA,c(nInd.train,9))
  for(j in 1:nInd.train){
    hit.recall[j,1]=ifelse(!length(x$alarm[[j]]$julian[x$alarm[[j]]$julian <(vitdrop[j]-1)]),0.001,length(unique(x$alarm[[j]]$julian[x$alarm[[j]]$julian <(vitdrop[j]-1)])))
    hit.recall[j,2]=(vitdrop[j]-1)-126+1
    hit.recall[j,3]=hit.recall[j,1]/hit.recall[j,2]
    hit.recall[j,4]=sum(x$alarm[[j]]$julian <=vitdrop[j]+1 & x$alarm[[j]]$julian >=vitdrop[j]-1)>0
    hit.recall[j,5]=sum(x$alarm[[j]]$julian <=vitdrop[j]+2 & x$alarm[[j]]$julian >=vitdrop[j]-2)>0
    hit.recall[j,6]=sum(x$alarm[[j]]$julian <=vitdrop[j]+1 & x$alarm[[j]]$julian >=vitdrop[j]-1)
    hit.recall[j,7]=sum(x$alarm[[j]]$julian <=vitdrop[j]+2 & x$alarm[[j]]$julian >=vitdrop[j]-2)
    hit.recall[j,8]=sum(x$alarm[[j]]$julian ==vitdrop[j])
    hit.recall[j,9]=sum(x$alarm[[j]]$julian <(vitdrop[j]-1))
  }
  hit.recall
}

output(fit.recall.v1)
hit.recall = lapply(fit.recall,output)


day1.results = cbind(hit.recall[[1]][,4],
hit.recall[[2]][,4],
hit.recall[[3]][,4],
hit.recall[[4]][,4])

day2.results = cbind(hit.recall[[1]][,5],
                     hit.recall[[2]][,5],
                     hit.recall[[3]][,5],
                     hit.recall[[4]][,5])

day2.results

results.variants = data.frame(cbind(c("<(vd-2)","1:nT","<pw",">=pw"),rbind(apply(hit.recall[[1]][,4:5],2,function(x){sum(x)/12}),
apply(hit.recall[[2]][,4:5],2,function(x){sum(x)/12}),
apply(hit.recall[[3]][,4:5],2,function(x){sum(x)/12}),
apply(hit.recall[[4]][,4:5],2,function(x){sum(x)/12}))),stringsAsFactors = FALSE)

names(results.variants) = c("Model","Prop 1 Day","Prop 2 Day")
results.variants

results.fp = data.frame(cbind(hit.recall[[1]][,1],hit.recall[[2]][,1],hit.recall[[3]][,1],hit.recall[[4]][,1]))
results.fp[results.fp==.001] = 0
names(results.fp)=c("<(vd-2)","1:nT","<pw",">=pw")
apply(results.fp,2,sum)


results.decide = data.frame(decide.versions[,2,])
names(results.decide)=c("<(vd-2)","1:nT","<pw",">=pw")
results.decide

# ##############################################################################################
# ###
# ### Multivariate Anomaly Detection Using F1 criteria with above epsilon
# ### changing for all featurs predictors
# ###
# ##############################################################################################

###precision final model
fit.prec = training(d.dos,eps=fit.loo$decide[,1],pw=part.window,vd=vitdrop)# fit model
eval.prec=evaluate(alarm=fit.prec$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=d.vit$juliandropped)

###recall final model
fit.recall = training(d.dos,eps=fit.loo$decide[,2],pw=part.window,vd=vitdrop)# fit model
eval.recall=evaluate(alarm=fit.recall$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=d.vit$juliandropped)

### F1 final model
fit.F1 = training(d.dos,eps=fit.loo$decide[,3],pw=part.window,vd=vitdrop)# fit model
eval.F1=evaluate(alarm=fit.F1$alarm,possible.hits=ph,nInd=nInd.train,vitdropday=d.vit$juliandropped)

###Plots for recall criteria final model
# 
# ###recall final model plots
# par(mfrow=c(4,1))
# for(j in 1:n.vit){
#   pdf(paste("outputplots_anomaly3b/",individs[j],"_recall_part1.pdf",sep=""),width=6,height=9)
#   par(mfrow=c(4,1))
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$mu.step[d.dos$lowtag==individs[j]],main=paste(individs[j],"mu.step"))
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.step,col=2)
#   abline(v=d.vit$juliandropped[j])
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$sig.step[d.dos$lowtag==individs[j]],main="sig.step")
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.step,col=2)
#   abline(v=d.vit$juliandropped[j])
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$mu.altitude[d.dos$lowtag==individs[j]],main="mu.altitude")
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.altitude,col=2)
#   abline(v=d.vit$juliandropped[j])
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$sig.altitude[d.dos$lowtag==individs[j]],main="sig.altitude")
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.altitude,col=2)
#   abline(v=d.vit$juliandropped[j])
#   dev.off()
# }
# for(j in 1:n.vit){
#   pdf(paste("outputplots_anomaly3b/",individs[j],"_recall_part2.pdf",sep=""),width=6,height=9)
#   par(mfrow=c(4,1))
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$mu.temp[d.dos$lowtag==individs[j]],main=paste(individs[j],"mu.temp"))
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.temp,col=2)
#   abline(v=d.vit$juliandropped[j])
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$sig.temp[d.dos$lowtag==individs[j]],main="sig.temp")
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.temp,col=2)
#   abline(v=d.vit$juliandropped[j])
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$mu.bearing[d.dos$lowtag==individs[j]],main="mu.bearing")
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.bearing,col=2)
#   abline(v=d.vit$juliandropped[j])
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$sig.bearing[d.dos$lowtag==individs[j]],main="sig.bearing")
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.bearing,col=2)
#   abline(v=d.vit$juliandropped[j])
#   dev.off()
# }
# 
# for(j in 1:n.vit){
#   pdf(paste("outputplots_anomaly3b/",individs[j],"_recall_part3.pdf",sep=""),width=6,height=9)
#   par(mfrow=c(4,1))
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$mu.rel.angle[d.dos$lowtag==individs[j]],main=paste(individs[j],"mu.rel.angle"))
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.rel.angle,col=2)
#   abline(v=d.vit$juliandropped[j])
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$sig.rel.angle[d.dos$lowtag==individs[j]],main="sig.rel.angle")
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.rel.angle,col=2)
#   abline(v=d.vit$juliandropped[j])
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$mu.R2n[d.dos$lowtag==individs[j]],main="mu.R2n")
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.R2n,col=2)
#   abline(v=d.vit$juliandropped[j])
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$sig.R2n[d.dos$lowtag==individs[j]],main="sig.R2n")
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.R2n,col=2)
#   abline(v=d.vit$juliandropped[j])
#   dev.off()
# }
# for(j in 1:n.vit){
#   pdf(paste("outputplots_anomaly3b/",individs[j],"_recall_part4.pdf",sep=""),width=6,height=9)
#   par(mfrow=c(4,1))
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$mu.dx[d.dos$lowtag==individs[j]],main=paste(individs[j],"mu.dx"))
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.dx,col=2)
#   abline(v=d.vit$juliandropped[j])
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$sig.dx[d.dos$lowtag==individs[j]],main="sig.dx")
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.dx,col=2)
#   abline(v=d.vit$juliandropped[j])
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$mu.dy[d.dos$lowtag==individs[j]],main="mu.dy")
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.dy,col=2)
#   abline(v=d.vit$juliandropped[j])
#   plot(d.dos$julian[d.dos$lowtag==individs[j]],d.dos$sig.dy[d.dos$lowtag==individs[j]],main="sig.dy")
#   points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.dy,col=2)
#   abline(v=d.vit$juliandropped[j])
#   dev.off()
# }
# 

###
### Evaluating predictions
###

#how many days are anomalies prior to vit drop day
#col1 is total number of days with hits prior to drop+1
#col2 is the total number of days prior to vit drop
#col3 is the ratio of hits to total possible days
#col4 is the logical, whether the the algorithm hit +/- 1 day of vit drop
#col5 is the logical, whether the the algorithm hit +/- 2 days of vit dropth
#col6 is the total number of hits within 1 days
#col7 is the total number of hits within 2 days
#col8 is the total number of hits on vit drop day
#col9 is the total number of hits prior to vit drop day
hit.dos.prec=matrix(NA,nr=nInd.train,nc = 9)
for(j in 1:nInd.train){
    hit.dos.prec[j,1]=length(unique(fit.prec$alarm[[j]]$julian[fit.prec$alarm[[j]]$julian <(d.vit$juliandropped[j]-1)]))
    hit.dos.prec[j,2]=(d.vit$juliandropped[j]-1)-126+1
    hit.dos.prec[j,3]=hit.dos.prec[j,1]/hit.dos.prec[j,2]
    hit.dos.prec[j,4]=sum(fit.prec$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.prec$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)>0
    hit.dos.prec[j,5]=sum(fit.prec$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.prec$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)>0
    hit.dos.prec[j,6]=sum(fit.prec$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.prec$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)
    hit.dos.prec[j,7]=sum(fit.prec$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.prec$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)
    hit.dos.prec[j,8]=sum(fit.prec$alarm[[j]]$julian ==d.vit$juliandropped[j])
    hit.dos.prec[j,9]=sum(fit.prec$alarm[[j]]$julian <(d.vit$juliandropped[j]-1))
    }
hit.dos.prec

# write.csv(hit.dos.prec,file="dos_summary_hitday_prec.csv")

hit.dos.recall=matrix(NA,nr=n.vit,nc = 9)
for(j in 1:n.vit){
    hit.dos.recall[j,1]=length(unique(fit.recall$alarm[[j]]$julian[fit.recall$alarm[[j]]$julian <(d.vit$juliandropped[j]-1)]))
    hit.dos.recall[j,2]=(d.vit$juliandropped[j]-1)-126+1
    hit.dos.recall[j,3]=hit.dos.recall[j,1]/hit.dos.recall[j,2]
    hit.dos.recall[j,4]=sum(fit.recall$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.recall$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)>0
    hit.dos.recall[j,5]=sum(fit.recall$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.recall$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)>0
    hit.dos.recall[j,6]=sum(fit.recall$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.recall$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)
    hit.dos.recall[j,7]=sum(fit.recall$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.recall$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)
    hit.dos.recall[j,8]=sum(fit.recall$alarm[[j]]$julian ==d.vit$juliandropped[j])
    hit.dos.recall[j,9]=sum(fit.recall$alarm[[j]]$julian <(d.vit$juliandropped[j]-1))
}
hit.dos.recall

# write.csv(hit.dos.recall,file="dos_summary_hitday_recall.csv")

fit.recall$alarm[[2]]


Metric=c("total number of days with hits prior to day before vit drop",
         "total number of days prior to vit drop",
         "ratio of hits to total possible days",
         "logical- whether the the algorithm hit +/- 1 day of vit drop",
         "logical- whether the the algorithm hit +/- 2 days of vit dropth",
         "total number of hits within 1 days",
         "total number of hits within 2 days",
         "total number of hits on vit drop day",
         "total number of hits prior to vit drop day")

hit.dos.recall[,3]=round(hit.dos.recall[,3],2)

out.recall=as.data.frame(cbind(Metric,t(hit.dos.recall)),stringsAsFactors=FALSE)
names(out.recall)=c("Metric",as.character(individs))
library(xtable)
xtable(out.recall,include.rownames=FALSE)

# write.csv(out.recall,file="dos_hitday_recall.csv",row.names=FALSE)



hit.dos.recall

hit.dos.F1=matrix(NA,nr=n.vit,nc = 9)
#col1 is total number of days with hits prior to day before vit drop
#col2 is the total number of days prior to vit drop
#col3 is the ratio of hits to total possible days
#col4 is the logical, whether the the algorithm hit +/- 1 day of vit drop
#col5 is the logical, whether the the algorithm hit +/- 3 days of vit dropth
#col6 is the total number of hits within 1 days
#col7 is the total number of hits within 3 days
#col8 is the total number of hits on vit drop day
#col9 is the total number of hits prior to vit drop day
for(j in 1:n.vit){
    hit.dos.F1[j,1]=length(unique(fit.F1$alarm[[j]]$julian[fit.F1$alarm[[j]]$julian <(d.vit$juliandropped[j]-1)]))
    hit.dos.F1[j,2]=(d.vit$juliandropped[j]-1)-126+1
    hit.dos.F1[j,3]=hit.dos.F1[j,1]/hit.dos.F1[j,2]
    hit.dos.F1[j,4]=sum(fit.F1$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.F1$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)>0
    hit.dos.F1[j,5]=sum(fit.F1$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.F1$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)>0
    hit.dos.F1[j,6]=sum(fit.F1$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.F1$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)
    hit.dos.F1[j,7]=sum(fit.F1$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.F1$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)
    hit.dos.F1[j,8]=sum(fit.F1$alarm[[j]]$julian ==d.vit$juliandropped[j])
    hit.dos.F1[j,9]=sum(fit.F1$alarm[[j]]$julian <(d.vit$juliandropped[j]-1))
}
hit.dos.F1

# write.csv(hit.dos.F1,file="dos_summary_hitday_F1.csv")

hit.recall = hit.dos.recall
hit.prec = hit.dos.prec
hit.F1 = hit.dos.F1




##############################################################################################
#SAVE working directory

save.image("anomaly_v3.Rdata")

##############################################################################################
detect.quant.test.v1=matrix(NA,nr=36,nc=16)
detect.quant.test.v2=matrix(NA,nr=36,nc=16)
detect.quant.test.v3=matrix(NA,nr=36,nc=16)
detect.quant.test.v4=matrix(NA,nr=36,nc=16)

for(k in 1:36){
  detect.quant.test.v1[k,]=fit.loo.v1$loo.fittest[[k]]$detect.quant
  detect.quant.test.v2[k,]=fit.loo.v2$loo.fittest[[k]]$detect.quant
  detect.quant.test.v3[k,]=fit.loo.v3$loo.fittest[[k]]$detect.quant
  detect.quant.test.v4[k,]=fit.loo.v4$loo.fittest[[k]]$detect.quant
}
detect.quant.test.v1
detect.quant.test.v2
detect.quant.test.v3
detect.quant.test.v4







length(fit.loo$loo.fittest[[1]]$detect.quant)

cbind(c(">=pw","<vd","nT","<=pw"),rbind(detect.quant.test[1,],detect.quant.test.v2[1,],detect.quant.test.v3[1,],detect.quant.test.v4[1,]))
