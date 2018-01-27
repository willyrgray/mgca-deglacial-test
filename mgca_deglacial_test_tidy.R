#packages 'pangaear', 'seacarb' and 'mgcv' need to be installed
require(pangaear)
require(seacarb)
require(mgcv)

#to calculate whole surface-ocean change in pH need atmospheric co2 data here I'm using the data from wais divide 
co2_dat<- read.table("ftp://ftp.ncdc.noaa.gov/pub/data/paleo/icecore/antarctica/wais2014co2.txt",header=TRUE)
co2_dat<- data.frame(age= co2_dat$age_calBP/1000, co2= co2_dat$CO2_ppm)

#calculate pH assuming an alkalinity of 2350 umol/kg, temperature of 28 celcius, and salinity of 34.5
alk <- 0.002350
c<-carb(flag=24,var1=co2_dat$co2,var2=alk,T=28,S=34.5,P=2.5,Pt=0.000002,Sit=0.00003,k1k2="m06",kf="dg",ks="d",pHscale="SWS",b="u74")
co2_dat$pH <- c$pH

#now model the pH data as a function of age using a GAM in 'mgcv' package
#later this will be used to predict the pH at different time intervals 
#there are some great posts on Gavin Simpson's blog about GAMs such as https://www.fromthebottomoftheheap.net/2011/06/12/additive-modelling-and-the-hadcrut3v-global-mean-temperature-series/
m_pH<-gam(pH~s(age,k=30),data=co2_dat)
p<- data.frame(age=seq(9,23,by=0.1))
p_pH<-predict(m_pH, newdata=p)

#plot up pH over deglaciaton
dev.new(width=6, height=4.5); par(mar=c(4,4,1,1),xpd=TRUE)
plot(p$age,p_pH,type='l', ylim=c(8.15,8.30), xlab='Age (ka)', ylab='pH (seawater scale)', col='black')

print(round(as.numeric(p_pH[match(20,p$age)]-p_pH[1]), digits=2))
#0.1

#convert to Dsst using -8% per 0.1 pH unit sensitivity and 6% per degree C temperature sensitivity
DpH<-p_pH-p_pH[1]
DT_pH <- -0.83*DpH/0.06

#plot apparent T change from pH
dev.new(width=6, height=4.5); par(mar=c(4,4,1,1),xpd=TRUE)
plot(p$age,DT_pH,type='l', ylim=c(-1.75,0.75), xlab='Age (ka)', ylab='apparent ∆T (°C)', col='red3')

#now salinity
#sealevel curve from sprat et al 2014
esl_dat<- read.table("https://www1.ncdc.noaa.gov/pub/data/paleo/contributions_by_author/spratt2016/spratt2016.txt",header=TRUE)
esl_dat<- data.frame(age= esl_dat$age_calkaBP, esl = esl_dat$SeaLev_shortPC1)
esl_dat<- subset(esl_dat, age < 24)
#model sea level as a function of age
m_esl<-gam(esl~s(age,k=10),data=esl_dat)
p_esl<-predict(m_esl, newdata=p)

#whole ocean salinity from adkins et al 2002
desl<- max(esl_dat$esl)-min(esl_dat$esl)
dsesl = 1.15
Ds<-  -p_esl*(dsesl/desl)

#plot apparent T change from S
DT_s<- 0.033*Ds/0.06
lines(p$age,DT_s,type='l', col='blue3')

#plot apparent T change from combined pH and S
DT_pH_s <-  DT_pH + DT_s
lines(p$age,DT_pH_s,type='l', col='grey17')

legend('bottomleft', lty=1, col=c('red3', 'blue3', 'grey17'), c('pH', 'salinity', 'pH + salinity'), bty='n')


#now to apply to real data
#first import mgca data from pagaea
#for help with pangaear see https://www.fromthebottomoftheheap.net/2016/12/16/pangaea-r-open-palaeo-data/
#mgca data from Mohtadi et al 2014
res<- pg_data(doi = "10.1594/PANGAEA.833322")
dat<- res[[1]]$data
dat<- dat[,1:3]
names(dat)<- c('age','d18o','sst')
dat<- data.frame(age=as.numeric(dat$age),d18o=as.numeric(dat$d18o),sst=as.numeric(dat$sst))
#data is listed as a sst rather than the raw mg/ca (not a fan of this) so we need to back out the mg/ca by solving the paleotemperature equation used in that study (anand 2003) for mg/ca
dat$mgca<- 0.38*exp(0.09*dat$sst)
#next we select only the deglacial data
dat<-subset(dat, age >9 & age <23)

#predict pH for each time interval using the GAM
dat$pH <- predict(m_pH, newdata= dat)

#predict sea level from GAM and convert to absolute salinity using modern value
modern_salinity<- 33.9
esl<-predict(m_esl, newdata=dat)
dat$salinity_sl<-  modern_salinity-esl*(dsesl/desl)

#now to plot up data
dev.new(width=6, height=4.5); par(mar=c(4,4,1,1))
plot(1,1,xlim=c(9,23),ylim=c(24,30),xaxp=c(10,22,6),type="p",pch=1,col=adjustcolor(col="grey99",alpha=0.01),xlab="Age (ka)",ylab="SST (°C)")

#t_mgca_dek_anand
dat$sst_dek_anand<- ((1/0.09)*log(dat$mgca/0.38))
#points(dat$age,dat$sst_dek_anand,type="p",pch=3,col=adjustcolor(col="grey27",alpha=0.3),cex=0.7)
m<-gam(sst_dek_anand~s(age,k=20),data=dat)
dat$sst_dek_anand_gam<-predict(m,dat)
points(dat$age,dat$sst_dek_anand_gam,type="l",col=adjustcolor("grey27",alpha=0.9),lwd=1.25,lty=2)

#t_mgca no pH, no S
dat$sst_nopH_noS<- (1/0.059759)*log(dat$mgca/exp(0.03313*modern_salinity-0.83066*(p_pH[1]-8)-1.06917))
#points(dat$age,dat$sst_nopH_noS,type="p",pch=1,col=adjustcolor(col="grey37",alpha=0.5),cex=1)
m<-gam(sst_nopH_noS~s(age,k=20),data=dat)
dat$sst_nopH_noS_gam<-predict(m,dat)
points(dat$age,dat$sst_nopH_noS_gam,type="l",col=adjustcolor("grey47",alpha=0.9),lwd=1.25)

#t_mgca no pH, S
dat$sst_nopH_S<- (1/0.059759)*log(dat$mgca/exp(0.03313*dat$salinity-0.83066*(p_pH[1]-8)-1.06917))
#points(dat$age,dat$sst_nopH_S,type="p",pch=1,col=adjustcolor(col="blue3",alpha=0.5),cex=1)
m<-gam(sst_nopH_S~s(age,k=20),data=dat)
dat$sst_nopH_S_gam<-predict(m,dat)
points(dat$age,dat$sst_nopH_S_gam,type="l",col=adjustcolor("blue3",alpha=0.9),lwd=1.25)

#t_mgca pH, no S
dat$sst_pH_noS<- (1/0.059759)*log(dat$mgca/exp(0.03313*modern_salinity-0.83066*(dat$pH-8)-1.06917))
#points(dat$age,dat$sst_pH_noS, type="p",pch=1,col=adjustcolor(col="red3",alpha=0.5))
m<-gam(sst_pH_noS~s(age,k=20),data=dat)
dat$sst_pH_noS_gam<-predict(m,dat)
points(dat$age,dat$sst_pH_noS_gam,type="l",col=adjustcolor("red3",alpha=0.9),lwd=1.25)


#t_mgca pH, S
dat$sst_pH_S<- (1/0.059759)*log(dat$mgca/exp(0.03313*dat$salinity-0.83066*(dat$pH-8)-1.06917))
#points(dat$age,dat$sst_pH_S, type="p",pch=1,col=adjustcolor(col="black",alpha=0.5))
m<-gam(sst_pH_S~s(age,k=20),data=dat)
dat$sst_pH_S_gam<-predict(m,dat)
points(dat$age,dat$sst_pH_S_gam,type="l",col=adjustcolor("black",alpha=0.9),lwd=1.25)

legend('topright', lty=c(2,1,1,1,1), col=c('grey27','grey47','blue3', 'red3', 'black' ), c('Dekens02/Anand03','no salinity + no pH', 'salinity + no pH', 'no salinity + pH', 'salinity + pH'), bty='n')
