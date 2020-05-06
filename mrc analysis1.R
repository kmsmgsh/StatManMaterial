## Real Rata (MRC) analysis
library(foreign)
library(nlme)
library(survival)
library(JM)
library(joineR)
library(dplyr) # near
library(statmod)
library(progress)
library(MASS)
library(mvtnorm)
library(tensor)
library(future.apply)

mrc.or=read.spss("mrc2_25visits_longitudinal.sav")
class(mrc.or)
attributes(mrc.or)
################################
mrc=data.frame(id=mrc.or$myid,trial.id=mrc.or$id,Time=mrc.or$time_censoring,delta=mrc.or$censoring,deathcode=mrc.or$deathcod,age=mrc.or$age,sex=mrc.or$sex,bmi=mrc.or$bmi,smo=mrc.or$smoke,tmt=mrc.or$tmt,
                 time=mrc.or$time,sbp=mrc.or$sbp_nm,month=mrc.or$month,whc=mrc.or$whco)
### creat common observe time (in week scale): mrc.y$week
mrc$week=rep(c(0,0.041,0.082,0.166,0.25,0.5,0.75,seq(1,5.5,by=0.25)),max(mrc$id))

class(mrc)
dim(mrc)
head(mrc)

# mrc.trial=read.csv("mrc_trials_2.csv")
# mrc.trial.T=data.frame(mrc.trial$id_trial,mrc.trial$death,mrc.trial$deathmi,mrc.trial$deathst,mrc.trial$deathcv)
# id.nondeathcv=mrc.trial.T$mrc.trial.id_trial[(mrc.trial.T$mrc.trial.death==1)&(mrc.trial.T$mrc.trial.deathcv==0)]
# unique(mrc$trial.id[(mrc$deathcode==4)|(mrc$deathcode==5)])
###### Delete rows with NA values 
mrc1=mrc[-(which((is.na(mrc$time))|(is.na(mrc$sbp))|(is.na(mrc$bmi)))),]
dim(mrc1)
mrc1=mrc1[-which((mrc1$id)%in%(unique(mrc1$id)[table(unlist(mrc1$id))<10])),]
dim(mrc1)
mrc1$delta=ifelse(mrc1$delta==1,1,0)
## mrc1$deathcode=4 or 5 indicate death due to noncardivascular events
mrc1$delta=ifelse(mrc1$deathcode>3,0,mrc1$delta)
internal_knots=c(0.25,1.5,3.5)
boundary_knots=c(0,5.8)
mrc.M=as.data.frame(mycubicbs(mrc1$time,internal_knots,boundary_knots)$mat)
names(mrc.M)=paste0("time",c(1:7))
mrc1=cbind(mrc1,mrc.M)


mrc1$tmt1=ifelse(mrc1$tmt==1,1,0)
mrc1$tmt2=ifelse(mrc1$tmt==2,1,0)
dim(mrc1)
head(mrc1,n=50)
######################### Longitudinal analysis ########
mrc1$sbp=mrc1$sbp/20
qqnorm(mrc1$sbp*20)
ctl <- lmeControl (msMaxIter=100)
## in this scale (/20 in sbp), random effect variance for whc is very small (0.056), almost not random.
lm2.mrc=lme(fixed=sbp~-1+time1+time2+time3+time4+time5+time6+time7+age+sex+smo+whc+
              tmt1:(time1+time2+time3+time4+time5+time6+time7)+
              tmt2:(time1+time2+time3+time4+time5+time6+time7),
            random =list(id=pdDiag(form=~-1+time2+time3+time4+time5+time6+time7)),data=mrc1,control=ctl) 
summary(lm2.mrc)

f=function(t){
  mycubicbs(t,internal_knots,boundary_knots)$mat
}

lm3.mrc=lme(fixed=sbp~-1+f(time)+age+sex+smo+whc,
            random=list(id=pdDiag(form=~-1+f(time))),data=mrc1,control=ctl)


plot(x=range(mrc1$time),y=c(6,12),type="n")

# predict random effects; same as lm2.mrc$coefficients[[2]]
raneff.pre=data.frame(id=unique(mrc1$id),ranef(lm2.mrc)) 
for (i in 1:10){
  
#points(x=mrc1$time[mrc1$id==i][-c(1:5)],y=mrc1$sbp[mrc1$id==i][-c(1:5)]-mrc1$whc[mrc1$id==i][-c(1:5)]*(lm2.mrc$coefficients[[1]][11]),type="b")
points(x=mrc1$time[mrc1$id==i][-c(1:5)],
       y=lm2.mrc$fitted[mrc1$id==i,2][-c(1:5)]-mrc1$whc[mrc1$id==i][-c(1:5)]*(lm2.mrc$coefficients[[1]][11]),
       col=mrc1$tmt[mrc1$id==i],pch=16,type="b")

}





beta=lm2.mrc$coefficients[[1]][1:7]
beta1=lm2.mrc$coefficients[[1]][1:7]+lm2.mrc$coefficients[[1]][12:18]
beta2=lm2.mrc$coefficients[[1]][1:7]+lm2.mrc$coefficients[[1]][19:25]
sigma2=lm2.mrc$sigma^2

D=diag(c(as.numeric(VarCorr(lm2.mrc)[1:6])))

####################### Cox model
mrc1.id=mrc1[!duplicated(mrc1$id),]
dim(mrc1.id)
Time=mrc1.id$Time
delta=mrc1.id$delta
Tmt1=mrc1.id$tmt1
Tmt2=mrc1.id$tmt2
cox.data=coxph(Surv(Time,delta)~tmt1+tmt2+sex,data=mrc1.id,x=TRUE)
summary(cox.data)# Tmt2 is not significant, we first delete tmt1 and tmt2 in survival submodel

#### Joint model##########
mrc1.joi=jointModel(lm2.mrc, cox.data,
           timeVar = "time1", method = "piecewise-PH-aGH")


##### Cox model with time dependent covariate
mrc1$start=mrc1$time
mrc1$stop=c(mrc1$time[-1],Time[length(Time)])
mrc1$event=numeric(dim(mrc1)[1])
cols=colnames(mrc1)
new_cols=c(cols[1],cols[3],cols[4],cols[11],cols[length(cols)-2],cols[length(cols)-1],cols[length(cols)],cols[12],cols[length(cols)-4],cols[length(cols)-3],cols[6:9],cols[13:22])
mrc1=mrc1[,new_cols]
mrc1$stop=ifelse(mrc1$stop<mrc1$start,mrc1$Time,mrc1$stop)
mrc1$event=ifelse((mrc1$stop==mrc1$Time)&(mrc1$delta==1),1,mrc1$event)
mrc1=mrc1[-which(mrc1$start==mrc1$stop),]
mrc1$sbp=mrc1$sbp/20
cox.extend=coxph(Surv(start,stop,event)~sbp,data=mrc1)
summary(cox.extend)



gamma=cox.data$coefficients[3] #effect of sex on hazard

alpha1=0
alpha2=0
cumbase=rbind(c(0,0),basehaz(cox.data, centered=FALSE))

#####################
### design matrix in Longitudinal data
knots=mycubicbs(0,internal_knots,boundary_knots)$knots
q=7
m=dim(mrc1.id)[1]
Q.2=matrix(0,nrow=q,ncol=q-2)
for(l in 1:(q-2)){
  Q.2[l,l]=6/((knots[l+4]-knots[l+2])*(knots[l+4]-knots[l+1]))
  Q.2[l+2,l]=6/((knots[l+4]-knots[l+2])*(knots[l+5]-knots[l+2]))
  Q.2[l+1,l]=-(Q.2[l,l]+Q.2[l+2,l])
}

des.Y=as.matrix(mrc.M)
des.T=mycubicbs(mrc1.id$Time[mrc1.id$delta==1],internal_knots,boundary_knots)$mat
N=nrow(des.Y)
l=as.data.frame(table(mrc1$id))$Freq
# diXZ=matrix(0,ncol=q*m,nrow=N) ##block diagnal version of XZ
# for(i in 1:m){
#   diXZ[c((cumsum(l)[i]-l[i]+1):(cumsum(l)[i])),c((q*(i-1)+1):(q*i))]=des.Y[c((cumsum(l)[i]-l[i]+1):(cumsum(l)[i])),]
# }
### collection of K(t)=Q.2%*%R.2(t)%*%t(Q.2) for t=data.id$Time[data.id$delta==1]
K.2=array(c(sapply(mrc1.id$Time[mrc1.id$delta==1],function(t) Q.2%*%R.2(t)%*%t(Q.2))),
          dim=c(q,q,length(mrc1.id$Time[mrc1.id$delta==1]))) 

####################
plot(x=c(0,5.8),y=c(140,165),"n")

for(t in seq(0,5.8,by=0.01)){
  points(t,mycubicbs(t,internal_knots=internal_knots,boundary_knots=boundary_knots)$mat%*%(beta),col=1)
}

for(t in seq(0,5.8,by=0.01)){
  points(t,mycubicbs(t,internal_knots=internal_knots,boundary_knots=boundary_knots)$mat%*%(beta),col=1)
}
plot(x=c(-0.5,20.5),y=c(0,1),"n")
for(t in seq(-3,20.5,by=0.01)){
  points(t,mycubicbs(t,internal_knots=internal_knots,boundary_knots=boundary_knots)$mat[,1],col=3)
  points(t,mycubicbs(t,internal_knots=internal_knots,boundary_knots=boundary_knots)$mat[,2],col=4)
  points(t,mycubicbs(t,internal_knots=internal_knots,boundary_knots=boundary_knots)$mat[,7],col=2)
  
  
}


plot(x=c(0,20),y=c(0,10),"n")

for(t in seq(0,5.8,by=0.01)){
  points(t,t(beta-100)%*%Q.2%*%R.2(t)%*%t(Q.2)%*%(beta-100),col=2)
}
