####### B-spline basis #########
## my cubic B-spline basis functions ####
mycubicbs=function(s,internal_knots,boundary_knots){
  df=length(internal_knots)+3+1
  r=numeric(df)
  out_knots=c(seq(from=min(boundary_knots),by=-0.1,length=4),seq(from=max(boundary_knots),by=0.1,length=4))
  knots=sort(c(internal_knots,out_knots))
  temp=function(x){
    my0bs=function(x,k){
      a=ifelse(x<knots[k+1]&x>=knots[k],1,0)
      return(a)
    }
    my1bs=function(x,k){
      (x-knots[k])/(knots[k+1]-knots[k])*my0bs(x,k)+
        (knots[k+2]-x)/(knots[k+2]-knots[k+1])*my0bs(x,k+1)
    }
    
    my2bs=function(x,k){
      (x-knots[k])/(knots[k+2]-knots[k])*my1bs(x,k)+
        (knots[k+3]-x)/(knots[k+3]-knots[k+1])*my1bs(x,k+1)
    }
    my3bs=function(x,k){
      (x-knots[k])/(knots[k+3]-knots[k])*my2bs(x,k)+
        (knots[k+4]-x)/(knots[k+4]-knots[k+1])*my2bs(x,k+1)
    }
    for(k in 1:df){
      r[k]=my3bs(x,k)
    }
    return(r)}
  return=list(mat=t(sapply(s,temp)),knots=knots)# a matrix with ncol=df
}

R.2=function(t){
  R.2=matrix(0,ncol=q-2,nrow=q-2)
  for (l in 1:(q-2)){
    if((t>knots[l+2])&(t<=knots[l+3])){
      R.2[l,l]=(t-knots[l+2])^3/(3*(knots[l+3]-knots[l+2])^2)
    }else{
      if((t>knots[l+3])&(t<knots[l+4])){
        R.2[l,l]=(knots[l+4]-knots[l+2])/3-(knots[l+4]-t)^3/(3*(knots[l+4]-knots[l+3])^2)
        R.2[l,l+1]=(-t^3/3+(knots[l+4]+knots[l+3])/2*t^2-knots[l+4]*knots[l+3]*t-(knots[l+3])^3/6+
                      (knots[l+3])^2*knots[l+4]/2)/(knots[l+4]-knots[l+3])^2
        R.2[l+1,l]= R.2[l,l+1]
      }else{
        if(t>=knots[l+4]){
          R.2[l,l]=(knots[l+4]-knots[l+2])/3
          R.2[l,l+1]=(knots[l+4]-knots[l+3])/6
          R.2[l+1,l]=(knots[l+4]-knots[l+3])/6
          
        }
      }
    }
  }
  R.2[1,1]=ifelse(t>=knots[5],(knots[5]-knots[4])/3,(knots[5]-knots[4])/3-(knots[5]-t)^3/(3*(knots[5]-knots[4])^2))
  return(R.2)
}

# true value of longitudinal biomarker (no measurement error)
longit.true=function(t,fix_eff, ran_eff){
  desmat=mycubicbs(t,internal_knots=internal_knots,boundary_knots=boundary_knots)$mat# design matrix
  desmat%*%(fix_eff+ran_eff)
}




## Estimation
lambda0=function(t){
  if(t %in% unique(data$Time[data$delta==1])){
    cumbase[which(cumbase[,2]==t),1]-cumbase[which(cumbase[,2]==t)-1,1]} ##warn: options(digits = 22) default 7
  else{0}
}

riskset=function(t){# individuals in the risk set
  unique(data$id[(data$Time>=t)])## "=" is important
}

