#transform design matrix to location matrix and location matrix to model matrix.


component_2_location<-function (x){
m=length(x[1,]);
n=length(x[,1]);
y=matrix(rep(0,n*m),nrow=n);
for(i in 1:n){
for(j in 1:m){
loc=x[i,j]
y[i,loc]=j
}}
return(y)
}


location_2_CPdesign<-function(x){
m=length(x[1,])
n=length(x[,1])
F=matrix(rep(0,n*((m-1)^2+1)),nrow=n)
F[,1]=1
for(k in 1:n){
for(i in 1:(m-1)){
for(j in 1:(m-1)){
F[k,(i-1)*(m-1)+j+1]=(x[k,i+1]==(j+1))
}}}
return(F)
}

sgn<-function(x){
  s=((2*(x>0)-1)+(-2*(x<0)+1))/2
  return(s)
}



location_2_PWOdesign<-function(x){
  m=length(x[1,])
  n=length(x[,1])
  F=matrix(rep(0,n*((m-1)*m/2+1)),nrow=n)
  F[,1]=1
  for(k in 1:n){
    count=2
    for(i in 1:(m-1)){
      for(j in (i+1):m){
        F[k,count]=sgn(x[k,j]-x[k,i])
        count=count+1
      }}}
  return(F)
}



location_2_firstorder<-function(x){
  m=length(x[1,])
  n=length(x[,1])
  F=matrix(rep(0,n*(m-1)+n),nrow=n)
  F[,1]=1
  for(i in 1:n){
    F[i,-1]=(x[i,1:(m-1)]-(m+1)/2) 
  }
  return(F)
}


location_2_quadratic<-function(x){
  m=length(x[1,])
  n=length(x[,1])
  F=matrix(rep(0,n*(m-1)*2+n),nrow=n)
  F[,1]=1
  for(i in 1:n){
    F[i,-1]=c((x[i,1:(m-1)]-(m+1)/2),(x[i,1:(m-1)]-(m+1)/2)^2-(m^2-1)/12)
  }
  return(F)
}


location_2_secondorder<-function(x){
  m=length(x[1,])
  n=length(x[,1])
  F=matrix(rep(0,n*(m-1)+n*(m-2)+n*(m-1)*(m-2)/2+n),nrow=n)
  F[,1]=1
  for(i in 1:n){
    F[i,2:(2*m-2)]=c((x[i,1:(m-1)]-(m+1)/2),(x[i,1:(m-2)]-(m+1)/2)^2-(m^2-1)/12)
    count=(2*m-2)+1
    for(j in 1:(m-2)){
      for(k in (j+1):(m-1)){
        F[i,count]=(x[i,j]-(m+1)/2)*(x[i,k]-(m+1)/2)
        count=count+1
      }
    }
  }
  return(F)
}

location_2_tPWOdesign<-function(x,taper=NULL){
  m=length(x[1,])
  n=length(x[,1])
  if(is.null(taper)){taper=1/(1:(m-1))}
  F=matrix(rep(0,n*((m-1)*m/2+1)),nrow=n)
  F[,1]=1
  for(k in 1:n){
    count=2
    for(i in 1:(m-1)){
      for(j in (i+1):m){
        F[k,count]=sgn(x[k,j]-x[k,i])*taper[abs(x[k,j]-x[k,i])]
        count=count+1
      }}}
  return(F)
}


location_2_tripletPWOdesign<-function(x){
  m=length(x[1,])
  n=length(x[,1])
  F=matrix(rep(0,n*((m-1)*m/2+1)),nrow=n)
  F[,1]=1
  for(k in 1:n){
    count=2
    for(i in 1:(m-1)){
      for(j in (i+1):m){
        F[k,count]=sgn(x[k,j]-x[k,i])
        count=count+1
      }}}
  
  for(i in 1:(m-2)){
    for(j in (i+1):(m-1)){
      for(k in (j+1):m){
        F=cbind(F,F[,1+(j-i)+(2*m-i)*(i-1)/2]*(F[,1+(k-i)+(2*m-i)*(i-1)/2]+F[,1+(k-j)+(2*m-j)*(j-1)/2]))
      }
    }
  }
  
  return(F)
}
