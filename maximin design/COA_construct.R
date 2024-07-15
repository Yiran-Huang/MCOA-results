#COA_construct2: core function. Generate COA from algorithm 2 in Yang et al. (2021).

cartesian<-function(list_all){
  output=cbind(list_all[[1]])
  s=output
  if(length(list_all)>1){
    for(i in 2:length(list_all)){
      v=rep(c(list_all[[i]][1]),length(s[,1]))
      for(j in 1:(length(list_all[[i]])-1)){
        output<-rbind(output,s)
        v<-c(v,rep(c(list_all[[i]][j+1]),length(s[,1])))
      }
      output<-cbind(output,v)
      s=output
    }}
  return(output)
}

check_pair<-function(x,m=NULL){
  if(is.null(m)){m=(max(x)+1)}
  if(min(x)<0){return(FALSE)}
  if(max(x)>(m-1)){return(FALSE)}
  s=table(paste(x[,1],x[,2]))
  if(length(s)!=m^2){return(FALSE)}else{{return(TRUE)}}
}

check_COA<-function(x){
  m=length(x[1,])
  if(min(x)<0){return(FALSE)}
  if(max(x)>(m-1)){return(FALSE)}
  lambda=length(x[,1])/m/(m-1)
  flag=TRUE
  for(i in 1:(m-1)){
    for (j in (i+1):m) {
      s=table(paste(x[,i],x[,j]))
      if(length(s)!=m*(m-1)){flag=FALSE;break}
      if(sum(s==lambda)!=m*(m-1)){flag=FALSE;break}
    }
    if(!flag){break}
  }
  if(flag){return(TRUE)}else{return(FALSE)}
}

check_COA_strenth3<-function(x){
  m=length(x[1,])
  if(min(x)<0){return(FALSE)}
  if(max(x)>(m-1)){return(FALSE)}
  t=3
  least=prod(m:(m-t+1))
  lambda=length(x[,1])/least
  flag=TRUE
  for(i in 1:(m-2)){
    for (j in (i+1):(m-1)) {
      for(k in (j+1):m){
        s=table(paste(x[,i],x[,j],x[,k]))
        if(length(s)!=least){flag=FALSE;break}
        if(sum(s==lambda)!=least){flag=FALSE;break}
      }
      if(!flag){break}
    }
    if(!flag){break}
  }
  if(flag){return(TRUE)}else{return(FALSE)}
}



is.prime<-function(x){
  upper=floor(x^0.5)
  if(upper<2){return(TRUE)}
  return(sum((x %% (2:upper))==0)==0)
}

multiple_poly<-function(a,b){
  if(length(a)==1 && length(b)==1){return(a*b)}
  if(length(a)>length(b)){
    b=c(b,rep(0,length(a)-length(b)))
  }else{a=c(a,rep(0,length(b)-length(a)))}
  p=length(a)
  n=2*(p-1)+1
  x=rep(0,n)
  for(i in 1:p){
    for(j in 1:i){
      x[i]=x[i]+a[j]*b[i+1-j]
    }
  }
  for(i in (p+1):n){
    s=n-i+1
    for(j in 1:s){
      x[i]=x[i]+a[p+1-j]*b[p-s+j]
    }
  }
  for(i in n:1){
    if(x[i]!=0){break}
  }
  x=x[1:i]
  return(x)
}

inver_element_in_prime_galois<-function(x,m){
  for(i in 1:(m-1)){
    if(((i*x)%%m)==1){
      return(i)}}}

divid_poly_in_galois<-function(a,b,m){
  for(i in length(a):1){
    if(a[i]!=0){break}
  }
  a=a[1:i]
  for(i in length(b):1){
    if(b[i]!=0){break}
  }
  b=b[1:i]
  #return a/b under module m
  leb=length(b)
  lea=length(a)
  s=lea-leb+1
  if(s<=0){
    return(a%%m)
  }
  a=a%%m
  b=b%%m
  inver=inver_element_in_prime_galois(b[leb],m)
  for(i in 1:s){
    loc=lea+1-i
    times=inver*a[loc]
    a[(loc-leb+1):loc]=a[(loc-leb+1):loc]-b*times
    a=a%%m
  }
  for(i in lea:1){if(a[i]!=0){return(a[1:i])}}
  return(a)
}

find_irreducible_poly_in_galois<-function(m,s){
  if(!is.prime(m)){print("***********m is not a prime number***********");return(NULL)}
  if(s==1){return(t(c(0,1)))}
  galois_field=cartesian(rep(list(0:(m-1)),s))
  candidate_list=cbind(galois_field,1)
  for(i in 1:m){
    times=c()
    for(j in 1:(s+1)){times=c(times,(i^{j-1})%%m)}
    a=colSums(t(candidate_list)*times)%%m
    if((sum(a!=0)==1)){candidate_list=t(as.matrix(candidate_list[a!=0,]));break}
    candidate_list=candidate_list[a!=0,]
  }
  for(i in 1:length(candidate_list[,1])){
    flag=TRUE
    for(j in (m+1):length(galois_field[,1])){
      module=divid_poly_in_galois(candidate_list[i,],galois_field[j,],m)
      if(sum(module!=0)==0){flag=FALSE;break}
    }
    if(flag==TRUE){return(candidate_list[i,])}
  }
}

galois_2_num<-function(x,m,s=NULL,ordering=NULL){
  if(is.null(s)){s=length(x[1,])}
  if(is.null(ordering)){ordering=1:s}
  n=m^s
  k=c()
  for(i in 1:s){
    k=c(k,m^(i-1))
  }
  k=k[ordering]
  k=t(t(k))
  times=length(x[1,])/s
  loc=matrix(nrow=length(x[,1]),ncol=times)
  for(i in 1:times){
    loc[,i]=x[,((i-1)*s+1):(i*s)]%*%k+1
  }
 return(loc) 
}

Latin_square<-function(m){
  output=c()
  for(i in 1:m){
    output=rbind(output,(1:m+(i-1))%%m)
  }
  output[which(output==0)]=m
  return(output)
}



COA_construct1<-function(x){
  n=length(x[,1])
  p=length(x[1,])
  m=p-1
  flag=TRUE
  if(n!=(m^2)){print("***********Input array is not an OA(m^2,m+1,2)***********");return(NULL)}
  #check whether it is a OA(m^2,m+1,2)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      flag=check_pair(x,m)
      if(!flag){print("***********Input array is not an OA(m^2,m+1,2)***********");return(NULL)}
    }
  }
  #step 1: recognize 0,....,m-1 in first column. Exchange the order of second column to 0,1,...,m-1
  loc=c()
  for(i in 1:m){
    loc=c(loc,which(x[,1]==(i-1)))
  }
  x=x[loc,]
  loc=c()
  for(i in 1:m){
    loc=c(loc,which(x[1:m,2]==(i-1)))
  }
  x[1:m,]=x[loc,]
  
  #step 2 to 3: Find the permutation in column 3 to m+1
  permutation=numeric(m)
  per_num=as.numeric(t(matrix(rep(1:m,m),m)))-1
  for(i in 3:(m+1)){
    permutation<-x[1:m,i]
    loc=c()
    for(j in 1:m){loc=c(loc,which(x[,i]==permutation[j]))}
    x[loc,i]=per_num
  }
  
  #step 4: Deleting first column and first m rows
  x=x[-(1:m),-1]
  return(x)
}

COA_construct2<-function(m,s,ordering=NULL){
  #Find the irreducible polynomial under m^s with module m
  ir_poly=find_irreducible_poly_in_galois(m,s)
  #step 1: construct the difference matrix, the matrix is ordered by column
  galois_field=cartesian(rep(list(0:(m-1)),s))
  D=matrix(ncol=m^s*s,nrow=(m^s-1))
  for(i in 2:(m^s)){
    for(j in 1:(m^s)){
      d_ij=multiple_poly(galois_field[i,],galois_field[j,])%%m
      d_ij=divid_poly_in_galois(d_ij,ir_poly,m)
      d_ij=c(d_ij,rep(0,s-length(d_ij)))
      D[(i-1),(1+s*(j-1)):(s*j)]=d_ij
    }
  }
  X=c()
  for(i in 1:(m^s)){
    X=cbind(X,((t(D)+galois_field[i,])%%m))
  }
  X=t(X)
  return(galois_2_num(X,m,s,ordering)-1)
}
