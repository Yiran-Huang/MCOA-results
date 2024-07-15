rm(list=ls())
script_path <- as.character(sys.frame(1)$ofile)
file_name= file.path(dirname(script_path))
source(paste(file_name,'/COA_construct.R',sep=""))
source(paste(file_name,'/algorithm_on_maximin.R',sep=""))
source(paste(file_name,'/functions.R',sep=""))

set.seed(1810029)
loop_times=10
m=16
source(paste(file_name,'/construct_approximate_PWO.R',sep=""))
load(paste(file_name,'/MCOA_vs_COA2_l1_3.RData',sep=""))
xx3=output16_3[[2]]
x=COA_construct2(2,4)+1
z=component_2_location(D+1)
shjesCOA3_16=c()
shjesMCOA3_16=c()
shjesD3_16=c()
shjesR3_16=c()
shjesRs3_16=c()

t1=Sys.time()
for(loop in 1:loop_times){
  print(paste("m=",m," ",loop,"/",loop_times,sep=""))
  loc31=sample(1:m,m)
  loc32=sample(1:m,m)
  loc33=sample(1:m,m)
  x1<-x
  x3<-rbind(x1,x1[,sample(1:m,m)],x1[,sample(1:m,m)])
  z3=z
  rept=15
  for(i in 1:(rept-1)){
    z3<-rbind(z3,z[,sample(1:16,16)])
  }
  rept2=20
  rs3=t(replicate(m*(m-1)*rept2*3,sample(1:m,m)))
  rept3=5
  rss3=t(replicate(m*(m-1)*rept3*3,sample(1:m,m)))
  w=sort(sample(1:sample(10:30,1),m,replace = TRUE))
  
  shj=c()
  for(j in 1:m){
    if(j==1){shj=c(shj,((w[1]-0)/(m)))
    }else{shj=c(shj,sum((w[1:j]-c(0,w[1:(j-1)]))/(m+1-1:j)))}
  }
  v<-function(x){max(x)}
  esti<-function(x,w){
    m=length(x)
    shjes=c()
    for(j in 1:m){
      loc=which(x<x[j])
      if(length(loc)==0){shjes=c(shjes,v(w[j]))}else{
        shjes=c(shjes,v(w[c(loc,j)])-v(w[c(loc)]))
      }
    }
    return(shjes)
  }
  temp<-function(x){esti(x,w)}
  shjesCOA3_16=rbind(shjesCOA3_16,(colMeans(t(apply(x3[,loc31],1,temp)))-shj)^2)
  shjesMCOA3_16=rbind(shjesMCOA3_16,(colMeans(t(apply(xx3[,loc32],1,temp)))-shj)^2)
  shjesD3_16=rbind(shjesD3_16,(colMeans(t(apply(z3[,loc33],1,temp)))-shj)^2)
  shjesR3_16=rbind(shjesR3_16,(colMeans(t(apply(rs3,1,temp)))-shj)^2)
  shjesRs3_16=rbind(shjesRs3_16,(colMeans(t(apply(rss3,1,temp)))-shj)^2)
}
aa=matrix(c(rbind(shjesMCOA3_16,shjesCOA3_16,shjesR3_16,shjesRs3_16,shjesD3_16)),ncol=m*5)
aa=aa[,(5*m-5*4+1):(5*m)]
pr_out=apply(aa,2,mean)
names(pr_out)=paste(rep(paste("Com ",(m-4+1):m,sep=""),each=5),rep(c("MCOA","COA","R","Rs","D"),4),sep=" ")

#m=23 ---------------------------------------------------------------------------
m=23
load(paste(file_name,'/MCOA_vs_COA5_l1_3.RData',sep=""))
xx3=output23_3[[2]]
x=COA_construct2(23,1)+1
shjesCOA3_23=c()
shjesMCOA3_23=c()
shjesR3_23=c()
shjesRs3_23=c()

for(loop in 1:loop_times){
  print(paste("m=",m," ",loop,"/",loop_times,sep=""))
  loc31=sample(1:m,m)
  loc32=sample(1:m,m)
  loc33=sample(1:m,m)
  x1<-x
  x3<-rbind(x1,x1[,sample(1:m,m)],x1[,sample(1:m,m)])
  rept2=20
  rs3=t(replicate(m*(m-1)*rept2*3,sample(1:m,m)))
  rept3=5
  rss3=t(replicate(m*(m-1)*rept3*3,sample(1:m,m)))
  w=sort(sample(1:sample(10:30,1),m,replace = TRUE))
  
  shj=c()
  for(j in 1:m){
    if(j==1){shj=c(shj,((w[1]-0)/(m)))
    }else{shj=c(shj,sum((w[1:j]-c(0,w[1:(j-1)]))/(m+1-1:j)))}
  }
  v<-function(x){max(x)}
  esti<-function(x,w){
    m=length(x)
    shjes=c()
    for(j in 1:m){
      loc=which(x<x[j])
      if(length(loc)==0){shjes=c(shjes,v(w[j]))}else{
        shjes=c(shjes,v(w[c(loc,j)])-v(w[c(loc)]))
      }
    }
    return(shjes)
  }
  temp<-function(x){esti(x,w)}
  shjesCOA3_23=rbind(shjesCOA3_23,(colMeans(t(apply(x3[,loc31],1,temp)))-shj)^2)
  shjesMCOA3_23=rbind(shjesMCOA3_23,(colMeans(t(apply(xx3[,loc32],1,temp)))-shj)^2)
  shjesR3_23=rbind(shjesR3_23,(colMeans(t(apply(rs3,1,temp)))-shj)^2)
  shjesRs3_23=rbind(shjesRs3_23,(colMeans(t(apply(rss3,1,temp)))-shj)^2)
}
aa=matrix(c(rbind(shjesMCOA3_23,shjesCOA3_23,shjesR3_23,shjesRs3_23)),ncol=m*4)
aa=aa[,(4*m-4*4+1):(4*m)]
pr_out=apply(aa,2,mean)
names(pr_out)=paste(rep(paste("Com ",(m-4+1):m,sep=""),each=4),rep(c("MCOA","COA","R","Rs"),4),sep=" ")


#m=21 ---------------------------------------------------------------------------
mfake=23
m=21

load(paste(file_name,'/MCOA_vs_COA5_l1_3.RData',sep=""))
xx3=output23_3[[2]]
x=COA_construct2(23,1)+1
shjesCOA3_21=c()
shjesMCOA3_21=c()
shjesR3_21=c()
shjesRs3_21=c()
for(loop in 1:loop_times){
  print(paste("m=",m," ",loop,"/",loop_times,sep=""))
  loc31=sample(1:mfake,mfake)
  loc32=sample(1:mfake,mfake)
  loc33=sample(1:mfake,mfake)
  x1<-x
  x3<-rbind(x1,x1[,sample(1:mfake,mfake)],x1[,sample(1:mfake,mfake)])
  rept2=20
  rs3=t(replicate(mfake*(mfake-1)*rept2*3,sample(1:m,m)))
  rept3=5
  rss3=t(replicate(mfake*(mfake-1)*rept3*3,sample(1:m,m)))
  w=sort(sample(1:sample(10:30,1),m,replace = TRUE))
  
  shj=c()
  for(j in 1:m){
    if(j==1){shj=c(shj,((w[1]-0)/(m)))
    }else{shj=c(shj,sum((w[1:j]-c(0,w[1:(j-1)]))/(m+1-1:j)))}
  }
  v<-function(x){max(x)}
  esti1<-function(x,w){
    m=length(w)
    shjes=c()
    for(j in 1:m){
      loc=which(x<x[j])
      if(length(loc)==0){shjes=c(shjes,v(w[j]))}else{
        shjes=c(shjes,v(w[c(loc,j)])-v(w[c(loc)]))
      }
    }
    return(shjes)
  }
  esti2<-function(x,w){
    m=length(w)
    shjes=c()
    for(j in 1:m){
      loc=which(x<x[j])
      loc=setdiff(loc,(m+1):mfake)
      if(length(loc)==0){shjes=c(shjes,v(w[j]))}else{
        shjes=c(shjes,v(w[c(loc,j)])-v(w[c(loc)]))
      }
    }
    return(shjes)
  }
  temp1<-function(x){esti1(x,w)}
  temp2<-function(x){esti2(x,w)}
  shjesCOA3_21=rbind(shjesCOA3_21,(colMeans(t(apply(x3[,loc31],1,temp2)))-shj)^2)
  shjesMCOA3_21=rbind(shjesMCOA3_21,(colMeans(t(apply(xx3[,loc32],1,temp2)))-shj)^2)
  shjesR3_21=rbind(shjesR3_21,(colMeans(t(apply(rs3,1,temp1)))-shj)^2)
  shjesRs3_21=rbind(shjesRs3_21,(colMeans(t(apply(rss3,1,temp1)))-shj)^2)
}
t2=Sys.time()
print(paste("Time: ",t2-t1))


result16=matrix(c(rbind(shjesMCOA3_16,shjesCOA3_16,shjesR3_16,shjesRs3_16,shjesD3_16)),ncol=16*5)
result23=matrix(c(rbind(shjesMCOA3_23,shjesCOA3_23,shjesR3_23,shjesRs3_23)),ncol=23*4)
result21=matrix(c(rbind(shjesMCOA3_21,shjesCOA3_21,shjesR3_21,shjesRs3_21)),ncol=21*4)
result=list(result16,result21,result23)

#Plot
m=16
aa=result[[1]]
aa=aa[,(5*m-5*4+1):(5*m)]
group=rep(paste("m",m," P",(m-4+1):m,sep=""),each=5*nrow(aa))
variable=rep(rep(c("MCOA","COA","R(20)","R(5)","D(15)"),each=nrow(aa)),4)
df1=data.frame(value=c(aa),group=group,variable=variable)
m=21
aa=result[[2]]
aa=aa[,(4*m-4*4+1):(4*m)]
group=rep(paste("m",m," P",(m-4+1):m,sep=""),each=4*nrow(aa))
variable=rep(rep(c("MCOA","COA","R(20)","R(5)"),each=nrow(aa)),4)
df2=data.frame(value=c(aa),group=group,variable=variable)
m=23
aa=result[[3]]
aa=aa[,(4*m-4*4+1):(4*m)]
group=rep(paste("m",m," P",(m-4+1):m,sep=""),each=4*nrow(aa))
variable=rep(rep(c("MCOA","COA","R(20)","R(5)"),each=nrow(aa)),4)
df3=data.frame(value=c(aa),group=group,variable=variable)
df=rbind(df1,df2,df3)

df$variable <- factor(df$variable, levels = c("MCOA","COA","R(20)","R(5)","D(15)"))
df$group <- factor(df$group, levels = c(paste("m",16," P",(16-4+1):16,sep=""),paste("m",21," P",(21-4+1):21,sep=""),
                                        paste("m",23," P",(23-4+1):23,sep="")))

colors <- c("MCOA" = "#F8766D", "COA" = "#A3A500","R(20)"="#00BF7D","R(5)"="#00B0F6","D(15)"="#E76BF3")
y_upper=mean(df$value)+0.9*sd(df$value)
if(mean(df$value>y_upper)==0){y_upper=max(df$value)}
y_lower <- min(df$value)
g=ggplot(df, aes(x = group, y = value, fill = variable)) +
  geom_boxplot(outlier.size = 0.1,width = 0.5, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme() +
  geom_vline(xintercept = (c(setdiff(1:11,c(4,8))*3+0.5)-3.5)/3+1.5, linetype = "dashed", color = "gray50",size = 1) +
  geom_vline(xintercept = (c(c(4,8)*3+0.5)-3.5)/3+1.5, linetype = "solid", color = "black",size = 1) +
  labs(y = "squared error", title = "Boxplots of squared errors for last four players under m=16, 21 and 23")+
  coord_cartesian(ylim = c(y_lower, y_upper))
grid.arrange(g,nrow=1)