rm(list=ls())
script_path <- as.character(sys.frame(1)$ofile)
file_name= file.path(dirname(script_path))
source(paste(file_name,'/COA_construct.R',sep=""))
source(paste(file_name,'/algorithm_on_maximin.R',sep=""))
source(paste(file_name,'/functions.R',sep=""))
source(paste(file_name,'/simulation_MCOA_function.R',sep=""))

set.seed(1810029)
loop_times=10


m=11

source(paste(file_name,'/construct_approximate_PWO.R',sep=""))
load(paste(file_name,'/MCOA_vs_COA1_l1.RData',sep=""))
xx=output11_3[[2]]-1


z=component_2_location(D+1)-1
count=0
fun1_re=c()
fun2_re=c()
fun3_re=c()
t1=Sys.time()
min_dis2=min(dist(z,method='minkowski',p=1))
cumuworktime<-function(x,worktime){
  x=matrix(x,ncol=m)
  n=nrow(x)
  temp<-function(y){
    loc=sort(y,index.return=TRUE)[[2]]
    output=rep(0,m)
    output[loc[1]]=worktime[loc[1]]
    for(i in 2:m){output[loc[i]]=output[loc[i-1]]+worktime[loc[i]]}
    return(output)
  }
  output=(apply(x,1,temp))
  return(output)
}
for(i in 1:loop_times){
  weight=rchisq(m,1)
  worktime=rchisq(m,1)
  #fun放入location matrix
  fun1<-function(x){c(colSums(cumuworktime(x,worktime)^2*weight))}
  fun2<-function(x){c(colSums(cumuworktime(x,worktime)^3*weight))}
  fun3<-function(x){c(colSums(cumuworktime(x,worktime)^4*weight))}
  x<-COA_construct2(11,1)
  x<-rbind(x,x[,sample(1:11,11)],x[,sample(1:11,11)])
  loc1=sample(1:m,m)
  loc2=sample(1:m,m)
  loc3=sample(1:m,m)
  fun1_re=rbind(fun1_re,simulation_function(x[,loc1],xx[,loc2],z[,loc3],fun=fun1))
  fun2_re=rbind(fun2_re,simulation_function(x[,loc1],xx[,loc2],z[,loc3],fun=fun2))
  fun3_re=rbind(fun3_re,simulation_function(x[,loc1],xx[,loc2],z[,loc3],fun=fun3))
  
  count=count+1
  print(paste("finish: ",count,", total: ",loop_times))
  if(count%%5==1){
  output<-list(fun1_re,fun2_re,fun3_re)
  names(output)<-c("fun1","fun2","fun3")
}}
t2=Sys.time()-t1
print(t2)
output<-list(fun1_re,fun2_re,fun3_re)
names(output)<-c("fun1","fun2","fun3")

library(ggplot2)
library(gridExtra)
aa=output[[1]][,1:21]
variable=rep(rep(c("MCOA","COA","D"),each=nrow(aa)),7)
group=rep(c("CP","PWO","First Order","Quadratic","Second Order","tPWO","Triplet"),each=3*nrow(aa))
df=data.frame(value=c(aa),group=group,variable=variable)
df$variable <- factor(df$variable, levels = c("MCOA","COA","D"))
df$group <- factor(df$group, levels = c("CP","First Order","Quadratic","Second Order","PWO","tPWO","Triplet"))
colors <- c("MCOA" = "#F8766D", "COA" = "#00BA38", "D" = "#619CFF")
y_upper=mean(df$value)+3.5*sd(df$value)
if(mean(df$value>y_upper)==0){y_upper=max(df$value)}
y_lower <- min(df$value)
g1=ggplot(df, aes(x = group, y = value, fill = variable)) +
  geom_boxplot() +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  geom_vline(xintercept = (c(3.5, 6.5, 9.5, 12.5, 15.5, 18.5)-3.5)/3+1.5, linetype = "dashed", color = "gray25") +
  labs(x = "Models", y="SMSPE",title="Quadratic Penalty")+
  coord_cartesian(ylim = c(y_lower, y_upper))

aa=output[[2]][,1:21]
variable=rep(rep(c("MCOA","COA","D"),each=nrow(aa)),7)
group=rep(c("CP","PWO","First Order","Quadratic","Second Order","tPWO","Triplet"),each=3*nrow(aa))
df=data.frame(value=c(aa),group=group,variable=variable)
df$variable <- factor(df$variable, levels = c("MCOA","COA","D"))
df$group <- factor(df$group, levels = c("CP","First Order","Quadratic","Second Order","PWO","tPWO","Triplet"))
colors <- c("MCOA" = "#F8766D", "COA" = "#00BA38", "D" = "#619CFF")
y_upper=mean(df$value)+2.5*sd(df$value)
if(mean(df$value>y_upper)==0){y_upper=max(df$value)}
y_lower <- min(df$value)
g2=ggplot(df, aes(x = group, y = value, fill = variable)) +
  geom_boxplot() +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  geom_vline(xintercept = (c(3.5, 6.5, 9.5, 12.5, 15.5, 18.5)-3.5)/3+1.5, linetype = "dashed", color = "gray25") +
  labs(x = "Models", y="SMSPE",title="Cubic Penalty")+
  coord_cartesian(ylim = c(y_lower, y_upper))

aa=output[[3]][,1:21]
variable=rep(rep(c("MCOA","COA","D"),each=nrow(aa)),7)
group=rep(c("CP","PWO","First Order","Quadratic","Second Order","tPWO","Triplet"),each=3*nrow(aa))
df=data.frame(value=c(aa),group=group,variable=variable)
df$variable <- factor(df$variable, levels = c("MCOA","COA","D"))
df$group <- factor(df$group, levels = c("CP","First Order","Quadratic","Second Order","PWO","tPWO","Triplet"))
colors <- c("MCOA" = "#F8766D", "COA" = "#00BA38", "D" = "#619CFF")
y_upper=mean(df$value)+2.5*sd(df$value)
if(mean(df$value>y_upper)==0){y_upper=max(df$value)}
y_lower <- min(df$value)
g3=ggplot(df, aes(x = group, y = value, fill = variable)) +
  geom_boxplot() +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  geom_vline(xintercept = (c(3.5, 6.5, 9.5, 12.5, 15.5, 18.5)-3.5)/3+1.5, linetype = "dashed", color = "gray25") +
  labs(x = "Models", y="SMSPE",title="Quartic Penalty")+
  coord_cartesian(ylim = c(y_lower, y_upper))

grid.arrange(g1,g2,g3, nrow = 3)
