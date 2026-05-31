##This program calculate the pvalue for the weighted cvm (not truncated)
##to test H0: X has Gumbel distribution for given data set x
##The weighted function is w(t)=1 as proposed in Csörgò and Szabó work (you 
##can change W in the definition weighted function located in the first line).
##Input x (data set)
##Input m (number of replications to calculate the pvalue of the test).

pvalue_gumbel_cvm=function(x,m=1000)
{
  n=length(x)
  W=function(x){x-x+1}
  #W=function(x){1/(x^2+1)}
  F=function(x){exp(-exp(-x))*W(x)}
  F1=function(x){(1-exp(-exp(-x)))^2*W(x)}
  F2=function(x){exp(-2*exp(-x))*W(x)}
  
  
  
  y=(pi/sqrt(6))*(x-mean(x))/sqrt(var(x))-digamma(1)# Gumbel standarisation.
  
  
  y=sort(y)
  C=n*as.numeric(integrate(F2,-Inf,y[n])[1])
  D=n*as.numeric(integrate(F1,y[n],Inf)[1])
  
  j=seq(1,(n-1))
  dif_y=y[2:n]-y[1:(n-1)]
  
  B=rep(NA,(n-1))
  for (jj in 1:(n-1) )
  {
    B[jj]=as.numeric(integrate(F,y[jj],y[(jj+1)])[1])
    
  }
  
  t=(1/n)*sum(j^2*dif_y)-2*sum(j*B)+C+D
  
  T=rep(NA,m)
  for (i in 1:m)
  {
    u=runif(n)
    x=-log(-log(u))
    
    y=(pi/sqrt(6))*(x-mean(x))/sqrt(var(x))-digamma(1)#standardise Gumbel.
    
    
    y=sort(y)
    C=n*as.numeric(integrate(F2,-Inf,y[n])[1])
    D=n*as.numeric(integrate(F1,y[n],Inf)[1])
    
    j=seq(1,(n-1))
    dif_y=y[2:n]-y[1:(n-1)]
    
    B=rep(NA,(n-1))
    for (jj in 1:(n-1) )
    {
      B[jj]=as.numeric(integrate(F,y[jj],y[(jj+1)])[1])
      
    }
    
    T[i]=(1/n)*sum(j^2*dif_y)-2*sum(j*B)+C+D
    #print(i)
    
  }
  pvalue=mean(T>t)
  pvalue
  paste("p-value=",pvalue)
  
  #return(pvalue)
  
}
