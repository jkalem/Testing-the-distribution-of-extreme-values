#En este archivo calculamos el p-valor para la prueba de hipótesis cuya hipótesis nula es que
#los datos distribuyen según la variable Gumbel.
#El test que aplicamos es del tipo de Cramér von Mises recortado (pero integrando entre a_n
#y b_n siendo a_n=F^(-1)(1/n) b_n=F^(-1)(1-1/n).

#Ingresar m (cantidad de réplicas del test) para calcular el p-valor como el porcentaje
#de veces que el valor del estadístico supera el observado)
#Ingresar x (datos a testear)

n=length(x)
F=function(x){exp(-exp(-x))}
F2=function(x){exp(-2*exp(-x))}

C=n*as.numeric(integrate(F2,-log(n),log(n))[1])

y=(pi/sqrt(6))*(x-mean(x))/sqrt(var(x))-digamma(1)#estandarizo a la Gumbel estandar
#como integro entre -log(n) y log(n), me quedo con los "y" en ese intervalo.
k=sum(y< -log(n))
y=y[y<log(n)]
y=y[y>-log(n)]

y=sort(y)
y=c(-log(n),y,log(n))
j=seq(0,(length(y)-2))
dif_y=y[2:length(y)]-y[1:(length(y)-1)]

B=rep(NA,(length(y)-1))
for (jj in 1:(length(y)-1) )
{
B[jj]=as.numeric(integrate(F,y[jj],y[(jj+1)])[1])

}

t=(1/n)*sum((j+k)^2*dif_y)-2*sum((j+k)*B)+C

T=rep(NA,m)
for (i in 1:m)
{
u=runif(n)
x=-log(-log(u))

y=(pi/sqrt(6))*(x-mean(x))/sqrt(var(x))-digamma(1)#estandarizo a la Gumbel estandar
#como integro entre -log(n) y log(n), me quedo con los "y" en ese intervalo.
k=sum(y< -log(n))
y=y[y<log(n)]
y=y[y>-log(n)]

y=sort(y)
y=c(-log(n),y,log(n))
j=seq(0,(length(y)-2))
dif_y=y[2:length(y)]-y[1:(length(y)-1)]

B=rep(NA,(length(y)-1))
for (jj in 1:(length(y)-1) )
{
B[jj]=as.numeric(integrate(F,y[jj],y[(jj+1)])[1])

}

T[i]=(1/n)*sum((j+k)^2*dif_y)-2*sum((j+k)*B)+C
print(i)

}
pvalue=mean(T>t)
pvalue


