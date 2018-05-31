######################fisher model A

T=2; t=0.25
C_t <-seq(-4 ,4 ,0.1)
betaN=4*exp(-0.085*C_t^2)

J= -exp(-1 +(-0.5*(1+betaN)-5*exp(-10*C_t^2))*T 
        -(-0.5*(1+betaN)-5*exp(-10*C_t^2))*t)

plot(C_t, J, type="l", ylim=c(-0.4,0))


T=2; t=0.5
C_t <-seq(-4 ,4 ,0.1)
betaN=4*exp(-0.085*C_t^2)

J1= -exp(-1 +(-0.5*(1+betaN)-5*exp(-10*C_t^2))*T 
        -(-0.5*(1+betaN)-5*exp(-10*C_t^2))*t)

plot(C_t, J1, type="l", ylim=c(-0.4,0))

T=2; t=0.75
C_t <-seq(-4 ,4 ,0.1)
betaN=2*exp(-0.09*C_t^2)

J2= -exp(-1 +(-0.5*(1+betaN)-5*exp(-10*C_t^2))*T 
        -(-0.5*(1+betaN)-5*exp(-10*C_t^2))*t)

plot(C_t, J2, type="l", ylim=c(-0.4,0))


T=2; t=1
C_t <-seq(-4 ,4 ,0.1)
betaN=2*exp(-0.09*C_t^2)

J3= -exp(-1 +(-0.5*(1+betaN)-5*exp(-10*C_t^2))*T 
        -(-0.5*(1+betaN)-5*exp(-10*C_t^2))*t)

plot(C_t, J3, type="l", ylim=c(-0.4,0))

###################bank model A

T=2; t=0.25
C_t <-seq(-4 ,4 ,0.1)
betaN=4*exp(-0.085*C_t^2)

J= -exp(-1 +(-0.5*(1+betaN))*T 
        -(-0.5*(1+betaN))*t)

plot(C_t, J, type="l", ylim=c(-0.4,0))


T=2; t=0.5
C_t <-seq(-4 ,4 ,0.1)
betaN=4*exp(-0.085*C_t^2)

J1= -exp(-1 +(-0.5*(1+betaN)-5*exp(-10*C_t^2))*T 
         -(-0.5*(1+betaN)-5*exp(-10*C_t^2))*t)

plot(C_t, J1, type="l", ylim=c(-0.4,0))

T=2; t=0.75
C_t <-seq(-4 ,4 ,0.1)
betaN=2*exp(-0.09*C_t^2)

J2= -exp(-1 +(-0.5*(1+betaN))*T 
         -(-0.5*(1+betaN))*t)

plot(C_t, J2, type="l", ylim=c(-0.4,0))


T=2; t=1
C_t <-seq(-4 ,4 ,0.1)
betaN=2*exp(-0.09*C_t^2)

J3= -exp(-1 +(-0.5*(1+betaN))*T 
         -(-0.5*(1+betaN))*t)

plot(C_t, J3, type="l", ylim=c(-0.4,0))