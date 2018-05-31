function Psi=boundary(cmax,cmin,ymax,tmax,NC,NY,M,thetaE,sigmas,sigmac,mus2,t)    %The function for the known boundary conditions

dy=ymax/NY;
dc=(cmax-cmin)/NC;
dt=tmax/M;

y= dy:dy:ymax;
c= (cmin+dc):dc:cmax;



Psi=zeros(NY*NC,1);
D=zeros(NY,NC);
E=zeros(NY,NC);

for j=1:NC
    for i=1:NY
      if (-c(j) > 0)
        muc2=0;
      else
        muc2= -c(j);
      end  
      D(i,j)=((dt*y(i)*mus2)/dy)+(y(i)*sigmas*dt)/(2*(dy)^2);
      E(i,j)=((dt*muc2)/dc)+(sigmac*dt)/(2*(dc)^2);   
    end
end

Psi(1,1)=D(1,1)*exp(thetaE*t) +E(1,1)*exp(cmin + (thetaE -0.5)*t);   %The first entry of the boundary condition 

%Psi(1,1)=D(1,1)*exp(cmin/sigmac) +E(1,1)*exp(cmin/sigmac); 

for j=2:NC
    %Psi(j,1)=D(1,j)*exp(-c(j)/sigmac);      %The boundary condition when Y is equal to zero
    Psi(j,1)=D(1,j)*exp(thetaE*t);
end

for i=1:NY-1
    %Psi(NC*i+1,1)=E(i+1,1)*exp(cmin/sigmac);       %The boundary condition when the Climate risk(C) is equal to zero
    Psi(NC*i+1,1)=E(i+1,1)*exp(cmin + (thetaE -0.5)*t);
end      
end






