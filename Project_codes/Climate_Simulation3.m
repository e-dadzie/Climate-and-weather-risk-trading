clc


%Initial value
NY=50;  NC=50;  M=50; ymax =4; cmax=5;cmin=-5;tmax=1.5; sigmas=1; mus1=1;mus2=0; sigmac=4.5; thetaE=1;


%Grid Initialization

dy=ymax/NY;
dc=(cmax-cmin)/NC;
dt=tmax/M;

% Discretised point
y= dy:dy:ymax;
c= (cmin+dc):dc:cmax;
t= dt:dt:tmax;




%The Coefficient of the five stencil points
for t=0:dt:tmax
for j=1:NC
  if (-8*(c(j)^3-c(j))-sin(2*pi*t) > 0)
    muc1=-8*(c(j)^3-c(j))-sin(2*pi*t); 
    muc2=0;
  else
    muc1=0 ;
    muc2=-8*(c(j)^3-c(j))-sin(2*pi);
  end  
  for i=1:NY   
    A(i,j)=1-(dt/dy)*y(i)*(mus1+mus2)-((y(i)*dt*sigmas)/(dy)^2)-(dt/dc)*(muc1+muc2)-((dt*sigmac)/(dc)^2);
    B(i,j)=((dt*y(i)*mus1)/dy)+((y(i)*sigmas*dt)/2*(dy)^2);
    C(i,j)=((dt*muc1)/dc)+(sigmac*dt)/(2*(dc)^2);
    D(i,j)=((dt*y(i)*mus2)/dy)+(y(i)*sigmas*dt)/(2*(dy)^2);
    E(i,j)=((dt*muc2)/dc)+(sigmac*dt)/(2*(dc)^2);   
  end
end


              % The matrix Omega

% This is the upper off Diagonal block 
Omega1up=[]; 
for i=1:NY-1
    Omegaup1=diag(B(i,1:NC),0);
    Omega1up=sparse(blkdiag(Omega1up,Omegaup1));
end


%This is the lower off-diagonal block 
Omega1down=[]; 
for i=2:NY
    Omegadown1= diag(D(i,1:NC),0);
    Omega1down=sparse(blkdiag(Omega1down,Omegadown1));
end


% making the matrix Omega1up same as Omega1down

Omegastarup=sparse([zeros((NY-1)*NC,NC), Omega1up;zeros(NC,NY*NC)]); 
Omegastardown=sparse([zeros(NC,NY*NC); Omega1down,zeros((NY-1)*NC,NC)]); 



       % The second matrix A2 is just a diagonal matrix with tridiagonal block
Omega2=[];
for i=1:NY
    Omegadiag=diag(A(i,1:NC),0)+diag(E(i,2:NC),-1)+ diag(C(i,1:NC-1),1);
    Omega2=sparse(blkdiag(Omega2,Omegadiag));
end



Omega=sparse(Omega2 + Omegastardown + Omegastarup);

%figure(1);
%spy(Omega);

%eigs(Omega);



                    %INITIAL CONDITION
f=zeros(NY*NC,1);
for i=1:NY
   for j=1:NC
       f((i-1)*NC+j)=exp(-c(j)/sigmac);      
   end
end



tic

                   %TIME DISCRETIZATION
      %Using explicit-upwind finite Difference method
      
%Psi=boundary(cmax,cmin,ymax,tmax,NC,NY,M,thetaE,sigmas,sigmac,mus2,0) % The boundary condition at t=0
 
     Psi=boundary(cmax,cmin,ymax,tmax,NC,NY,M,thetaE,sigmas,sigmac,mus2,t); %The boundary condition at any time
     fstar=Omega*f + Psi; 
     f=fstar;  
 end

toc

F=zeros(NY*NC,1);
for i=1:NY
   for j=1:NC
       F((i-1)*NC+j)=1 + dt*(thetaE-0.5+(c(j)/sigmac));  
   end
end



fstar1=zeros(NY*NC,1);
for i=1:NY*NC
  fstar1(i,1)=f(i,1)/F(i,1);
end  


%fstar2=max(fstar1,0);
explicit_upwind=reshape(fstar1,NY,NC);         
                               %The surface graph of the expectation of insurance asset
figure(1)
surf(y,c,explicit_upwind);
xlabel('Stock Process');
ylabel('climate Risk');
zlabel('Moment of Insurance Stock');
 
ee = cputime-time             %This update the computational time
 




