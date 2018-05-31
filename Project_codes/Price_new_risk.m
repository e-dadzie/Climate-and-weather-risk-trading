clc

%Initial value
NY=2;  NC=5;  Mt=50; ymax =4; cmax=3;cmin=-3;tmax=1.5; sigmas=1; mus1=1;mus2=0; sigmac1=0.5; sigmac2=0;


%Grid Initialization

dy=ymax/NY;
dc=(cmax-cmin)/NC;
dt=tmax/Mt;

% Discretised point
y= dy:dy:ymax;
c= (cmin+dc):dc:cmax;
t= dt:dt:tmax;




                                        %Terminal condition                
u=zeros(NY*NC,1);

                                        % Boundary Conditions
Psi=zeros(NY*NC,1);


%The Coefficient of the five stencil points

for j=1:NC
  if (-c(j) > 0)
    muc1=-c(j); 
    muc2=0;
  else
    muc1=0 ;
    muc2=-c(j);
  end  
  for i=1:NY   
    M(i,j)=1-(dt*(muc1-muc2)/dc)-(dt*(y(i)^2)*(sigmas^2)/dy^2)-(dt*(sigmac1^2)/dc^2)-(dt*u((i-1)*NC+j)*0.5*(sigmac1^2-sigmac2^2)/dc);
    P(i,j)=(dt*muc1/dc)+ dt*(sigmac1^2)/(2*(dc^2))+ (dt*0.5*(sigmac1^2)*u((i-1)*NC+j)/dc);
    Q(i,j)=(dt*muc2/dc)+(dt*sigmac1^2/(2*(dc)^2))+(dt*0.5*(sigmac2^2)*u((i-1)*NC+j)/dc);
    N(i,j)=(dt*y(i)^2*sigmas^2)/(2*(dy^2));
    O(i,j)=(dt*y(i)^2*sigmas^2)/(2*(dy^2));
end
end

              % The matrix Omega

% This is the upper off Diagonal block 
Omega1up=[]; 
for i=1:NY-1
    Omegaup1=diag(N(i,1:NC),0);
    Omega1up=sparse(blkdiag(Omega1up,Omegaup1));
end


%This is the lower off-diagonal block 
Omega1down=[]; 
for i=2:NY
    Omegadown1= diag(O(i,1:NC),0);
    Omega1down=sparse(blkdiag(Omega1down,Omegadown1));
end


% making the matrix Omega1up same as Omega1down

Omegastarup=sparse([zeros((NY-1)*NC,NC), Omega1up;zeros(NC,NY*NC)]); 
Omegastardown=sparse([zeros(NC,NY*NC); Omega1down,zeros((NY-1)*NC,NC)]); 



       % The second matrix A2 is just a diagonal matrix with tridiagonal block
Omega2=[];
for i=1:NY
    Omegadiag=diag(M(i,1:NC),0)+diag(Q(i,2:NC),-1)+ diag(P(i,1:NC-1),1);
    Omega2=sparse(blkdiag(Omega2,Omegadiag));
end



Omega=sparse(Omega2 + Omegastardown + Omegastarup);

%figure(1);
%spy(Omega);

%eigs(Omega);


   %TIME DISCRETIZATION
      %Using explicit-upwind finite Difference method
      
%Psi=boundary(cmax,cmin,ymax,tmax,NC,NY,M,thetaE,sigmas,sigmac,mus2,0) % The boundary condition at t=0
 for t=0:dt:tmax
     fstar=Omega*u + Psi; 
     u=fstar;  
 end
u

