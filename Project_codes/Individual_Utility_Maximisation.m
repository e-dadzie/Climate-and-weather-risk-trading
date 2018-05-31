clc


%Initial value
NY=5;  NC=5;  M=5;ymas =4; cmax=3;cmin=-3;tmax=1.5; sigmas=1; mus1=1;mus2=0; sigmac=0.5; betaN=1;


%Grid Initialization

dy=ymas/NY;
dc=(cmax-cmin)/NC;
dt=tmax/M;

% Discretised point
y= dy:dy:ymas;
c= (cmin+dc):dc:cmax;
t= dt:dt:tmax;


%The Coefficient of the five stencil points

for j=1:NC
  if(-c(j)> 0.5)
    muc1=-c(j); betaN1=betaN; muc2=0; betaN2=0;
  else
    muc1=0; betaN1=0; muc2=-c(j); betaN2=betaN;
  end
  for i=1:NY  
    G(i,j)=1-(dt/dc)*(muc1-sigmac*betaN1)-(dt/dc)*(muc2-sigmac*betaN2)-(dt/(dy)^2)*y(i)^2*sigmas^2-(dt/(dc)^2)*sigmac^2;
    H(i,j)=(dt*y(i)^2*sigmas^2)/(2*(dy)^2);
    I(i,j)=(dt/dc)*(muc1-sigmac*betaN1)+(dt*sigmac^2)/(2*(dc)^2);
    J(i,j)=(dt*y(i)^2*sigmas^2)/(2*(dy)^2);
    K(i,j)=(dt/dc)*(muc2-sigmac*betaN2)+(dt*sigmac^2)/(2*(dc)^2);
  end
end

                  % The matrix Omega
                    

% This is the upper off Diagonal block 
Omega1up=[]; 
for i=1:NY-1
    Omegaup1=diag(H(i,1:NC),0);
    Omega1up=sparse(blkdiag(Omega1up,Omegaup1));
end


%This is the lower off-diagonal block 
Omega1down=[]; 
for i=2:NY
    Omegadown1= diag(J(i,1:NC),0);
    Omega1down=sparse(blkdiag(Omega1down,Omegadown1));
end


% making the matrix Omega1up same as Omega1down

Omegastarup=sparse([zeros((NY-1)*NC,NC), Omega1up;zeros(NC,NY*NC)]); 
Omegastardown=sparse([zeros(NC,NY*NC); Omega1down,zeros((NY-1)*NC,NC)]); 


       % The second matrix A2 is just a diagonal matrix with tridiagonal block
Omega2=[];
for i=1:NY
    Omegadiag=diag(G(i,1:NY),0)+diag(K(i,2:NY),-1)+ diag(I(i,1:NY-1),1);
    Omega2=sparse(blkdiag(Omega2,Omegadiag));
end

Omega=sparse(Omega2 + Omegastardown + Omegastarup);

figure(1);
spy(Omega);

eigs(Omega);


                    %Terminal  Condition
f=zeros(NY*NC,1);
for i=1:NY
   for j=1:NC
       f((i-1)*NC+j)=(1+betaN-(1/2)-(-c(j)/sigmac))*exp(-c(j)/sigmac);      
   end
end

