

  function [Escattcal,Z,I,s,C,ETotal,M,N,x,y] = solver(Reconstructed) 
    
%% implementation of method of Moment
    
f=3e9;  % Frequency of wave (in meters) 
c=3e8; % speed of light (in meters) 
lam=c/f*100; % wavelength - remember to multiply it by 100 to convert the units in CM
k=2*pi/lam; %propagation constant
kext=2*pi/(lam/sqrt(76));


  

% defination of the geometry

L1=12;  % length of the square (in centimeters)
L2=12;  % height of the square (in centimeters)

M=11;  % no. of division of square(N=M)

N=11;

x=linspace(-6,6,M+1);

y=linspace(-6,6,N+1);  
    
                                % finding the center point of the cell.

for i = 1:M
    x_mid(i)=(x(i+1)+x(i))/2;
end

for i = 1:N
    y_mid(i)=(y(i+1)+y(i))/2;
end

                            % making the center points for the grid (square geometry)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
xcc1=repmat(x_mid,[1,M]);

ycc1=repmat(y_mid,[N,1]);

ycc1=ycc1(:)';

water = 76+14.4i;

% Integral Equation formulation
Er =Reconstructed(:) ;      % permittivity of the object
s = diag(Er-water) ;                  % constrast in relation to the background(Air)
a=sqrt((L1*L2/(M*N))/pi);   % radius of the equivalent circular cell

theta=0:45:315 ;            % position of the 8 Transmitters for sending plane waves
 
 for m = 1:length(theta)
     
     
  for i = 1:M*N
      
Ei(i)=exp(1j*kext*(xcc1(i)*cosd(theta(m))+ycc1(i)*sind(theta(m))));  

   end
Einc(:,m)= Ei(:);

%% Distance between the center of the cells : pmn

 for i = 1:length(xcc1)
     for j = 1:length(ycc1)
         
         rho_mn(i,j) = sqrt((xcc1(i)-xcc1(j)).^2 +(ycc1(i)-ycc1(j)).^2);
         
     end
 end

%% calculation of the coefficient matrix C(m,n)
C= zeros(M*N,M*N);
I = eye(M*N);
for ii = 1:length(xcc1)
    for jj= 1:length(ycc1)
        
        if ii==jj
            C(ii,jj)=(1j/2)*(pi*kext*a*besselh(1,2,kext*a)-2*1j);
       else
            C(ii,jj)=(1i*pi*kext*a/2)*besselj(1,k*a)*besselh(0,2,kext.*rho_mn(ii,jj));
     
       end
    end
end
 
% cond(C)

 
%  The total field inside the cells   
                                                                                                                                                                                                      
 ETotal(:,m) = inv(I+C*s)*Einc(:,m) ;
 
                                   % calculation of scattered field Es(x,y) 

 ang=0:5.6250:354.3750; % Reciever placed on 64 points on on the circle of Radius of 3.5lambda=3.41cm  
 
 R=7;
 
 X=R*cosd(ang); 
 Y=R*sind(ang);
 

 for kk=1:length(X)
 for jj=1:length(xcc1)
            Const = 1i*(pi*kext/2)*a*besselj(1,kext*a);
            rhon(kk,jj)=sqrt((X(kk)-xcc1(jj)).^2 + (Y(kk)-ycc1(jj)).^2);%distance between reciever points and the object cells
         
            Z(kk,jj) = Const*besselh(0,2,kext*rhon(kk,jj));       
 end
 end
 
 Escattcal(:,m) =- Z*diag(ETotal(:,m))*diag(s);
  
 end 
  end
 
 
  