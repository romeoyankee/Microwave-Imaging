
%% difference between the measured and calculated scattered field
for i= 1:10
[Escattcal,Z,I,s,C,ETotal,M,N,x,y] = solver(Reconstructed);
 load('Escattmeasured');

 delta_Es = Escattmeasured-Escattcal ;
 
 %% Least square sense
 
Q=norm(delta_Es).^2;
 
 
 %% Calculating the Jacobian Matrix 'D' for each Transmitter
n=8;
data= cell(n,1);

 for m= 1:8
  data{m}=-Z*inv(I+s*C)*diag(ETotal(:,m));
 end
 
 D = cellfun(@(col) vertcat(col{:}), num2cell(data, 1), 'UniformOutput', false); % CONCATENATE THE CELL ARRAY OF 'DATA' TO GET THE MATRIX 'D'
 fnm = sprintf('Data%d.mat');
 save(fnm,'D');
 
 load('Data.mat')

%% calculating the change in Epsr :The G matrix from LM method with Tikhonov Regularization
beta=0.01;
alpha =beta*(trace(ctranspose(D{1})*D{1})/(M*N))*(norm(delta_Es).^2)/(norm(Escattmeasured).^2);
    
G = inv(ctranspose(D{1})*D{1}+alpha*eye(N*N))*ctranspose(D{1});
 
 
%% Change in contrast
 
   delta_Er= G*delta_Es(:);
   

%% new epsr
   water = 76+14.4i;
  Er_new = reshape((diag(s)+delta_Er)+ water,[M,N])  ;
  
  

%% prior information

% object boundary
Er_new(1,1:4)= water; Er_new(1,8:11)= water; 
Er_new(2,1:2)= water; Er_new(2,10:11)= water; 
Er_new(3,1)= water; Er_new(3,11)= water;
Er_new(4,1)= water; Er_new(4,11)= water;
Er_new(8,1)= water; Er_new(8,11)= water;
Er_new(9,1)= water; Er_new(9,11)= water;
Er_new(10,1:2)= water; Er_new(10,10:11)= water;
Er_new(11,1:4)= water; Er_new(11,8:11)= water;

% permittivity levels


%% image reconstruction 
 
  Reconstructed= Er_new;
    
%%  Error
error = sqrt(norm(Q).^2/norm(Escattmeasured).^2)*100

% errorpermittivity=(norm(object-Reconstructed))/(norm(object));
 cnd = cond(ctranspose(D{1})*D{1});
 figure()
subplot(1,2,1), imagesc(x,y,real(Reconstructed)),colorbar,xlabel('12 cm'),ylabel('12cm'),title('Reconstructed - Real part of permittivity')
title(['Estimated $Re(\epsilon)$ with Q=',num2str(Q)],'interpreter', 'latex', 'fontsize', 14)

subplot(1,2,2),imagesc(x,y,imag(Reconstructed)),colorbar,xlabel('12 cm'),ylabel('12cm'),title('Reconstructed - imaginary part of permittivity')
title(['Estimated $Imag(\epsilon)$ with Q=',num2str(Q)],'interpreter', 'latex', 'fontsize', 14)
end