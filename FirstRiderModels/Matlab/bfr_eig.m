% File: bfr_eig.m
% Date: Feb 24, 2007
% Author: Jason Moore
% Description: 
clear
close all
global C1 C2 C3 D1 D2 D3 G;
global IC11 IC22 IC33;
global ID11 ID12 ID22 ID23 ID31 ID33;
global IF11 IF12 IF22 IF23 IF31 IF33;
global IG11 IG22 IG33;
global IH11 IH12 IH22 IH23 IH31 IH33;
global IJ11 IJ12 IJ22 IJ23 IJ31 IJ33;
global K1 K2 K3;
global L1 L10 L2 L3 L4 L5 L6 L7 L8 L9 LAMBDA;
global MC MD MF MG MH MJ RF RR;
w                               =   1.02   ;                  %m
c                               =   0.08   ;                 %m
lambda                          =   pi/10  ;                  %rad
g                               =   9.81   ;                  %N/kg
v                               =   4.6    ;                  %m/s
rR                              =   0.3    ;                  %m
mR                              =   2      ;                  %kg
IRxx                            =   0.0603 ;                  %kg*m^2
IRyy                            =   0.12   ;                  %kg*m^2
xB                              =   0.3    ;                  %m
zB                              = - 0.9    ;                  %m
mB                              =  85      ;                  %kg
IBxx                            =   9.2    ;                 %kg*m^2
IBxz                            =   2.4    ;                 %kg*m^2
IByy                            =  11      ;                 %kg*m^2
IBzz                            =   2.8    ;                 %kg*m^2
xH                              =   0.9    ;                 %m
zH                              = - 0.7    ;                 %m
mH                              =   4      ;                 %kg
IHxx                            =   0.05892;                 %kg*m^2
IHxz                            = - 0.00756;                 %kg*m^2
IHyy                            =   0.06   ;                %kg*m^2
IHzz                            =   0.00708;                 %kg*m^2
rF                              =   0.35   ;                 %m
mF                              =   3      ;                 %kg
IFxx                            =   0.1405 ;                 %kg*m^2
IFyy                            =   0.28   ;                 %kg*m^2
IH                              = [IHxx 0    IHxz
                                   0    IHyy 0   
                                   IHxz 0    IHzz];
cf                              = [cos(lambda) 0 -sin(lambda)
                                   0           1  0          
                                   sin(lambda) 0  cos(lambda)];
IHrot                           = cf*IH*cf';                             
%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
C1                              =  0;                       % UNITS               Crotch
C2                              =  0;                    % UNITS               Lean
C3                              =  0;                    % UNITS               Twist
D1                              =  w*cos(lambda)+c*cos(lambda)-sin(lambda)*rR;                    % UNITS               Constant
D2                              =  w*sin(lambda)+rR*cos(lambda)-rF*cos(lambda);                    % UNITS               Constant
D3                              =  -c*cos(lambda)+rF*sin(lambda);                   % UNITS               Constant
G                               =  g;                      % UNITS               Constant
IC11                            =  IRxx;                   % UNITS               Constant
IC22                            =  IRyy;                   % UNITS               Constant
IC33                            =  IRxx;                   % UNITS               Constant
ID11                            =  0.1383;                   % UNITS               Constant
ID12                            =  0.0;                    % UNITS               Constant
ID22                            =  0.3132;                    % UNITS               Constant
ID23                            =  0.0;                    % UNITS               Constant
ID31                            =  0.0380;                   % UNITS               Constant
ID33                            =  0.1781;                   % UNITS               Constant
IF11                            =  IHrot(1,1);             % UNITS               Constant
IF12                            =  IHrot(1,2);             % UNITS               Constant
IF22                            =  IHrot(2,2);             % UNITS               Constant
IF23                            =  IHrot(2,3);             % UNITS               Constant
IF31                            =  IHrot(3,1);             % UNITS               Constant
IF33                            =  IHrot(3,3);             % UNITS               Constant
IG11                            =  IFxx;                   % UNITS               Constant
IG22                            =  IFyy;                   % UNITS               Constant
IG33                            =  IFxx;                   % UNITS               Constant
IH11                            =  1.7592;                   % UNITS               Constant
IH12                            =  0.0;                    % UNITS               Constant
IH22                            =  1.6478;                    % UNITS               Constant
IH23                            =  0.0;                    % UNITS               Constant
IH31                            =  -0.3849;                   % UNITS               Constant
IH33                            =  0.6885;                   % UNITS               Constant
IJ11                            =  2.0261;                   % UNITS               Constant
IJ12                            =  0.0;                    % UNITS               Constant
IJ22                            =  1.6851;                    % UNITS               Constant
IJ23                            =  0.0;                    % UNITS               Constant
IJ31                            =  0.08750;                   % UNITS               Constant
IJ33                            =  1.11030;                   % UNITS               Constant
K1                              =  876;                    % UNITS               Crotch
K2                              =  127.8;                    % UNITS               Lean
K3                              =  20.5;                     % UNITS               Twist
L1                              =  w/2;                    % UNITS               Constant
L2                              =  -rR;                 % UNITS               Constant
L3                              =  0.02;                   % UNITS               Constant
L4                              =  w/2;                    % UNITS               Constant
L5                              =  -0.45;                  % UNITS               Constant
L6                              =  -.77;                   % UNITS               Constant
L7                              =  w/2;                 % UNITS               Constant
L8                              =  -0.3;                   % UNITS               Constant
L9                              =  0;                   % UNITS               Constant
L10                             =  0;                   % UNITS               Constant
LAMBDA                          =  lambda;                  % UNITS               Constant
MC                              =  mR;                     % UNITS               Constant
MD                              =  15;                     % UNITS               Constant
MF                              =  mH;                     % UNITS               Constant
MG                              =  mF;                     % UNITS               Constant
MH                              =  30;                     % UNITS               Constant
MJ                              =  40;                     % UNITS               Constant
RF                              =  rF;                     % UNITS               Constant
RR                              =  rR;                     % UNITS               Constant
%-------------------Construct the Stability Matrix and Calculate Eigenvalues for
%-------------------Various Velocities
delta = 1e-11; % perturbance value
vmax = 10; % max foward velocity of the rear wheel to be calculated
n = 1000; % number of iterations
for i=1:n  
    v(i) = (i-1)/n*vmax; % ith velocity
    NU5(i) = -v(i)/RR; % ith angular velocity of the rear wheel 
    % "bfr_evalprimes" computes the derivatives of the state variables, the
    % equations of motion were generated in Autolev
    % compute nominal solution for the ith velocity
    nominal(:,i) = bfr_evalprimes([0;0;0;0;0;0;0;0;0;0;0;NU5(i);0;0;0;0]);
    % build the stability matrix by numerically calculating the partial
    % derivatives of each differential equation with respect to each state
    % variable
    for j=1:16;
        perturb1 = [0;0;0;0;0;0;0;0;0;0;0;NU5(i);0;0;0;0]; %initialize function input
        
        perturb2 = [0;0;0;0;0;0;0;0;0;0;0;NU5(i);0;0;0;0]; %initialize function input
        perturb1(j) = perturb1(j) + delta; %perturb the jth variable
        perturb2(j) = perturb2(j) - delta;
        % solve differential equations for perturbance
        prime1 = bfr_evalprimes(perturb1);
        prime2 = bfr_evalprimes(perturb2);
        m(:,j) = (prime1-prime2)./2./delta;  % compute partial derivative
    end
    % reduce stability matrix to 10 x 10 matrix
    N = [4 6 8 9 10 11 13 14 15 16]';
    stab=zeros(10);
    for k1=1:length(N)
        for k2 = 1:length(N)
            stab(k1,k2)=m(N(k1),N(k2));
        end
    end
    A(:,:,i)=stab;
    % calculate the eigenvalues for the reduced stability matrix
    % calculate the eigenvalues for stability matrix
    [V,D]=eig(stab);
    eigval(1:length(diag(D)),i)=diag(D);
    eigvec(:,:,i)=V;
end
%-------------------Plot the Eigenvalues
figure(1)
hold on
for i=1:length(eigval(:,1))
    plot(v,real(eigval( i,1:length(v))),'.k')
%     plot(v,imag(eigval( i,1:length(v))),'.r')
end
plot(v,zeros(length(v),1),'k')
hold off

% Create color gradient vector
color = colormap(cool(length(v)));
% Plot eigenvalues as a function of speed
figure(2)
hold on
plot([min(min(real(eigval)));max(max(real(eigval)))],[0;0],'k',[0;0],[min(min(imag(eigval)));max(max(imag(eigval)))],'k')
for i=1:length(v)
    plot(real(eigval(:,i)),imag(eigval(:,i)),'.','markeredgecolor',color(i,:),'markersize',5)
end
title('Eigenvalue Loci as a Function of Speed')
xlabel('Re(\lambda) [1/s]')
ylabel('Imag(\lambda) [1/s]')
axis image
box on
hold off
colormap cool
caxis([0 vmax])
colorbar('YTickLabel',...
    {'0 m/s','1','2','3',...
     '4','5','6','7','8','9','10 m/s'})

%-------------------Organize Eigenvalue Matrix
eigorg(:,1)=eigval(:,1);
eigvecorg(:,:,1)=eigvec(:,:,1);
eigred=eigval(:,2);
eigvecred=eigvec(:,:,2);
% Step through each velocity
for i = 2:length(v),
%     Now step through each eigenvalue
    for j = 1:length(eigorg(:,i-1))
        x1 = real(eigorg(j,i-1));
        y1 = imag(eigorg(j,i-1));
        for k = 1:length(eigred)
            x2 = real(eigred(k));
            y2 = imag(eigred(k));
            dist(k) = abs(sqrt((x2-x1)^2+(y2-y1)^2));
        end
        [mindist indice]=min(dist);
        eigorg(j,i)=eigred(indice);
        eigvecorg(:,j,i)=eigvecred(:,indice);
        eigred(indice)=[];
        eigvecred(:,indice)=[];
        clear dist
    end
    if i<length(v),
        eigred=eigval(:,i+1);
        eigvecred=eigvec(:,:,i+1);
    end
end

for i=1:length(eigval(:,1)),
    for j=1:length(eigval(:,1)),
        for k = 2:length(v),
            if (real(eigvecorg(i,j,k)/eigvecorg(i,j,k-1))+1)<0.1,
                eigvecorg(i,j,k)=-eigvecorg(i,j,k);
            end
        end
    end
end

%-------------------Plot the Real Parts of the Eigenvalues
figure(3)
hold on
color = colormap(jet(length(eigval(:,1))));
for i=1:length(eigval(:,1))
    plot(v,real(eigorg(i,1:length(eigorg))),'.','markeredgecolor',color(i,:),'markersize',10)
end
plot(v,zeros(length(v),1),'k')
xlabel('Magnitude of the velocity of the center of the rear wheel [m/s]')
ylabel('Real part of the eigenvalues')
legend()
hold off


% color = colormap(jet(length(eigval(:,1))/2));
% for j = 1:length(eigval(:,1))
%     figure(j+3)
%     hold on
%     for i=1:length(eigval(:,1))/2
%         plot(squeeze(real(eigvecorg(i,j,:))),squeeze(imag(eigvecorg(i,j,:))),'.','markeredgecolor',color(i,:),'markersize',10)
%     end
% %     t=min(min((real(eigvecorg(1:5,j,:))))):0.1:max(max((real(eigvecorg(1:5,j,:)))));
% %     plot(t,t,t,-t)
%     hold off
% end
