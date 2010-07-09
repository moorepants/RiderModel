% File: brr_eig_05.m
% Date: January 22, 2008
% Author: Jason Moore
clear all
close all
clc
global   C G IBXX IBXZ IBYY IBZZ IFXX IFYY IHXX IHXZ IHYY IHZZ IRXX IRYY LAMBDA M_B M_F M_H M_R RF RR W XB XH ZB ZH;
%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
C                               =  0.08;                   % M                   Constant
G                               =  9.81;                   % N/KG                Constant
IBXX                            =  9.2;                    % KG*M^2              Constant
IBXZ                            =  2.4;                    % KG*M^2              Constant
IBYY                            =  11;                     % KG*M^2              Constant
IBZZ                            =  2.8;                    % KG*M^2              Constant
IFXX                            =  0.1405;                 % KG*M^2              Constant
IFYY                            =  0.28;                   % KG*M^2              Constant
IHXX                            =  0.05892;                % KG*M^2              Constant
IHXZ                            = -0.00756;                % KG*M^2              Constant
IHYY                            =  0.06;                   % KG*M^2              Constant
IHZZ                            =  0.00708;                % KG*M^2              Constant
IRXX                            =  0.0603;                 % KG*M^2              Constant
IRYY                            =  0.12;                   % KG*M^2              Constant
LAMBDA                          =  0.3141592653589793;     % RAD                 Constant
M_B                             =  85;                     % KG                  Constant
M_F                             =  3;                      % KG                  Constant
M_H                             =  4;                      % KG                  Constant
M_R                             =  2;                      % KG                  Constant
RF                              =  0.35;                   % M                   Constant
RR                              =  0.3;                    % M                   Constant
W                               =  1.02;                   % M                   Constant
XB                              =  0.3;                    % M                   Constant
XH                              =  0.9;                    % M                   Constant
ZB                              = -0.9;                    % M                   Constant
ZH                              = -0.7;                    % M                   Constant

%-------------------Construct the Stability Matrix and Calculate Eigenvalues for
%-------------------Various Velocities
delta = 1e-11; % perturbance value
vmax = 10; % max foward velocity of the rear wheel to be calculated
n = 100; % number of iterations
for i=1:n  
    v(i) = (i-1)/n*vmax; % ith velocity
    NU5(i) = -v(i)/RR; % ith angular velocity of the rear wheel 
    % "bfr_evalprimes" computes the derivatives of the state variables, the
    % equations of motion were generated in Autolev
    % compute nominal solution for the ith velocity
    nominal(:,i) = brr_ep_05([zeros(9,1);NU5(i);0]);
    % build the stability matrix by numerically calculating the partial
    % derivatives of each differential equation with respect to each state
    % variable
    for j=1:11;
        perturb1 = [zeros(9,1);NU5(i);0]; %initialize function input
        perturb2 = [zeros(9,1);NU5(i);0]; %initialize function input
        perturb1(j) = perturb1(j) + delta; %perturb the jth variable
        perturb2(j) = perturb2(j) - delta;
        % solve differential equations for perturbance
        prime1 = brr_ep_05(perturb1);
        prime2 = brr_ep_05(perturb2);
        m(:,j) = (prime1-prime2)./2./delta;  % compute partial derivative
    end
    % reduce stability matrix to 10 x 10 matrix
    N = [4 7 9 11]';
    stab=zeros(4);
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
    plot(v,imag(eigval( i,1:length(v))),'.r')
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


color = colormap(jet(length(eigval(:,1))/2));
for j = 1:length(eigval(:,1))
    figure(j+3)
    hold on
    for i=1:length(eigval(:,1))/2
        plot(squeeze(real(eigvecorg(i,j,:))),squeeze(imag(eigvecorg(i,j,:))),'.','markeredgecolor',color(i,:),'markersize',10)
    end
%     t=min(min((real(eigvecorg(1:5,j,:))))):0.1:max(max((real(eigvecorg(1:5,j,:)))));
%     plot(t,t,t,-t)
    hold off
end
