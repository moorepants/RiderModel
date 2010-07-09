% File: whipple.m
% Date: June 15, 2007
% Author: Jason Moore
clear
close all
clc
%---------Define Schwab's Parameters
%----General
w     = 1.02;            % wheelbase [m]
t     = 0.08;            % trail [m]
alpha = atan(3);         % head tube angle [m]
g     = 9.81;            % acceleration due to gravity [m/s^2]
%----Rear Wheel
Rrw   = 0.3;             % radius [m]
xrw   = 0;               % mass center x coordinate [m]
yrw   = 0;               % mass center y coordinate [m]
zrw   = -Rrw;            % mass center z coordinate [m]
cmrw  = [xrw;yrw;zrw];   % mass center vector [m]
mrw   = 2;               % mass [kg]
Axx   = 0.06;            % moment of inertia about x axis [kg*m^2]
Ayy   = 0.12;            % moment of inertia about y axis [kg*m^2]
Azz   = 0.06;            % moment of inertia about z axis [kg*m^2]
Ai    = [Axx,   0,   0;
           0, Ayy,   0;
           0,   0, Azz]; % inertia tensor
%----Rear Frame
xrf   = 0.3;             % mass center x coordinate [m]
yrf   = 0;               % mass center y coordinate [m]
zrf   = -0.9;            % mass center z coordinate [m]
cmrf  = [xrf;yrf;zrf];   % mass center vector [m]
mrf   = 85;              % mass [kg]
Bxx   = 9.2;             % moment of inertia about x axis [kg*m^2]
Bxy   = 0;               % product of inertia [kg*m^2]
Bxz   = 2.4;             % product of inertia [kg*m^2]
Byy   = 11;              % moment of inertia about y axis [kg*m^2]
Byz   = 0;               % product of inertia [kg*m^2]
Bzz   = 2.8;             % moment of inertia about z axis [kg*m^2]
Bi    = [Bxx, Bxy, Bxz;
         Bxy, Byy, Byz;
         Bxz, Byz, Bzz]; % inertia tensor
%----Front Frame
xff   = 0.9;             % mass center x coordinate [m]
yff   = 0;               % mass center y coordinate [m]
zff   = -0.7;            % mass center z coordinate [m]
cmff  = [xff;yff;zff];   % mass center vector [m]
mff   = 4;               % mass [kg]
Cxx   = 0.0546;          % moment of inertia about x axis [kg*m^2]
Cxy   = 0;               % product of inertia [kg*m^2]
Cxz   = -0.0162;         % product of inertia [kg*m^2]
Cyy   = 0.06;            % moment of inertia about y axis [kg*m^2]
Cyz   = 0;               % product of inertia [kg*m^2]
Czz   = 0.0114;          % moment of inertia about z axis [kg*m^2]
Ci    = [Cxx, Cxy, Cxz;
         Cxy, Cyy, Cyz;
         Cxz, Cyz, Czz]; % inertia tensor
%----Front Wheel
Rfw   = 0.35;            % radius [m]
xfw   = w;               % mass center x coordinate [m]
yfw   = 0;               % mass center y coordinate [m]
zfw   = -Rfw;            % mass center z coordinate [m]
cmfw  = [xfw;yfw;zfw];   % mass center vector [m]
mfw   = 3;               % mass [kg]
Dxx   = 0.14;            % moment of inertia about x axis [kg*m^2]
Dyy   = 0.28;            % moment of inertia about y axis [kg*m^2]
Dzz   = 0.14;            % moment of inertia about z axis [kg*m^2]
Di    = [Dxx,   0,   0;
           0, Dyy,   0;
           0,   0, Dzz]; % inertia tensor
%---------Convert Schwab's Parameters To Moore's Parameters
wb = w;         % wheelbase [m]
tr = t;         % trail [m]
lambda = alpha; % head tube angle [m]
rr = Rrw;       % rear wheel radius [m]
rf = Rfw;       % front wheel radius [m]
dr = rr*2;      % rear wheel diamter [m]
df = rf*2;      % front wheel diameter [m]
fo = rf*cos(lambda)-tr*sin(lambda); % fork offset [m]
%----rotate to moore's global reference frame
sCm = [1,  0,  0;
       0, -1,  0;
       0,  0, -1]; % direction cosine matrix (Moore relative to Schwab)
rwheelcenter = sCm*cmrw;
framecenter = sCm*cmrf;
forkcenter = sCm*cmff;
fwheelcenter = sCm*cmfw;
rwheelinertia = sCm*Ai*sCm';
frameinertia = sCm*Bi*sCm';
forkinertia = sCm*Ci*sCm';
fwheelinertia = sCm*Di*sCm';
%----rotate fork reference frame through head tube angle
theta = pi/2 - lambda; % complement of the head tube angle [rad]
Cfork = [ cos(theta), 0, sin(theta);
          0,          1,          0;
         -sin(theta), 0, cos(theta) ];
forkinertia = Cfork*forkinertia*Cfork';
%---------Define Autolev Constants
G = g; 
IC11 = rwheelinertia(1,1);
IC22 = rwheelinertia(2,2);
IC33 = rwheelinertia(3,3);
ID11 = frameinertia(1,1);
ID12 = frameinertia(1,2);
ID22 = frameinertia(2,2);
ID23 = frameinertia(2,3);
ID31 = frameinertia(3,1);
ID33 = frameinertia(3,3);
IF11 = forkinertia(1,1);
IF12 = forkinertia(1,2);
IF22 = forkinertia(2,2);
IF23 = forkinertia(2,3);
IF31 = forkinertia(3,1);
IF33 = forkinertia(3,3);
IG11 = fwheelinertia(1,1);
IG22 = fwheelinertia(2,2);
IG33 = fwheelinertia(3,3);
L1 = wb-fo*sin(lambda)-((dr-df)/2+fo*cos(lambda))/tan(lambda);
L2 = ((dr-df)/2+fo*cos(lambda))/sin(lambda);
L3 = fo;
L4 = framecenter(1) - rwheelcenter(1);
L5 = framecenter(3) - rwheelcenter(3);
temp1 = forkcenter(1) - fwheelcenter(1);
temp2 = forkcenter(3) - fwheelcenter(3);
temp3 = temp1^2 + temp2^2;
temp4 = temp1/temp2;
L6 = sqrt((temp3)/(1+(tan(pi/2-theta-atan(temp4)))^2));
L7 = sqrt(temp3-L6^2);
MC = mrw;
MD = mrf;
MF = mff;
MG = mfw;
RF = rf;
RR = rr;
THETA = theta;
autolev_constants = [G;
                     IC11;IC22;IC33;
                     ID11;ID12;ID22;ID23;ID31;ID33;
                     IF11;IF12;IF22;IF23;IF31;IF33;
                     IG11;IG22;IG33;
                     L1;L2;L3;L4;L5;L6;L7;
                     MC;MD;MF;MG;
                     RF;RR;THETA];
%-------------------Construct the Stability Matrix and Calculate
%-------------------Eigenvalues for Various Velocities
delta = 1e-11; % perturbance value
vmax = 10; % max foward velocity of the rear wheel to be calculated
n = 1000; % number of iterations
for i=1:n  
    v(i) = (i-1)/n*vmax; % ith velocity
    NU5(i) = v(i)/RR; % ith angular velocity of the rear wheel 
    % "evalprimes" computes the derivatives of the state variables, the
    % equations of motion were generated in Autolev
    % compute nominal solution for the ith velocity
    nominal(:,i) = evalprimes([0;0;0;0;0;0;0;0;NU5(i);0],autolev_constants);
    % build the stability matrix by numerically calculating the partial
    % derivatives of each differential equation with respect to each state
    % variable
    for j=1:10;
        perturb1 = [0;0;0;0;0;0;0;0;NU5(i);0]; %initialize function input
        perturb2 = [0;0;0;0;0;0;0;0;NU5(i);0]; %initialize function input
        perturb1(j) = perturb1(j) + delta; %perturb the jth variable
        perturb2(j) = perturb2(j) - delta;
        %solve differential equations for perturbance
        prime1 = evalprimes(perturb1,autolev_constants);
        prime2 = evalprimes(perturb2,autolev_constants);
        %compute partial derivative
        matrix(:,j) = (prime1-prime2)./2./delta;
    end
    % reduce stability matrix to 4 x 4 matrix for steer and roll
    stab(1,1:4) = [matrix(4,4)  matrix(4,6)  matrix(4,8)  matrix(4,10) ];
    stab(2,1:4) = [matrix(6,4)  matrix(6,6)  matrix(6,8)  matrix(6,10) ];
    stab(3,1:4) = [matrix(8,4)  matrix(8,6)  matrix(8,8)  matrix(8,10) ];
    stab(4,1:4) = [matrix(10,4) matrix(10,6) matrix(10,8) matrix(10,10)];
    A(:,:,i)=stab;
    % calculate the eigenvalues for the reduced stability matrix
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
