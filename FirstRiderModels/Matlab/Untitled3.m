% File: bfr_eig.m
% Date: Feb 24, 2007
% Author: Jason Moore
% Description: Determines the critical velocites for the weave and capsize
% eigenmodes of a simplified bicycle model made of four rigid bodies.
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
%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
C1                              =  60.4;                   % UNITS               Constant
C2                              =  0.0;                    % UNITS               Constant
C3                              =  0.0;                    % UNITS               Constant
D1                              =  0.8;                    % UNITS               Constant
D2                              =  0.2;                    % UNITS               Constant
D3                              =  0.05;                   % UNITS               Constant
G                               =  9.81;                   % UNITS               Constant
IC11                            =  0.0603;                 % UNITS               Constant
IC22                            =  0.12;                   % UNITS               Constant
IC33                            =  0.0603;                 % UNITS               Constant
ID11                            =  1.84;                   % UNITS               Constant
ID12                            =  0.0;                    % UNITS               Constant
ID22                            =  2.2;                    % UNITS               Constant
ID23                            =  0.0;                    % UNITS               Constant
ID31                            =  0.48;                   % UNITS               Constant
ID33                            =  0.56;                   % UNITS               Constant
IF11                            =  0.0584;                 % UNITS               Constant
IF12                            =  0.0;                    % UNITS               Constant
IF22                            =  0.06;                   % UNITS               Constant
IF23                            =  0.0;                    % UNITS               Constant
IF31                            =  0.0091;                 % UNITS               Constant
IF33                            =  0.0076;                 % UNITS               Constant
IG11                            =  0.1405;                 % UNITS               Constant
IG22                            =  0.28;                   % UNITS               Constant
IG33                            =  0.1405;                 % UNITS               Constant
IH11                            =  3.68;                   % UNITS               Constant
IH12                            =  0.0;                    % UNITS               Constant
IH22                            =  4.4;                    % UNITS               Constant
IH23                            =  0.0;                    % UNITS               Constant
IH31                            =  0.96;                   % UNITS               Constant
IH33                            =  1.12;                   % UNITS               Constant
IJ11                            =  3.68;                   % UNITS               Constant
IJ12                            =  0.0;                    % UNITS               Constant
IJ22                            =  4.4;                    % UNITS               Constant
IJ23                            =  0.0;                    % UNITS               Constant
IJ31                            =  0.96;                   % UNITS               Constant
IJ33                            =  1.12;                   % UNITS               Constant
K1                              =  876;                    % UNITS               Constant
K2                              =  128;                    % UNITS               Constant
K3                              =  20;                     % UNITS               Constant
L1                              =  0.4;                    % UNITS               Constant
L2                              =  -0.275;                 % UNITS               Constant
L3                              =  0.45;                   % UNITS               Constant
L4                              =  0.1;                    % UNITS               Constant
L5                              =  -0.05;                  % UNITS               Constant
L6                              =  -.45;                   % UNITS               Constant
L7                              =  -0.225;                 % UNITS               Constant
L8                              =  -0.8;                   % UNITS               Constant
L9                              =  0.25;                   % UNITS               Constant
L10                             =  -0.2;                   % UNITS               Constant
L11                             =  0.0;                    % UNITS               Constant
L12                             =  -0.1;                   % UNITS               Constant
LAMBDA                          =  pi/10;                  % UNITS               Constant
MC                              =  2;                      % UNITS               Constant
MD                              =  15;                     % UNITS               Constant
MF                              =  4;                      % UNITS               Constant
MG                              =  3;                      % UNITS               Constant
MH                              =  30;                     % UNITS               Constant
MJ                              =  38;                     % UNITS               Constant
RF                              =  0.35;                   % UNITS               Constant
RR                              =  0.3;                    % UNITS               Constant
%-------------------Construct the Stability Matrix and Calculate Eigenvalues for
%-------------------Various Velocities
delta = 1e-11; % perturbance value
vmax = 20; % max foward velocity of the rear wheel to be calculated
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
        perturb1 = [0;0;0;0;0;0;0;0;NU5(i);0]; %initialize function input
        perturb2 = [0;0;0;0;0;0;0;0;NU5(i);0]; %initialize function input
        perturb1(j) = perturb1(j) + delta; %perturb the jth variable
        perturb2(j) = perturb2(j) - delta;
        %solve differential equations for perturbance
        prime1 = bfr_evalprimes(perturb1);
        prime2 = bfr_evalprimes(perturb2);
        m(:,j) = (prime1-prime2)./2./delta;  %compute partial derivative
    end
    % reduce stability matrix to 10 x 10 matrix
    N = [4 6 8 9 10 11 13 14 15 16]';
    for k1=1:length(N)
        for k2 = 1:length(N)
            stab(k1,k2)=m(N(k1),N(k2));
        end
    end
    % calculate the eigenvalues for the reduced stability matrix
    % calculate the eigenvalues for stability matrix
    eigval(1:10,i)=eig(stab);
end
%-------------------Plot the Real Parts of the Eigenvalues
eigvalreal = real(eigval); % extract real parts of eigenvalues
figure(2)
hold on
plot(v,eigvalreal(1,1:length(eigvalreal)),'.r')
plot(v,eigvalreal(2,1:length(eigvalreal)),'.b')
plot(v,eigvalreal(3,1:length(eigvalreal)),'.g')
plot(v,eigvalreal(4,1:length(eigvalreal)),'.m')
plot(v,eigvalreal(1,1:length(eigvalreal)),'.y')
plot(v,eigvalreal(2,1:length(eigvalreal)),'+b')
plot(v,eigvalreal(3,1:length(eigvalreal)),'+g')
plot(v,eigvalreal(4,1:length(eigvalreal)),'+m')
plot(v,eigvalreal(3,1:length(eigvalreal)),'+y')
plot(v,eigvalreal(4,1:length(eigvalreal)),'+r')
plot(v,zeros(length(v),1),'k')
hold off
% %-------------------Organize Eigenvalue Matrix
% % set first column in organized matrix to the first column of the original matrix
% eigorg(:,1)=eigval(:,1);
% % rearrange the eigenvalue columns by calcuating the absolute value of the
% % difference between values of the successicve column and the preceeding
% % column
% for i = 1:length(v)-1
%     for j = 1:4
%         first = eigorg(j,i);
%         for k =1:4
%             second = eigval(k,i+1);
%             diff(k) = abs(second - first);
%         end
%         [mindiff indice]=min(diff);
%         eigorg(j,i+1)=eigval(indice,i+1);
%     end
% end
% %-------------------Extract Weave, Capsize, and Caster Eigenmodes
% found = 0;
% for i = 1:size(eigorg,1)
%     if isreal(eigorg(i,:)) == 0 & found==0
%         weave(1,:)=real(eigorg(i,:));
%         found = 1;
%     elseif isreal(eigorg(i,:)) == 0 & found==1
%         weave(2,:)=real(eigorg(i,:));
%     elseif isreal(eigorg(i,:)) == 1 & (eigorg(i,size(eigorg,2))-eigorg(i,1))>0
%         capsize=real(eigorg(i,:));
%     else
%         caster=real(eigorg(i,:));
%     end
% end
% %-------------------Plot Eigenmodes vs Velocity
% figure(3)
% hold on
% axis([0 10 -10 10])
% title('Eigenmodes')
% xlabel('Velocity [m/s]')
% ylabel('Eigenvalues (Real) [1/s]')
% plot(v,weave(1,:),'.b')
% plot(v,capsize,'.y')
% plot(v,caster,'.r')
% legend('Weave','Capsize','Caster')
% plot(v,weave(2,:),'.b')
% plot(v,zeros(length(v),1),'k') %plot horizontal line at zero
% hold off
% %-------------------Calculate Critical Velocities
% % find the velocity at which the weave mode becomes stable
% for i=1:length(weave(1,:))
%     if weave(1,i)<=0
%         index=i;
%         break
%     end
% end
% vw=mean([v(index-1),v(index)]) % weave critical velocity
% % find the velocity at which the capsize mode becomes unstable
% for i=1:length(capsize)
%     if capsize(i)>=0
%         index=i;
%         break
%     end
% end
% vc=mean([v(index-1),v(index)]) % capsize critical velocity
% %-------------------Plot the Critical Velocities
% hold on
% plot(vw,0,'ok',zeros(1,20)+vw,-19:1:0,'k',vc,0,'ok',zeros(1,20)+vc,-19:1:0,'k')
% hold off