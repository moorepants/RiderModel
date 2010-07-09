% File: bfr_eig.m
% Date: Feb 24, 2007
% Author: Jason Moore
% Description: 
clear
close all
%-------------------Construct the Stability Matrix and Calculate Eigenvalues for
%-------------------Various Velocities
delta = 1e-11; % perturbance value
vmax = 20; % max foward velocity of the rear wheel to be calculated
n = 100; % number of iterations
for i=1:n  
    v(i) = (i-1)/n*vmax; % ith velocity
    nu2(i) = v(i); % ith angular velocity of the rear wheel 
    % "bfr_evalprimes" computes the derivatives of the state variables, the
    % equations of motion were generated in Autolev
    % compute nominal solution for the ith velocity
    nominal(:,i) = brendan_bfr_evalprimes([0;0;0;0;0;0;0;nu2(i);0;0]);
    % build the stability matrix by numerically calculating the partial
    % derivatives of each differential equation with respect to each state
    % variable
    for j=1:10;
        perturb1 = [0;0;0;0;0;0;0;nu2(i);0;0]; %initialize function input
        
        perturb2 = [0;0;0;0;0;0;0;nu2(i);0;0]; %initialize function input
        perturb1(j) = perturb1(j) + delta; %perturb the jth variable
        perturb2(j) = perturb2(j) - delta;
        % solve differential equations for perturbance
        prime1 = brendan_bfr_evalprimes(perturb1);
        prime2 = brendan_bfr_evalprimes(perturb2);
        m(:,j) = (prime1-prime2)./2./delta;  % compute partial derivative
    end
    % reduce stability matrix to 10 x 10 matrix
    N = [4 5 9 10]';
    stab=zeros(10);
    for k1=1:length(N)
        for k2 = 1:length(N)
            stab(k1,k2)=m(N(k1),N(k2));
        end
    end
    % calculate the eigenvalues for the reduced stability matrix
    % calculate the eigenvalues for stability matrix
    [V,D]=eig(stab);
    eigval(1:length(diag(D)),i)=diag(D);
    eigvec(i)={V};
end
%-------------------Plot the Eigenvalues
eigvalreal = real(eigval); % extract real parts of eigenvalues
% eigvalimag = imag(eigval); % extract the imaginary parts
figure(2)
hold on
for i=1:10
    plot(v,eigvalreal( i,1:length(v)),'.k')
%     plot(v,eigvalimag( i,1:length(v)),'.r')
end
plot(v,zeros(length(v),1),'k')
hold off
% %-------------------Organize Eigenvalue Matrix
% % set first column in organized matrix to the first column of the original matrix
% eigorg(:,1)=eigvalreal(:,1);
% % rearrange the eigenvalue columns by calcuating the absolute value of the
% % difference between values of the successicve column and the preceeding
% % column
% for i = 1:length(v)-1
%     for j = 1:10
%         first = eigorg(j,i);
%         for k =1:10
%             second = eigvalreal(k,i+1);
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