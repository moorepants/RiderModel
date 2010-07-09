% Old code when using L1 L2 L3 to describe the bicycle frame
% syms w c rF rR lambda L1 L2 L3
% 
% % w=L1 + L2*sin(lambda) + L3*cos(lambda)
% % c=rR*tan(lambda) - L2*sin(lambda) - L3*cos(lambda)
% % rF=rR - L2*cos(lambda) + L3*sin(lambda)
% 
% 
% s=solve('w=L1 + L2*sin(lambda) + L3*cos(lambda)','c=rR*tan(lambda) - L2*sin(lambda) - L3*cos(lambda)','rF=rR - L2*cos(lambda) + L3*sin(lambda)','L1,L2,L3')
% 
% L1=simple(s.L1)
% L2=simple(s.L2)
% L3=simple(s.L3)
% 
clear all;
close all;
clc;
syms w c d1 d2 d3 RR RF lambda

soln=solve( 'tan(lambda)=d2/(d1+d3+(RR-RF)/(sin(lambda)))',...
            'w=(d1+d3)/cos(lambda)+(RR-RF)*tan(lambda)',...
            'c=(RF-d3/sin(lambda))*tan(lambda)','d1,d2,d3');
d1=simple(soln.d1)
d2=simple(soln.d2)
d3=simple(soln.d3)

% 
% soln=solve( 'tan(lambda)=d2/(d1+d3+(RR-RF)/(sin(lambda)))',...
%             'w=(d1+d3)/cos(lambda)+(RR-RF)*tan(lambda)',...
%             'c=(RF-d3/sin(lambda))*tan(lambda)','w,c,lambda');
% w=simple(soln.w)
% c=simple(soln.c)
% lambda=simple(soln.lambda)


soln=solve( 'tan(lambda)=d2/(d1+d3)',...
            'w=(d1+d3)/cos(lambda)',...
            'c=(RF-d3/sin(lambda))*tan(lambda)','w,c,lambda');
w=simple(soln.w)
c=simple(soln.c)
lambda=simple(soln.lambda)