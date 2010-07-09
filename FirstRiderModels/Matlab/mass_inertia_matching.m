clear all;
close all;
clc;
syms L1 L2 L3 L4 L5 L6 L7 L8 RR RF xb zb
syms MD MH MJ mb
syms IDCxx IDCyy IDCzz IDCxz IHCxx IHCyy IHCzz IHCxz IJCxx IJCyy IJCzz IJCxz
syms IDxx IDyy IDzz IDxz IHxx IHyy IHzz IHxz IJxx IJyy IJzz IJxz 
syms IBxx IByy IBzz IBxz

IDC=[IDCxx 0 IDCxz; 0 IDCyy 0; IDCxz 0 IDCzz];
IHC=[IHCxx 0 IHCxz; 0 IHCyy 0; IHCxz 0 IHCzz];
IJC=[IJCxx 0 IJCxz; 0 IJCyy 0; IJCxz 0 IJCzz];

a=L1;
b=0;
c=L2-RR;
ID=parallelaxis(IDC, MD, a, b, c)
IDxx = ID(1,1);
IDyy = ID(2,2);
IDzz = ID(3,3);
IDxz = ID(1,3);

a=L4;
b=0;
c=L3+L5-RR;
IH=parallelaxis(IHC, MH, a, b, c)
IHxx = IH(1,1);
IHyy = IH(2,2);
IHzz = IH(3,3);
IHxz = IH(1,3);

a=L7;
b=0;
c=L3+L6+L8-RR;
IJ=parallelaxis(IJC, MJ, a, b, c)
IJxx = IJ(1,1);
IJyy = IJ(2,2);
IJzz = IJ(3,3);
IJxz = IJ(1,3);



ID=[IDxx 0 IDxz; 0 IDyy 0; IDxz 0 IDzz];
IH=[IHxx 0 IHxz; 0 IHyy 0; IHxz 0 IHzz];
IJ=[IJxx 0 IJxz; 0 IJyy 0; IJxz 0 IJzz];


system = cell(7,1);
system{1}=MD + MH + MJ - mb;
system{2}=MD*L1 + MH*L4 + MJ*L7 - mb*xb;
system{3}=MD*L2 + MH*(L3 + L5) + MJ*(L3 + L6 + L8) - mb*(zb + RR);
system{4}=IDxx + IHxx + IJxx - IBxx;
system{5}=IDyy + IHyy + IJyy - IByy;
system{6}=IDzz + IHzz + IJzz - IBzz;
system{7}=IDxz + IHxz + IJxz - IBxz;

% From Meijaard et. al
mb=85;
xb=0.3;
zb=-0.9;
IBxx=9.2;
IByy=11;
IBzz=2.8;
IBxz=2.4;
RR=0.3;
RF=0.35;

% Estimates of mass
MD = 15;
MH = 30;
MJ = 40;
% Estimates of pivot locations
L3=RR-0.28;
L6=-1.05+RR-L3;
% 12 Inertia Components from
% IDCxx IDCyy IDCzz IDCxz IHCxx IHCyy IHCzz IHCxz IJCxx IJCyy IJCzz IJCxz
IDCxx =  0.1383;
IDCyy =  0.3132;
IDCzz =  0.1781;
IDCxz =  0.0380;
IHCxx =  1.7592;
IHCyy =  1.6478;
IHCzz =  0.6885;
IHCxz = -0.3849;
IJCxx =  2.0261;
IJCyy =  1.6851;
IJCzz =  1.1103;
IJCxz =  0.0875;
% Estimates of the 4 inertia components for each of the three bodies

L1=0.41003;
L2 = RR-0.58801;

p1=eval(system{1})
p2=eval(system{2})
p3=eval(system{3})
p4=eval(system{4})
p5=eval(system{5})
p6=eval(system{6})
p7=eval(system{7})

solve(p2,p3,p4,p5,p6,p7,'L4,L5,L7,L8')

% solve('15*L1+30*L4+40*L7-51/2',...
%       '15*L2+108/5+30*L5+40*L8',...
%       '-13191/2500+15*(L2-3/10)^2+30*(-7/25+L5)^2+40*(-21/20+L8)^2',...
%       '-73539/10000+15*(L2-3/10)^2+15*L1^2+30*(-7/25+L5)^2+30*L4^2+40*(-21/20+L8)^2+40*L7^2',...
%       '-8231/10000+15*L1^2+30*L4^2+40*L7^2',...
%       '-13297/5000-15*L1*(L2-3/10)-30*L4*(-7/25+L5)-40*L7*(-21/20+L8)',...
%       'L1,L2,L4,L5,L7,L8')