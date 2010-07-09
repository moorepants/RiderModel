%Author: Jason Moore
%File Name: pcon_syms.m
%Date: November 29, 2007
%Description: takes the output from PCON.AL and simplifies the expressions
clear all
close all
clc

syms lambda rR rF
syms w c xB zB xH zH
syms D1 D2 D3 L1 L2 L3 L4

D1 = cos(lambda)*(c+w-rR*tan(lambda))
D2 = -cos(lambda)*(rF-rR-w*tan(lambda))
D3 = -cos(lambda)*(c-rF*tan(lambda))
L1 = xB
L2 = rR + zB
L3 = cos(lambda)*(xH-(c+w)*cos(lambda)^2) - sin(lambda)*(zH+c*sin...
     (lambda)*cos(lambda)+w*sin(lambda)*cos(lambda))
L4 = rR*cos(lambda) + xH*sin(lambda) + zH*cos(lambda)

D1=simple(D1)
D2=simple(D2)
D3=simple(D3)
L1=simple(L1)
L2=simple(L2)
L3=simple(L3)
L4=simple(L4)