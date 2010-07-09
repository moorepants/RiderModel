%Author: Jason Moore
%File Name: pcon.m
%Date: November 29, 2007
%Description: Converts the benchmark parameters to our parameter set
clear all
close all
clc
format long

w      =   1.02   ;                  % wheelbase          [m]
c      =   0.08   ;                  % trail              [m]
lambda =   pi/10  ;                  % steer axis tilt    [rad]
rR     =   0.3    ;                  % rear wheel radius  [m]
rF     =   0.35   ;                  % front wheel radius [m]
xB     =   0.3    ;                  % frame COM          [m]
zB     = - 0.9    ;                  % frame COM          [m]
xH     =   0.9    ;                  % fork COM           [m]
zH     = - 0.7    ;                  % fork COM           [m]

D1 = cos(lambda)*(c+w-rR*tan(lambda))
D2 = -cos(lambda)*(rF-rR-w*tan(lambda))
D3 = -c*cos(lambda)+rF*sin(lambda)
L1 = xB
L2 = rR + zB
L3 = cos(lambda)*xH-sin(lambda)*zH-c*cos(lambda)-w*cos(lambda)
L4 = rR*cos(lambda) + xH*sin(lambda) + zH*cos(lambda)

clear w c

w = D1*cos(lambda)+D2*sin(lambda)+D3*cos(lambda)
c = rF*tan(lambda)-D3/cos(lambda)