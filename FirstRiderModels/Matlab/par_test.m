% parameters from Meijaard, et al 2007
w      =   1.02   ;                  % wheelbase         [m]
c      =   0.08   ;                  % trail             [m]
lambda =   pi/10  ;                  % steer axis tilt   [rad]
g      =   9.81   ;                  % gravity           [N/kg]
v      =   4.60   ;                  % foward speed      [m/s]
rR     =   0.30   ;                  % rear wheel radius [m]
mR     =   2.00   ;                  % rear wheel mass   [kg] 
IRxx   =   0.0603 ;                  %                   [kg*m^2]
IRyy   =   0.12   ;                  %                   [kg*m^2]
xB     =   0.30   ;                  %                   [m]
zB     = - 0.90   ;                  %                   [m]
mB     =  85.00   ;                  %                   [kg]
IBxx   =   9.20   ;                  %kg*m^2
IBxz   =   2.40   ;                  %kg*m^2
IByy   =  11.00   ;                  %kg*m^2
IBzz   =   2.80   ;                  %kg*m^2
xH     =   0.90   ;                  %m
zH     = - 0.70   ;                  %m
mH     =   4.00   ;                  %kg
IHxx   =   0.05892;                  %kg*m^2
IHxz   = - 0.00756;                  %kg*m^2
IHyy   =   0.06   ;                  %kg*m^2
IHzz   =   0.00708;                  %kg*m^2
rF     =   0.35   ;                  %m
mF     =   3.00   ;                  %kg
IFxx   =   0.1405 ;                  %kg*m^2
IFyy   =   0.28   ;                  %kg*m^2
IH     = [IHxx 0    IHxz
          0    IHyy 0   
          IHxz 0    IHzz];
cf     = [cos(lambda) 0 -sin(lambda)
          0           1  0          
          sin(lambda) 0  cos(lambda)];
IHrot  = cf*IH*cf';                             

%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
D1                              =  w*cos(lambda)+c*cos(lambda)-sin(lambda)*rR;                    % UNITS               Constant
D2                              =  w*sin(lambda)+rR*cos(lambda)-rF*cos(lambda);                    % UNITS               Constant
D3                              =  -c*cos(lambda)+rF*sin(lambda);                   % UNITS               Constant
G                               =  g;                      % UNITS               Constant
IC11                            =  IRxx;                   % UNITS               Constant
IC22                            =  IRyy;                   % UNITS               Constant
IC33                            =  IRxx;                   % UNITS               Constant
ID11                            =  IBxx;                    % UNITS               Constant
ID12                            =  0.0;                    % UNITS               Constant
ID22                            =  IByy;                    % UNITS               Constant
ID23                            =  0.0;                    % UNITS               Constant
ID31                            =  IBxz;                    % UNITS               Constant
ID33                            =  IBzz;                    % UNITS               Constant
IF11                            =  IHrot(1,1);             % UNITS               Constant
IF12                            =  IHrot(1,2);             % UNITS               Constant
IF22                            =  IHrot(2,2);             % UNITS               Constant
IF23                            =  IHrot(2,3);             % UNITS               Constant
IF31                            =  IHrot(3,1);             % UNITS               Constant
IF33                            =  IHrot(3,3);             % UNITS               Constant
IG11                            =  IFxx;                   % UNITS               Constant
IG22                            =  IFyy;                   % UNITS               Constant
IG33                            =  IFxx;                   % UNITS               Constant
L1                              =  0.3;                    % UNITS               Constant
L2                              = -0.9;                 % UNITS               Constant
L3                              =  0.026;                    % UNITS               Constant
L4                              = -0.102;                    % UNITS               Constant
LAMBDA                          =  lambda;                  % UNITS               Constant
MC                              =  mR;                     % UNITS               Constant
MD                              =  mB;                    % UNITS               Constant
MF                              =  mH;                     % UNITS               Constant
MG                              =  mF;                     % UNITS               Constant
RF                              =  rF;                     % UNITS               Constant
RR                              =  rR;                     % UNITS               Constant

Q1                              =  0.0;                    % UNITS               Initial Value
Q2                              =  0.0;                    % UNITS               Initial Value
Q3                              =  0.0;                    % UNITS               Initial Value
Q4                              =  0.0;                    % UNITS               Initial Value
Q5                              =  0.0;                    % UNITS               Initial Value
Q6                              =  0.0;                    % UNITS               Initial Value
Q7                              =  0.0;                    % UNITS               Initial Value
Q8                              =  0.0;                    % UNITS               Initial Value
U4                              =  0.5;                    % UNITS               Initial Value
U5                              = -v/RR;                    % UNITS               Initial Value
U7                              =  0.0;                    % UNITS               Initial Value

TINITIAL                        =  0.0;                    % sec                 Initial Time
TFINAL                          =  5.0;                    % sec                 Final Time
INTEGSTP                        =  0.01;                    % sec                 Integration Step
PRINTINT                        =  1;                      % Positive Integer    Print-Integer
ABSERR                          =  1.0E-08;                %                     Absolute Error
RELERR                          =  1.0E-07 ;               %                     Relative Error

% someName = load( 'BRR_04.2' );
% figure(1)
% time = someName(:,1);
% rollrate = someName(:,5);
% steerrate = someName(:,8);
% velocity = -someName(:,6).*0.3;
% [AX,H1,H2] = plotyy(time,rollrate,time,velocity);
% set(get(AX(1),'Ylabel'),'String','Angular Rate [rad/s]')
% set(get(AX(2),'Ylabel'),'String','Forward Velocity [m/s]')
% hold on
% plot(time,steerrate,'r')
% xlabel('Time [s]')
