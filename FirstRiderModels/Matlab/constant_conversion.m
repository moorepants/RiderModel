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
G                               =  g;                      % UNITS               Constant
IC11                            =  IRxx;                   % UNITS               Constant
IC22                            =  IRyy;                   % UNITS               Constant
IC33                            =  IRxx;                   % UNITS               Constant
ID11                            =  IBxx;                   % UNITS               Constant
ID12                            =  0.0;                    % UNITS               Constant
ID22                            =  IByy;                   % UNITS               Constant
ID23                            =  0.0;                    % UNITS               Constant
ID31                            =  IBxz;                   % UNITS               Constant
ID33                            =  IBzz;                   % UNITS               Constant
IF11                            =  IHrot(1,1);             % UNITS               Constant
IF12                            =  IHrot(1,2);             % UNITS               Constant
IF22                            =  IHrot(2,2);             % UNITS               Constant
IF23                            =  IHrot(2,3);             % UNITS               Constant
IF31                            =  IHrot(3,1);             % UNITS               Constant
IF33                            =  IHrot(3,3);             % UNITS               Constant
IG11                            =  IFxx;                   % UNITS               Constant
IG22                            =  IFyy;                   % UNITS               Constant
IG33                            =  IFxx;                   % UNITS               Constant
L1                              =  w+c-rR*tan(lambda);     % UNITS               Constant
L2                              =  (rR-rF)/cos(lambda)+...
                                   sin(lambda)*(rF*tan...
                                   (lambda)-c);            % UNITS               Constant
L3                              =  rF*sin(lambda)-c*cos...
                                   (lambda);               % UNITS               Constant
L4                              =  xB;                     % UNITS               Constant
L5                              =  zB+rR;                  % UNITS               Constant
L6                              =  (xH-w)*cos(lambda)-...
                                   (zH+rF)*sin(lambda);    % UNITS               Constant
L7                              =  (xH-w)*sin(lambda)+...
                                   (zH+rF)*cos(lambda);    % UNITS               Constant
MC                              =  mR;                     % UNITS               Constant
MD                              =  mB;                     % UNITS               Constant
MF                              =  mH;                     % UNITS               Constant
MG                              =  mF;                     % UNITS               Constant
RF                              =  rF;                     % UNITS               Constant
RR                              =  rR;                     % UNITS               Constant
THETA                           =  lambda;                 % UNITS               Constant

Q1                              =  0.0;                    % UNITS               Initial Value
Q2                              =  0.0;                    % UNITS               Initial Value
Q3                              =  0.0;                    % UNITS               Initial Value
Q4                              =  0.0;                    % UNITS               Initial Value
Q5                              =  0.0;                    % UNITS               Initial Value
Q7                              =  0.0;                    % UNITS               Initial Value
Q8                              =  0.0;                    % UNITS               Initial Value
U4                              =  0.5;                    % UNITS               Initial Value
U5                              =  -v/rR;                    % UNITS               Initial Value
U7                              =  0.0;                    % UNITS               Initial Value

TINITIAL                        =  0.0;                    % sec                 Initial Time
TFINAL                          =  5.0;                    % sec                 Final Time
INTEGSTP                        =  0.01;                  % sec                 Integration Step
PRINTINT                        =  1;                      % Positive Integer    Print-Integer
ABSERR                          =  1.0E-08;                %                     Absolute Error
RELERR                          =  1.0E-07 ;               %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------