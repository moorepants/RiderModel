%---------------------------------------------------------------------%
%    FILE: LeanTorque.AL
%    DATE: JULY 6, 2010
%  AUTHOR: JASON MOORE
% PROBLEM: GENERATES THE LINEARIZED EQUATIONS OF MOTION FOR A LEANING
%          RIDER NO-HANDS BICYCLE MADE UP OF FIVE RIGID BODIES (I.E.
%          THE WHIPPLE MODEL PLUS RIDER LEAN).
%   NOTES: OUTPUTS THE A AND B MATRICES AND THE CANONICAL FORM
%---------------------------------------------------------------------%
%         DEFAULT SETTINGS
%---------------------------------------------------------------------%

AUTOZ ON
AUTORHS OFF
OVERWRITE ALL
BEEPSOUND OFF

%---------------------------------------------------------------------%
%         NEWTONIAN, BODIES, FRAMES, PARTICLES, POINTS
%---------------------------------------------------------------------%

% DECLARE THE INERTIAL REFERENCE FRAME

NEWTONIAN N

% DECLARE TWO INTERMEDIATE FRAMES
% A: YAW FRAME
% B: ROLL FRAME
% H: RIDER HINGE FRAME

FRAMES A,B,H

% DECLARE SIX BODIES
% C: BICYCLE FRAME
% D: REAR WHEEL
% E: FORK/HANDLEBAR
% F: FRONT WHEEL
% G: RIDER UPPER BODY

BODIES C,D,E,F,G

% DECLARE FOUR POINTS
% ND: REAR CONTACT POINT ON GROUND
% DN: REAR CONTACT POINT ON WHEEL
% NF: FRONT CONTACT POINT ON GROUND
% FN: FRONT CONTACT POINT ON WHEEL

POINTS ND,DN,NF,FN

%---------------------------------------------------------------------%
%         CONSTANTS AND VARIABLES
%---------------------------------------------------------------------%
% DEFINE THE BENCHMARK PARAMETERS FROM MEIJAARD ET. AL, 2007

% W:        WHEELBASE                          [M]
% C:        TRAIL                              [M]
% LAMBDA:   STEER AXIS TILT                    [RAD]
% G:        GRAVITY                            [N/KG]
% V:        FORWARD SPEED                      [M/S]
% RR:       REAR WHEEL RADIUS                  [M]
% M_R:      REAR WHEEL MASS                    [KG]
% IRXX:     REAR WHEEL MASS MOMENT OF INERTIA  [KG*M^2]
% IRYY:     REAR WHEEL MASS MOMENT OF INERTIA  [KG*M^2]
% XB:       REAR BODY CENTER OF MASS LOCATION  [M]
% ZB:       REAR BODY CENTER OF MASS LOCATION  [M]
% M_B:      REAR BODY MASS                     [KG]
% IBXX:     REAR BODY MASS MOMENT OF INERTIA   [KG*M^2]
% IBYY:     REAR BODY MASS MOMENT OF INERTIA   [KG*M^2]
% IBZZ:     REAR BODY MASS MOMENT OF INERTIA   [KG*M^2]
% IBXZ:     REAR BODY MASS MOMENT OF INERTIA   [KG*M^2]
% XH:       FORK CENTER OF MASS LOCATION       [M]
% ZH:       FORK CENTER OF MASS LOCATION       [M]
% M_H:      FORK MASS                          [KG]
% IHXX:     FORK MASS MOMENT OF INERTIA        [KG*M^2]
% IHYY:     FORK MASS MOMENT OF INERTIA        [KG*M^2]
% IHZZ:     FORK MASS MOMENT OF INERTIA        [KG*M^2]
% IHXZ:     FORK MASS MOMENT OF INERTIA        [KG*M^2]
% RF:       FRONT WHEEL RADIUS                 [M]
% M_F:      FRONT WHEEL MASS                   [KG]
% IFXX:     FRONT WHEEL MASS MOMENT OF INERTIA [KG*M^2]
% IFYY:     FRONT WHEEL MASS MOMENT OF INERTIA [KG*M^2]
% TPHI:     LEAN TORQUE                        [N*M]
% T_DELTA:  STEER TORQUE                       [N*M]
% T_THETAR: REAR WHEEL TORQUE                  [N*M]

CONSTANTS W,C,LAMBDA,G,V
CONSTANTS RR,M_R,IRXX,IRYY
CONSTANTS XB,ZB,M_B,IBXX,IBYY,IBZZ,IBXZ
CONSTANTS XH,ZH,M_H,IHXX,IHYY,IHZZ,IHXZ
CONSTANTS RF,M_F,IFXX,IFYY
SPECIFIED T_PHI,T_DELTA,T_THETAR

% ADD SOME NEW PARAMETERS FOR THE UPPER BODY LEAN
% T_PHI_R: RIDER LEAN TORQUE [N*M]
% HH: HINGE HEIGHT [M]
% IBRXX
% IBRYY
% IBRZZ
% IBRXZ
% M_BR
% XBR
% ZBR

SPECIFIED T_PHI_R
CONSTANTS HH
CONSTANTS XBR,ZBR,M_BR,IBRXX,IBRYY,IBRZZ,IBRXZ

% CONVERT THE BENCHMARK CONSTANTS TO MOORE'S CONSTANTS
% RF: RADIUS OF FRONT WHEEL
% RR: RADIUS OF REAR WHEEL
% D1: THE PERPENDICULAR DISTANCE FROM THE HEAD TUBE AXIS TO THE CENTER
%     OF THE REAR WHEEL
% D3: THE PERPENDICULAR DISTANCE FROM THE HEAD TUBE AXIS TO THE CENTER
%     OF THE FRONT WHEEL (FORK OFFSET)
% L1: THE DISTANCE IN THE D1> DIRECTION FROM THE CENTER OF THE REAR
%     WHEEL TO THE FRAME CENTER OF MASS
% L2: THE DISTANCE IN THE D3> DIRECTION FROM THE CENTER OF THE REAR
%     WHEEL TO THE FRAME CENTER OF MASS
% L3: THE DISTANCE IN THE F1> DIRECTION FROM THE STEER POINT TO THE
%     CENTER OF MASS OF THE FORK
% L4: THE DISTANCE IN THE F3> DIRECTION FROM THE STEER POINT TO THE
%     CENTER OF MASS OF THE FORK
% L5: UPPER BODY CG LOCATION
% L6: UPPER BODY CG LOCATION

D1    =  COS(LAMBDA)*(C+W-RR*TAN(LAMBDA))
D3    = -COS(LAMBDA)*(C-RF*TAN(LAMBDA))
D2    = (RR+D1*SIN(LAMBDA)-RF+D3*SIN(LAMBDA))/COS(LAMBDA)
% REAR WHEEL INERTIA
ID11  =  IRXX
ID22  =  IRYY
ID33  =  IRXX
CF    =  [COS(LAMBDA),0,-SIN(LAMBDA);0,1,0;SIN(LAMBDA),0,COS(LAMBDA)]
% ROTATE BICYCLE FRAME INERTIA THROUGH LAMBDA
IB    =  [IBXX,0,IBXZ;0,IBYY,0;IBXZ,0,IBZZ]
IBROT =  CF*IB*TRANSPOSE(CF)
% BICYCLE FRAME INERTIA
IC11  =  IBROT[1,1]
IC12  =  IBROT[1,2]
IC22  =  IBROT[2,2]
IC23  =  IBROT[2,3]
IC31  =  IBROT[3,1]
IC33  =  IBROT[3,3]
% ROTATE FORK INERTIA MATRIX THROUGH LAMBDA
IH    =  [IHXX,0,IHXZ;0,IHYY,0;IHXZ,0,IHZZ]
IHROT =  CF*IH*TRANSPOSE(CF)
% FORK/HANDLEBAR INERTIA
IE11  =  IHROT[1,1]
IE12  =  IHROT[1,2]
IE22  =  IHROT[2,2]
IE23  =  IHROT[2,3]
IE31  =  IHROT[3,1]
IE33  =  IHROT[3,3]
% FRONT WHEEL INERTIA
IF11  =  IFXX
IF22  =  IFYY
IF33  =  IFXX
% UPPER BODY INERTIA
IG11  =  IBRXX
IG12  =  0.0
IG22  =  IBRYY
IG23  =  0.0
IG31  =  IBRXZ
IG33  =  IBRZZ
% MASS CENTER LOCATIONS
L1    =  XB*COS(LAMBDA)-ZB*SIN(LAMBDA)-RR*SIN(LAMBDA)
L2    =  XB*SIN(LAMBDA)+ZB*COS(LAMBDA)+RR*COS(LAMBDA)
L3    =  COS(LAMBDA)*XH-SIN(LAMBDA)*ZH-C*COS(LAMBDA)-W*COS(LAMBDA)
L4    =  RR*COS(LAMBDA)+XH*SIN(LAMBDA)+ZH*COS(LAMBDA)
L5    =  XBR
L6    =  ZBR-HH
% MASSES
MC    =  M_B
MD    =  M_R
ME    =  M_H
MF    =  M_F
MG    =  M_BR
% INPUT TORQUES
T4    =  T_PHI
T6    =  T_THETAR
T7    =  T_DELTA
T9    =  T_PHI_R

% DECLARE THE GENERALIZED COORDINATES
% Q1:  PERPENDICULAR DISTANCE FROM THE N2> AXIS TO THE REAR CONTACT
%      POINT IN THE GROUND PLANE
% Q2:  PERPENDICULAR DISTANCE FROM THE N1> AXIS TO THE REAR CONTACT
%      POINT IN THE GROUND PLANE
% Q3:  FRAME YAW ANGLE
% Q4:  FRAME ROLL ANGLE
% Q5:  FRAME PITCH ANGLE
% Q6:  REAR WHEEL ROTATION ANGLE
% Q7:  STEERING ROTATION ANGLE
% Q8:  FRONT WHEEL ROTATION ANGLE
% Q9:  RIDER LEAN ANGLE

VARIABLES Q{9}'

%---------------------------------------------------------------------%
%         GENERALIZED SPEEDS
%---------------------------------------------------------------------%

MOTIONVARIABLES' U{9}'

%---------------------------------------------------------------------%
%         MASS AND INERTIA PROPERTIES
%---------------------------------------------------------------------%

MASS C=MC,D=MD,E=ME,F=MF,G=MG
INERTIA C,IC11,IC22,IC33,IC12,IC23,IC31
INERTIA D,ID11,ID22,ID33
INERTIA E,IE11,IE22,IE33,IE12,IE23,IE31
INERTIA F,IF11,IF22,IF33
INERTIA G,IG11,IG22,IG33,IG12,IG23,IG31

%---------------------------------------------------------------------%
%         ANGULAR RELATIONSHIPS                                       %
%---------------------------------------------------------------------%

% FRAME YAW

SIMPROT(N,A,3,Q3)

% FRAME ROLL

SIMPROT(A,B,1,Q4)

% FRAME PITCH

SIMPROT(B,C,2,Q5+LAMBDA)

% REAR WHEEL ROTATION

SIMPROT(C,D,2,Q6)

% STEERING ANGLE

SIMPROT(C,E,3,Q7)

% FRONT WHEEL ROTATION

SIMPROT(E,F,2,Q8)

% HINGE FRAME

SIMPROT(C,H,2,-LAMBDA)

% RIDER LEAN

SIMPROT(H,G,1,Q9)

%---------------------------------------------------------------------%
%         POSITION VECTORS
%---------------------------------------------------------------------%

% LOCATE THE CENTER OF MASS FOR EACH BODY

P_NO_DO>=Q1*N1>+Q2*N2>-RR*B3> % NEWTONIAN ORIGIN TO REAR WHEEL CENTER

P_DO_CO>=L1*C1>+L2*C3> % REAR WHEEL CENTER TO BICYCLE FRAME CENTER

P_DO_EO>=D1*C1>+L3*E1>+L4*E3> % REAR WHEEL CENTER TO FORK/HANDLEBAR CENTER

% REAR WHEEL CENTER TO THE FRONT WHEEL CENTER

P_DO_FO>=D1*C1>+D2*E3>+D3*E1>

% REAR WHEEL CENTER TO THE UPPER BODY CENTER

P_DO_GO>=(RR+HH)*H3>+L5*G1>+L6*G3>

% LOCATE THE GROUND CONTACT POINTS

P_DO_DN>=RR*B3>

P_DN_ND>=0>

P_FO_FN>=RF*UNITVEC(N3>-DOT(E2>,N3>)*E2>)

P_FN_NF>=0>

%---------------------------------------------------------------------%
%         DEFINE THE GENERALIZED SPEEDS
%---------------------------------------------------------------------%

Q1'=U1
Q2'=U2
Q3'=U3
Q4'=U4
Q5'=U5
Q6'=U6
Q7'=U7
Q8'=U8
Q9'=U9

%---------------------------------------------------------------------%
%         ANGULAR VELOCITIES
%---------------------------------------------------------------------%

ANGVEL(N,A)
ANGVEL(N,B)
ANGVEL(N,C)
ANGVEL(N,D)
ANGVEL(N,E)
ANGVEL(N,F)
ANGVEL(N,G)
ANGVEL(N,H)

%---------------------------------------------------------------------%
%         VELOCITIES
%---------------------------------------------------------------------%

V_CO_N>=DT(P_NO_CO>,N)
V_DO_N>=DT(P_NO_DO>,N)
V_EO_N>=DT(P_NO_EO>,N)
V_FO_N>=DT(P_NO_FO>,N)
V_GO_N>=DT(P_NO_GO>,N)

V2PTS(N,D,DO,DN)
V2PTS(N,F,FO,FN)

%---------------------------------------------------------------------%
%         DEFINE THE PITCH CONFIGURATION CONSTRAINT
%---------------------------------------------------------------------%

% SET THE N3> COMPONENT OF P_ND_NF> EQUAL TO ZERO
PZERO=DOT(P_ND_NF>,N3>)

%---------------------------------------------------------------------%
%         MOTION CONSTRAINTS
%---------------------------------------------------------------------%

% DUE TO THE ASSUMPTIONS OF NO SIDE SLIP AND NO SLIP ROLLING THE
% VELOCITIES OF THE FRONT AND REAR WHEEL CONTACT POINTS, CN AND GN,
% CANNOT HAVE COMPONENTS OF VELOCITY IN THE GROUND PLANE

DEPENDENT[1]=DOT(V_DN_N>,A1>)
DEPENDENT[2]=DOT(V_DN_N>,A2>)
DEPENDENT[3]=DOT(V_FN_N>,A1>)
DEPENDENT[4]=DOT(V_FN_N>,A2>)
DEPENDENT[5]=DT(PZERO)

% THE REAR WHEEL ANGULAR SPEED, U6, THE ROLL RATE, U4,THE
% STEERING RATE, U7, AND THE RIDER LEAN RATE, U9 ARE TAKEN TO BE THE INDEPENDENT
% GENERALIZED SPEEDS

CONSTRAIN(DEPENDENT[U1,U2,U3,U5,U8])

%---------------------------------------------------------------------%
%         ANGULAR ACCELERATIONS
%---------------------------------------------------------------------%

ALF_C_N>=DT(W_C_N>,N)

ALF_D_N>=DT(W_D_N>,N)

ALF_E_N>=DT(W_E_N>,N)

ALF_F_N>=DT(W_F_N>,N)

ALF_G_N>=DT(W_G_N>,N)

%---------------------------------------------------------------------%
%         ACCELERATIONS
%---------------------------------------------------------------------%

A_CO_N>=DT(V_CO_N>,N)
A_DO_N>=DT(V_DO_N>,N)
A_EO_N>=DT(V_EO_N>,N)
A_FO_N>=DT(V_FO_N>,N)
A_GO_N>=DT(V_GO_N>,N)

%---------------------------------------------------------------------%
%         FORCES AND TORQUES
%---------------------------------------------------------------------%

GRAVITY(G*N3>,C,D,E,F,G)
TORQUE(A/B,T4*A1>) % ROLL TORQUE
TORQUE(C/D,T6*C2>) % REAR WHEEL TORQUE
TORQUE(C/E,T7*C3>) % STEER TORQUE
TORQUE(C/G,T9*G1>) % RIDER LEAN TORQUE

%---------------------------------------------------------------------%
%         EQUATIONS OF MOTION
%---------------------------------------------------------------------%

ZERO=FR()+FRSTAR()
SOLVE(ZERO,U4',U6',U7',U9')

%---------------------------------------------------------------------%
%         LINEARIZATION
%---------------------------------------------------------------------%
% LINEARIZES THE EQUATIONS OF MOTION ABOUT THE UPRIGHT CONFIGURATION
% AND CONSTANT FORWARD SPEED

A[1,1]=D(Q4',Q4)
A[1,2]=D(Q4',Q7)
A[1,3]=D(Q4',Q9)
A[1,4]=D(Q4',U4)
A[1,5]=D(Q4',U7)
A[1,6]=D(Q4',U9)

A[2,1]=D(Q7',Q4)
A[2,2]=D(Q7',Q7)
A[2,3]=D(Q7',Q9)
A[2,4]=D(Q7',U4)
A[2,5]=D(Q7',U7)
A[2,6]=D(Q7',U9)

A[3,1]=D(Q9',Q4)
A[3,2]=D(Q9',Q7)
A[3,3]=D(Q9',Q9)
A[3,4]=D(Q9',U4)
A[3,5]=D(Q9',U7)
A[3,6]=D(Q9',U9)

A[4,1]=D(U4',Q4)
A[4,2]=D(U4',Q7)
A[4,3]=D(U4',Q9)
A[4,4]=D(U4',U4)
A[4,5]=D(U4',U7)
A[4,6]=D(U4',U9)

A[5,1]=D(U7',Q4)
A[5,2]=D(U7',Q7)
A[5,3]=D(U7',Q9)
A[5,4]=D(U7',U4)
A[5,5]=D(U7',U7)
A[5,6]=D(U7',U9)

A[6,1]=D(U9',Q4)
A[6,2]=D(U9',Q7)
A[6,3]=D(U9',Q9)
A[6,4]=D(U9',U4)
A[6,5]=D(U9',U7)
A[6,6]=D(U9',U9)

B[1,1]=D(Q4',T4)
B[1,2]=D(Q4',T7)
B[1,3]=D(Q4',T9)

B[2,1]=D(Q7',T4)
B[2,2]=D(Q7',T7)
B[2,3]=D(Q7',T9)

B[3,1]=D(Q9',T4)
B[3,2]=D(Q9',T7)
B[3,3]=D(Q9',T9)

B[4,1]=D(U4',T4)
B[4,2]=D(U4',T7)
B[4,3]=D(U4',T9)

B[5,1]=D(U7',T4)
B[5,2]=D(U7',T7)
B[5,3]=D(U7',T9)

B[6,1]=D(U9',T4)
B[6,2]=D(U9',T7)
B[6,3]=D(U9',T9)

% PUT THE SYSTEM IN CANONICAL FORM: M*Q''+V*C1*Q'+(G*K0+V^2*K2)*Q=T

M=INV(ROWS(B,4:6))
CMAT=-M*COLS(ROWS(A,4:6),4:6)
K=-M*COLS(ROWS(A,4:6),1:3)

%C1=CMAT/V
%
%VSQ=V^2
%
%K0[1,1]=COEF(K[1,1],G)
%K0[1,2]=COEF(K[1,2],G)
%K0[1,3]=COEF(K[1,3],G)
%K0[2,1]=COEF(K[2,1],G)
%K0[2,2]=COEF(K[2,2],G)
%K0[2,3]=COEF(K[2,3],G)
%K0[3,1]=COEF(K[3,1],G)
%K0[3,2]=COEF(K[3,2],G)
%K0[3,3]=COEF(K[3,3],G)
%
%K2[1,1]=COEF(REPLACE(EXPLICIT(K[1,1]),V^2=VSQ),VSQ)
%K2[1,2]=COEF(REPLACE(EXPLICIT(K[1,2]),V^2=VSQ),VSQ)
%K2[1,3]=COEF(REPLACE(EXPLICIT(K[1,3]),V^2=VSQ),VSQ)
%K2[2,1]=COEF(REPLACE(EXPLICIT(K[2,1]),V^2=VSQ),VSQ)
%K2[2,2]=COEF(REPLACE(EXPLICIT(K[2,2]),V^2=VSQ),VSQ)
%K2[2,3]=COEF(REPLACE(EXPLICIT(K[2,3]),V^2=VSQ),VSQ)
%K2[3,1]=COEF(REPLACE(EXPLICIT(K[3,1]),V^2=VSQ),VSQ)
%K2[3,2]=COEF(REPLACE(EXPLICIT(K[3,2]),V^2=VSQ),VSQ)
%K2[3,3]=COEF(REPLACE(EXPLICIT(K[3,3]),V^2=VSQ),VSQ)
%
%ZEE(K0)
%ZEE(K2)

UNITSYSTEM  KG,M,S

%INPUT W=1.02 M,C=0.08 M,LAMBDA=3.141592653589793/10 RAD,G=9.81 KG,&
      %V=4.6 M/S
%
%INPUT RR=0.3 M,M_R=2 KG,IRXX=0.0603 KG*M^2,IRYY=0.12 KG*M^2
%
%INPUT XB=0.3 M,ZB=-0.9 M,M_B=85 KG
%
%INPUT IBXX=9.2 KG*M^2,IBYY=11 KG*M^2,IBZZ=2.8 KG*M^2,IBXZ=2.4 KG*M^2
%
%INPUT XH=0.9 M,ZH=-0.7 M,M_H=4 KG
%
%INPUT IHXX=0.05892 KG*M^2,IHYY=0.06 KG*M^2,IHZZ=0.00708 KG*M^2
%
%INPUT IHXZ=-0.00756 KG*M^2
%
%INPUT RF=0.35 M,M_F=3 KG,IFXX=0.1405 KG*M^2,IFYY=0.28 KG*M^2

%OUTPUT A,B,M,C1,K0,K2

ENCODE A,B,M,CMAT,K

CODE ALGEBRAIC() LeanTorque.m

%---------------------------------------------------------------------%
%         SAVE OUTPUT
%---------------------------------------------------------------------%

SAVE LeanTorque.ALL

%---------------------------------------------------------------------%
