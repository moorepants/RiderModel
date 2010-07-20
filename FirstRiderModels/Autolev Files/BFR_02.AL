%---------------------------------------------------------------------%
%    FILE: BFR_02.AL
%    DATE: November 1, 2007
%  AUTHOR: JASON MOORE
% PROBLEM: GENERATES THE EQUATIONS OF MOTION FOR A FLEXIBLE RIDER 
%          NO-HANDS BICYCLE MADE UP OF SIX RIGID BODIES.
%   NOTES: IN REVISION 2, THE PITCH LINEARIZATION IS REMOVED AND THE
%          MODEL IS RECONSTRUCTED WITH AN ADDITIONAL DIFFERENTIAL
%          EQUATION FOR THE PITCH "DEGREE OF FREEDOM"
%---------------------------------------------------------------------%
%         DEFAULT SETTINGS
%---------------------------------------------------------------------%

AUTOZ ON
AUTORHS OFF
OVERWRITE ALL

%---------------------------------------------------------------------%
%         NEWTONIAN, BODIES, FRAMES, PARTICLES, POINTS
%---------------------------------------------------------------------%

% DECLARE THE INERTIAL REFERENCE FRAME

NEWTONIAN N

% DECLARE FOUR INTERMEDIATE FRAMES
% A: YAW FRAME
% B: ROLL FRAME
% E: HEAD TUBE ANGLE FRAME
% I: UPPER BODY LEAN FRAME

FRAMES A,B,E,I

% DECLARE SIX BODIES
% C: REAR WHEEL
% D: FRAME
% F: FORK
% G: FRONT WHEEL
% H: LOWER BODY
% J: UPPER BODY

BODIES C,D,F,G,H,J

% DECLARE SEVEN POINTS
% NC: REAR CONTACT POINT ON GROUND
% CN: REAR CONTACT POINT ON WHEEL
% NG: FRONT CONTACT POINT ON GROUND
% GN: FRONT CONTACT POINT ON WHEEL
% EF: STEER PIVOT POINT
% DH: LOWER BODY PIVOT POINT
% HI: UPPER BODY PIVOT POINT
% DS: CROTCH SPRING CONNECTION POINT

POINTS NC,CN,NG,GN,EF,DH,HI,DS

%---------------------------------------------------------------------%
%         CONSTANTS AND VARIABLES
%---------------------------------------------------------------------%

% RF: RADIUS OF FRONT WHEEL
% RR: RADIUS OF REAR WHEEL
% D1: THE PERPENDICULAR DISTANCE FROM THE HEAD TUBE AXIS TO THE CENTER
%     OF THE REAR WHEEL
% D2: THE DISTANCE ALONG THE HEAD TUBE AXIS
% D3: THE PERPENDICULAR DISTANCE FROM THE HEAD TUBE AXIS TO THE CENTER
%     OF THE FRONT WHEEL
% L1: THE DISTANCE IN THE D1> DIRECTION FROM THE CENTER OF THE REAR
%     WHEEL TO THE FRAME CENTER OF MASS
% L2: THE DISTANCE IN THE D3> DIRECTION FROM THE CENTER OF THE REAR
%     WHEEL TO THE FRAME CENTER OF MASS
% L3: THE DISTANCE IN THE D3> DIRECTION FROM THE CENTER OF THE REAR
%     WHEEL TO THE LOCATION OF THE LOWER BODY PIVOT POINT
% L4: THE DISTANCE IN THE H1> DIRECTION FROM THE LOWER BODY PIVOT POINT
%     TO THE CENTER OF MASS OF THE LOWER BODY
% L5: THE DISTANCE IN THE H3> DIRECTION FROM THE LOWER BODY PIVOT POINT
%     TO THE CENTER OF MASS OF THE LOWER BODY
% L6: THE DISTANCE IN THE H3> DIRECTION FROM THE LOWER BODY PIVOT POINT
%     TO THE UPPER BODY PIVOT POINT
% L7: THE DISTANCE IN THE J1> DIRECTION FROM THE UPPER BODY PIVOT POINT
%     TO THE CENTER OF MASS OF THE UPPER BODY
% L8: THE DISTANCE IN THE J3> DIRECTION FROM THE UPPER BODY PIVOT POINT
%     TO THE CENTER OF MASS OF THE UPPER BODY
% L9: THE DISTANCE IN THE F1> DIRECTION FROM THE STEER POINT TO THE
%     CENTER OF MASS OF THE FORK
% L10:THE DISTANCE IN THE F3> DIRECTION FROM THE STEER POINT TO THE
%     CENTER OF MASS OF THE FORK
% LAMBDA: COMPLEMENT TO THE HEAD TUBE ANGLE
% K1 : CROTCH STIFFNESS
% K2 : LEAN STIFFNESS
% K3 : TWIST STIFFNESS
% C1 : CROTCH DAMPING
% C2 : LEAN DAMPING
% C3 : TWIST DAMPING
% G:  ACCELERATION DUE TO GRAVITY

CONSTANTS RF,RR,D{3},L{10},LAMBDA,K{3},C{3},G

% DECLARE THE GENERALIZED COORDINATES
% Q1:  PERPENDICULAR DISTANCE FROM THE N2> AXIS TO THE REAR CONTACT
%      POINT IN THE GROUND PLANE
% Q2:  PERPENDICULAR DISTANCE FROM THE N1> AXIS TO THE REAR CONTACT
%      POINT IN THE GROUND PLANE
% Q3:  FRAME YAW ANGLE
% Q4:  FRAME ROLL ANGLE
% Q5:  REAR WHEEL ROTATION ANGLE
% Q6:  FRAME PITCH ANGLE
% Q7:  STEERING ROTATION ANGLE
% Q8:  FRONT WHEEL ROTATION ANGLE
% Q9:  LOWER BODY LEAN ANGLE
% Q10: UPPER BODY LEAN ANGLE
% Q11: UPPER BODY TWIST

VARIABLES Q{11}'

%---------------------------------------------------------------------%
%         GENERALIZED SPEEDS
%---------------------------------------------------------------------%

MOTIONVARIABLES' U{5}',U{6:11}'

%---------------------------------------------------------------------%
%         MASS AND INERTIA PROPERTIES
%---------------------------------------------------------------------%

MASS C=MC,D=MD,F=MF,G=MG,H=MH,J=MJ
INERTIA C,IC11,IC22,IC33
INERTIA D,ID11,ID22,ID33,ID12,ID23,ID31
INERTIA F,IF11,IF22,IF33,IF12,IF23,IF31
INERTIA G,IG11,IG22,IG33
INERTIA H,IH11,IH22,IH33,IH12,IH23,IH31
INERTIA J,IJ11,IJ22,IJ33,IJ12,IJ23,IJ31

%---------------------------------------------------------------------%
%         ANGULAR RELATIONSHIPS                                       %
%---------------------------------------------------------------------%

% FRAME YAW

SIMPROT(N,A,3,Q3)

% FRAME ROLL

SIMPROT(A,B,1,Q4)

% REAR WHEEL ROTATION

SIMPROT(B,C,2,Q5)

% FRAME PITCH

SIMPROT(B,D,2,Q6)

% HEAD TUBE ANGLE
SIMPROT(D,E,2,LAMBDA)

% STEERING ANGLE

SIMPROT(E,F,3,Q7)

% FRONT WHEEL ROTATION

SIMPROT(F,G,2,Q8)

% LOWER BODY

SIMPROT(D,H,1,Q9)

% UPPER BODY LEAN ANGLE

SIMPROT(H,I,1,Q10)

% UPPER BODY TWIST ANGLE

SIMPROT(I,J,3,Q11)

%---------------------------------------------------------------------%
%         POSITION VECTORS
%---------------------------------------------------------------------%

P_NO_NC>=Q1*N1>+Q2*N2>

P_NC_CN>=0>

P_CN_CO>=-RR*B3>

P_CO_EF>=D1*E1>

P_EF_GO>=D2*F3>+D3*F1>

P_GO_GN>=RF*UNITVEC(N3>-DOT(F2>,N3>)*F2>)

P_GN_NG>=0>

P_CO_DO>=L1*D1>+L2*D3>

P_CO_DH>=L3*D3>

P_DH_HO>=L4*H1>+L5*H3>

P_DH_HI>=L6*H3>

P_HI_JO>=L7*J1>+L8*J3>

P_EF_FO>=L9*F1>+L10*F3>

P_DH_DS>=DOT(P_DH_HI>,D3>)*D3>

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
Q10'=U10
Q11'=U11

%---------------------------------------------------------------------%
%         ANGULAR VELOCITIES
%---------------------------------------------------------------------%

W_A_N>=U3*A3>

W_B_A>=U4*B1>

W_C_B>=U5*C2>

W_D_B>=U6*D2>

W_E_D>=0*E2>

W_F_E>=U7*F3>

W_G_F>=U8*G2>

W_H_D>=U9*H1>

W_I_H>=U10*I1>

W_J_I>=U11*J3>

%---------------------------------------------------------------------%
%         VELOCITIES
%---------------------------------------------------------------------%

V_CO_N>=DT(P_NO_CO>,N)

V2PTS(N,C,CO,CN)

V_GO_N>=DT(P_NO_GO>,N)

V2PTS(N,G,GO,GN)

V2PTS(N,D,CO,DO)

V2PTS(N,F,GO,FO)

V2PTS(N,D,CO,DH)

V2PTS(N,H,DH,HO)

V2PTS(N,I,HO,JO)

V2PTS(N,D,CO,DS)

V2PTS(N,H,DH,HI)

%---------------------------------------------------------------------%
%         MOTION CONSTRAINTS
%---------------------------------------------------------------------%

% DUE TO THE ASSUMPTIONS OF NO SIDE SLIP AND NO SLIP ROLLING THE
% VELOCITIES OF THE FRONT AND REAR WHEEL CONTACT POINTS, CN AND GN,
% CANNOT HAVE COMPONENTS OF VELOCITY IN THE GROUND PLANE

DEPENDENT[1]=DOT(V_CN_N>,N1>)
DEPENDENT[2]=DOT(V_CN_N>,N2>)
DEPENDENT[3]=DOT(V_GN_N>,N1>)
DEPENDENT[4]=DOT(V_GN_N>,N2>)

% U4,U5,U6,U7,U9,U10,U11 ARE TAKEN TO BE THE INDEPENDENT GENERALIZED
% SPEEDS

CONSTRAIN(DEPENDENT[U1,U2,U3,U8])

%---------------------------------------------------------------------%
%         ANGULAR ACCELERATIONS
%---------------------------------------------------------------------%

ALF_A_N>=DT(W_A_N>,N)

ALF_B_A>=DT(W_B_A>,N)

ALF_C_B>=DT(W_C_B>,N)

ALF_D_B>=DT(W_D_B>,N)

ALF_E_D>=DT(W_E_D>,N)

ALF_F_E>=DT(W_F_E>,N)

ALF_G_F>=DT(W_G_F>,N)

ALF_H_D>=DT(W_H_D>,N)

ALF_I_H>=DT(W_I_H>,N)

ALF_J_I>=DT(W_J_I>,N)

%---------------------------------------------------------------------%
%         ACCELERATIONS
%---------------------------------------------------------------------%

A_CO_N>=DT(V_CO_N>,N)

A_GO_N>=DT(V_GO_N>,N)

A2PTS(N,D,CO,DO)

A2PTS(N,F,GO,FO)

A2PTS(N,D,CO,DH)

A2PTS(N,H,DH,HO)

A2PTS(N,I,HO,JO)

A2PTS(N,D,CO,DS)

A2PTS(N,H,DH,HI)

%---------------------------------------------------------------------%
%         FORCES
%---------------------------------------------------------------------%

GRAVITY(G*N3>,C,D,F,G,H,J)

FORCE(DS/HI,(-K1*DOT(P_DH_HI>,D2>)-C1*DT(DOT(P_DH_HI>,D2>)))*D2>)

TORQUE(H/I,(-K2*Q10-C2*U10)*H1>)

TORQUE(I/J,(-K3*Q11-C3*U11)*I3>)

%---------------------------------------------------------------------%
%         EQUATIONS OF MOTION
%---------------------------------------------------------------------%

ZERO=FR()+FRSTAR()
KANE()

%---------------------------------------------------------------------%
%         SIMULATION
%---------------------------------------------------------------------%

UnitSystem  kg,meter,sec

OUTPUT T S, Q1 M,  Q2 M,  Q3 RAD,  Q4 RAD,  Q5 RAD,  Q7 RAD,  Q8 RAD,&
Q9 RAD, Q10 RAD, Q11 RAD

OUTPUT T S, U1 M/S,U2 M/S,U3 RAD/S,U4 RAD/S,U5 RAD/S,U7 RAD/S,&
U8 RAD/S, U9 RAD/S, U10 RAD/S, U11 RAD/S

CODE DYNAMICS() BFR_02.M

%---------------------------------------------------------------------%
%         SAVE OUTPUT
%---------------------------------------------------------------------%

SAVE BFR_02.ALL
