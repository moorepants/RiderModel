function [Q6,U6] = pitch(Q4,Q7,U4,U7,D1,D2,D3,LAMBDA,RF,RR)
D1=D1;
D2=D2;
D3=D3;
LAMBDA=LAMBDA;
RF=RF;
RR=RR;
Q4=Q4;
Q7=Q7;
U4=U4;
U7=U7;

Q6 = fzero(@(x) pzero(x,Q4,Q7,D1,D2,D3,LAMBDA,RF,RR),0);

U6 = ((sin(Q4)*cos(Q7)+sin(Q7)*cos(Q4)*sin(LAMBDA+Q6))*(D3+RF*(sin(Q4)*...
    sin(Q7)-cos(Q4)*cos(Q7)*sin(LAMBDA+Q6))/(1-(sin(Q4)*cos(Q7)+sin(Q7)*...
     cos(Q4)*sin(LAMBDA+Q6))^2)^0.5)*U7+(RR*sin(Q4)+D1*sin(Q4)*sin(...
     LAMBDA+Q6)+D3*(sin(Q7)*cos(Q4)+sin(Q4)*cos(Q7)*sin(LAMBDA+Q6))-D2*...
     sin(Q4)*cos(LAMBDA+Q6)-RF*(sin(Q4)*cos(Q7)+sin(Q7)*cos(Q4)*sin(...
     LAMBDA+Q6))*(cos(Q4)*cos(Q7)-sin(Q4)*sin(Q7)*sin(LAMBDA+Q6))/(1-(...
     sin(Q4)*cos(Q7)+sin(Q7)*cos(Q4)*sin(LAMBDA+Q6))^2)^0.5)*U4)/(cos(...
     Q4)*(D1*cos(LAMBDA+Q6)+D2*sin(LAMBDA+Q6)+D3*cos(Q7)*cos(LAMBDA+Q6)+...
     RF*sin(Q7)*cos(LAMBDA+Q6)*(sin(Q4)*cos(Q7)+sin(Q7)*cos(Q4)*sin(...
     LAMBDA+Q6))/(1-(sin(Q4)*cos(Q7)+sin(Q7)*cos(Q4)*sin(LAMBDA+Q6))^2)...
     ^0.5));
