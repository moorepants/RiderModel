function f = pzero(x,Q4,Q7,D1,D2,D3,LAMBDA,RF,RR)

f = D2*cos(Q4)*cos(LAMBDA+x) + D3*(sin(Q4)*sin(Q7)-cos(Q4)*cos(Q7)*...
        sin(LAMBDA+x)) + RF*(1-(sin(Q4)*cos(Q7)+sin(Q7)*cos(Q4)*sin(...
        LAMBDA+x))^2)^0.5 - RR*cos(Q4) - D1*cos(Q4)*sin(LAMBDA+x);