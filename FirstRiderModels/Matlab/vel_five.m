close all
clear all
%whipple
eval_w = -0.7634 + 4.3251i;
evec_w = [  0.0016 + 0.1388i       % roll
            0.0301 + 0.1706i ];    % steer
%extended
eval_e = 1.1135 +3.4323i;
evec_e = [  0.0103 - 0.1272i       % roll
           -0.0323 - 0.1022i       % steer
            0.0329 - 0.0185i       % lower body lean 
            0.1206 + 0.0579i       % upper body lean
           -0.0480 + 0.1480i ];    % upper body twist
evec_e = evec_e*-1;
%plot
subplot(121)
hold on
plot([0;real(evec_w(1))],[0;imag(evec_w(1))],'r',...
     [0;real(evec_w(2))],[0;imag(evec_w(2))],'b')
legend('Roll','Steer')
plot(real(evec_w(1)),imag(evec_w(1)),'or',...
     real(evec_w(2)),imag(evec_w(2)),'ob')
axis([-0.2 0.2 -0.2 0.2])
axis square
box on
title({'Whipple Model Stable Weave Mode';'\lambda = (-0.76, 4.33 j), v = 5 m/s'})
xlabel('Imaginary parts of the eigenvector components [1/s]')
ylabel('Real parts of the eigenvector components [1/s]')
hold off
subplot(122)
hold on
plot([0;real(evec_e(1))],[0;imag(evec_e(1))],'r',...
     [0;real(evec_e(2))],[0;imag(evec_e(2))],'b',...
     [0;real(evec_e(3))],[0;imag(evec_e(3))],'g',...
     [0;real(evec_e(4))],[0;imag(evec_e(4))],'c',...
     [0;real(evec_e(5))],[0;imag(evec_e(5))],'m')
legend('Roll',...
       'Steer',...
       'Lower Body Lean',...
       'Upper Body Lean',...
       'Upper Body Twist')
plot(real(evec_e(1)),imag(evec_e(1)),'or',...
     real(evec_e(2)),imag(evec_e(2)),'ob',...
     real(evec_e(3)),imag(evec_e(3)),'og',...
     real(evec_e(4)),imag(evec_e(4)),'oc',...
     real(evec_e(5)),imag(evec_e(5)),'om')
axis([-0.2 0.2 -0.2 0.2])
axis square
box on
title({'Extended Model Unstable Weave Mode';'\lambda = (1.11, 3.43 j), v = 5 m/s'})
xlabel('Imaginary parts of the eigenvector components [1/s]')
ylabel('Real parts of the eigenvector components [1/s]')
hold off