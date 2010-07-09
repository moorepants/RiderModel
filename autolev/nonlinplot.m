clear all
close all
clc

% filename = input('What is the filename? (no extensions)','s')
% 
% time       = load([filename '.1']);
% qs         = load([filename '.2']);
% us         = load([filename '.3']);
% energy     = load([filename '.4']);
% con        = load([filename '.5']);

time       = load('WhipplePlusRiderLeanDyn.1');
qs         = load('WhipplePlusRiderLeanDyn.2');
us         = load('WhipplePlusRiderLeanDyn.3');
energy   =0 % = load('BRR_05.4');
con      =0%  = load('BRR_05.5');

%plot the q's
figure(1)
subplot(2,2,1)
%plot(qs(:,1),qs(:,2),'-',con(:,2),con(:,3),':')
title('Path of the wheel contacts')
xlabel('Distance [m]')
ylabel('Distance [m]')
legend('rear wheel','front wheel')
subplot(2,2,2)
plot(time,qs(:,1),'-',time,qs(:,2),':')
xlabel('Time [s]')
ylabel('Distance [m]')
legend('q1','q2')
subplot(2,2,3)
plot(time,qs(:,3),'-',time,qs(:,4),':',time,qs(:,6),'--',time,qs(:,7),'-.')
xlabel('Time [s]')
ylabel('Angle [rad]')
legend('q3 [yaw]','q4 [roll]','q6 [pitch]','q7 [steer]')
subplot(2,2,4)
plot(time,qs(:,5),'-',time,qs(:,8),':')
xlabel('Time [s]')
ylabel('Angle [rad]')
legend('q5 [rear wheel rotation]','q8 [front wheel rotation]')

%plot the u's
figure(2)
subplot(3,1,1)
plot(time,us(:,1),'-',time,us(:,2),':')
xlabel('Time [s]')
ylabel('Speed [m/s]')
legend('u1','u2')
subplot(3,1,2)
plot(time,us(:,3),'-',time,us(:,4),':',time,us(:,6),'--',time,us(:,7),'-.')
xlabel('Time [s]')
ylabel('Rotation Rate [rad]')
legend('u3 [yaw]','u4 [roll]','u6 [pitch]','u7 [steer]')
subplot(3,1,3)
plot(time,us(:,5),'-',time,us(:,8),':')
xlabel('Time [s]')
ylabel('Rotation Rate [rad/s]')
legend('u5 [rear wheel]','u8 [front wheel]')

%plot the benchmark graph
figure(3)
rollrate = us(:,4);
steerrate = us(:,7);
velocity = -us(:,5).*0.3;
[AX,H1,H2] = plotyy(time,rollrate,time,velocity);
set(get(AX(1),'Ylabel'),'String','Angular Rate [rad/s]')
set(AX(1),'YLim',[-0.5 1])
set(AX(1),'YTick',[-0.5;0;0.5;1])
set(AX(1),'XTick',[0;1;2;3;4;5])
set(get(AX(2),'Ylabel'),'String','Forward Velocity [m/s]')
set(AX(2),'YLim',[4.55 4.70])
set(AX(2),'YTick',[4.55;4.60;4.65;4.70])
set(AX(2),'XTick',[0;1;2;3;4;5])
hold on
plot(time,steerrate,'r')
xlabel('Time [s]')
grid on

%plot the system energy
figure(4)
subplot(3,1,1)
plot(time,energy(:,1))
xlabel('Time [s]')
ylabel('Kinetic Energy [J]')
subplot(3,1,2)
plot(time,energy(:,2))
xlabel('Time [s]')
ylabel('Potential Energy [J]')
subplot(3,1,3)
tot_en = energy(:,1)+energy(:,2);
plot(time,(energy(:,1)+energy(:,2)))
xlabel('Time [s]')
ylabel('Kinetic + Potential [J]')

%plot
figure(5)
plot(time,con(:,1))
xlabel('Time [s]')
ylabel('N3 component of the contact point vector')