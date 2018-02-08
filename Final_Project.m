%% Project Name- BallBot by Daksh Dhingra
clc;
close all;
clear all;
%% Define Parameters
J0= 4.5e-8;
J1= 6.7e-3;
J2= 0.9375;
R0=0.025;
R1= 0.124;
M1= 0.65;
M2= 30;
L= 0.5;
be=1.85e-3;
g=9.81;
h1 = J1 +J0*(R1^2/R0^2) + (M1+M2)*R1^2;
h2 = M2*L*R1;
h3 = be*(R1^2/R0^2);
h4 = R1/R0;
h5 = J2+ M2* L^2;
h6 = -M2*g*L;
%% 2.1 Generation of continuous non linear model.
% The system vaiables to define the function.
syms('theta1','theta2','theta1_d','theta2_d','theta1_dd','theta2_dd','tau')

% The value of FR0, FR1, FH, Fv are computed mathematically and plugged in
% the equations below
fh = M2*R1*theta1_dd-M2*L*theta2_d^2*sin(theta2)+M2*L*theta2_dd*cos(theta2);
Fv = -M2*L*theta2_dd*sin(theta2)-M2*L*theta2_d*cos(theta2)+M2*g;
FR1 = M1*R1*theta1_dd+M2*R1*theta1_dd-M2*L*theta2_d^2*sin(theta2)+M2*L*theta2_dd*cos(theta2);
FR0 = (-J0*(R1)*theta1_dd+tau*R0-be*R1*theta1_d)/R0^2;

% functions of theta1_dd
func1 = (R1*FR0-R1*FR1)/J1; 
% funtion of theta2_dd
func2 = (L*Fv*sin(theta2)-L*fh*cos(theta2)+tau)/J2;

eq1 = solve(func1 == theta1_dd,theta1_dd); % equation of theta1_dd with theta2_dd in it. 
eq2 = solve(func2 == theta2_dd,theta2_dd); % equation of theta2_dd with theta1_dd in it.
eq_theta1_dd = subs(eq1,theta2_dd,eq2);  
sol_theta1_dd = solve(eq_theta1_dd == theta1_dd,theta1_dd); % final equation of theta 1 double dot

eq_theta2_dd = subs(eq2,theta1_dd,eq1);
sol_theta2_dd = solve(eq_theta2_dd == theta2_dd,theta2_dd); % final equation of theta 2 double dot
% The equations obtained above are fed to the continuous Non linear model
% function CNL_func
%% Part 3- Continuous Linear Time Invariant system
% State Matrices of Linearized model
A=[0 1 0 0;
    0 -(h3*h5)/((h1*h5)-h2^2) (h6*h2)/(h1*h5-h2^2) 0;
    0 0 0 1;
    0 (h2/h5)*(h3*h5)/((h1*h5-h2^2)) -(h6/h5)-((h6*h2^2)/(h5*(h1*h5-h2^2))) 0];
%A= [2,1,0,0;0,-0.669812166454980,-477.394460398384,0;2,0,2,1;0,0.147656370916298,122.678956603377,1];
B= [0; (-h2+(h4*h5))/((h1*h5)-h2^2); 0; ((h2^2-h4*h5*h2)/(h1*h5^2-h2^2*h5))+1/h5];
C= [1 0 0 0; 0 0 1 0];
D= [0;0];
%Transfer function of CLTI system
sys = ss(A,B,C,D);
tf(sys);
% 3.3 Stability of Linear system
[vec,Lam]=eig(A); % Eigen values are [0,-11.3696,-0.0952,10.7949] 
% System is unstable due to positive eigen value here.
% 3.4 Controllability and Observability
n= length(A);
a= rank(ctrb(sys));
if a==n
    disp('Continuos system is controllable');
end
% The system comes out to be controllable
c= rank(obsv(sys));
if c==n
    disp('Continuos system is Observable');
end
% The system is observable
%% 4- System Control
% 4.1- Simulating open loop model for the system

time = 0.1;
timestep = 0.0001;
dt = 0:timestep:time;
x0 = [0; 0; 0.01; 0];% the inital condition of the robot
x1 = [0; 0; 0; 0];% the final condition of the robot
fun = @(t) expm(A*(t))*(B*B')*expm(A'*(t)); % function for controlability Gramian.
Wc = integral(fun, 0,time, 'ArrayValued', true); % integration for controlability Gramian.
count=1;
for n=0:timestep:time
    ut(count)= -B'*expm(A'*(time-n))*inv(Wc)*((expm(A*time))*x0-x1);
    count =count+1;
end

t= 0:timestep:time;
ug=[t' ut'];
figure
plot(t,ut)
hold on
title("Input plotted against the time")
xlabel('Time (t)')
ylabel('Minimal Input Energy (u)')
[st,m,n]=lsim(sys,ut,t,x0);
figure
plot(m,st)
legend('Output 1','Output 2')
xlabel('Time(t)')
ylabel('Amplitude')
title('Output Vs Time for an Open Loop System')

%% 4.2 Closed Loop state feedback system
time= 2;
timestep= 0.01
pole= [-17.003 -10.0708 -10 -3.1217]; % stable poles of the system.
KvalueCT= place(A,B,pole); % k value after placing the poles
ACT=A-B*KvalueCT;% New A matrix of the system
sys_cl= ss(ACT,B,C,D); % new sytsem 
fun= @(t) expm(ACT.*t)*(B*B')*expm(ACT'.*t);
q= integral(fun,0,time,'ArrayValued',true);
count=1;
for n=0:timestep:time
    ut2(count)= -B'*expm(ACT'*(time-n))*inv(q)*((expm(ACT*time))*x0-x1);
    count =count+1;
end
sh= 0:timestep:time;
u2g=[sh' ut2'];
figure
plot(0:timestep:time,ut2)
hold on
title("Input plotted against the time")
xlabel('Time (t)')
ylabel('Minimal Input Energy (u)')
figure
[st3,m,n]=lsim(sys_cl,ut2,0:timestep:time,x0);
plot(m,st3)
legend('Output 1','Output 2')
xlabel('Time(t)')
ylabel('Amplitude')
title('Output Vs Time for a Closed Loop System')

%% 4.3 Open Loop Estimator with no feedback
time = 0.1;
timestep=0.0001;
fun = @(t) expm(A*(t))*(B*B')*expm(A'*(t)); % function for controlability Gramian.
Wc = integral(fun, 0,time, 'ArrayValued', true); % integration for controlability Gramian.
count=1;
for n=0:timestep:time
    ut3(count)= -B'*expm(A'*(time-n))*inv(Wc)*((expm(A*time))*x0-x1);
    count =count+1;
end
figure
plot(0:timestep:time,ut3)
hold on
title("Input plotted against the time")
xlabel('Time (t)')
ylabel('Minimal Input Energy (u)')
syses= ss(A,B,eye(4),[0;0;0;0]);
[st1,m,n]=lsim(syses,ut3,0:timestep:time,x0);
figure
plot(m,st1)
hold on
plot(m,st,'--')
legend('Estimated Output 1','Estimated Output 2','Estimated Output 3','Estimated Output 4','Output 1','Output 2')
xlabel('Time(t)')
ylabel('Amplitude')
title('Output Vs Time for a Closed Loop Estimator  without Feedback')

%% 4.4 Closed Loop System Estimator no feedback
pole= [-17.003 -10.0708 -10 -3.1217];
Lpole=5*pole;
LvalueCT= place(A',C',Lpole)'; %- See Simulink
%% 4.5 See simulink
%% 5.1 Continuous Non Linear system- In Simulink
%% 5.2 COntinuous Non Linear system with Noise- In Simulink