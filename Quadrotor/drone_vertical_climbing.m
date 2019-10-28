% This code runs a simulation for vertical climbing of the drone and also
% computes the matrices for the linearization:
clear all; close all; clc;
% Parameters of the given system:
m    = 1;       % Kg
I    = 0.01;    % Kg*m^2
g    = 9.81;    % m/s^2
b    = 0.1;     % Kg/s
T1   = g;       % Thrust vector for right-side propellers.
T2   = g;       % Thrust vector for right-side propellers.
% Conditions of the non-equilibrium trajectory:
x10     = 0; % Orientation of the drone in space. Given by qtheta. Rad
x1_t0   = 0; % Rotational velocity and acceleration of the drone.  Rad
x1_tt0  = 0;
x20     = 5; % Initial horizontal position of drone in space.
x2_t0   = 0; % Initial spatial velocity in Hz direction.
x2_tt0  = 0;
x30     = 5; 
x3_t0   = 0; % The drone starts climbing from rest position.
% Intial vertical position of the drone in space.
% Since the drone is climbing in accelerated manner the
% accleration in the vertical direction is non-zero.

% Let us define out simulation time characteristics:
T   = 2;           % Total simulation time
dt  = 0.01;         % Step-size
ns  = round(T/dt);  % Number of steps
t   = 0:dt:T-dt;       % Time series vector

% Let us initialize our time-series state vectors:
x1    = zeros(1,ns);
x2    = zeros(1,ns);
x3    = zeros(1,ns);
x1_t  = zeros(1,ns);
x2_t  = zeros(1,ns);
x3_t  = zeros(1,ns);
x1_tt = zeros(1,ns);
x2_tt = zeros(1,ns);
x3_tt = zeros(1,ns);
TR1   = zeros(1,ns);
TL1   = zeros(1,ns);
% Simulation:
for i = 1:ns
    if i == 1
        x1(i)   = x10;
        x2(i)   = x20;
        x3(i)   = x30;
        x1_t(i) = x1_t0;
        x2_t(i) = x2_t0;
        x3_t(i) = x3_t0;
        x1_tt(i)= x1_tt0;
        x2_tt(i)= x2_tt0;
        x3_tt(i)= -g + ((T1+T2)*cos(x1(i)))/m - (b*x3_t(i))/m;
        TR1(i)   = T1;
        TL1(i)   = T1;
    else
        x1(i)   = x1(i-1) + x1_t(i-1)*dt + x1_tt(i-1)*(dt^2)/2;
        x2(i)   = x2(i-1) + x2_t(i-1)*dt + x2_tt(i-1)*(dt^2)/2;
        x3(i)   = x3(i-1) + x3_t(i-1)*dt + x3_tt(i-1)*(dt^2)/2;
        x1_t(i) = x1_t(i-1) + x1_tt(i-1)*dt;
        x2_t(i) = x2_t(i-1) + x2_tt(i-1)*dt;
        x3_t(i) = x3_t(i-1) + x3_tt(i-1)*dt;
        x1_tt(i)= (T1-T2)/I;
        x2_tt(i)= ((T1+T2)*sin(x1(i)))/m - (b*x2_t(i))/m;
        x3_tt(i)= -g + ((T1+T2)*cos(x1(i)))/m - (b*x3_t(i))/m;
        TR1(i)   = T1;
        TL1(i)   = T1;
    end
end

figure()
plot(t,x3,'LineWidth',1.2);
xlabel('Time(s)');
ylabel('Vertical position');
title('Vertical position Vs time');
legend('Vertical Position of Drone');
figure()
hold on; grid on;
plot(t,x2,'LineWidth',1.2);
plot(t,x1,'LineWidth',1.2);
plot(t,TR1,'LineWidth',1.2);
xlabel('Time(s)');
ylabel('Horizontal position and Rotation angle');
title('Position and Orientation plot for Vertically climbing Drone:');
legend('Horizontal Position of Drone','Orientation of Drone','Thrust');
% The vertical position first increases exponentially and then increases in
% a linear fashion much later. This is because of the damping term in
% acceleration which makes it slowly drop. After acceleration becomes zero,
% the velocity is constant and hence the position change becomes linear.

% Here the trajectory is chosen in such a way that the resulting
% linearization is still Time-invariant.

% Plot of spatial trajectory of the drone:
figure()
plot(x2,x3,'LineWidth',1.2);
xlabel('Horizontal Position');
ylabel('Vertical Position');
title('Trajectory taken by drone:');

% Interesting trajectory 2:
% Now what if we make our torques sinusoidal and time-varying. Let us try
% that:

% Simulation 2:
for i = 1:ns
    if i == 1
        x1(i)   = x10;
        x2(i)   = x20;
        x3(i)   = x30;
        x1_t(i) = x1_t0;
        x2_t(i) = x2_t0;
        x3_t(i) = x3_t0;
        x1_tt(i)= x1_tt0;
        x2_tt(i)= x2_tt0;
        x3_tt(i)= -g + ((T1*sin(i*dt)+T2*cos(i*dt))*cos(x1(i)))/m - (b*x3_t(i))/m;
        TR(i)   = T1*sin(i*dt);
        TL(i)   = T2*cos(i*dt);
    else
        x1(i)   = x1(i-1) + x1_t(i-1)*dt + x1_tt(i-1)*(dt^2)/2;
        x2(i)   = x2(i-1) + x2_t(i-1)*dt + x2_tt(i-1)*(dt^2)/2;
        x3(i)   = x3(i-1) + x3_t(i-1)*dt + x3_tt(i-1)*(dt^2)/2;
        x1_t(i) = x1_t(i-1) + x1_tt(i-1)*dt;
        x2_t(i) = x2_t(i-1) + x2_tt(i-1)*dt;
        x3_t(i) = x3_t(i-1) + x3_tt(i-1)*dt;
        x1_tt(i)= (T1*sin(i*dt)-T2*cos(i*dt))/I;
        x2_tt(i)= ((T1*sin(i*dt)+T2*cos(i*dt))*sin(x1(i)))/m - (b*x2_t(i))/m;
        x3_tt(i)= -g + ((T1*sin(i*dt)+T2*cos(i*dt))*cos(x1(i)))/m - (b*x3_t(i))/m;
        TR(i)   = T1*sin(i*dt);
        TL(i)   = T2*cos(i*dt);
    end
end
figure()
plot(t,x3,'LineWidth',1.2);
hold on; grid on;
plot(t,x2,'LineWidth',1.2);
xlabel('Time(s)');
ylabel('Vertical position');
title('Vertical position Vs time');
legend('Vertical Position of Drone','Horizontal Position of Drone');
figure()
hold on; grid on;
plot(t,x1,'LineWidth',1.2);
plot(t,TR,'LineWidth',1.2);
plot(t,TL,'LineWidth',1.2);
xlabel('Time(s)');
ylabel('Rotation angle and thrusts');
title('Orientation and thrusts plot for Vertically climbing Drone:');
legend('Orientation of Drone','Right Thrust','Left Thrust');

% Plot of spatial trajectory of the drone:
figure()
plot(x2,x3,'LineWidth',1.2);
xlabel('Horizontal Position');
ylabel('Vertical Position');
title('Trajectory taken by drone:');

% In this case, the thrust from the right-side propellers starts at 0 as
% it is defined as a Sine function. Hence, the drone banks to the left and
% falls off!! Then if we extend the plot, it will fall off in a sinusoidal
% manner to the left of its satrting position. If we swap the functions for
% the thrusts, the same thing will happen on the other side of the drone.

% Interesting trajectory 3:
% Now what if we make our torques square-sinusoidal and time-varying. Let us try
% that:
% Simulation 3:
for i = 1:ns
    if i == 1
        x1(i)   = x10;
        x2(i)   = x20;
        x3(i)   = x30;
        x1_t(i) = x1_t0;
        x2_t(i) = x2_t0;
        x3_t(i) = x3_t0;
        x1_tt(i)= x1_tt0;
        x2_tt(i)= x2_tt0;
        x3_tt(i)= -g + ((T1*(sin(i*dt)^2)+T2*(cos(i*dt))^2)*cos(x1(i)))/m - (b*x3_t(i))/m;
        TR(i)   = T1*(sin(i*dt)^2);
        TL(i)   = T2*(cos(i*dt)^2);
    else
        x1(i)   = x1(i-1) + x1_t(i-1)*dt + x1_tt(i-1)*(dt^2)/2;
        x2(i)   = x2(i-1) + x2_t(i-1)*dt + x2_tt(i-1)*(dt^2)/2;
        x3(i)   = x3(i-1) + x3_t(i-1)*dt + x3_tt(i-1)*(dt^2)/2;
        x1_t(i) = x1_t(i-1) + x1_tt(i-1)*dt;
        x2_t(i) = x2_t(i-1) + x2_tt(i-1)*dt;
        x3_t(i) = x3_t(i-1) + x3_tt(i-1)*dt;
        x1_tt(i)= (T1*(sin(i*dt)^2)-T2*(cos(i*dt)^2))/I;
        x2_tt(i)= ((T1*(sin(i*dt)^2)+T2*(cos(i*dt)^2))*sin(x1(i)))/m - (b*x2_t(i))/m;
        x3_tt(i)= -g + ((T1*(sin(i*dt)^2)+T2*(cos(i*dt)^2))*cos(x1(i)))/m - (b*x3_t(i))/m;
        TR(i)   = T1*(sin(i*dt)^2);
        TL(i)   = T2*(cos(i*dt)^2);
    end
end
figure()
plot(t,x3,'LineWidth',1.2);
hold on; grid on;
plot(t,x2,'LineWidth',1.2);
xlabel('Time(s)');
ylabel('Vertical position');
title('Vertical position Vs time');
legend('Vertical Position of Drone','Horizontal Position of Drone');
figure()
hold on; grid on;
plot(t,x1,'LineWidth',1.2);
plot(t,TR,'LineWidth',1.2);
plot(t,TL,'LineWidth',1.2);
xlabel('Time(s)');
ylabel('Rotation angle and thrusts');
title('Orientation and thrusts plot for Vertically climbing Drone:');
legend('Orientation of Drone','Right Thrust','Left Thrust');

% Plot of spatial trajectory of the drone:
figure()
plot(x2,x3,'LineWidth',1.2);
xlabel('Horizontal Position');
ylabel('Vertical Position');
title('Trajectory taken by drone:');

%% Computational routine to generate the state-space for a given time t:
% For the given drone system and chosen initial conditions:
prompt = '\nEnter a time value to linearize the "vertical climbing" system:\n';
%t = input(prompt); this can be used to ask for inputs.
t = 2;
A = [zeros(3,3),eye(3,3);
    0,0,0,0,0,0;
    (T1+T2)/m,0,0,0,-b/m,0;
    0,0,0,0,0,-b/m];
B = [zeros(3,3);
    1/I,-1/I,0;
    0,0,0;
    1/m,1/m,-1;];
C = [eye(3,3),zeros(3,3)];
D = 0;
fprintf('\nVertically climbing system:\n');
fprintf('\nThe A matrix at the given time t:\n');
A
fprintf('\nThe B matrix at the given time t:\n');
B
fprintf('\nThe C matrix at the given time t:\n');
C
fprintf('\nThe D matrix at the given time t:\n');
D

% For the given drone system and chosen initial conditions:
prompt = '\nEnter a time value to linearize the "sinusoidally falling" system:\n';
%t = input(prompt); this can be used to ask for inputs.
t = 2;
A = [zeros(3,3),eye(3,3);
    0,0,0,0,0,0;
    (T1*sin(t)+T2*cos(t))/m,0,0,0,-b/m,0;
    0,0,0,0,0,-b/m];
B = [zeros(3,3);
    1/I,-1/I,0;
    0,0,0;
    1/m,1/m,-1;];
C = [eye(3,3),zeros(3,3)];
D = 0;
fprintf('\nSinusoidally falling system:\n');
fprintf('\nThe A matrix at the given time t:\n');
A
fprintf('\nThe B matrix at the given time t:\n');
B
fprintf('\nThe C matrix at the given time t:\n');
C
fprintf('\nThe D matrix at the given time t:\n');
D

% For the given drone system and chosen initial conditions:
prompt = '\nEnter a time value to linearize the "square sinusoidally falling" system:\n';
%t = input(prompt); this can be used to ask for inputs.
t = 2;
A = [zeros(3,3),eye(3,3);
    0,0,0,0,0,0;
    (T1*(sin(t)^2)+T2*(cos(t)^2))/m,0,0,0,-b/m,0;
    0,0,0,0,0,-b/m];
B = [zeros(3,3);
    1/I,-1/I,0;
    0,0,0;
    1/m,1/m,-1;];
C = [eye(3,3),zeros(3,3)];
D = 0;
fprintf('\nsquare-ssinusoidally falling system:\n');
fprintf('\nThe A matrix at the given time t:\n');
A
fprintf('\nThe B matrix at the given time t:\n');
B
fprintf('\nThe C matrix at the given time t:\n');
C
fprintf('\nThe D matrix at the given time t:\n');
D