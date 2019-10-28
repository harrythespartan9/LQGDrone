% This code runs a simulation for vertical climbing of the drone which is
% an equilibrium trajectory. We would also be computing an LQR based cost
% at each instant.
clear all; close all; clc;

% Cost Function:

% There isn't a fixed way to choose my Q and R matrices. They
% (as in their corresponding Eigenvalues, preferably a diagonal matrix) are
% chosen for the Vertical altitude adjustment of Drone. Here some states are
% important, whereas some states don't change with time. Hence, the matrix
% weights (or Eigenvalues) are chosen accordingly.

% Let jstar be my quadratic cost at that instant. The jstar calculated will
% be a time series vector which will then be integrated in each simulation
% until the latest value to obtain the actual cost.

% Cost matrix for states Q(6x6):
Q = [0.1*eye(4,4),zeros(4,2);zeros(2,4),0.2*eye(2,2)]; % This is the cost matrix
% for the states. It is 6x6 and higher priority is given to the last two
% states of vertical position and velocity which is the most related to our
% application.
fprintf('\nCost matrix for the states Q(6x6) is given by:\n');
Q

% Cost matrix for control inputs R(3x3):
R = [0.2,0,0;0,0.2,0;0,0,0.1]; % The third value is basically an uncontrollable
% input of gravitational acceleration. The first two include the
% controllable inputs or thrusts from the right and left motors
% respectively. Hence, it is given a higher cost and priority.
fprintf('\nCost matrix for the states R(3x3) is given by:\n');
R

% Perturbation:

% So here, we are gonna choose 10 trajectories to perturb the drone from
% hovering condition which is our equilibrium. Since we have chosen
% vertical altitude adjustment, I will be doing 5 trajectories with
% climbing and 5 with descent.
% Also I've chosen to characterize my perturbation as the amount of thrust
% acting on the drone.
a1 = [];
a1(1) = 0; % Perturbation Magnitude for the first trajectory.

% We have taken 10 separate trajectories based off of the perturbation
% magnitude in the positive and negative direction.
% Let us introduce a variable which contains the final cost at the
% end of each simulation:
JJ = [];

% Parameters of the given system:
m    = 1;       % Kg
I    = 0.01;    % Kg*m^2
g    = 9.81;    % m/s^2
b    = 0.1;     % Kg/s
T1   = g/2;       % Thrust vector for right-side propellers.
T2   = g/2;       % Thrust vector for right-side propellers.
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
T   = 1;           % Total simulation time
dt  = 0.01;         % Step-size
ns  = round(T/dt);  % Number of steps
t   = 0:dt:T-dt;       % Time series vector

% Let us initialize our time-series process variables:
x1    = [];
x2    = [];
x3    = [];
x1_t  = [];
x2_t  = [];
x3_t  = [];
x1_tt = [];
x2_tt = [];
x3_tt = [];
TR1   = [];
TL1   = [];
X     = [];
U     = [];
jstar = [];
J     = [];

%%%%% Trajectory 1: (Hovering with no perturbation - Cost will be zero)
%%%%% %%%%%
% This is basically a test simulation to check if our cost structure works.
% This is the first Trajectory with a perturbation magnitude of 0;
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
        x3_tt(i)= -g + ((T1 + a1(1)*g + T2 + a1(1)*g)*cos(x1(i)))/m - (b*x3_t(i))/m;
        TR1(i)  = T1 + a1(1)*g;
        TL1(i)  = T2 + a1(1)*g;
        X(:,i)    = [x1(i); x1_t(i); x2(i); x2_t(i); x3(i); x3_t(i)];
        U(:,i)    = [TR1(i); TL1(i); g];
        jstar(i)= 0;
        J(i)    = 0;
    else
        x1(i)   = x1(i-1) + x1_t(i-1)*dt + x1_tt(i-1)*(dt^2)/2;
        x2(i)   = x2(i-1) + x2_t(i-1)*dt + x2_tt(i-1)*(dt^2)/2;
        x3(i)   = x3(i-1) + x3_t(i-1)*dt + x3_tt(i-1)*(dt^2)/2;
        x1_t(i) = x1_t(i-1) + x1_tt(i-1)*dt;
        x2_t(i) = x2_t(i-1) + x2_tt(i-1)*dt;
        x3_t(i) = x3_t(i-1) + x3_tt(i-1)*dt;
        x1_tt(i)= (T1-T2)/I;
        x2_tt(i)= ((T1 + a1(1)*g + T2 + a1(1)*g)*sin(x1(i)))/m - (b*x2_t(i))/m;
        x3_tt(i)= -g + ((T1 + a1(1)*g + T2 + a1(1)*g)*cos(x1(i)))/m - (b*x3_t(i))/m;
        TR1(i)   = T1 + a1(1)*g;
        TL1(i)   = T2 + a1(1)*g;
        X(:,i)    = [x1(i); x1_t(i); x2(i); x2_t(i); x3(i); x3_t(i)];
        U(:,i)    = [TR1(i); TL1(i); g];
        jstar(i)= (X(:,i)-X(:,1))'*Q*(X(:,i)-X(:,1)) + (U(:,i)-U(:,1))'*R*(U(:,i)-U(:,1)); 
        J(i)    = trapz(dt,jstar); % Calculates instantaneous cost value.
    end
end

% This keeps track of the cost accumulated at the end of each trajectory.
JJ(1) = J(end);

% The following plots are shown for the first two trajectories:
figure()
plot(t,x3,'LineWidth',1.2);
xlabel('Time(s)');
ylabel('Vertical position');
title('Vertical position Vs time: T1 Hovering');
legend('Vertical Position of Drone');
figure()
hold on; grid on;
plot(t,x2,'LineWidth',1.2);
plot(t,x1,'LineWidth',1.2);
plot(t,TR1,'LineWidth',1.2);
xlabel('Time(s)');
ylabel('Horizontal position and Rotation angle');
title('Position and Orientation plot: T1 Hovering');
legend('Horizontal Position of Drone','Orientation of Drone','Thrust');

% Plot of spatial trajectory of the drone:
figure()
plot(x2,x3,'kx','LineWidth',1.2); % Small x's in the plot denote positions.
xlabel('Horizontal Position');
ylabel('Vertical Position');
title('Trajectory taken by drone: T1 Hovering');

%%%%% Trajectory 2: %%%%%

a1(2) = 1; % Second trajectory perturbation magnitude.

% Let us initialize our time-series process variables:
x1    = [];
x2    = [];
x3    = [];
x1_t  = [];
x2_t  = [];
x3_t  = [];
x1_tt = [];
x2_tt = [];
x3_tt = [];
TR1   = [];
TL1   = [];
X     = [];
U     = [];
jstar = [];
J     = [];

%%% Simulation: %%%

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
        x3_tt(i)= -g + ((T1 + a1(2)*g + T2 + a1(2)*g)*cos(x1(i)))/m - (b*x3_t(i))/m;
        TR1(i)  = T1 + a1(2)*g;
        TL1(i)  = T2 + a1(2)*g;
        X(:,i)    = [x1(i); x1_t(i); x2(i); x2_t(i); x3(i); x3_t(i)];
        U(:,i)    = [TR1(i); TL1(i); g];
        jstar(i)= 0;
        J(i)    = 0;
    else
        x1(i)   = x1(i-1) + x1_t(i-1)*dt + x1_tt(i-1)*(dt^2)/2;
        x2(i)   = x2(i-1) + x2_t(i-1)*dt + x2_tt(i-1)*(dt^2)/2;
        x3(i)   = x3(i-1) + x3_t(i-1)*dt + x3_tt(i-1)*(dt^2)/2;
        x1_t(i) = x1_t(i-1) + x1_tt(i-1)*dt;
        x2_t(i) = x2_t(i-1) + x2_tt(i-1)*dt;
        x3_t(i) = x3_t(i-1) + x3_tt(i-1)*dt;
        x1_tt(i)= (T1-T2)/I;
        x2_tt(i)= ((T1 + a1(2)*g + T2 + a1(2)*g)*sin(x1(i)))/m - (b*x2_t(i))/m;
        x3_tt(i)= -g + ((T1 + a1(2)*g + T2 + a1(2)*g)*cos(x1(i)))/m - (b*x3_t(i))/m;
        TR1(i)   = T1 + a1(2)*g;
        TL1(i)   = T2 + a1(2)*g;
        X(:,i)    = [x1(i); x1_t(i); x2(i); x2_t(i); x3(i); x3_t(i)];
        U(:,i)    = [TR1(i); TL1(i); g];
        jstar(i)= (X(:,i)-X(:,1))'*Q*(X(:,i)-X(:,1)) + (U(:,i)-U(:,1))'*R*(U(:,i)-U(:,1)); 
        J(i)    = trapz(dt,jstar); % Calculates instantaneous cost value.
    end
end
JJ(2) = J(end);
figure()
plot(t,x3,'LineWidth',1.2);
xlabel('Time(s)');
ylabel('Vertical position');
title('Vertical position Vs time: T2 Climbing');
legend('Vertical Position of Drone');
figure()
hold on; grid on;
plot(t,x2,'LineWidth',1.2);
plot(t,x1,'LineWidth',1.2);
plot(t,TR1,'LineWidth',1.2);
xlabel('Time(s)');
ylabel('Horizontal position and Rotation angle');
title('Position and Orientation plot: T2 Climbing');
legend('Horizontal Position of Drone','Orientation of Drone','Thrust');

% Plot of spatial trajectory of the drone:
figure()
plot(x2,x3,'kx','LineWidth',1.2); % Small x's in the plot denote positions.
xlabel('Horizontal Position');
ylabel('Vertical Position');
title('Trajectory taken by drone: T2 Climbing');

% The other trajectories will be run iteratively:

for k = 3:10
    a1(k) = a1(k-1) + 1;
    % Let us initialize our time-series process variables:
    x1    = [];
    x2    = [];
    x3    = [];
    x1_t  = [];
    x2_t  = [];
    x3_t  = [];
    x1_tt = [];
    x2_tt = [];
    x3_tt = [];
    TR1   = [];
    TL1   = [];
    X     = [];
    U     = [];
    jstar = [];
    J     = [];

    %%% Simulation: %%%

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
            x3_tt(i)= -g + ((T1 + a1(k)*g + T2 + a1(k)*g)*cos(x1(i)))/m - (b*x3_t(i))/m;
            TR1(i)  = T1 + a1(k)*g;
            TL1(i)  = T2 + a1(k)*g;
            X(:,i)    = [x1(i); x1_t(i); x2(i); x2_t(i); x3(i); x3_t(i)];
            U(:,i)    = [TR1(i); TL1(i); g];
            jstar(i)= 0;
            J(i)    = 0;
        else
            x1(i)   = x1(i-1) + x1_t(i-1)*dt + x1_tt(i-1)*(dt^2)/2;
            x2(i)   = x2(i-1) + x2_t(i-1)*dt + x2_tt(i-1)*(dt^2)/2;
            x3(i)   = x3(i-1) + x3_t(i-1)*dt + x3_tt(i-1)*(dt^2)/2;
            x1_t(i) = x1_t(i-1) + x1_tt(i-1)*dt;
            x2_t(i) = x2_t(i-1) + x2_tt(i-1)*dt;
            x3_t(i) = x3_t(i-1) + x3_tt(i-1)*dt;
            x1_tt(i)= (T1-T2)/I;
            x2_tt(i)= ((T1 + a1(k)*g + T2 + a1(k)*g)*sin(x1(i)))/m - (b*x2_t(i))/m;
            x3_tt(i)= -g + ((T1 + a1(k)*g + T2 + a1(k)*g)*cos(x1(i)))/m - (b*x3_t(i))/m;
            TR1(i)   = T1 + a1(k)*g;
            TL1(i)   = T2 + a1(k)*g;
            X(:,i)    = [x1(i); x1_t(i); x2(i); x2_t(i); x3(i); x3_t(i)];
            U(:,i)    = [TR1(i); TL1(i); g];
            jstar(i)= (X(:,i)-X(:,1))'*Q*(X(:,i)-X(:,1)) + (U(:,i)-U(:,1))'*R*(U(:,i)-U(:,1)); 
            J(i)    = trapz(dt,jstar); % Calculates instantaneous cost value.
        end
    end
    JJ(k) = J(end);
end

% Plot of cost vs perturbation magnitude:
% Perturbation magnitude is given by our a1 vector.
% Cost is given by our JJ vector.
figure()
xlabel('Perturbation Magnitude');
ylabel('Cost for each trajectory');
title('Cost vs PertMag - For 10 trajectories!')
hold on; grid on;
for l = 1:10
    plot(a1(l),JJ(l),'bx');
end