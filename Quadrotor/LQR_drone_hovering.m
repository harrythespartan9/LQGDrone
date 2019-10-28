%% AUTHOR : HARI KRISHNA HARI PRASAD        STUDENT ID: 1869297

% This script runs a numerical algorithm solving the LQR control problem
% for the drone project system about the hovering condition:
% This code will take a different state-space kinda approach compared to
% the previous drone cost evaluation problem:
% clear all;
close all; clc;

%% STATE-SPACE:
% The given system is linearized about the hovering fixed
% point/equilibrium. Then a simple discretization process is carried out
% and the new matrices are denoted by Abar and Bbar. This also
% distinguishes them from their CT counterparts.

% Other properties:
% mass m = 1; I = 0.01; T1 = g/2; T2 = g/2;
% We can see that T1 and T2 together balance the drone perfectly if there
% are no external disturbances.

Abar = [eye(3,3),0.01*eye(3,3); % Discretized system matrix.
       0,0,0,1,0,0;
       0.0981,0,0,0,0.999,0
       0,0,0,0,0,0.999];
fprintf('\nDiscretized system matrix Abar:\n');
disp(Abar);
% The determinant of Abar is non-zero. Hence invertible.

fprintf('\nDiscretized input matrix Bbar:\n');
Bbar = [zeros(3,3);
       1,-1,0;
       zeros(1,3);
       0.01,0.01,-0.01];
disp(Bbar);

%% LQR BASED COST:
% Here we will describe the LQR associated cost matrices and they will
% remain constant throught the control time horizon. Meaning our control
% objective will not change during the simulation.

q = (1)*eye(6,6); % Cost associated with state deviation;
r = 0.001*eye(3,3);    % Cost associated with input usage;
% Basically, I am choosing to penalize deviation from optimal state 10^4
% more than input usage.

% Control objective:
% I will transition the system from hovering at a certain horizontal
% position to a different horizontal position. Hence, the final states will
% be cost accordingly.

p =(1)*eye(6,6);% Let use first try to use this cost description and see
                   % the behaviour. Same as q.
                   
%% SIMULATION: IC1

% Simulation Parameters:
t         = 5;          % Length of simulation.
dt        = 0.01;         % Step size.
N         = t/dt + 1;     % Number of steps.
T         = 0:dt:t;       % This is our time series vector.

% State space definition/initialization:
n = 6;            % Number of states.
m = 3;            % Number of inputs.
x = zeros(n,1,N);   % Time-series of states.
u = zeros(m,1,N);   % Time-series of applied input.
A = zeros(n,n,N);   % Time-series of system matrix.
B = zeros(n,m,N);   % Time-series of input matrix.

A(:,:,1) = Abar; % State-space is constant in time.
B(:,:,1) = Bbar;

x(:,1,1) = [20,10,10,2,5,3]'; % Let us attempt without considering the final
                            % state for now.
                            % Literally dropping my drone from a height of
                            % 10 meters.
u(:,1,1) = [0,0,9.81];     % Intial input.                            
% xf       = [0,0,15,0,0,0]'; % Final state.

K = zeros(m,n,N); % Coefficient of Optimal input.
P = zeros(n,n,N); % Coefficient of Optimal cost.
J = zeros(1,N);   % Our scalar cost at each time.
Q = zeros(n,n,N)*dt; % State-Cost time series.
R = zeros(m,m,N)*dt; % Input-Cost time series.

% Iteratively find our value for P and K matrix backwards:
for i = N:-1:1
    if i == N
        P(:,:,i) = p; % Cost associated with final state.
        K(:,:,i) = zeros(m,n); % Optimal input coefficient at final time.
    else
        Q(:,:,i) = q*dt; % Discretized Q and R.
        R(:,:,i) = r*dt;
        A(:,:,i) = A(:,:,1);
        B(:,:,i) = B(:,:,1);
        
        K(:,:,i) = (R(:,:,i) + transpose(B(:,:,i))*P(:,:,i+1)*B(:,:,i))\...
        transpose(B(:,:,i))*P(:,:,i+1)*A(:,:,i);
        P(:,:,i) = transpose(A(:,:,i) - B(:,:,i)*K(:,:,i))*P(:,:,i+1)*(A(:,:,i)...
        - B(:,:,i)*K(:,:,i)) + transpose(K(:,:,i))*R(:,:,i)*K(:,:,i) + Q(:,:,i);
    end
end
% We know that our third input is uncontrollable. Hence we make the bottom
% row of K always zero. The third input is the ACCELERATION DUE TO GRAVITY!

% Start Simulation:

for i = 1:N
    if i == 1
        % Initial state is already defined.
        u(:,1,i) = -K(:,:,i)*x(:,1,i); % Calculating intial optimal input.
        if x(3,1,i) < 10^-16 && x(3,1,i) < -1*(10^-16)
            u(m,1,i) = 0;
            x(3,1,i) = 0;
        else
            u(m,1,i) = 9.81;
        end
        if x(3,1,i) < 0
            x(3,1,i) = 0;
        end
        J(1,i)   = 0.5*transpose(x(:,1,i))*P(:,:,i)*x(:,1,i);
    else
        x(:,1,i) = A(:,:,i-1)*x(:,1,i-1) + B(:,:,i-1)*u(:,1,i-1); %DT system.
        u(:,1,i) = -K(:,:,i)*x(:,1,i);
        if x(3,1,i) < 10^-16 && x(3,1,i) < -1*(10^-16)
            u(m,1,i) = 0;
            x(3,1,i) = 0;
        else
            u(m,1,i) = 9.81;
        end
        if x(3,1,i) < 0
            x(3,1,i) = 0;
        end
        J(1,i)   = 0.5*transpose(x(:,1,i))*P(:,:,i)*x(:,1,i);
    end
end

%% PLOTTING RESULTS: PART (B): IC1:
% Now let us plot the results we just obtained.
theta = zeros(1,N);
H = zeros(1,N);
V = zeros(1,N);
T1= zeros(1,N); % Input 1 - Thrust 1.
T2= zeros(1,N); % Input 2 - Thrust 2.

for j = 1:N
    theta(1,j) = x(1,1,j); % Angular position fo drone.
    H(1,j) = x(2,1,j); % Horizontal position of drone.
    V(1,j) = x(3,1,j); % Vertical position of drone.
    T1(1,j)= u(1,1,j); % Total drone Left thrust.
    T2(1,j)= u(2,1,j); % Total drone Right thrust.
end

% Progression of states:
figure()
plot(T,theta,'LineWidth',1.2);
hold on; grid on;
plot(T,V,'LineWidth',1.2);
plot(T,H,'LineWidth',1.2);
xlabel('Time');
ylabel('Magnitude');
title('State Progression: IC1');
legend('Angular Position of Drone. \theta','Vertical Position of Drone. H','Horizontal Position of Drone. V');
% THE ANGULAR POSITION IS APPROXIMATELY 0. HENCE NOT VISIBLE IN THE PLOT.

% Trajectory taken by the drone:
figure()
plot(H(1),V(1),'rx','LineWidth',1.5);
hold on; grid on;
plot(H,V,'b','LineWidth',1.5);
xlabel('Horizontal position. H');
ylabel('Vertical position. V');
title('Drone trajectory: IC1');
legend('Starting position','Trajectory taken');

% Progression of inputs:
figure()
plot(T,T1,'LineWidth',1.2);
xlabel('Time');
ylabel('Magnitude');
title('Inputs T1&T2 Progression: IC1');

% Progression of Cost:
figure()
plot(T,J,'LineWidth',1.2);
xlabel('Time');
ylabel('Cost');
title('Progression of Cost: IC1');
% Successfully shows exponential decay!! :))

%% SIMULATION: IC2

% Initial condition 2:
x(:,1,1) = [-20,-20,10,-10,-5,-3]';

% Iteratively find our value for P and K matrix backwards:
for i = N:-1:1
    if i == N
        P(:,:,i) = p; % Cost associated with final state.
        K(:,:,i) = zeros(m,n); % Optimal input coefficient at final time.
    else
        Q(:,:,i) = q*dt;
        R(:,:,i) = r*dt;
        A(:,:,i) = A(:,:,1);
        B(:,:,i) = B(:,:,1);
        
        K(:,:,i) = (R(:,:,i) + transpose(B(:,:,i))*P(:,:,i+1)*B(:,:,i))\...
        transpose(B(:,:,i))*P(:,:,i+1)*A(:,:,i);
        P(:,:,i) = transpose(A(:,:,i) - B(:,:,i)*K(:,:,i))*P(:,:,i+1)*(A(:,:,i)...
        - B(:,:,i)*K(:,:,i)) + transpose(K(:,:,i))*R(:,:,i)*K(:,:,i) + Q(:,:,i);
    end
end
% We know that our third input is uncontrollable. Hence we make the bottom
% row of K always zero. The third input is the ACCELERATION DUE TO GRAVITY!

% Start Simulation:

for i = 1:N
    if i == 1
        % Initial state is already defined.
        u(:,1,i) = -K(:,:,i)*x(:,1,i); % Calculating intial optimal input.
        if x(3,1,i) < 10^-16 && x(3,1,i) < -1*(10^-16)
            u(m,1,i) = 0;
            x(3,1,i) = 0;
        else
            u(m,1,i) = 9.81;
        end
        if x(3,1,i) < 0
            x(3,1,i) = 0;
        end
        J(1,i)   = 0.5*transpose(x(:,1,i))*P(:,:,i)*x(:,1,i);
    else
        x(:,1,i) = A(:,:,i-1)*x(:,1,i-1) + B(:,:,i-1)*u(:,1,i-1); %DT system.
        u(:,1,i) = -K(:,:,i)*x(:,1,i);
        if x(3,1,i) < 10^-16 && x(3,1,i) < -1*(10^-16)
            u(m,1,i) = 0;
            x(3,1,i) = 0;
        else
            u(m,1,i) = 9.81;
        end
        if x(3,1,i) < 0
            x(3,1,i) = 0;
        end
        J(1,i)   = 0.5*transpose(x(:,1,i))*P(:,:,i)*x(:,1,i);
    end
end

%% PLOTTING RESULTS: PART (B): IC2:
% Now let us plot the results we just obtained.
theta = zeros(1,N);
H = zeros(1,N);
V = zeros(1,N);
T1= zeros(1,N); % Input 1 - Thrust 1.
T2= zeros(1,N); % Input 2 - Thrust 2.

for j = 1:N
    theta(1,j) = x(1,1,j); % Angular position fo drone.
    H(1,j) = x(2,1,j); % Horizontal position of drone.
    V(1,j) = x(3,1,j); % Vertical position of drone.
    T1(1,j)= u(1,1,j); % Total drone Left thrust.
    T2(1,j)= u(2,1,j); % Total drone Right thrust.
end

% Progression of states:
figure()
plot(T,theta,'LineWidth',1.2);
hold on; grid on;
plot(T,V,'LineWidth',1.2);
plot(T,H,'LineWidth',1.2);
xlabel('Time');
ylabel('Magnitude');
title('State Progression: IC2');
legend('Angular Position of Drone. \theta','Vertical Position of Drone. H','Horizontal Position of Drone. V');
% THE ANGULAR POSITION IS APPROXIMATELY 0. HENCE NOT VISIBLE IN THE PLOT.

% Trajectory taken by the drone:
figure()
plot(H(1),V(1),'rx','LineWidth',1.5);
hold on; grid on;
plot(H,V,'b','LineWidth',1.5);
xlabel('Horizontal position. H');
ylabel('Vertical position. V');
title('Drone trajectory: IC2');
legend('Starting position','Trajectory taken');

% Progression of inputs:
figure()
plot(T,T1,'LineWidth',1.2);
xlabel('Time');
ylabel('Magnitude');
title('Inputs T1&T2 Progression: IC2');

% Progression of Cost:
figure()
plot(T,J,'LineWidth',1.2);
xlabel('Time');
ylabel('Cost');
title('Progression of Cost: IC2');
% Successfully shows exponential decay!! :))

%% NUMERICAL OPTIMIZATION: STEEPEST DESCENT APPROACH:
% Now we will approach the same problem with the Hamiltonian approach.
% Here we will be solving constrained optimization problems at each step.
% Then the co-state/Lagrange multipliers will be calculated and the inputs
% will be iteratively updated.

% Simulation Parameters:
t         = 10;          % Length of simulation.
dt        = 0.01;         % Step size.
N         = t/dt + 1;     % Number of steps.
T         = 0:dt:t;       % This is our time series vector.
s         = 1:N;          % Steps vector.

% State space definition/initialization:
n = 6;            % Number of states.
m = 3;            % Number of inputs.
x2 = zeros(n,1,N);   % Time-series of states.
l  = zeros(1,m,N);   % Co-state/Lagrange multiplier.
J2  = zeros(1,N);     % Re-defining our cost.

x2(:,1,1) = [-20,-20,10,-10,-5,-3]'; % Matching the initial condition with
                            % the IC of previous method.
                            % Let us attempt without considering the final
                            % state for now.
                            % Literally dropping my drone from a height of
                            % 10 meters. 
                            
u2 = rand(m,1,N);   % Time-series of randomly applied inputs.

% No initial condition for the optimal input. We will use the analytical
% relationship derived from the Hamiltonian.

lamb= zeros(N,n); % Co-state vector initialization.

% The system and input matrices are the same as last case.
A = Abar;
B = Bbar;

% Cost definition/initialization:
Q = (1)*eye(6,6)*dt; % Cost associated with state deviation;
R = 0.1*eye(3,3)*dt;    % Cost associated with input usage;

% Input perturbation magnitude:
a = 0.01;

% Descent iteration variable:
iter = 1;

% Start simulation:
for i = 1:N
    if i == 1
        % Iteratively finding the optimal input:
        DuJstar = (J(i+1)-J(i))/(u(:,1,i+1)-u(:,1,i));
        while iter < 3*(10^3) || DuJ(1) < DuJstar(1) || DuJ(2) < DuJstar(2)
            DuJ = Jacobiu(x2,u2,i,Q,R);
            u2(:,1,i) = u2(:,1,i) - a*transpose(DuJ);
            iter = iter + 1;
        end
        iter = 1; % Resetting the descent counter.
        % Checking if our copter's vertical position is close to zero:
        if x2(3,1,i) < 10^-16 && x2(3,1,i) > -1*(10^-16)
            u2(m,1,i) = 0;
            x2(3,1,i) = 0; % If zero making gravity zero.
        else
            u2(m,1,i) = 9.81;
        end
        if x2(3,1,i) < 0
            x2(3,1,i) = 0;
        end
        J2(1,i) = 0.5*transpose(x2(:,1,i))*Q*x2(:,1,i)... 
            + 0.5*transpose(u2(:,1,i))*R*u2(:,1,i);
    elseif i == N
        x2(:,1,i) = A*x2(:,1,i-1) + B*u2(:,1,i-1);
        u2(:,1,i) = 0;
        if x2(3,1,i) < 10^-16 && x2(3,1,i) > -1*(10^-16)
            u2(m,1,i) = 0;
            x2(3,1,i) = 0;
        else
            u2(m,1,i) = 9.81;
        end
        if x2(3,1,i) < 0
            x2(3,1,i) = 0;
        end
        J2(1,i) = 0.5*transpose(x2(:,1,i))*Q*x2(:,1,i)... 
            + 0.5*transpose(u2(:,1,i))*R*u2(:,1,i) + sum(J2);
    else
        x2(:,1,i) = A*x2(:,1,i-1) + B*u2(:,1,i-1);
        % Iteratively finding the optimal input:
        DuJstar = (J(i+1)-J(i))/(u(:,1,i+1)-u(:,1,i));
        while iter < 3*(10^3) || DuJ(1) < DuJstar(1) || DuJ(2) < DuJstar(2)
            DuJ = Jacobiu(x2,u2,i,Q,R);
            u2(:,1,i) = u2(:,1,i) - a*transpose(DuJ);
            iter = iter + 1;
        end
        iter = 1; % Resetting the descent counter.
        if x2(3,1,i) < 10^-16 && x2(3,1,i) > -1*(10^-16)
            u2(m,1,i) = 0;
            x2(3,1,i) = 0;
        else
            u2(m,1,i) = 9.81;
        end
        if x2(3,1,i) < 0
            x2(3,1,i) = 0;
        end
        J2(1,i) = 0.5*transpose(x2(:,1,i))*Q*x2(:,1,i)... 
            + 0.5*transpose(u2(:,1,i))*R*u2(:,1,i) + sum(J2);
    end
end

%% PLOTTING RESULTS: STEEPEST DESCENT:
% Now let us plot the results we just obtained.
theta1 = zeros(1,N);
H1 = zeros(1,N);
V1 = zeros(1,N);
T1= zeros(1,N); % Input 1 - Thrust 1.
T2= zeros(1,N); % Input 2 - Thrust 2.

for j = 1:N
    theta1(1,j) = x2(1,1,j); % Angular position of drone.
    H1(1,j) = x2(2,1,j); % Horizontal position of drone.
    V1(1,j) = x2(3,1,j); % Vertical position of drone.
    T1(1,j)= u2(1,1,j); % Total drone Left thrust.
    T2(1,j)= u2(2,1,j); % Total drone Right thrust.
end

% Progression of states:
figure()
plot(T,theta1,'LineWidth',1.2);
hold on; grid on;
plot(T,V1,'LineWidth',1.2);
plot(T,H1,'LineWidth',1.2);
xlabel('Time');
ylabel('Magnitude');
title('State Progression: STEEPEST DESCENT');
legend('Angular Position of Drone. \theta','Vertical Position of Drone. H','Horizontal Position of Drone. V');
% THE ANGULAR POSITION IS APPROXIMATELY 0. HENCE NOT VISIBLE IN THE PLOT.

% Trajectory taken by the drone:
figure()
plot(H1(1),V1(1),'rx','LineWidth',1.5);
hold on; grid on;
plot(H1,V1,'b','LineWidth',1.5);
xlabel('Horizontal position. H');
ylabel('Vertical position. V');
title('Drone trajectory: STEEPEST DESCENT');
legend('Starting position','Trajectory taken');

% Progression of inputs:
figure()
plot(T,T1,'LineWidth',1.2);
xlabel('Time');
ylabel('Magnitude');
title('Inputs T1&T2 Progression: STEEPEST DESCENT');

% Progression of Cost:
figure()
plot(T,J2,'LineWidth',1.2);
xlabel('Time');
ylabel('Cost');
title('Progression of Cost: STEEPEST DESCENT');

%% Calculating the difference between the inputs: SD:

U = zeros(1,N); % To store the norm of the difference between the inputs. 

for c = 1:2
    for b = 1:N
        U(c,b) = abs(u(c,1,b) - u2(c,1,b));
    end
end

figure()
plot(T,U(1,:),'r','LineWidth',1.5);
grid on; hold on;
plot(T,U(2,:),'b','LineWidth',1.5);
xlabel('Time');
ylabel('Magnitude of difference');
title('Norm of the difference in inputs: U1 and U2');
legend('Input 1 difference','Input 2 difference');

%% Start simulation: STEEPEST DESCENT WITH CO-STATE:
for i = 1:N
    if i == 1
        % Iteratively finding the optimal input:
        while iter < 3*(10^3)
            DxJ = Jacobix(x2,u2,i,Q,R,A,B);
            u2(:,1,i) = u2(:,1,i) - a*transpose((-DxJ/A)*B);
            iter = iter + 1;
        end
        iter = 1; % Resetting the descent counter.
        % Checking if our copter's vertical position is close to zero:
        if x2(3,1,i) < 10^-16 && x2(3,1,i) > -1*(10^-16)
            u2(m,1,i) = 0;
            x2(3,1,i) = 0; % If zero making gravity zero.
        else
            u2(m,1,i) = 9.81;
        end
        if x2(3,1,i) < 0
            x2(3,1,i) = 0;
        end
        J2(1,i) = 0.5*transpose(x2(:,1,i))*Q*x2(:,1,i)... 
            + 0.5*transpose(u2(:,1,i))*R*u2(:,1,i);
    elseif i == N
        x2(:,1,i) = A*x2(:,1,i-1) + B*u2(:,1,i-1);
        u2(:,1,i) = 0;
        if x2(3,1,i) < 10^-16 && x2(3,1,i) > -1*(10^-16)
            u2(m,1,i) = 0;
            x2(3,1,i) = 0;
        else
            u2(m,1,i) = 9.81;
        end
        if x2(3,1,i) < 0
            x2(3,1,i) = 0;
        end
        J2(1,i) = 0.5*transpose(x2(:,1,i))*Q*x2(:,1,i)... 
            + 0.5*transpose(u2(:,1,i))*R*u2(:,1,i) + sum(J2);
    else
        x2(:,1,i) = A*x2(:,1,i-1) + B*u2(:,1,i-1);
        % Iteratively finding the optimal input:
        while iter < 3*(10^3)
            DxJ = Jacobix(x2,u2,i,Q,R,A,B);
            u2(:,1,i) = u2(:,1,i) - a*transpose((-DxJ/A)*B);
            iter = iter + 1;
        end
        iter = 1; % Resetting the descent counter.
        if x2(3,1,i) < 10^-16 && x2(3,1,i) > -1*(10^-16)
            u2(m,1,i) = 0;
            x2(3,1,i) = 0;
        else
            u2(m,1,i) = 9.81;
        end
        if x2(3,1,i) < 0
            x2(3,1,i) = 0;
        end
        J2(1,i) = 0.5*transpose(x2(:,1,i))*Q*x2(:,1,i)... 
            + 0.5*transpose(u2(:,1,i))*R*u2(:,1,i) + sum(J2);
    end
end

%% PLOTTING RESULTS: STEEPEST DESCENT:
% Now let us plot the results we just obtained.
theta1 = zeros(1,N);
H1 = zeros(1,N);
V1 = zeros(1,N);
T1= zeros(1,N); % Input 1 - Thrust 1.
T2= zeros(1,N); % Input 2 - Thrust 2.

for j = 1:N
    theta1(1,j) = x2(1,1,j); % Angular position of drone.
    H1(1,j) = x2(2,1,j); % Horizontal position of drone.
    V1(1,j) = x2(3,1,j); % Vertical position of drone.
    T1(1,j)= u2(1,1,j); % Total drone Left thrust.
    T2(1,j)= u2(2,1,j); % Total drone Right thrust.
end

% Progression of states:
figure()
plot(T,theta1,'LineWidth',1.2);
hold on; grid on;
plot(T,V1,'LineWidth',1.2);
plot(T,H1,'LineWidth',1.2);
xlabel('Time');
ylabel('Magnitude');
title('State Progression: STEEPEST DESCENT');
legend('Angular Position of Drone. \theta','Vertical Position of Drone. H','Horizontal Position of Drone. V');
% THE ANGULAR POSITION IS APPROXIMATELY 0. HENCE NOT VISIBLE IN THE PLOT.

% Trajectory taken by the drone:
figure()
plot(H1(1),V1(1),'rx','LineWidth',1.5);
hold on; grid on;
plot(H1,V1,'b','LineWidth',1.5);
xlabel('Horizontal position. H');
ylabel('Vertical position. V');
title('Drone trajectory: STEEPEST DESCENT');
legend('Starting position','Trajectory taken');

% Progression of inputs:
figure()
plot(T,T1,'LineWidth',1.2);
xlabel('Time');
ylabel('Magnitude');
title('Inputs T1&T2 Progression: STEEPEST DESCENT');

% Progression of Cost:
figure()
plot(T,J2,'LineWidth',1.2);
xlabel('Time');
ylabel('Cost');
title('Progression of Cost: STEEPEST DESCENT');

%% Calculating the difference between the inputs: SD with CO-STATE!

U = zeros(1,N); % To store the norm of the difference between the inputs. 

for c = 1:2
    for b = 1:N
        U(c,b) = abs(u(c,1,b) - u2(c,1,b));
    end
end

figure()
plot(T,U(1,:),'g','LineWidth',1.5);
grid on; hold on;
plot(T,U(2,:),'c','LineWidth',1.5);
xlabel('Time');
ylabel('Magnitude of difference');
title('Norm of the difference in inputs: U1 and U2');
legend('Input 1 difference','Input 2 difference');

%% EXTERNAL FUNCTION DEFINITION:
% We will now compute the Jacobian of cost with respect to the two inputs:
%%% WRONG NEED TO CHANGE!!! %%%
function DuJ = Jacobiu(x,u,i,Q,R)
DuJ = zeros(1,3); % Initializing the Jacobian.
% Cost at this time-step for current input:
Ji = 0.5*transpose(x(:,1,i))*Q*x(:,1,i)... 
            + 0.5*transpose(u(:,1,i))*R*u(:,1,i);
% Difference in cost with respect to the next input:
Jii = 0.5*transpose(x(:,1,i))*Q*x(:,1,i)... 
            + 0.5*transpose(u(:,1,i+1))*R*u(:,1,i+1);
for j = 1:3
    DuJ(j) = (Jii - Ji)/(u(j,1,i+1)-u(j,1,i));
    if isnan(DuJ(j)) || DuJ(j) == Inf || DuJ(j) == -Inf
        DuJ(j) = 0;
    end
end
% Now the rectified DuJ will be returned.
end

% We will now compute the Jacobian of cost with respect to the two states:
function DxJ = Jacobix(x,u,i,Q,R,A,B)
DxJ = zeros(1,6); % Initializing the Jacobian.
% Cost at this time-step for current input:
Ji = 0.5*transpose(x(:,1,i))*Q*x(:,1,i)... 
            + 0.5*transpose(u(:,1,i))*R*u(:,1,i);
x(:,1,i+1) = A*x(:,1,i) + B*u(:,1,i);
% Difference in cost with respect to the next input:
Jii = 0.5*transpose(x(:,1,i+1))*Q*x(:,1,i+1)... 
    + 0.5*transpose(u(:,1,i))*R*u(:,1,i);
for j = 1:6
    DxJ(j) = (Jii - Ji)/(x(j,1,i+1)-x(j,1,i));
    if isnan(DxJ(j)) || DxJ(j) == Inf || DxJ(j) == -Inf
        DxJ(j) = 0;
    end
end
% Now the rectified DxJ will be returned.
% Lambda or costate is given by -DxJ*inv(A); This is what we will use to do
% a more specific descent.
end