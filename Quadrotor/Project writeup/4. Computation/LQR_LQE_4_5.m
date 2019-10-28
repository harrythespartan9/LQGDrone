% This script runs an LQR based controller with LQ State Estimation.
clear all; close all; clc;

%% PARAMTERS:

% Simulation Parameters:
t         = 5;          % Length of simulation.
dt        = 0.01;         % Step size.
N         = t/dt + 1;     % Number of steps.
T         = 0:dt:t;       % This is our time series vector.

% State Disturbance Formulation:
% Here we want to create a state disturbance vector which is as long as the
% number of steps.
mu1        = zeros(6,1); % Gaussian/"white".
sigma1     = zeros(6,6,N); % Covariance of disturbance through time.

for i = 1:N
    sigma1(:,:,i)     = diag(diag(rand(6))); % Randomized uncorrelated covariance
                                         % at every time step.
end

% Observation Noise Formulation:
% Similar to the previous one.
mu2        = zeros(2,1); % Gaussian/"white".
sigma2     = zeros(2,2,N); % Covariance of disturbance through time.

for i = 1:N
    sigma2(:,:,i)     = diag(diag(rand(2))); % Randomized uncorrelated covariance
                                         % at every time step.
end

%% STATE-SPACE:
m = 3; % Number of inputs.
n = 6; % Number of states.
q = 2; % Number of observations.

% Linearization and Discretization done for Hovering:
Abar = [zeros(3,3),eye(3,3); % Discretized system matrix.
       zeros(1,6);
       9.81,0,0,0,-0.1,0;
       0,0,0,0,0,-0.1];

Abar = eye(6,6) + Abar*dt; % 1st order approximation of Ad.

Bbar = [zeros(3,3);
       100,-100,0;
       zeros(1,3);
       1,1,-1];

Bbar = Bbar*dt; % 1st order approximation of Bd.

% Let us say the linearization is done at (10,5) where the drone is
% hovering:
H = 10; V =5; % Horizontal and Vertical position of drone.
Cbar = [0, V^2/((H^2+V^2)^1.5), -H*V/((H^2+V^2)^1.5), 0, 0, 0;
        0, -H*V/((H^2+V^2)^1.5), H^2/((H^2+V^2)^1.5), 0, 0, 0];

% Let us initialize our variables and matrices:
Ot = zeros(6,6,N);
O  = zeros(6,6,N);
L = zeros(6,2,N); % Kalman gain.

A = zeros(6,6,N);
B = zeros(6,3,N);
B62 = zeros(n,m-1,N); % Time-series of truncated input matrix.
% Only this matrix will be used for the input computation.
C = zeros(2,6,N);

% Initialization of variables:
% Let our Linearization start exactly on the optimal trajectory:
x   = zeros(n,1,N); % Actual state.
xhat= zeros(n,1,N); % Estimated state.
y   = zeros(q,1,N); % Actual Observation.
yhat= zeros(q,1,N); % Predicted observation.
delta = zeros(n,1,N); % State disturbance.
eta = zeros(q,1,N); % Observation noise.
u = zeros(m,1,N);   % Time-series of applied input.

% States initial condition:
x(:,1,1) = 20*rand(6,1);

K = zeros(m-1,n,N); % Coefficient of Optimal input. Third state unc!! 2x6
P = zeros(n,n,N); % Coefficient of Optimal cost.
J = zeros(1,N);   % Our scalar cost at each time.
Q = zeros(n,n,N)*dt; % State-Cost time series.
R = zeros(m,m,N)*dt; % Input-Cost time series.
R22 = zeros(m-1,m-1,N)*dt; % Input-comp-Cost time series. 2x2

%% LQR BASED COST:
% Here we will describe the LQR associated cost matrices and they will
% remain constant throught the control time horizon. Meaning our control
% objective will not change during the simulation.

q = (10)*eye(n,n); % Cost associated with state deviation;
r = 0.001*eye(m,m);    % Cost associated with input usage;
                       % Third input is uncontrollable.
r(3,3) = 0;
% Basically, I am choosing to penalize deviation from optimal state 10^4
% more than input usage.

p =(10)*eye(n,n); % Final/Terminal state cost.

%% INPUT AND COST COEFFICIENT COMPUTATION:
% Iteratively find our value for P and K matrix backwards:
for i = N:-1:1
    if i == N
        P(:,:,i) = p; % Cost associated with final state.
        K(:,:,i) = zeros(m-1,n); % Optimal input coefficient at final time.
    else
        Q(:,:,i) = q*dt; % Discretized Q and R.
        R(:,:,i) = r*dt;
        R22(:,:,i) = R(1:m-1,1:m-1,i); % Computing recursive block.
        A(:,:,i) = Abar;
        B(:,:,i) = Bbar;
        B62(:,:,i) = Bbar(:,1:m-1,1);
        
        K(:,:,i) = (R22(:,:,i) + transpose(B62(:,:,i))*P(:,:,i+1)*B62(:,:,i))\...
        transpose(B62(:,:,i))*P(:,:,i+1)*A(:,:,i);
        P(:,:,i) = transpose(A(:,:,i) - B62(:,:,i)*K(:,:,i))*P(:,:,i+1)*(A(:,:,i)...
        - B62(:,:,i)*K(:,:,i)) + transpose(K(:,:,i))*R22(:,:,i)*K(:,:,i) + Q(:,:,i);
    end
end

%% KALMAN GAIN AND COVARIANCE COMPUTATION:
% Iteratively find our value for L and O matrix forward in time:
for i = 1:N
    if i == 1
        % Initially:
        C(:,:,i) = Cbar;
        Rbar = inv(sigma1(:,:,i));
        Ot(:,:,i) = Rbar; % Covariance of initial state disturbance.
        Qbar = inv(sigma2(:,:,i)); % Covariance of observation noise.
        % Find the Kalman gain:
        L(:,:,i)  = Ot(:,:,i)*transpose(C(:,:,i))/...
            (Qbar + C(:,:,i)*Ot(:,:,i)*transpose(C(:,:,i)));
        % Update state Covariance:
        O(:,:,i) = (eye(6,6)-L(:,:,i)*C(:,:,i))*Ot(:,:,i)*...
            transpose((eye(6,6)-L(:,:,i)*C(:,:,i))) + L(:,:,i)*Qbar*...
            transpose(L(:,:,i));
    else
        C(:,:,i) = Cbar;
        Rbar = inv(sigma1(:,:,i));
        Ot(:,:,i) = A(:,:,i-1)*O(:,:,i-1)*transpose(A(:,:,i-1))...
            + Rbar;
        Qbar = inv(sigma2(:,:,i));
        L(:,:,i)  = Ot(:,:,i)*transpose(C(:,:,i))/...
            (Qbar + C(:,:,i)*Ot(:,:,i)*transpose(C(:,:,i)));
        O(:,:,i) = (eye(6,6)-L(:,:,i)*C(:,:,i))*Ot(:,:,i)*...
            transpose((eye(6,6)-L(:,:,i)*C(:,:,i))) + L(:,:,i)*Qbar*...
            transpose(L(:,:,i));
    end
end

% Now we have both our Kalman Gain and Input coefficient for the whole time
% horizon.

%% SIMULATION:
for i = 1:N
    if i == 1
        % Intially we will have our observation. Calculating noise:
        eta(:,1,i) = 0.01*mvnrnd(mu2,sigma2(:,:,i)); % Obs Noise.
        delta(:,1,i) = 0.01*mvnrnd(mu1,sigma1(:,:,i)); % State Disturbance.
        % Deterministic system dynamics:
        y(:,1,i) = C(:,:,i)*x(:,1,i) + eta(:,1,i); % Actual Observation.
        % Estimation system dynamics:
        xhat(:,1,i) = x(:,1,i); % Intial state we know.
        % Input computation from estimated state:
        u(1:m-1,1,i) = -K(:,:,i)*x(:,1,i); % Saying I perfectly know my intial state.
        u(m,1,i)     = 9.81; % Uncontrolled g! 3rd input.
        J(1,i)   = 0.5*transpose(x(:,1,i))*P(:,:,i)*x(:,1,i);
        if ((x(3,1,i) > -10^(-2)) && (x(3,1,i) < 10^(-2))) && (x(2,1,i) > -10^(-2)) && (x(2,1,i) < 10^(-2))
            u(m,1,i) = 0;
        end
    else
        % From second time-step we have to do our shit!
        eta(:,1,i) = 0.01*mvnrnd(mu2,sigma2(:,:,i)); % Obs Noise.
        delta(:,1,i) = 0.01*mvnrnd(mu1,sigma1(:,:,i)); % State Disturbance.
        % Deterministic system dynamics:
        x(:,1,i) = A(:,:,i-1)*x(:,1,i-1) + B(:,:,i-1)*u(:,:,i-1) + delta(:,1,i); % Actual State with disturbance.
        y(:,1,i) = C(:,:,i)*x(:,1,i) + eta(:,1,i); % Actual Observation.
        % Estimation system dynamics:
        yhat(:,1,i) = C(:,:,i)*(A(:,:,i-1)*xhat(:,1,i-1) + B(:,:,i-1)*u(:,:,i-1)); % Observation prediction.
        xhat(:,1,i) = A(:,:,i-1)*xhat(:,1,i-1) + B(:,:,i-1)*u(:,:,i-1) +...
            L(:,:,i)*(y(:,1,i)-yhat(:,1,i)); 
        % State prediction from previous state estimate.
        % Input computation from estimated state:
        u(1:m-1,1,i) = -K(:,:,i)*xhat(:,1,i); % Calculating intial optimal input.
        u(m,1,i)     = 9.81; % Uncontrolled g! 3rd input.
        J(1,i)   = 0.5*transpose(xhat(:,1,i))*P(:,:,i)*xhat(:,1,i);
        if ((xhat(3,1,i) > -10^(-2)) && (xhat(3,1,i) < 10^(-2))) && (xhat(2,1,i) > -10^(-2)) && (xhat(2,1,i) < 10^(-2))
            u(m,1,i) = 0;
        end
    end
end

%% PLOTTING THE RESULTS:

H = zeros(1,N);
V = zeros(1,N);
He = zeros(1,N);
Ve = zeros(1,N);
T1= zeros(1,N); % Input 1 - Thrust 1.
T2= zeros(1,N); % Input 2 - Thrust 2.

for j = 1:N
    H(1,j) = x(2,1,j); % Horizontal position of drone.
    V(1,j) = x(3,1,j); % Vertical position of drone.
    He(1,j) = xhat(2,1,j); % Horizontal position of drone.
    Ve(1,j) = xhat(3,1,j); % Vertical position of drone.
    T1(1,j)= u(1,1,j); % Total drone Left thrust.
    T2(1,j)= u(2,1,j); % Total drone Right thrust.
end

% Actual state vs Estimate state:
figure()
plot(T,H,'r','LineWidth',1.2);
grid on; hold on;
plot(T,He,'g--','LineWidth',2);
plot(T,V,'b','LineWidth',1.2);
plot(T,Ve,'y--','LineWidth',2);
xlabel('Time (s)');
ylabel('Magnitude');
title('Estimated Position vs Actual Position:');
legend('Horizontal Pos','Estimated Horizontal Pos','Vertical Pos','Estimated Vertical Pos');

% Trajectory taken by the drone:
figure()
plot(H,V,'r','LineWidth',1.2);
grid on; hold on;
plot(H(1),V(1),'x','LineWidth',1.2);
plot(He,Ve,'b--','LineWidth',1.2);
xlabel('Horizontal Position');
ylabel('Vertical Position');
title('Actual and Estimated Trajectory:');
legend('Actual Trajectory','Starting point','Estimated Trajectory');

% Variation of Inputs:
figure()
plot(T,T1,'r','LineWidth',1.2);
grid on; hold on;
plot(T,T2,'b','LineWidth',1.2);
xlabel('Time (s)');
ylabel('Magnitude');
title('Input Thrust:');
legend('Thrust 1','Thrust 2');

% Variation of Cost with the Estimated states:
figure()
plot(T,J,'r','LineWidth',1.2);
grid on;
xlabel('Time (s)');
ylabel('Magnitude');
title('Cost based on estimated states:')