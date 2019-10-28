% This script will run an LQE procedure for a DLTV system which is a linear
% approximation of a CNL. Here we will compare the estimated states for an
% arbitrary set of inputs.

%% Linear Quadratic Regulation LQE:

% Simulation Parameters:
t         = 5;          % Length of simulation.
dt        = 0.01;         % Step size.
N         = t/dt + 1;     % Number of steps.
T         = 0:dt:t;       % This is our time series vector.

g         = 9.81;        % accln due to gravity.

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

% From the solution of a 5th order polynomial we have obtained our inputs
% as a function of time for our non-equilibrium trajectory:
u = 9.81*ones(3,1,N);
u(1,1,:) = 0.0048*T.^4 + 0.144*T.^3 - 1.32*T.^2 + 2.4*T + g/2;
u(2,1,:) = 0.0048*T.^4 + 0.144*T.^3 - 1.32*T.^2 + 2.4*T + g/2;
% This input raises the drone which is initially at rest (1,0). Then brings
% it to 10m vt pos in 5 seconds.
% Now we have inputs of the required format.

%% STATE-SPACE:
% This is a DLTV system whose states vary in time.
% Linearization and Discretization done for various instants in time about
% a nominal non-equilibrium climbing trajectory:

% Initializing our system matrices:
Abar = zeros(6,6,N);
Bbar = zeros(6,3,N);
Cbar = zeros(2,6,N);


for i = 1:N
    % Finding the state space representation at each time instant:
    s = dt*i; % Current time instant.
    
    Abar(:,:,i) = [zeros(3,3),eye(3,3); % Discretized system matrix.
       zeros(1,6);
       u(1,1,i)+u(2,1,i),0,0,0,-0.1,0;
       0,0,0,0,0,-0.1];

    Abar(:,:,i) = eye(6,6) + Abar(:,:,i)*dt; % 1st order approximation of Ad.

    Bbar(:,:,i) = [zeros(3,3);
            100,-100,0;
            zeros(1,3);
            1,1,-1];

    Bbar(:,:,i) = Bbar(:,:,i)*dt; % 1st order approximation of Bd.
    
    % Horizontal position at time s:
    H = 1; % The drone is climbing in the vertical direction. In the ideal
           % case there won't be deviation in the non-equi traj.
    V = 0.8*s^3 - 0.24*s^4 + 0.0192*s^5;

    Cbar(:,:,i) = [0, V^2/((H^2+V^2)^1.5), -H*V/((H^2+V^2)^1.5), 0, 0, 0;
        0, -H*V/((H^2+V^2)^1.5), H^2/((H^2+V^2)^1.5), 0, 0, 0];
end

% Let us initialize our variables and matrices:
Ot = zeros(6,6,N);
O  = zeros(6,6,N);
L = zeros(6,2,N); % Kalman gain.

A = zeros(6,6,N);
B = zeros(6,3,N);
C = zeros(2,6,N);

% Initial condition:
% Let our Linearization start exactly on the optimal trajectory:
x   = zeros(6,1,N); % Actual state.
xhat= zeros(6,1,N); % Estimated state.
y   = zeros(2,1,N); % Actual Observation.
yhat= zeros(2,1,N); % Predicted observation.
delta = zeros(6,1,N); % State disturbance.
eta = zeros(2,1,N); % Observation noise.

x(:,1,1) = [0,1,0,0,0,0]'; % Intial rest position of drone.

% We have everything we need now. Let us go to the simulation.
% Iteratively find our value for L and O matrix forward in time:
for i = 1:N
    if i == 1
        % Initially:
        A(:,:,i) = Abar(:,:,i);
        B(:,:,i) = Bbar(:,:,i);
        C(:,:,i) = Cbar(:,:,i);
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
        A(:,:,i) = Abar(:,:,i);
        B(:,:,i) = Bbar(:,:,i);
        C(:,:,i) = Cbar(:,:,i);
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

% Now we have our Kalman Gain values for the whole time horizon. Let is
% simulate our system and see how good our estimation is:

% Start Simulation:
for i = 1:N
    if i == 1
        % Intially we will have our observation. Calculating noise:
        eta(:,1,i) = 0.01*mvnrnd(mu2,sigma2(:,:,i)); % Obs Noise.
        delta(:,1,i) = 0.01*mvnrnd(mu1,sigma1(:,:,i)); % State Disturbance.
        % Deterministic system dynamics:
        y(:,1,i) = C(:,:,i)*x(:,1,i) + eta(:,1,i); % Actual Observation.
        % Estimation system dynamics:
        xhat(:,1,i) = x(:,1,i); % Intial state is known.
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
plot(T,Ve,'k--','LineWidth',2);
xlabel('Time (s)');
ylabel('Magnitude');
title('Estimated Position vs Actual Position:');
legend('Horizontal Pos','Estimated Horizontal Pos','Vertical Pos','Estimated Vertical Pos');

% Actual Trajectory vs Estimated Trajectory:
figure()
plot(H,V,'r','LineWidth',1.2);
grid on; hold on;
plot(H(1),V(1),'x','LineWidth',2);
plot(He,Ve,'g--','LineWidth',2);
xlabel('THorizontal Position (m)');
ylabel('Vertical Position (m)');
title('Estimated Trajectory vs Actual Trajectory:');
legend('Actual Trajectory','Starting Point','Estimated Trajectory');


% Variation of Inputs:
figure()
plot(T,T1,'r','LineWidth',1.2);
grid on; hold on;
plot(T,T2,'b','LineWidth',1.2);
xlabel('Time (s)');
ylabel('Magnitude');
title('Input Thrust:');
legend('Thrust 1','Thrust 2');
