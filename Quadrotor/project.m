close all, clear all, clc
% problem 3:
% 3.2:
m = 1;
I = 0.01;
g = 9.81;
b = 0.1;
u1_lin = m*g/2;
u2_lin = m*g/2;
qt_lin = 0;
A_lin = dynamics_mat(u1_lin,u2_lin,qt_lin,m,b);
B_lin = input_mat(I,qt_lin,m);
u_lin = [m*g/2; m*g/2; m*g];
D_lin = 0;
qv_lin = 40;
qh_lin = 0;
C_lin = output_mat(qv_lin,qh_lin);

%-----------------------------------------------------------

%3.3:
p0 = 0;
p1 = 0;
p2 = g/2;
tf = 10;
traj_coeff = traj_mat(p2,tf);
dt = 0.0001;
t = (0:dt:tf);
qv_nlin = zeros(size(t));
qv_dot_nlin = zeros(size(t));
qv_double_dot_nlin = zeros(size(t));
qv_nlin(1) = 0;
qv_dot_nlin(1) = 0;
qv_double_dot_nlin(1) = g;
traj_coeff = traj_mat(p2,tf);

for i=2:length(t)
    qv_vals = qv_traj(traj_coeff,t(i));
    qv_nlin(i) = qv_vals(1);
    qv_dot_nlin(i) = qv_vals(2);
    qv_double_dot_nlin(i) = qv_vals(3);
end

thrusts = 0.5.*((m.*qv_double_dot_nlin)+(m*g)+(b*qv_dot_nlin));
u1_nlin = thrusts;
u2_nlin = thrusts;
u3_nlin = m*g;
qt_nlin = 0;
states_nlin = zeros(6,length(t));
states_nlin(3,1) = 3;
B_nlin = input_mat(I,qt_nlin,m);

for i=1:1:(length(t)-1)
    A_nlin = dynamics_mat(u1_nlin(i),u2_nlin(i),qt_nlin,m,b);
    Ad_nlin = eye(6)+(dt.*A_nlin);
    Bd_nlin = dt.*B_nlin;
    states_nlin(:,i+1) = Ad_nlin*states_nlin(:,i)+(Bd_nlin*[u1_nlin(i);u2_nlin(i);u3_nlin]);         
end  
%-----------------------------------------------
%3.4:
% figure(1)
% subplot(6,1,1)
% plot(t,states_nlin(1,:));
% title('qtheta vs Time')
% xlabel('Time')
% ylabel('qtheta')
% subplot(6,1,2)
% plot(t,states_nlin(2,:));
% title('qtheta dot vs Time')
% xlabel('Time')
% ylabel('qtheta dot')
% subplot(6,1,3)
% plot(t,states_nlin(3,:));
% title('qH vs Time')
% xlabel('Time')
% ylabel('qH')
% subplot(6,1,4)
% plot(t,states_nlin(4,:));
% title('qH dot vs Time')
% xlabel('Time')
% ylabel('qH dot')
% subplot(6,1,5)
% plot(t,states_nlin(5,:));
% title('qV vs Time')
% xlabel('Time')
% ylabel('qV')
% subplot(6,1,6)
% plot(t,states_nlin(6,:));
% title('qV dot vs Time')
% xlabel('Time')
% ylabel('qV dot')

%------------------------------------------------------
%Q4:
%4.1:
time_step = 0.01;
Ad_lin = eye(6)+(time_step*A_lin);
Bd_lin = (time_step*B_lin);
Cd_lin = output_mat(1,0.2);
Dd_lin = 0;
t_horizon = (0:time_step:10);

%-------------------------------------------------------
%4.2:
J = zeros(size(t_horizon));
P = zeros(6,6,length(t_horizon));
K = zeros(2,6,length(t_horizon));
Ud_lin = zeros(3,length(t_horizon));

Q = [10 0 0 0 0 0; 0 10 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 10 0; 0 0 0 0 0 10]*time_step;
R = 0.001*eye(3)*time_step;

R(3,3) = 0;
P(:,:,end) = Q;
K(:,:,end) = zeros(2,6);

for i=(length(t_horizon)-1):-1:1
     K_p1 = (transpose(Bd_lin(:,1:2))*P(:,:,i+1)*(Bd_lin(:,1:2))+(R(1:2,1:2)));
     K_p2 = transpose(Bd_lin(:,1:2))*P(:,:,i+1)*Ad_lin;
     K(:,:,i) = inv(K_p1)*K_p2;
     P_p1 = Ad_lin - (Bd_lin(:,1:2)*K(:,:,i));
     P_p2 = P(:,:,i+1);
     P_p3 = Ad_lin - (Bd_lin(:,1:2)*K(:,:,i));
     P_p4 = (transpose(K(:,:,i))*R(1:2,1:2)*K(:,:,i)) + Q;
     P(:,:,i) = (transpose(P_p1)*P_p2*P_p3)+P_p4;
 end

Xd_lin = zeros(6,length(t_horizon));
X0 = [2;1;2;3.3;4.4;1.2];
Xd_lin(:,1) = X0;

for i=2:1:(length(t_horizon))
    if ((Xd_lin(5,i-1) < 10^(-2))&& (Xd_lin(5,i-1)>-10^(-2)))
        Ud_lin(3,i-1) = 0;
    else
        Ud_lin(3,i-1) = g;
    end
    
    Ud_lin(1:2,i-1) = -K(:,:,i-1)*Xd_lin(:,i-1);
    Xd_lin(:,i) = (Ad_lin*Xd_lin(:,i-1)) + (Bd_lin*Ud_lin(:,i-1));
        
end

for i=1:1:length(t_horizon)
    J(i) = 0.5.*(transpose(Xd_lin(:,i))*P(:,:,i)*Xd_lin(:,i));
end

% figure(2)
% plot(t_horizon,J)
% title("Cost Function vs Time")
% xlabel("Time")
% ylabel("Cost Function (J)")
% 
% figure(3)
% subplot(6,1,1)
% plot(t_horizon,Xd_lin(1,:));
% title('qtheta vs Time to Minimize J')
% xlabel('Time')
% ylabel('qtheta')
% subplot(6,1,2)
% plot(t_horizon,Xd_lin(2,:));
% title('qtheta dot vs Time to Minimize J')
% xlabel('Time')
% ylabel('qtheta dot')
% subplot(6,1,3)
% plot(t_horizon,Xd_lin(3,:));
% title('qH vs Time to Minimize J')
% xlabel('Time')
% ylabel('qH')
% subplot(6,1,4)
% plot(t_horizon,Xd_lin(4,:));
% title('qH dot vs Time to Minimize J')
% xlabel('Time')
% ylabel('qH dot')
% subplot(6,1,5)
% plot(t_horizon,Xd_lin(5,:));
% title('qV vs Time to Minimize J')
% xlabel('Time')
% ylabel('qV')
% subplot(6,1,6)
% plot(t_horizon,Xd_lin(6,:));
% title('qV dot vs Time to Minimize J')
% xlabel('Time')
% ylabel('qV dot')

%-------------------------------------------------------------------

%4.4:
P_k_tilda = zeros(6,6,length(t));
P_k = zeros(6,6,length(t));
K_k = zeros(6,2,length(t));
cov_SD = diag(diag(rand(6)));
P_k_tilda(:,:,1) = cov_SD;
A_k = zeros(6,6,length(t));
B_k = B_nlin*dt;
C_k = zeros(2,6,length(t));
D_k = 0;
Xactual = states_nlin;
Ydist = zeros(2,1,length(t));


for i=1:1:length(t)
    cov_OD = diag(diag(rand(2)));
    if(i>=2)
        cov_SD = diag(diag(rand(6)));
    end
    A_k(:,:,i) = dynamics_mat(u1_nlin(i),u2_nlin(i),0,m,b)*dt + eye(6);
    
    C_k(:,:,i) = output_mat(states_nlin(5,i),states_nlin(3,i));
    
    if (i>=2)
        P_k_tilda(:,:,i) = (A_k(:,:,i-1)*P_k(:,:,i-1)*transpose(A_k(:,:,i-1)))+cov_SD;
    end
    K_kp1 = P_k_tilda(:,:,i)*transpose(C_k(:,:,i));
    K_kp2 = cov_OD+(C_k(:,:,i)*P_k_tilda(:,:,i)*transpose(C_k(:,:,i)));
    K_k(:,:,i) = K_kp1*inv(K_kp2);
    P_kp1 = (eye(6)-(K_k(:,:,i)*C_k(:,:,i)));
    P_kp2 = P_k_tilda(:,:,i);
    P_kp3 = (K_k(:,:,i)*cov_OD*transpose(K_k(:,:,i)));
    P_k(:,:,i) = ((P_kp1)*(P_kp2)*transpose(P_kp1))+P_kp3;
    Ydist(:,:,i) = output_mat(Xactual(5,i),Xactual(3,i))*Xactual(:,i);
    dist = 0.05.*transpose(mvnrnd([0;0],cov_OD));
    Ydist(:,:,i) = Ydist(:,:,i)+(dist);
end

Xestimate = zeros(6,length(t));
Xestimate(:,1) = Xactual(:,1);
Xtilda = zeros(6,1);

for i=2:1:length(t)
    Xtp1 = ((dynamics_mat(u1_nlin(i-1),u2_nlin(i-1),0,m,b))*dt)+eye(6);
    Xtp2 = Xestimate(:,i-1);
    Xtp3 = B_k;
    
    if (Xestimate(5,i-1)>-10^(-5) && Xestimate(5,i-1)<10^(-5))
        input_now = [u1_nlin(i-1);u2_nlin(i-1);0];
    else
        input_now = [u1_nlin(i-1);u2_nlin(i-1);9.81];
    end
    Xtilda = (Xtp1*Xtp2)+(B_k*input_now);
    diff = Ydist(:,:,i)-(output_mat(Xtilda(5),Xtilda(3))*Xtilda);
    Xestimate(:,i) = Xtilda+(K_k(:,:,i)*diff);
    
end

figure(4)
plot(t,Xestimate(5,:))
hold on
plot(t,Xactual(5,:))
title("Comparison between estimated and actual Vertical position qV (DLTV system)")
xlabel("Time")
ylabel("qV")
legend("Estimated","Actual")

%---------------------------------------------------------------------
%Building estimator for DLTI from 4.4 and then doing 4.5:

P_s_tilda = zeros(6,6,length(t_horizon));
P_s = zeros(6,6,length(t_horizon));
K_s = zeros(6,2,length(t_horizon));
cov_SD_lti = diag(diag(rand(6)));
P_s_tilda(:,:,1) = cov_SD_lti;
A_s = Ad_lin;
B_s = Bd_lin;
C_s = Cd_lin;
D_s = 0;
Xactual_lti = Xd_lin;
Ydist_lti = zeros(2,1,length(t_horizon));

for i=1:1:length(t_horizon)
    cov_OD_lti = diag(diag(rand(2)));
    if (i>=2)
        cov_SD_lti = diag(diag(rand(6)));
        P_s_tilda(:,:,i) = (A_s*P_s(:,:,i-1)*transpose(A_s))+cov_SD_lti;
    end
    K_kp1 = P_s_tilda(:,:,i)*transpose(C_s);
    K_kp2 = cov_OD_lti+(C_s*P_s_tilda(:,:,i)*transpose(C_s));
    K_s(:,:,i) = K_kp1*inv(K_kp2);
    P_kp1 = (eye(6)-(K_s(:,:,i)*C_s));
    P_kp2 = P_s_tilda(:,:,i);
    P_kp3 = (K_s(:,:,i)*cov_OD_lti*transpose(K_s(:,:,i)));
    P_s(:,:,i) = ((P_kp1)*(P_kp2)*transpose(P_kp1))+P_kp3;
end


Xestimate_lti = zeros(6,length(t_horizon));
Xestimate_lti(:,1) = Xactual_lti(:,1);
Xtilda_lti = zeros(6,1);

for i=2:1:length(t_horizon)
    Xtp1 = (A_s);
    Xtp2 = Xestimate_lti(:,i-1);
    Xtp3 = B_s;
    cov_OD_lti = diag(diag(rand(2)));
    if ((Xestimate_lti(5,i-1)>-10^(-2) && Xestimate_lti(5,i-1)<10^(-2)))
        input_now1 = -K(:,:,i-1)*Xestimate_lti(:,i-1);
        input_now2 = 0;
        input_now_lti = [input_now1;input_now2];
    else
        input_now1 = -K(:,:,i-1)*Xestimate_lti(:,i-1);
        input_now2 = g;
        input_now_lti = [input_now1;input_now2];
    end
    if(i==2)
        Ydist_lti(:,:,i-1) = C_s*Xestimate_lti(:,i-1);
        dist = 0.02.*transpose(mvnrnd([0;0],cov_OD_lti));
        Ydist_lti(:,:,i-1) = Ydist_lti(:,:,i-1)+(dist);
    end
    Xtilda_lti = (Xtp1*Xtp2)+(Xtp3*input_now_lti);
    Ydist_lti(:,:,i) = C_s*Xtilda_lti;
    dist = 0.17.*transpose(mvnrnd([0;0],cov_OD_lti));
    Ydist_lti(:,:,i) = Ydist_lti(:,:,i)+(dist);
    diff_lti = Ydist_lti(:,:,i)-(C_s*Xtilda_lti);
    Xestimate_lti(:,i) = Xtilda_lti+(K_s(:,:,i)*diff_lti);
end

% figure(5)
% plot(t_horizon,Xestimate_lti(5,:))
% hold on
% plot(t_horizon, Xd_lin(5,:))
%------------------------------------------------------------------------------

%part 5:
%5.1:
lqr_time = length(t_horizon);
tf_cnl = 1000;
t_cnl = linspace(0,tf_cnl,1000.*lqr_time);
X_cnl_perf = [0;0;0;0;40;0];
u_cnl_perf = [g/2;g/2;g];
X_cnl_dist = zeros(6,length(t_cnl));
cov_SD_cnl = diag(diag(rand(6)));
X_clti_stab = zeros(6,lqr_time);

Xprev = zeros(6,1);
for i=1:lqr_time:length(t_cnl)
    dist = 4.*transpose(mvnrnd([0;0;0;0;0;0],cov_SD_cnl));
    if(i==1)
        X_cnl_dist(:,1) = X_cnl_perf+dist;
    else
        X_cnl_dist(:,i) = Xprev+dist;
    
    end
    X_clti_stab(:,1) = X_cnl_perf - X_cnl_dist(:,i);
    Ydist_lti = zeros(2,1,lqr_time);
    X_stab_estimate = zeros(6,lqr_time);
    X_stab_estimate(:,1) = X_clti_stab(:,1);
    X_stab_tilda = zeros(6,1);
    
    for j=2:1:lqr_time
        if(j==2)
            cov_OD_lti = diag(diag(rand(2)));
            Ydist_lti(:,:,j-1) = C_s*X_stab_estimate(:,j-1);
            dist_est = 0.0007.*transpose(mvnrnd([0;0],cov_OD_lti));
            Ydist_lti(:,:,j-1) = Ydist_lti(:,:,j-1)+(dist_est);
        end
            
        if ((X_stab_estimate(5,j-1)>-10^(-2) && X_stab_estimate(5,j-1)<10^(-2)))
            input_now1 = -K(:,:,j-1)*X_stab_estimate(:,j-1);
            input_now2 = 0;
            input_stab = [input_now1;input_now2];
        else
            input_now1 = -K(:,:,j-1)*X_stab_estimate(:,j-1);
            input_now2 = g;
            input_stab = [input_now1;input_now2];
        end
        Xtp1 = (A_s);
        Xtp2 = X_stab_estimate(:,j-1);
        Xtp3 = B_s;
        X_stab_tilda = (Xtp1*Xtp2)+(Xtp3*input_stab);
        cov_OD_lti = diag(diag(rand(2)));
        Ydist_lti(:,:,j) = C_s*X_stab_tilda;
        dist_est = 0.007.*transpose(mvnrnd([0;0],cov_OD_lti));
        Ydist_lti(:,:,j) = Ydist_lti(:,:,j)+(dist_est);        
        diff_lti = Ydist_lti(:,:,j)-(C_s*X_stab_tilda);
        X_stab_estimate(:,j) = X_stab_tilda+(K_s(:,:,j)*diff_lti);
        X_cnl_dist(:,i+j-1) = X_cnl_perf+X_stab_estimate(:,j);
    end
         
    Xprev = X_cnl_dist(:,i+(lqr_time-1));
end
plot(t_cnl,X_cnl_dist(5,:));


%-------------------Functions:---------------------------------------------------
function A = dynamics_mat(u1,u2,qt,m,b)
A1 = [0 1 0 0 0 0];
A2 = [0 0 0 0 0 0];
A3 = [0 0 0 1 0 0];
A4 = [((u1+u2)*cos(qt))/m 0 0 -b/m 0 0];
A5 = [0 0 0 0 0 1];
A6 = [-(u1+u2)*sin(qt)/m 0 0 0 0 -b/m];
A = [A1;A2;A3;A4;A5;A6];
end

function B = input_mat(I,qt,m)
B1 = [0 0 0];
B2 = [inv(I) -inv(I) 0];
B3 = [0 0 0];
B4 = [sin(qt)/m sin(qt)/m 0];
B5 = [0 0 0];
B6 = [cos(qt)/m cos(qt/m) -1/m];
B = [B1;B2;B3;B4;B5;B6];
end

function C = output_mat(qv,qh)
q = (qv^2+qh^2);
C1 = [0 0 qv^2/q^(3/2) 0 -qh*qv/(q)^(3/2) 0];
C2 = [0 0 -qh*qv/(q)^(3/2) 0 qh^2/q^(3/2) 0];
C = [C1;C2];
end

function p = traj_mat(p2,tf)
r1 = [tf^3 tf^4 tf^5];
r2 = [3*tf^2 4*tf^3 5*tf^4];
r3 = [6*tf 12*tf^2 20*tf^3];
rhs = [r1;r2;r3];
p1 = 0;
p0 = 0;
lhs = [40-(p2*tf^2);-2*p2*tf;-2*p2];
sol = inv(rhs)*lhs;
p = [p0;p1;p2;sol(1);sol(2);sol(3)];
end

function qv_vals = qv_traj(p,t)
qv = p(1)+(p(2)*t)+(p(3)*t^2)+(p(4)*t^3)+(p(5)*t^4)+(p(6)*t^5);
qv_dot = p(2)+(2*p(3)*t)+(3*p(4)*t^2)+(4*p(5)*t^3)+(5*p(6)*t^4);
qv_double_dot = (2*p(3))+(6*p(4)*t)+(12*p(5)*t^2)+(20*p(6)*t^3);
qv_vals = [qv;qv_dot;qv_double_dot];
end
