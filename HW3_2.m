clc; clear; close all;
rng('default') % seed RNG

%% Part D

% Initialization
mu = 3.986004418*10^5; % earth gravitational parameter (km^3/s^2)
x0 = [6753.137 0 0 0.00113975]'; % initial state [r, rdot, theta, thetadot]
state_dim = length(x0); % store number of states

a = -mu/(2*((x0(1)^2*x0(4)^2 + x0(2)^2)/2 - mu/x0(1))); % semi-major axis (km)
T = 2*pi*sqrt(a^3/mu); % orbit period (s)

dt = 30; % measurement time step (s)
t_meas = 0:dt:T; % measurement times (s)
t_noise = linspace(0,T,10^4); % noise time vector (s)

G = [0 0; 1 0; 0 0; 0 1]; % noise state mapping matrix

std_r = 10^-4; % radius process noise standard deviation
std_theta = 10^-8; % theta process noise standard deviation
std_meas = 10^-2; % measurement noise standard deviation

w_r = std_r*randn(1,length(t_noise)); % radius process noise
w_theta = std_theta*randn(1,length(t_noise)); % theta process noise
w = [w_r; w_theta]; % process noise matrix
R = std_meas^2; % measurement noise covariance
R_chol = chol(R,'lower'); % measurement noise covariance cholesky factor

H_k = [0 0 1 0]; % state to measurement transformation
z_k = NaN(1,length(t_meas));

% Calculations
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[tx, x] = ode45(@ODEC, t_meas, x0, opts, w, mu, G, t_noise);
x = x.'; % transpose state matrix
num_time_steps = length(tx.'); % number of time steps

for i = 2:num_time_steps
    z_k(:,i) = H_k*x(:,i) + R_chol*randn;
end

%% Part E

% Initialization
xkm_hist = NaN(state_dim,num_time_steps); % initialize priori state estimate
xkp_hist = NaN(state_dim,num_time_steps); % initialize posteriori state estimate

xkm_hist(:,1) = [6.75681367*10^3; -1.37015864*10^-1; 4.45127583*10^-2; 1.08593276*10^-3]; % initial state estimate
xkp_hist(:,1) = [6.75681367*10^3; -1.37015864*10^-1; 4.45127583*10^-2; 1.08593276*10^-3]; % initial state estimate

Pkm_hist = NaN(state_dim,state_dim,num_time_steps); % initialize priori covariance matrix
Pkp_hist = NaN(state_dim,state_dim,num_time_steps); % initialize posteriori covariance matrix

Pkm_hist(:,:,1) = diag([40 10^-2 10^-3 10^-8]); % initial covariance matrix
Pkp_hist(:,:,1) = diag([40 10^-2 10^-3 10^-8]); % initial covariance matrix

% initialize filter recursion
xkm1 = xkm_hist(:,1);
Pkm1 = Pkm_hist(:,:,1);

tkm1 = 0; % initialize integration time (s)
Qs = diag([std_r^2, std_theta^2]); % process noise covariance

for i = 2:num_time_steps
    % get current time and measurement
    tk = t_meas(i); % integration final time (s)
    tspan = [tkm1 tk]; % integration time interval (s)
    zk = z_k(:,i);
    
    % predict
    [~, xhat] = ode45(@ODED, tspan, [xkm1; reshape(Pkm1,16,1)], opts, mu, G, Qs); % propogate state estimate
    xkm = xhat(end,1:state_dim).'; % propogated states
    Pkm = reshape(xhat(end,5:end),state_dim,state_dim); % propagated covariance

    % update
    zhat_k = H_k*xkm; % update measurement estimate
    W_k = H_k*Pkm*H_k.' + R; % innovation covariance gain
    C_k = Pkm*H_k.'; % cross covariance gain
    K_k = C_k/W_k; % kalman gain

    xkp = xkm + K_k*(zk - zhat_k); % update state estimate
    Pkp = Pkm - C_k*K_k.' - K_k*C_k.' + K_k*W_k*K_k.'; % update covariance
    Pkp = (Pkp + Pkp')/2; % enforce symmetry

    % store values
    xkm_hist(:,i) = xkm;
    xkp_hist(:,i) = xkp;
    Pkm_hist(:,:,i) = Pkm;
    Pkp_hist(:,:,i) = Pkp;

    % recurse
    xkm1 = xkp;
    Pkm1 = Pkp;
    tkm1 = tk;
end

x_pos = x(1,:).*cos(x(3,:)); % true x position
y_pos = x(1,:).*sin(x(3,:)); % true y position

xhat_pos = xkp_hist(1,:).*cos(xkp_hist(3,:)); % estimated x position (km)
yhat_pos = xkp_hist(1,:).*sin(xkp_hist(3,:)); % estimated y position (km)

sigp(1, :) = squeeze(sqrt(Pkp_hist(1,1,:)))'; % r posteriori standard deviation bound
sigp(2, :) = squeeze(sqrt(Pkp_hist(2,2,:)))'; % r dot posteriori standard deviation bound
sigp(3, :) = squeeze(sqrt(Pkp_hist(3,3,:)))'; % theta posteriori standard deviation bound
sigp(4, :) = squeeze(sqrt(Pkp_hist(4,4,:)))'; % theta dot posteriori standard deviation bound

sigm(1, :) = squeeze(sqrt(Pkm_hist(1,1,:)))'; % r priori standard deviation bound
sigm(2, :) = squeeze(sqrt(Pkm_hist(2,2,:)))'; % r dot priori standard deviation bound
sigm(3, :) = squeeze(sqrt(Pkm_hist(3,3,:)))'; % theta priori standard deviation bound
sigm(4, :) = squeeze(sqrt(Pkm_hist(4,4,:)))'; % theta priori dot standard deviation bound

% initialize combined arrays for plotting
time_comb = reshape([t_meas;t_meas],[],1);
sig_comb = NaN(state_dim,2*num_time_steps);
error_comb = NaN(state_dim,2*num_time_steps);

for i = 1:num_time_steps
    % standard deviation combined array
    sig_comb(:,2*i - 1) = sigm(:,i);
    sig_comb(:,2*i) = sigp(:,i);

    % state estimation error combined array
    error_comb(:,2*i - 1) = x(:,i) - xkm_hist(:,i);
    error_comb(:,2*i) = x(:,i) - xkp_hist(:,i);
end

% Plotted Figures and Displays
figure
plot(x_pos,y_pos,xhat_pos,yhat_pos)
title('True and Estimated Position - Sean Bohne')
xlabel('X Position (km)')
ylabel('Y Position (km)')
legend('True Position','Estimated Position','Location','best')
axis equal
grid on

figure
plot(t_meas,x(1,:) - xkp_hist(1,:),t_meas,3*sigp(1,:),'r', t_meas,-3*sigp(1,:),'r')
title('Radius Estimation Error - Sean Bohne')
xlabel('Time (s)')
ylabel('Radial Position Error (km)')
legend('Radius Error','\pm3\sigma')
grid on

figure
plot(t_meas,x(2,:) - xkp_hist(2,:),t_meas, 3*sigp(2,:),'r',t_meas,-3*sigp(2,:),'r')
title('Velocity Estimation Error - Sean Bohne')
xlabel('Time (s)')
ylabel('Velocity Error (km/s)')
legend('Velocity Error','\pm3\sigma')
grid on

figure
plot(t_meas,x(3,:) - xkp_hist(3,:),t_meas, 3*sigp(3,:),'r',t_meas,-3*sigp(3,:),'r')
title('Angular Position Estimation Error - Sean Bohne')
xlabel('Time (s)')
ylabel('Angular Position Error (rad)')
legend('Angular Position Error','\pm3\sigma')
grid on

figure
plot(t_meas,x(4,:) - xkp_hist(4,:),t_meas,3*sigp(4,:),'r', t_meas,-3*sigp(4,:),'r')
title('Angular Velocity Estimation Error - Sean Bohne')
xlabel('Time (s)')
ylabel('Angular Velocity Error (rad/s)')
legend('Angular Velocity Error','\pm3\sigma')
grid on

figure
plot(time_comb,error_comb(1,:),time_comb,3*sig_comb(1,:),'r',time_comb,-3*sig_comb(1,:),'r')
title('Radius Priori and Posterior Estimation Error - Sean Bohne','FontSize',8)
xlabel('Time (s)')
ylabel('Radial Position Error (km)')
legend('Radius Error','\pm3\sigma')
grid on

figure
plot(time_comb,error_comb(2,:),time_comb,3*sig_comb(2,:),'r',time_comb,-3*sig_comb(2,:),'r')
title('Velocity Prior and Posterior Estimation Error - Sean Bohne','FontSize',8)
xlabel('Time (s)')
ylabel('Velocity Error (km/s)')
legend('Velocity Error','\pm3\sigma')
grid on

figure
plot(time_comb,error_comb(3,:),time_comb,3*sig_comb(3,:),'r',time_comb,-3*sig_comb(3,:),'r')
title('Angular Position Prior and Posterior Estimation Error - Sean Bohne','FontSize',8)
xlabel('Time (s)')
ylabel('Angular Position Error (rad)')
legend('Angular Position Error','\pm3\sigma')
grid on

figure
plot(time_comb,error_comb(4,:),time_comb,3*sig_comb(4,:),'r',time_comb,-3*sig_comb(4,:),'r')
title('Angular Velocity Prior and Posterior Estimation Error - Sean Bohne','FontSize',8)
xlabel('Time (s)')
ylabel('Angular Velocity Error (rad/s)')
legend('Angular Velocity Error','\pm3\sigma')
grid on

%% Function Definitions
function dxdt = ODEC(t, x, w, mu, G, t_noise)
idx = find(t_noise <= t, 1, 'last');
% store states to meaningful variables
r = x(1);
rdot = x(2);
theta = x(3);
thetadot = x(4);

dxdt = [rdot; r*thetadot^2 - mu/r^2; thetadot; -2*thetadot*rdot/r] + G*w(:,idx);
end

function J = Jacobian(x,mu)
x = x(:);
% store states to meaningful variables
r = x(1);
rdot = x(2);
theta = x(3);
thetadot = x(4);

J = [0, 1, 0, 0; thetadot^2 + 2*mu/r^3, 0, 0, 2*r*thetadot; 0,... 
    0, 0, 1; 2*thetadot*rdot/r^2, -2*thetadot/r, 0, -2*rdot/r];
end

function dPdt = ODED(~, x, mu, G, Qs)
% store states to meaningful variables
r = x(1);
rdot = x(2);
theta = x(3);
thetadot = x(4);

dxdt = [rdot; r*thetadot^2 - mu/r^2; thetadot; -2*thetadot*rdot/r];

F = Jacobian(x(1:4), mu);
P = reshape(x(5:end), 4, 4);
Pdot = F*P + P*F' + G*Qs*G';
dPdt = [dxdt;reshape(Pdot,16,1)];
end