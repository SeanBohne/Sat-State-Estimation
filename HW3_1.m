clc; clear all; close all;
%% Part A

% Initialization
seed = 1; % initialize seed
rng(seed) % seed RNG
month = 36; % simulation time (months)
t = 0:month; % time vector (months)
deorbit = 0.1; % monthly satellite deorbit rate (%)
sat_cost = 1; % cost per satellite (millions of dollars)
std_b = 1.5; % monthly budget standard deviation (millions of dollars)
std_ep = 7; % satellite # error standard deviation
gamma_k = [0;1]; % process noise mapping matrix
x = NaN(2,month + 1); % initialize state vector
p_k = NaN(2,month + 1); % initialize satellite # vector
b_k = NaN(2,month + 1); % initialize budget vector
x(:,1) = [200; 50]; % initial state: [# satellites; budget]
phi_k = [(1 - deorbit), 1/sat_cost; 0, 1]; % state transition matrix

% Calculations
for i = 2:length(x)
    w_k = std_b*randn;
    x(:,i) = phi_k*x(:,i - 1) + gamma_k*w_k;
end

figure(1)
plot(t,x(1,:))
title('Number of Satellites Over 36 Months - Sean Bohne')
xlabel('Time (Months)')
ylabel('Number of Satellites')
grid on

figure(2)
plot(t,x(2,:))
title('Satellite Budget Over 36 Months - Sean Bohne')
xlabel('Time (Months)')
ylabel('Monthly Budget (Millions of Dollars)')
grid on

%% Part C

% Initialization
xhat_k = NaN(2,month + 1); % initialize state esimate
xhat_k(:,1) = [170; 80]; % initial state estimate
std_pe = 50; % satellite number error standard deviation
std_be = 20; % budget error standard deviation (millions of dollars)
std_z = 7; % measurement noise standard deviation
R_k = std_z^2; % measurement noise covariance
P_k = NaN(2,2,month + 1); % initial error covariance
P_k(:,:,1) = diag([std_pe^2, std_be^2]); % initial error covariance matrix
H_k = [1 0]; % state to measurement mapping
z_k = H_k*x + std_z*randn(1,month + 1); % satellite # noisy measurements
Q_k = std_b^2; % process noise covariance

for i = 2:length(x)
    xhat_k(:,i) = phi_k*xhat_k(:,i - 1); % propogate state estimate
    P_k(:,:,i) = phi_k*P_k(:,:,i - 1)*phi_k' + gamma_k*Q_k*gamma_k'; % propogate covariance

    W_k = H_k*P_k(:,:,i)*H_k' + R_k; % innovation covariance gain
    C_k = P_k(:,:,i)*H_k'; % cross covariance gain
    K_k = C_k/W_k; % kalman gain

    zhat_k = H_k*xhat_k(:,i); % update measurement estimate
    xhat_k(:,i) = xhat_k(:,i) + K_k*(z_k(:,i) - zhat_k); % update state estimate
    P_k(:,:,i) = P_k(:,:,i) - C_k*K_k' - K_k*C_k' + K_k*W_k*K_k'; % update covariance
end

sig_bound(1, :) = squeeze(sqrt(P_k(1,1,:)))'; % constellation standard deviation bound
sig_bound(2, :) = squeeze(sqrt(P_k(2,2,:)))'; % budget standard deviation bound

figure(3)
plot(t,xhat_k(1,:),t,x(1,:))
title('True and Estimated Constellation Size - Sean Bohne')
xlabel('Time (Months)')
ylabel('Number of Satellites')
legend('Estimated Satellites','True Satellites','Location','best')
grid on

figure(4)
plot(t,xhat_k(2,:),t,x(2,:))
title('True and Estimated Budget - Sean Bohne')
xlabel('Time (Months)')
ylabel('Monthly Budget (Millions of Dollars)')
legend('Estimated Budget','True Satellites','Location','best')
grid on

figure(5)
plot(t,x(1,:)-xhat_k(1,:),t,sig_bound(1,:),'r--',t,-sig_bound(1,:),'r--')
title('Satellite Estimation Error with \pm1\sigma Bounds - Sean Bohne')
xlabel('Time (Months)')
ylabel('Satellite Number Estimation Error')
legend('Estimation Error', '\pm1\sigma Upper Bound', 'Location', 'best')
grid on

figure(6)
plot(t,x(2,:)-xhat_k(2,:),t,sig_bound(2,:),'r--',t,-sig_bound(2,:),'r--')
title('Budget Error with \pm1\sigma Bounds - Sean Bohne')
xlabel('Time (Months)')
ylabel('Budget Estimation Error (Millions of Dollars)')
legend('Estimation Error', '\pm1\sigma Upper Bound', 'Location', 'best')
grid on