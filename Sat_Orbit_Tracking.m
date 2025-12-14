clc; clear; close all;
rng(1)
% setting meas_model to 1 uses optical measurements (angles only) and 2
% uses radar (range and range-rate). If set to 3, the script runs with
% optical and radar measurement models.
meas_model = 2;

% select number of Monte Carlo trials
MC_trials = 100;

%% Initialization
% Initialize keplerian orbital elements
param.mu_Earth = 3.986*10^5; % Earth gravitational parameter (km^3/s^2)
r_Earth = 6378; % radius of Earth (km)

h_alt = 1000; % satellite orbit altitude at periapsis (km)
a =  r_Earth + h_alt; % semi major axis (km)
e = 0.01; % eccentricity
incl = 45*pi/180; % inclination (rad)
right_ascen = 0*pi/180; % longitude of the ascending node (rad)
arg_perigee = 0*pi/180; % argument of perigee (rad)
true_anom = 0*pi/180; % true anomoly (rad)

r_sat = a*(1 - e^2)/(1 + e*cos(true_anom)); % satellite initial radial position from center of Earth (km)
r_sat_pf = [r_sat*cos(true_anom); r_sat*sin(true_anom); 0]; % satellite perifocal position (km)

% transform satellite perifocal position to ECI (km)
r_sat_ECI = DCM3(right_ascen)*DCM1(incl)*DCM3(arg_perigee)*r_sat_pf;

% determine satellite perifocal velocity (km/s)
h = sqrt(param.mu_Earth*a*(1 - e^2)); % specific angular momentum (km^2/s)
rdot_sat = h*e*sin(true_anom)/(a*(1 - e^2)); % radial velocity (km/s)
thetadot_sat = h/r_sat; % angular velocity (km/s)
v_sat = [rdot_sat*cos(true_anom) - thetadot_sat*sin(true_anom); rdot_sat*sin(true_anom) + thetadot_sat*cos(true_anom); 0]; % satellite velocity in perifocal frame (km/s)

% transform satellite perifocal velocity to ECI (km)
v_sat_ECI = DCM3(right_ascen)*DCM1(incl)*DCM3(arg_perigee)*v_sat;

% initial satellite states
x_sat = r_sat_ECI(1); % satellite initial x position (km)
y_sat = r_sat_ECI(2); % satellite initial y position (km)
z_sat = r_sat_ECI(3); % satellite initial z position (km)

vx_sat = v_sat_ECI(1); % satellite initial x velocity (km/s)
vy_sat = v_sat_ECI(2); % satellite initial x velocity (km/s)
vz_sat = v_sat_ECI(3); % satellite initial x velocity (km/s)

num_orbits = 5; % number of orbit period simulations
v0 = sqrt(2*param.mu_Earth/r_sat - param.mu_Earth/a); % satellite initial velocity (km/s)
T = num_orbits*2*pi*sqrt(a^3/param.mu_Earth); % orbit period (s)
x0_sat = [x_sat; vx_sat; y_sat; vy_sat; z_sat; vz_sat]; % satellite initial state in ECI [x; xdot; y; ydot; z; zdot]

% initialize ground station location
r_gs = r_Earth; % ground station distance from center of Earth (km)
lat_gs = 0*pi/180; % ground station latitude position (rad)
lon_gs = 0*pi/180; % ground station longitude position (rad)

% ground station measurement limits
el_lim_low = 0*pi/180; % elevation limit lower bound (rad)
el_lim_high = 90*pi/180; % elevation limit upper bound (rad)

% ground station position in ECEF coordinates
r_gs_ECEF = [r_gs*cos(lat_gs)*cos(lon_gs);
                 r_gs*cos(lat_gs)*sin(lon_gs);
                 r_gs*sin(lat_gs)];

dt = 10; % measurement time step (s)
t_meas = 0:dt:T; % measurement times (s)
t_comb = reshape([t_meas/60/60/24; t_meas/60/60/24],[],1); % plot time in days (days) for combined array
t_ASNEES = t_meas/60/60/24; % plot time in days for ASNEES (days)
num_time_steps = length(t_meas); % number of measurement time steps

n_sat = length(x0_sat); % state dimension

% determine measurement dimension
if meas_model == 1
    n_meas = 2;
elseif meas_model == 2
    n_meas = 2;
else
    n_meas = 4;
end

% process noise parameters
sig_vx = 1.5e-6; % x velocity process noise standard deviation (km/s)
sig_vy = 1e-6; % y velocity process noise standard deviation (km/s)
sig_vz = 1e-6; % z velocity process noise standard deviation (km/s)
G = [0 0 0; 1 0 0; 0 0 0; 0 1 0; 0 0 0; 0 0 1]; % process noise mapping matrix
Q = dt*G*diag([sig_vx^2, sig_vy^2, sig_vz^2])*G';

% measurement noise parameters
sig_range = 0.05; % range measurement noise standard deviation (km)
sig_range_rate = 5e-5; % range-rate measurement noise standard deviation (km/s)
sig_el = 10*pi/(180*3600); % elevation measurement noise standard deviation (rad)
sig_az = 10*pi/(180*3600); % azimuth measurement noise standard deviation (rad)

% choose measurement noise covariance based on measurement model
if meas_model == 1
    R = diag([sig_el^2, sig_az^2]);
elseif meas_model == 2
    R = diag([sig_range^2, sig_range_rate^2]);
else
    R = diag([sig_range^2, sig_range_rate^2, sig_el^2, sig_az^2]);
end

% satellite physical parameters based on 6U CubeSat
A = 0.085; % 1U cubesat frontal area (m^2)
m = 10; % 1U cubesat mass (kg)
param.area_mass = A/m; % 1U cubesat area to mass ratio for drag and radiation forces

% drag and solar pressure force parameters
param.rho_ref = 1.225; % reference density (kg/m^3)
param.dens_scale = 1/7.8; % km density scaling parameter
param.h0 = r_Earth; % reference altitude (km)
param.CD = 2.13; % 1U cubesat drag coefficient
param.P = 4.56*10^-6; % momentum flux from the sun from Tapley
param.nu = 1; % eclipse factor (assume satellite is always exposed to sun)
param.CR = 1; % reflectivity coefficient approximately 1 according to Tapley

%% True Orbit and Measurement Simulation
% simulate dynamics
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);

% preallocate measurement variable memory
z_k_range = NaN(1, num_time_steps);
z_k_range_rate = NaN(1, num_time_steps);
z_k_el = NaN(1, num_time_steps);
z_k_az = NaN(1, num_time_steps);

% preallocate measurement to ECI conversion
x_meas = NaN(1, num_time_steps);
y_meas = NaN(1, num_time_steps);
z_meas = NaN(1, num_time_steps);
r_rel_ECEF = NaN(3, num_time_steps);

w_Earth = 7.292115*10^-5; % Earth angular velocity (rad/s)
T_t = ECEF2Topo(lat_gs, lon_gs); % ground station position in topocentric frame

z = NaN(n_meas, num_time_steps); % preallocate combined measurements
Z_UKF = NaN(n_meas, num_time_steps); % preallocate UKF measurements
Xkm = NaN(n_sat, 2*n_sat + 1); % preallocate sigma point propogation

%% EKF and UKF loop
% initialize SUT weighting parameters
alpha = 1e-3;
beta = 2;
kappa = 0;

% parameters to determine output vectors in ODE_Low_Fidelity function
% EKF = 1 and UKF = 2
EKF = 1;
UKF = 2;

x_high_fidel = NaN(n_sat, num_time_steps,MC_trials); % preallocate simulated dynamics

P0 = diag([0.05, 10^-7, 0.02, 10^-8, 0.02, 10^-8]); % initial covariance matrix
Pkm_hist_UKF = NaN(n_sat,n_sat,num_time_steps,MC_trials); % initialize priori covariance matrix
Pkp_hist_UKF = NaN(n_sat,n_sat,num_time_steps,MC_trials); % initialize posterior covariance matrix
Pkm_hist_EKF = NaN(n_sat,n_sat,num_time_steps,MC_trials); % initialize priori covariance matrix
Pkp_hist_EKF = NaN(n_sat,n_sat,num_time_steps,MC_trials); % initialize posterior covariance matrix

% initialize UKF a priori and posterior states and covariances
xkm_hist_UKF = NaN(n_sat,num_time_steps,MC_trials); % initialize priori state estimate
xkp_hist_UKF = NaN(n_sat,num_time_steps,MC_trials); % initialize posterior state estimate
xkm_hist_EKF = NaN(n_sat,num_time_steps,MC_trials); % initialize priori state estimate
xkp_hist_EKF = NaN(n_sat,num_time_steps,MC_trials); % initialize posterior state estimate

% initialize runtime analysis variables
EKF_tot_time = 0;
UKF_tot_time = 0;

% initialize innovations variable
inno = NaN(n_meas,num_time_steps,MC_trials);

% start monte carlo trials
for i = 1:MC_trials
    Pkm_hist_EKF(:,:,1,i) = P0; % initial covariance matrix
    Pkp_hist_EKF(:,:,1,i) = P0; % initial covariance matrix
    Pkm_hist_UKF(:,:,1,i) = P0; % initial covariance matrix
    Pkp_hist_UKF(:,:,1,i) = P0; % initial covariance matrix

    sig_initial = chol(Pkm_hist_UKF(:,:,1,i), 'lower'); % initial state standard deviation
    
    x0 = x0_sat + sig_initial*randn(n_sat,1);
    fprintf('Trial %d: Initial state error = [%.3f, %.3f, %.3f] km\n', ...
        i, x0(1)-x0_sat(1), x0(3)-x0_sat(3), x0(5)-x0_sat(5));
    xkm_hist_EKF(:,1,i) = x0; % initial state estimate
    xkp_hist_EKF(:,1,i) = x0; % initial state estimate
    xkm_hist_UKF(:,1,i) = x0; % initial state estimate
    xkp_hist_UKF(:,1,i) = x0; % initial state estimate

    % initialize filter recursion
    xkm1_EKF = xkm_hist_EKF(:,1,i);
    Pkm1_EKF = Pkm_hist_EKF(:,:,1,i);
    tkm1_EKF = 0;
    xkm1_UKF = xkm_hist_UKF(:,1,i);
    Pkm1_UKF = Pkm_hist_UKF(:,:,1,i);
    tkm1_UKF = 0;

    % sample measurement and process noise
    w_vx_i = sig_vx*randn(1,num_time_steps);
    w_vy_i = sig_vy*randn(1,num_time_steps);
    w_vz_i = sig_vz*randn(1,num_time_steps);
    
    v_range_i = sig_range*randn(1,num_time_steps);
    v_range_rate_i = sig_range_rate*randn(1,num_time_steps);
    v_el_i = sig_el*randn(1,num_time_steps);
    v_az_i = sig_az*randn(1,num_time_steps);
    
    % combine process noise and measurement noise
    w = [w_vx_i; w_vy_i; w_vz_i];
    v = [v_range_i; v_range_rate_i; v_el_i; v_az_i];

    x_high_fidel(:,1,i) = x0; % initial high fidelity state

    % simulate true dynamics for current monte carlo trial
    [t, x_intermed] = ode45(@ODE_High_Fidelity, t_meas, x_high_fidel(:,1,i), opts, param, w_Earth);
    x_high_fidel(:,:,i) = x_intermed.';

    for j = 1:num_time_steps
        z(:,j) = meas_func(x_high_fidel(:,j,i), r_gs_ECEF, t_meas(j), T_t, w_Earth, meas_model, el_lim_low, el_lim_high, v(:,j));
    end

    % run EKF loop for current monte carlo trial
    tic
    for k = 2:num_time_steps        
        % check for positive definiteness
        if any(eig(Pkm1_EKF) < 1e-12) || ~isequal(Pkm1_EKF, Pkm1_EKF')
            fprintf(['Warning: Covariance matrix not positive ' ...
                     'definite! Stopping the program.\n'])
            fprintf('Failure after %i iterations.\n',k - 1)
            break
        end

        % get current time and measurement
        tk = t_meas(k);
        tspan = [tkm1_EKF tk];
        zk = z(:,k);

        % propogate step
        [~, xhat] = ode45(@ODE_Low_Fidelity, tspan, [xkm1_EKF; reshape(Pkm1_EKF,n_sat^2,1)], opts, param.mu_Earth, Q, n_sat, EKF);
        xkm_EKF = xhat(end,1:n_sat)';
        Pkm_EKF = reshape(xhat(end,n_sat + 1:end),n_sat,n_sat);

        % skip update step if measurement unavailable
        if any(isnan(zk))
            % store values
            xkm_hist_EKF(:,k,i) = xkm_EKF;
            xkp_hist_EKF(:,k,i) = xkm_EKF;
            Pkm_hist_EKF(:,:,k,i) = Pkm_EKF;
            Pkp_hist_EKF(:,:,k,i) = Pkm_EKF;

            % Update recursion variables
            xkm1_EKF = xkm_EKF;
            Pkm1_EKF = (Pkm_EKF + Pkm_EKF')/2;
            tkm1_EKF = tk;
            
            continue
        end

        % update step
        H_k = Meas_Jacobian(xkm_EKF, r_gs_ECEF, tk, w_Earth, T_t, el_lim_low, el_lim_high, meas_model);
        zhat_k = meas_func(xkm_EKF, r_gs_ECEF, tk, T_t, w_Earth, meas_model, el_lim_low, el_lim_high);
        W_k = H_k*Pkm_EKF*H_k.' + R; % innovation covariance gain
        C_k = Pkm_EKF*H_k.'; % cross covariance gain
        K_k = C_k/W_k; % kalman gain

        xkp_EKF = xkm_EKF + K_k*(zk - zhat_k); % update state estimate
        Pkp_EKF = Pkm_EKF - C_k*K_k.' - K_k*C_k.' + K_k*W_k*K_k.'; % update covariance
        Pkp_EKF = (Pkp_EKF + Pkp_EKF')/2; % enforce symmetry

        % store values
        xkm_hist_EKF(:,k,i) = xkm_EKF;
        xkp_hist_EKF(:,k,i) = xkp_EKF;
        Pkm_hist_EKF(:,:,k,i) = Pkm_EKF;
        Pkp_hist_EKF(:,:,k,i) = Pkp_EKF;
    
        % recurse
        xkm1_EKF = xkp_EKF;
        Pkm1_EKF = Pkp_EKF;
        tkm1_EKF = tk;
    end

    EKF_time = toc;

    % check if measurement is available
    if ~any(isnan(zk))
        % display runtimes for each trial and determine running total
        EKF_tot_time = EKF_tot_time + EKF_time;
        disp(['EKF elapsed time is ' num2str(EKF_time) ' seconds for Monte Carlo trial ' num2str(i)]);
    end

    % run UKF loop for current monte carlo trial
    tic
    for k = 2:num_time_steps
        % check for positive definiteness
        if any(eig(Pkm1_UKF) < 1e-12) || ~isequal(Pkm1_UKF, Pkm1_UKF')
            fprintf(['Warning: Covariance matrix not positive ' ...
                     'definite! Stopping the program.\n'])
            fprintf('Failure after %i iterations.\n',k - 1)
            break
        end
        
        % get current time and measurement
        tk = t_meas(k);
        tspan = [tkm1_UKF tk];
        zk = z(:,k);
    
        % form sigma points and calculate weights
        [sig_points, wm, wc] = SigPoints(Pkm1_UKF, xkm1_UKF, alpha, beta, kappa, n_sat);
    
        % propogate sigma points
        for j = 1:(2*n_sat + 1)
            [~, Xhat] = ode45(@ODE_Low_Fidelity, tspan, sig_points(:,j), opts, param.mu_Earth, Q, n_sat, UKF);
            Xkm(:,j) = Xhat(end,:).';
        end
        
        % calculate SUT propogated mean estimate
        xkm_UKF = wm(1)*Xkm(:,1);
    
        for j = 2:(2*n_sat + 1)
            xkm_UKF = xkm_UKF + wm(j)*Xkm(:,j);
        end
    
        % calculate SUT propogated covariance estimate
        Pkm_UKF = wc(1)*(Xkm(:,1) - xkm_UKF)*(Xkm(:,1) - xkm_UKF)' + Q;
        
        for j = 2:(2*n_sat + 1)
            Pkm_UKF = Pkm_UKF + wc(j)*(Xkm(:,j) - xkm_UKF)*(Xkm(:,j) - xkm_UKF)';
        end
    
        % skip update step if measurement unavailable
        if any(isnan(zk))
            % store values
            xkm_hist_UKF(:,k,i) = xkm_UKF;
            xkp_hist_UKF(:,k,i) = xkm_UKF;
            Pkm_hist_UKF(:,:,k,i) = (Pkm_UKF + Pkm_UKF')/2;
            Pkp_hist_UKF(:,:,k,i) = (Pkm_UKF + Pkm_UKF')/2;

            % Update recursion variables
            xkm1_UKF = xkm_UKF;
            Pkm1_UKF = (Pkm_UKF + Pkm_UKF')/2;
            tkm1_UKF = tk;

            continue
        end
    
        % update sigma points and weights
        [sig_points, wm, wc] = SigPoints(Pkm_UKF, xkm_UKF, alpha, beta, kappa, n_sat);

        for l = 1:(2*n_sat + 1)
            Z_UKF(:,l) = meas_func(sig_points(:,l), r_gs_ECEF, tk, T_t, w_Earth, meas_model, el_lim_low, el_lim_high);
        end
    
        % calculate predicted measurement mean
        zhat_k_UKF = wm(1)*Z_UKF(:,1);
    
        for l = 2:(2*n_sat + 1)
            zhat_k_UKF = zhat_k_UKF + wm(l)*Z_UKF(:,l);
        end
    
        % calculate predicted measurement covariance and cross covariance
        Pzzk_UKF = wc(1)*(Z_UKF(:,1) - zhat_k_UKF)*(Z_UKF(:,1) - zhat_k_UKF)' + R;
        Pxzk_UKF = wc(1)*(sig_points(:,1) - xkm_UKF)*(Z_UKF(:,1) - zhat_k_UKF)';
    
        for l = 2:(2*n_sat + 1)
            Pzzk_UKF = Pzzk_UKF + wc(l)*(Z_UKF(:,l) - zhat_k_UKF)*(Z_UKF(:,l) - zhat_k_UKF)';
            Pxzk_UKF = Pxzk_UKF + wc(l)*(sig_points(:,l) - xkm_UKF)*(Z_UKF(:,l) - zhat_k_UKF)';
        end
    
        % compute Kalman gain and update state mean and covariance
        K_k_UKF = Pxzk_UKF/Pzzk_UKF;
        xkp_UKF = xkm_UKF + K_k_UKF*(zk - zhat_k_UKF);
        Pkp_UKF = Pkm_UKF - Pxzk_UKF*K_k_UKF' - K_k_UKF*Pxzk_UKF' + K_k_UKF*Pzzk_UKF*K_k_UKF';
    
        Pkp_UKF = (Pkp_UKF + Pkp_UKF')/2; % enforce symmetry
    
        % store values
        xkm_hist_UKF(:,k,i) = xkm_UKF;
        xkp_hist_UKF(:,k,i) = xkp_UKF;
        Pkm_hist_UKF(:,:,k,i) = Pkm_UKF;
        Pkp_hist_UKF(:,:,k,i) = Pkp_UKF;
        inno(:,k,i) = zk - zhat_k_UKF;
    
        % recurse
        xkm1_UKF = xkp_UKF;
        Pkm1_UKF = Pkp_UKF;
        tkm1_UKF = tk;
    end

    UKF_time = toc;

    % check if measurement is available
    if ~any(isnan(zk))
        % display runtimes for each trial and determine running total
        UKF_tot_time = UKF_tot_time + UKF_time;
        disp(['UKF elapsed time is ' num2str(UKF_time) ' seconds for Monte Carlo trial ' num2str(i)]);
    end
end

% displayed total time elapsed for each filter
disp('-------------------------------------------------------------------------')
disp(['EKF total run time is ' num2str(EKF_tot_time) ' seconds across ' num2str(MC_trials) ' Monte Carlo trial']);
disp(['UKF total elapsed time is ' num2str(UKF_tot_time) ' seconds across ' num2str(MC_trials) ' Monte Carlo trial']);

% display average time elapsed for each Monte Carlo trial
disp('-------------------------------------------------------------------------')
disp(['EKF average run time is ' num2str(EKF_tot_time/MC_trials) ' seconds across ' num2str(MC_trials) ' Monte Carlo trial']);
disp(['UKF average run time is ' num2str(UKF_tot_time/MC_trials) ' seconds across ' num2str(MC_trials) ' Monte Carlo trial']);
disp('-------------------------------------------------------------------------')

% display how much faster EKF is compared to UKF on average
disp(['EKF is ' num2str(UKF_tot_time/EKF_tot_time) ' times faster on average compared to UKF over ' num2str(MC_trials) ' Monte Carlo trials']);
disp('-------------------------------------------------------------------------')

% initialize EKF and UKF combined error and sigma arrays for plotting
x_pos_error_comb_EKF = NaN(MC_trials,2*num_time_steps);
y_pos_error_comb_EKF = NaN(MC_trials,2*num_time_steps);
z_pos_error_comb_EKF = NaN(MC_trials,2*num_time_steps);
x_vel_error_comb_EKF = NaN(MC_trials,2*num_time_steps);
y_vel_error_comb_EKF = NaN(MC_trials,2*num_time_steps);
z_vel_error_comb_EKF = NaN(MC_trials,2*num_time_steps);

x_pos_sig_comb_EKF = NaN(MC_trials,2*num_time_steps);
y_pos_sig_comb_EKF = NaN(MC_trials,2*num_time_steps);
z_pos_sig_comb_EKF = NaN(MC_trials,2*num_time_steps);
x_vel_sig_comb_EKF = NaN(MC_trials,2*num_time_steps);
y_vel_sig_comb_EKF = NaN(MC_trials,2*num_time_steps);
z_vel_sig_comb_EKF = NaN(MC_trials,2*num_time_steps);

x_pos_error_comb_UKF = NaN(MC_trials,2*num_time_steps);
y_pos_error_comb_UKF = NaN(MC_trials,2*num_time_steps);
z_pos_error_comb_UKF = NaN(MC_trials,2*num_time_steps);
x_vel_error_comb_UKF = NaN(MC_trials,2*num_time_steps);
y_vel_error_comb_UKF = NaN(MC_trials,2*num_time_steps);
z_vel_error_comb_UKF = NaN(MC_trials,2*num_time_steps);

x_pos_sig_comb_UKF = NaN(MC_trials,2*num_time_steps);
y_pos_sig_comb_UKF = NaN(MC_trials,2*num_time_steps);
z_pos_sig_comb_UKF = NaN(MC_trials,2*num_time_steps);
x_vel_sig_comb_UKF = NaN(MC_trials,2*num_time_steps);
y_vel_sig_comb_UKF = NaN(MC_trials,2*num_time_steps);
z_vel_sig_comb_UKF = NaN(MC_trials,2*num_time_steps);

% calculate ASNEES values and plot estimation errors with 3 sigma bounds
for i = 1:MC_trials
    % EKF and UKF position priori and posterior errors
    xkp_pos_error_EKF = squeeze(x_high_fidel(1,:,i) - xkp_hist_EKF(1,:,i));
    ykp_pos_error_EKF = squeeze(x_high_fidel(3,:,i) - xkp_hist_EKF(3,:,i));
    zkp_pos_error_EKF = squeeze(x_high_fidel(5,:,i) - xkp_hist_EKF(5,:,i));
    xkm_pos_error_EKF = squeeze(x_high_fidel(1,:,i) - xkm_hist_EKF(1,:,i));
    ykm_pos_error_EKF = squeeze(x_high_fidel(3,:,i) - xkm_hist_EKF(3,:,i));
    zkm_pos_error_EKF = squeeze(x_high_fidel(5,:,i) - xkm_hist_EKF(5,:,i));

    xkp_pos_error_UKF = squeeze(x_high_fidel(1,:,i) - xkp_hist_UKF(1,:,i));
    ykp_pos_error_UKF = squeeze(x_high_fidel(3,:,i) - xkp_hist_UKF(3,:,i));
    zkp_pos_error_UKF = squeeze(x_high_fidel(5,:,i) - xkp_hist_UKF(5,:,i));
    xkm_pos_error_UKF = squeeze(x_high_fidel(1,:,i) - xkm_hist_UKF(1,:,i));
    ykm_pos_error_UKF = squeeze(x_high_fidel(3,:,i) - xkm_hist_UKF(3,:,i));
    zkm_pos_error_UKF = squeeze(x_high_fidel(5,:,i) - xkm_hist_UKF(5,:,i));
    
    % EKF and UKF velocity priori and posterior errors
    xkp_vel_error_EKF = squeeze(x_high_fidel(2,:,i) - xkp_hist_EKF(2,:,i));
    ykp_vel_error_EKF = squeeze(x_high_fidel(4,:,i) - xkp_hist_EKF(4,:,i));
    zkp_vel_error_EKF = squeeze(x_high_fidel(6,:,i) - xkp_hist_EKF(6,:,i));
    xkm_vel_error_EKF = squeeze(x_high_fidel(2,:,i) - xkm_hist_EKF(2,:,i));
    ykm_vel_error_EKF = squeeze(x_high_fidel(4,:,i) - xkm_hist_EKF(4,:,i));
    zkm_vel_error_EKF = squeeze(x_high_fidel(6,:,i) - xkm_hist_EKF(6,:,i));

    xkp_vel_error_UKF = squeeze(x_high_fidel(2,:,i) - xkp_hist_UKF(2,:,i));
    ykp_vel_error_UKF = squeeze(x_high_fidel(4,:,i) - xkp_hist_UKF(4,:,i));
    zkp_vel_error_UKF = squeeze(x_high_fidel(6,:,i) - xkp_hist_UKF(6,:,i));
    xkm_vel_error_UKF = squeeze(x_high_fidel(2,:,i) - xkm_hist_UKF(2,:,i));
    ykm_vel_error_UKF = squeeze(x_high_fidel(4,:,i) - xkm_hist_UKF(4,:,i));
    zkm_vel_error_UKF = squeeze(x_high_fidel(6,:,i) - xkm_hist_UKF(6,:,i));
    
    % EKF and UKF position priori and posterior standard deviation bounds
    xkp_pos_sig_EKF = squeeze(sqrt(Pkp_hist_EKF(1,1,:,i)))';
    ykp_pos_sig_EKF = squeeze(sqrt(Pkp_hist_EKF(3,3,:,i)))';
    zkp_pos_sig_EKF = squeeze(sqrt(Pkp_hist_EKF(5,5,:,i)))';
    xkm_pos_sig_EKF = squeeze(sqrt(Pkm_hist_EKF(1,1,:,i)))';
    ykm_pos_sig_EKF = squeeze(sqrt(Pkm_hist_EKF(3,3,:,i)))';
    zkm_pos_sig_EKF = squeeze(sqrt(Pkm_hist_EKF(5,5,:,i)))';
    
    xkp_pos_sig_UKF = squeeze(sqrt(Pkp_hist_UKF(1,1,:,i)))';
    ykp_pos_sig_UKF = squeeze(sqrt(Pkp_hist_UKF(3,3,:,i)))';
    zkp_pos_sig_UKF = squeeze(sqrt(Pkp_hist_UKF(5,5,:,i)))';
    xkm_pos_sig_UKF = squeeze(sqrt(Pkm_hist_UKF(1,1,:,i)))';
    ykm_pos_sig_UKF = squeeze(sqrt(Pkm_hist_UKF(3,3,:,i)))';
    zkm_pos_sig_UKF = squeeze(sqrt(Pkm_hist_UKF(5,5,:,i)))';
    
    % EKF and UKF velocity priori and posterior standard deviation bounds
    xkp_vel_sig_EKF = squeeze(sqrt(Pkp_hist_EKF(2,2,:,i)))';
    ykp_vel_sig_EKF = squeeze(sqrt(Pkp_hist_EKF(4,4,:,i)))';
    zkp_vel_sig_EKF = squeeze(sqrt(Pkp_hist_EKF(6,6,:,i)))';
    xkm_vel_sig_EKF = squeeze(sqrt(Pkm_hist_EKF(2,2,:,i)))';
    ykm_vel_sig_EKF = squeeze(sqrt(Pkm_hist_EKF(4,4,:,i)))';
    zkm_vel_sig_EKF = squeeze(sqrt(Pkm_hist_EKF(6,6,:,i)))';
    
    xkp_vel_sig_UKF = squeeze(sqrt(Pkp_hist_UKF(2,2,:,i)))';
    ykp_vel_sig_UKF = squeeze(sqrt(Pkp_hist_UKF(4,4,:,i)))';
    zkp_vel_sig_UKF = squeeze(sqrt(Pkp_hist_UKF(6,6,:,i)))';
    xkm_vel_sig_UKF = squeeze(sqrt(Pkm_hist_UKF(2,2,:,i)))';
    ykm_vel_sig_UKF = squeeze(sqrt(Pkm_hist_UKF(4,4,:,i)))';
    zkm_vel_sig_UKF = squeeze(sqrt(Pkm_hist_UKF(6,6,:,i)))';
    
    for j = 1:num_time_steps
        % standard deviation combined arrays
        x_pos_sig_comb_EKF(i,2*j - 1) = xkm_pos_sig_EKF(:,j);
        x_pos_sig_comb_EKF(i,2*j) = xkp_pos_sig_EKF(:,j);
        y_pos_sig_comb_EKF(i,2*j - 1) = ykm_pos_sig_EKF(:,j);
        y_pos_sig_comb_EKF(i,2*j) = ykp_pos_sig_EKF(:,j);
        z_pos_sig_comb_EKF(i,2*j - 1) = zkm_pos_sig_EKF(:,j);
        z_pos_sig_comb_EKF(i,2*j) = zkp_pos_sig_EKF(:,j);
        x_vel_sig_comb_EKF(i,2*j - 1) = xkm_vel_sig_EKF(:,j);
        x_vel_sig_comb_EKF(i,2*j) = xkp_vel_sig_EKF(:,j);
        y_vel_sig_comb_EKF(i,2*j - 1) = ykm_vel_sig_EKF(:,j);
        y_vel_sig_comb_EKF(i,2*j) = ykp_vel_sig_EKF(:,j);
        z_vel_sig_comb_EKF(i,2*j - 1) = zkm_vel_sig_EKF(:,j);
        z_vel_sig_comb_EKF(i,2*j) = zkp_vel_sig_EKF(:,j);

        x_pos_sig_comb_UKF(i,2*j - 1) = xkm_pos_sig_UKF(:,j);
        x_pos_sig_comb_UKF(i,2*j) = xkp_pos_sig_UKF(:,j);
        y_pos_sig_comb_UKF(i,2*j - 1) = ykm_pos_sig_UKF(:,j);
        y_pos_sig_comb_UKF(i,2*j) = ykp_pos_sig_UKF(:,j);
        z_pos_sig_comb_UKF(i,2*j - 1) = zkm_pos_sig_UKF(:,j);
        z_pos_sig_comb_UKF(i,2*j) = zkp_pos_sig_UKF(:,j);
        x_vel_sig_comb_UKF(i,2*j - 1) = xkm_vel_sig_UKF(:,j);
        x_vel_sig_comb_UKF(i,2*j) = xkp_vel_sig_UKF(:,j);
        y_vel_sig_comb_UKF(i,2*j - 1) = ykm_vel_sig_UKF(:,j);
        y_vel_sig_comb_UKF(i,2*j) = ykp_vel_sig_UKF(:,j);
        z_vel_sig_comb_UKF(i,2*j - 1) = zkm_vel_sig_UKF(:,j);
        z_vel_sig_comb_UKF(i,2*j) = zkp_vel_sig_UKF(:,j);
    
        % state estimation error combined arrays
        x_pos_error_comb_EKF(i,2*j - 1) = xkm_pos_error_EKF(:,j);
        x_pos_error_comb_EKF(i,2*j) = xkp_pos_error_EKF(:,j);
        y_pos_error_comb_EKF(i,2*j - 1) = ykm_pos_error_EKF(:,j);
        y_pos_error_comb_EKF(i,2*j) = ykp_pos_error_EKF(:,j);
        z_pos_error_comb_EKF(i,2*j - 1) = zkm_pos_error_EKF(:,j);
        z_pos_error_comb_EKF(i,2*j) = zkp_pos_error_EKF(:,j);
        x_vel_error_comb_EKF(i,2*j - 1) = xkm_vel_error_EKF(:,j);
        x_vel_error_comb_EKF(i,2*j) = xkp_vel_error_EKF(:,j);
        y_vel_error_comb_EKF(i,2*j - 1) = ykm_vel_error_EKF(:,j);
        y_vel_error_comb_EKF(i,2*j) = ykp_vel_error_EKF(:,j);
        z_vel_error_comb_EKF(i,2*j - 1) = zkm_vel_error_EKF(:,j);
        z_vel_error_comb_EKF(i,2*j) = zkp_vel_error_EKF(:,j);

        x_pos_error_comb_UKF(i,2*j - 1) = xkm_pos_error_UKF(:,j);
        x_pos_error_comb_UKF(i,2*j) = xkp_pos_error_UKF(:,j);
        y_pos_error_comb_UKF(i,2*j - 1) = ykm_pos_error_UKF(:,j);
        y_pos_error_comb_UKF(i,2*j) = ykp_pos_error_UKF(:,j);
        z_pos_error_comb_UKF(i,2*j - 1) = zkm_pos_error_UKF(:,j);
        z_pos_error_comb_UKF(i,2*j) = zkp_pos_error_UKF(:,j);
        x_vel_error_comb_UKF(i,2*j - 1) = xkm_vel_error_UKF(:,j);
        x_vel_error_comb_UKF(i,2*j) = xkp_vel_error_UKF(:,j);
        y_vel_error_comb_UKF(i,2*j - 1) = ykm_vel_error_UKF(:,j);
        y_vel_error_comb_UKF(i,2*j) = ykp_vel_error_UKF(:,j);
        z_vel_error_comb_UKF(i,2*j - 1) = zkm_vel_error_UKF(:,j);
        z_vel_error_comb_UKF(i,2*j) = zkp_vel_error_UKF(:,j);
    end
end

% color arrays for plotting
color_EKF_error = [0.9, 0.3, 0.1, 0.1];
color_UKF_error = [0.0, 0.4, 0.7, 0.4];
color_EKF_sigma = [1.0, 0.0, 1.0, 0.1];
color_UKF_sigma = [0.5, 0.7, 0.2, 0.1];

% color arrays for legend
color_EKF_error_legend = [0.9, 0.3, 0.1, 1];
color_UKF_error_legend = [0.0, 0.4, 0.7, 1];
color_EKF_sigma_legend = [1.0, 0.0, 1.0, 1];
color_UKF_sigma_legend = [0.5, 0.7, 0.2, 1];

% plot position errors with 3 sigma bounds
figure(1)
for i = 1:MC_trials
    plot(t_comb,x_pos_error_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_error)
    hold on
    plot(t_comb,3*x_pos_sig_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_sigma)
    plot(t_comb,-3*x_pos_sig_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_sigma)
    plot(t_comb,x_pos_error_comb_EKF(i,:),'LineWidth',0.1,'Color',color_EKF_error)
    plot(t_comb,3*x_pos_sig_comb_EKF(i,:),'LineWidth',0.1,'color',color_EKF_sigma)
    plot(t_comb,-3*x_pos_sig_comb_EKF(i,:),'LineWidth',0.1,'color',color_EKF_sigma)
end
title('EKF and UKF x Position Error with 3\sigma Bounds in ECI Frame')
xlabel('Time (days)')
ylabel('x Position Error (km)')
h1 = plot(NaN,NaN,'Color',color_UKF_error_legend,'LineWidth',1.5);
h2 = plot(NaN,NaN,'Color',color_UKF_sigma_legend,'LineWidth',1.5);
h3 = plot(NaN,NaN,'Color',color_EKF_error_legend,'LineWidth',1.5);
h4 = plot(NaN,NaN,'Color',color_EKF_sigma_legend,'LineWidth',1.5);
legend([h1 h2 h3 h4],'UKF Error','UKF 3\sigma','EKF Error','EKF 3\sigma','Location','northeast')
grid on
    
figure(2)
for i = 1:MC_trials
    plot(t_comb,y_pos_error_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_error)
    hold on
    plot(t_comb,3*y_pos_sig_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_sigma)
    plot(t_comb,-3*y_pos_sig_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_sigma)
    plot(t_comb,y_pos_error_comb_EKF(i,:),'LineWidth',0.1,'Color',color_EKF_error)
    plot(t_comb,3*y_pos_sig_comb_EKF(i,:),'LineWidth',0.1,'color',color_EKF_sigma)
    plot(t_comb,-3*y_pos_sig_comb_EKF(i,:),'LineWidth',0.1,'color',color_EKF_sigma)
end
title('EKF and UKF y Position Error with 3\sigma Bounds in ECI Frame')
xlabel('Time (days)')
ylabel('y Position Error (km)')
h1 = plot(NaN,NaN,'Color',color_UKF_error_legend,'LineWidth',1.5);
h2 = plot(NaN,NaN,'Color',color_UKF_sigma_legend,'LineWidth',1.5);
h3 = plot(NaN,NaN,'Color',color_EKF_error_legend,'LineWidth',1.5);
h4 = plot(NaN,NaN,'Color',color_EKF_sigma_legend,'LineWidth',1.5);
legend([h1 h2 h3 h4],'UKF Error','UKF 3\sigma','EKF Error','EKF 3\sigma','Location','northeast')
grid on
    
figure(3)
for i = 1:MC_trials
    plot(t_comb,z_pos_error_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_error)
    hold on
    plot(t_comb,3*z_pos_sig_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_sigma)
    plot(t_comb,-3*z_pos_sig_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_sigma)
    plot(t_comb,z_pos_error_comb_EKF(i,:),'LineWidth',0.1,'Color',color_EKF_error)
    plot(t_comb,3*z_pos_sig_comb_EKF(i,:),'LineWidth',0.1,'color',color_EKF_sigma)
    plot(t_comb,-3*z_pos_sig_comb_EKF(i,:),'LineWidth',0.1,'color',color_EKF_sigma)
end
title('EKF and UKF z Position Error with 3\sigma Bounds in ECI Frame')
xlabel('Time (days)')
ylabel('z Position Error (km)')
h1 = plot(NaN,NaN,'Color',color_UKF_error_legend,'LineWidth',1.5);
h2 = plot(NaN,NaN,'Color',color_UKF_sigma_legend,'LineWidth',1.5);
h3 = plot(NaN,NaN,'Color',color_EKF_error_legend,'LineWidth',1.5);
h4 = plot(NaN,NaN,'Color',color_EKF_sigma_legend,'LineWidth',1.5);
legend([h1 h2 h3 h4],'UKF Error','UKF 3\sigma','EKF Error','EKF 3\sigma','Location','northeast')
grid on
    
% plot velocity errors and 3 sigma bounds
figure(4)
for i = 1:MC_trials
    plot(t_comb,x_vel_error_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_error)
    hold on
    plot(t_comb,3*x_vel_sig_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_sigma)
    plot(t_comb,-3*x_vel_sig_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_sigma)
    plot(t_comb,x_vel_error_comb_EKF(i,:),'LineWidth',0.1,'Color',color_EKF_error)
    plot(t_comb,3*x_vel_sig_comb_EKF(i,:),'LineWidth',0.1,'color',color_EKF_sigma)
    plot(t_comb,-3*x_vel_sig_comb_EKF(i,:),'LineWidth',0.1,'color',color_EKF_sigma)
end
title('EKF and UKF x Velocity Error with 3\sigma Bounds in ECI Frame')
xlabel('Time (days)')
ylabel('x Velocity Error (km/s)')
h1 = plot(NaN,NaN,'Color',color_UKF_error_legend,'LineWidth',1.5);
h2 = plot(NaN,NaN,'Color',color_UKF_sigma_legend,'LineWidth',1.5);
h3 = plot(NaN,NaN,'Color',color_EKF_error_legend,'LineWidth',1.5);
h4 = plot(NaN,NaN,'Color',color_EKF_sigma_legend,'LineWidth',1.5);
legend([h1 h2 h3 h4],'UKF Error','UKF 3\sigma','EKF Error','EKF 3\sigma','Location','northeast')
grid on
    
figure(5)
for i = 1:MC_trials
    plot(t_comb,y_vel_error_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_error)
    hold on
    plot(t_comb,3*y_vel_sig_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_sigma)
    plot(t_comb,-3*y_vel_sig_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_sigma)
    plot(t_comb,y_vel_error_comb_EKF(i,:),'LineWidth',0.1,'Color',color_EKF_error)
    plot(t_comb,3*y_vel_sig_comb_EKF(i,:),'LineWidth',0.1,'color',color_EKF_sigma)
    plot(t_comb,-3*y_vel_sig_comb_EKF(i,:),'LineWidth',0.1,'color',color_EKF_sigma)
end
title('EKF and UKF y Velocity Error with 3\sigma Bounds in ECI Frame')
xlabel('Time (days)')
ylabel('y Velocity Error (km/s)')
h1 = plot(NaN,NaN,'Color',color_UKF_error_legend,'LineWidth',1.5);
h2 = plot(NaN,NaN,'Color',color_UKF_sigma_legend,'LineWidth',1.5);
h3 = plot(NaN,NaN,'Color',color_EKF_error_legend,'LineWidth',1.5);
h4 = plot(NaN,NaN,'Color',color_EKF_sigma_legend,'LineWidth',1.5);
legend([h1 h2 h3 h4],'UKF Error','UKF 3\sigma','EKF Error','EKF 3\sigma','Location','northeast')
grid on
    
figure(6)
for i = 1:MC_trials
    plot(t_comb,z_vel_error_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_error)
    hold on
    plot(t_comb,3*z_vel_sig_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_sigma)
    plot(t_comb,-3*z_vel_sig_comb_UKF(i,:),'LineWidth',0.1,'Color',color_UKF_sigma)
    plot(t_comb,z_vel_error_comb_EKF(i,:),'LineWidth',0.1,'Color',color_EKF_error)
    plot(t_comb,3*z_vel_sig_comb_EKF(i,:),'LineWidth',0.1,'color',color_EKF_sigma)
    plot(t_comb,-3*z_vel_sig_comb_EKF(i,:),'LineWidth',0.1,'color',color_EKF_sigma)
end
title('EKF and UKF z Velocity Error with 3\sigma Bounds in ECI Frame')
xlabel('Time (days)')
ylabel('z Velocity Error (km/s)')
h1 = plot(NaN,NaN,'Color',color_UKF_error_legend,'LineWidth',1.5);
h2 = plot(NaN,NaN,'Color',color_UKF_sigma_legend,'LineWidth',1.5);
h3 = plot(NaN,NaN,'Color',color_EKF_error_legend,'LineWidth',1.5);
h4 = plot(NaN,NaN,'Color',color_EKF_sigma_legend,'LineWidth',1.5);
legend([h1 h2 h3 h4],'UKF Error','UKF 3\sigma','EKF Error','EKF 3\sigma','Location','northeast')
grid on

% preallocate NEES and ASNEES values
NEES_EKF = NaN(1,MC_trials);
NEES_UKF = NaN(1,MC_trials);
ASNEES_EKF = NaN(1,num_time_steps);
ASNEES_UKF = NaN(1,num_time_steps);

% calculate NEES and ASNEES values for consistency analysis
for j = 1:num_time_steps
    for k = 1:MC_trials
        Pk_UKF = squeeze(Pkp_hist_UKF(:,:,j,k));
        ek_UKF = x_high_fidel(:,j,k) - xkp_hist_UKF(:,j,k);

        Pk_EKF = squeeze(Pkp_hist_EKF(:,:,j,k));
        ek_EKF = x_high_fidel(:,j,k) - xkp_hist_EKF(:,j,k);

        NEES_EKF(k) = ek_EKF'/Pk_EKF*ek_EKF;
        NEES_UKF(k) = ek_UKF'/Pk_UKF*ek_UKF;
    end
    ASNEES_EKF(j) = 1/(n_sat*MC_trials)*sum(NEES_EKF);
    ASNEES_UKF(j) = 1/(n_sat*MC_trials)*sum(NEES_UKF);
end

% plot ASNEES over time
figure(7)
plot(t_ASNEES,ASNEES_UKF,'b')
hold on
plot(t_ASNEES,ASNEES_EKF,'r')
title('EKF and UKF ASNEES Over Time')
xlabel('Time (days)')
ylabel('ASNEES')
legend('UKF', 'EKF', 'Location', 'northeast')
grid on

% plot measurements, high fidelity simulation, and UKF posterior position
% estimates
figure(8)
plot3(x_meas, y_meas, z_meas)
hold on
plot3(x_high_fidel(1,:), x_high_fidel(3,:), x_high_fidel(5,:))
plot3(xkp_hist_UKF(1,:), xkp_hist_UKF(3,:), xkp_hist_UKF(5,:))
title('Satellite Position and Measurements in ECI Frame')
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
legend('Measurements','True Orbit','UKF Estimate')
axis equal
grid on

%% Function Definitions
% high fidelity dynamics function
function dxdt = ODE_High_Fidelity(t, x, param, w_Earth)
    % extract parameters from param struct
    mu = param.mu_Earth;
    A_m = param.area_mass;
    dens_scale = param.dens_scale;
    h0 = param.h0;
    rho_ref = param.rho_ref;
    CD = param.CD;
    P = param.P;
    nu = param.nu;
    CR = param.CR;

    % store states to meaningful variables
    x_pos = x(1);
    xdot = x(2);
    y_pos = x(3);
    ydot = x(4);
    z_pos = x(5);
    zdot = x(6);

    % determine radial distance and relative velocity to atmosphere
    r = sqrt(x_pos^2 + y_pos^2 + z_pos^2);
    v_rel = [xdot; ydot; zdot] - cross([0; 0; w_Earth], [x_pos; y_pos; z_pos]);
    v_mag = sqrt(v_rel(1)^2 + v_rel(2)^2 + v_rel(3)^2);

    % calculate drag and solar pressure forces
    dens = rho_ref*exp(-dens_scale*(r - h0));
    drag = 0.5*dens*CD*A_m*v_mag*v_rel*1000; % multiply by 1000 to convert from km to m
    solar = -P*nu*A_m*CR*simple_sun_vec(t);

    dxdt(1) = xdot;
    dxdt(2) = -mu*x_pos/r^3 - drag(1) - solar(1);
    dxdt(3) = ydot;
    dxdt(4) = -mu*y_pos/r^3 - drag(2) - solar(2);
    dxdt(5) = zdot;
    dxdt(6) = -mu*z_pos/r^3 - drag(3) - solar(3);
    dxdt = dxdt';
end

function sun_vec = simple_sun_vec(t)
    years2days = 365.25; % days in a year
    omega_sun = 2*pi / (years2days*24*3600); % rotation around sun (rad/s)

    theta = omega_sun*t; % rotation around sun (rad)
    sun_vec = [cos(theta); sin(theta); 0];
end

% low fidelity dynamics function for model mismatch
function dxdt = ODE_Low_Fidelity(t, x, mu, Q, n_sat, EKF_UKF)
    % store states to meaningful variables
    x_pos = x(1);
    xdot = x(2);
    y_pos = x(3);
    ydot = x(4);
    z_pos = x(5);
    zdot = x(6);

    r = sqrt(x_pos^2 + y_pos^2 + z_pos^2);
    dxdt(1) = xdot;
    dxdt(2) = -mu*x_pos/r^3;
    dxdt(3) = ydot;
    dxdt(4) = -mu*y_pos/r^3;
    dxdt(5) = zdot;
    dxdt(6) = -mu*z_pos/r^3;
    dxdt = dxdt';

    if EKF_UKF == 1
        F = Dynamics_Jacobian(x(1:n_sat),mu);
        P = reshape(x(n_sat + 1:end),n_sat,n_sat);
        Pdot = F*P + P*F' + Q;
        dxdt = [dxdt; reshape(Pdot,n_sat^2,1)];
    end
end

% dynamics jacobian matrix calculation
function F = Dynamics_Jacobian(state, mu)
    state = state(:);
    % store states to meaningful variables
    x = state(1);
    xdot = state(2);
    y = state(3);
    ydot = state(4);
    z = state(5);
    zdot = state(6);
    
    r = sqrt(x^2 + y^2 + z^2);
    
    F = [            0             1               0             0               0             0;
         -mu*(1/r^3 - 3*x^2/r^5)   0         mu*3*x*y/r^5        0         mu*3*x*z/r^5        0;
                     0             0               0             1               0             0;
               mu*3*x*y/r^5        0   -mu*(1/r^3 - 3*y^2/r^5)   0         mu*3*y*z/r^5        0;
                     0             0               0             0               0             1;
               mu*3*x*z/r^5        0         mu*3*y*z/r^5        0   -mu*(1/r^3 - 3*z^2/r^5)   0];
end

function H = Meas_Jacobian(x_sat, r_gs_ECEF, t, w_Earth, T_t, el_low, el_high, meas_model)
    % extract satellite position and velocity
    r_sat = [x_sat(1); x_sat(3); x_sat(5)];
    v_sat = [x_sat(2); x_sat(4); x_sat(6)];

    % ground station position and velocity in ECI
    r_gs_ECI = ECEF2ECI(r_gs_ECEF, t, w_Earth);
    v_gs_ECI = cross([0; 0; w_Earth], r_gs_ECI);

    % relative position and velocity in ECI
    r_rel = r_sat - r_gs_ECI;
    v_rel = v_sat - v_gs_ECI;

    r_rel_ECEF = ECI2ECEF(r_rel, t, w_Earth);

    rho = norm(r_rel);
    rho_dot = r_rel'*v_rel/rho;

    % convert relative ECEF position to topocentric frame [east, north, up]
    r_t_vec = T_t*r_rel_ECEF;

    x_t = r_t_vec(1); % east
    y_t = r_t_vec(2); % north
    z_t = r_t_vec(3); % up
    r_t = norm(r_t_vec);
    r_xy = sqrt(x_t^2 + y_t^2);

    if meas_model == 1
        % Use numerical differentiation for angular measurements
        eps = 1e-9;
        n_state = length(x_sat);
        n_meas = 2;
        H = zeros(n_meas, n_state);

        % calculate derivate using central differencing
        for i = 1:n_state
            x_pert_plus = x_sat;
            x_pert_minus = x_sat;
            x_pert_plus(i) = x_pert_plus(i) + eps;
            x_pert_minus(i) = x_pert_minus(i) - eps;
            z_pert_plus = meas_func(x_pert_plus, r_gs_ECEF, t, T_t, w_Earth, 1, el_low, el_high);
            z_pert_minus = meas_func(x_pert_minus, r_gs_ECEF, t, T_t, w_Earth, 1, el_low, el_high);

            % Handle angle wrapping for azimuth
            dz = z_pert_plus - z_pert_minus;
            if abs(dz(2)) > pi
                dz(2) = dz(2) - sign(dz(2))*2*pi;
            end

            H(:,i) = dz/(2*eps);
        end
    elseif meas_model == 2
        H = [               r_rel(1)/rho                     0                       r_rel(2)/rho                     0                       r_rel(3)/rho                     0;
             (v_rel(1) - rho_dot*(r_rel(1)/rho))/rho   r_rel(1)/rho   (v_rel(2) - rho_dot*(r_rel(2)/rho))/rho   r_rel(2)/rho   (v_rel(3) - rho_dot*(r_rel(3)/rho))/rho   r_rel(3)/rho];
    else
        % Use numerical differentiation for angular measurements
        eps = 1e-9;
        n_state = length(x_sat);
        n_meas = 2;
        H_angles = zeros(n_meas, n_state);

        % calculate derivate using central differencing
        for i = 1:n_state
            x_pert_plus = x_sat;
            x_pert_minus = x_sat;
            x_pert_plus(i) = x_pert_plus(i) + eps;
            x_pert_minus(i) = x_pert_minus(i) - eps;
            z_pert_plus = meas_func(x_pert_plus, r_gs_ECEF, t, T_t, w_Earth, 1, el_low, el_high);
            z_pert_minus = meas_func(x_pert_minus, r_gs_ECEF, t, T_t, w_Earth, 1, el_low, el_high);

            % Handle angle wrapping for azimuth
            dz = z_pert_plus - z_pert_minus;
            if abs(dz(2)) > pi
                dz(2) = dz(2) - sign(dz(2))*2*pi;
            end

            H_angles(:,i) = dz/(2*eps);
        end
        
        H_range_rate = [               r_rel(1)/rho                     0                       r_rel(2)/rho                     0                       r_rel(3)/rho                     0;
                        (v_rel(1) - rho_dot*(r_rel(1)/rho))/rho   r_rel(1)/rho   (v_rel(2) - rho_dot*(r_rel(2)/rho))/rho   r_rel(2)/rho   (v_rel(3) - rho_dot*(r_rel(3)/rho))/rho   r_rel(3)/rho];
        H = [H_range_rate; H_angles];
    end
end

% measurement model with logic to select chosen model
function z = meas_func(x_sat, r_antenna, t, T_t, w_Earth, meas_model, low, high, v)
    % handle if noise not given
    if nargin < 9
        v_range_i = 0;
        v_range_rate_i = 0;
        v_el_i = 0;
        v_az_i = 0;
    else
        % extract measurement noise terms
        v_range_i = v(1);
        v_range_rate_i = v(2);
        v_el_i = v(3);
        v_az_i = v(4);
    end

    % convert satellite position from ECI to ECEF
    r_sat_ECI = [x_sat(1); x_sat(3); x_sat(5)];
    r_sat_ECEF = ECI2ECEF(r_sat_ECI, t, w_Earth);

    v_sat_ECI = [x_sat(2); x_sat(4); x_sat(6)];

    % calculate relative velocity in ECI
    r_antenna_ECI = ECEF2ECI(r_antenna, t, w_Earth);
    v_antenna_ECI = cross([0; 0; w_Earth], r_antenna_ECI);

    % calculate relative positon and velocity of satellite to ground station
    r_rel_ECEF = r_sat_ECEF - r_antenna;
    r_rel_ECI = r_sat_ECI - r_antenna_ECI;
    r_rel = norm(r_rel_ECEF);
    v_rel_ECI = v_sat_ECI - v_antenna_ECI;
    
    % convert relative ECEF position to topocentric frame [east, north, up]
    r_t_vec = T_t*r_rel_ECEF;
    
    x_t = r_t_vec(1); % east
    y_t = r_t_vec(2); % north
    z_t = r_t_vec(3); % up
    r_t = norm(r_t_vec);
    
    % calculate range and angle measurements
    z_k_el = asin(z_t/r_t);
    
    % set measurement limits. Ran into issues when implementing
    % measurement limits. Covariance matrix would result in all "NaN"
    % entries and could not be resolved. Uncomment the code below to
    % implement measurement limits
    % if z_k_el < low || z_k_el > high
    %     z_k_range = NaN;
    %     z_k_range_rate = NaN;
    %     z_k_el = NaN;
    %     z_k_az = NaN;
    % else
        z_k_range = r_rel + v_range_i;
        z_k_range_rate = r_rel_ECI'*v_rel_ECI/r_rel + v_range_rate_i;
        z_k_el = asin(z_t/r_t) + v_el_i;
        z_k_az = wrapTo2Pi(atan2(x_t,y_t) + v_az_i);
    % end

    % combine measurements based on chosen measurement model
    if meas_model == 1
        z = [z_k_el; z_k_az];
    elseif meas_model == 2
        z = [z_k_range; z_k_range_rate];
    else
        z = [z_k_range; z_k_range_rate; z_k_el; z_k_az];
    end
end

% calculates sigma points and mean and covariance weights for SUT
function [sig_points, wm, wc] = SigPoints(P, m, alpha, beta, kappa, n)
    % SUT weighting parameters
    lambda = alpha^2*(n + kappa) - n; % weighting parameter
    S = chol(P,'lower'); % square root factor of the covariance matrix
    
    % calculate sigma points
    sig_points = zeros(n,2*n + 1);
    sig_points(:,1) = m;
    
    for i = 2:(n + 1)
        sig_points(:,i) = m + sqrt(n + lambda)*S(:,i - 1);
        sig_points(:,i + n) = m - sqrt(n + lambda)*S(:,i - 1);
    end
    
    % determine weighting terms
    w_0m = lambda/(n + lambda);
    w_0c = lambda/(n + lambda) + (1 - alpha^2 + beta);
    w_im = 1/(2*(n + lambda));
    w_ic = 1/(2*(n + lambda));

    % combine weighting terms into corresponding arrays
    wm = [w_0m, repmat(w_im, 1, 2*n)];
    wc = [w_0c, repmat(w_ic, 1, 2*n)];
end

% ECEF to topocentric frame [x_t = east, y_t = north, z_t = up]
function T = ECEF2Topo(phi, lambda)
    C = @(angle) cos(angle);
    S = @(angle) sin(angle);
    T = [    -S(lambda),          C(lambda),              0;
         -S(phi)*C(lambda),   -S(phi)*S(lambda),     C(phi);
          C(phi)*C(lambda),     C(phi)*S(lambda),    S(phi)];
end

% simple ECEF to ECI coordinate transformation assuming Earth spins at
% constant angular velocity
function r_ECEF = ECI2ECEF(r_ECI, t, w)
    theta = w*t; % simple Earth rotation
    Rz = [cos(theta),   -sin(theta),   0;
          sin(theta),    cos(theta),   0;
              0,             0,        1];
    r_ECEF = Rz*r_ECI;
end

function r_ECI = ECEF2ECI(r_ECEF, t, omegaE)
    theta = omegaE*t;
    Rz = [cos(theta),   -sin(theta),   0;
          sin(theta),    cos(theta),   0;
              0,             0,        1];
    r_ECI = Rz'*r_ECEF;
end

% rotation matrix about axis 1 for 3-1-3 sequence
function R = DCM1(angle)
    R = [1       0              0;
         0   cos(angle)   -sin(angle);
         0   sin(angle)    cos(angle)];
end

% rotation matrix about axis 3 for 3-1-3 sequence
function R = DCM3(angle)
    R = [cos(angle)   -sin(angle)   0;
         sin(angle)    cos(angle)   0;
             0             0        1];
end
