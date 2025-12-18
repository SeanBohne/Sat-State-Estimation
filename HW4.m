clc; clear; close all;

%% Part A
% Initialization
Pxx = [4 0.5; 0.5 1]; % original covariance matrix
mx = [7; 0]; % state mean [r, theta]
points = 500; % determines number of (r,theta) pairs
r_range = linspace(0,15,points); % set range of r values
theta_range = linspace(-pi,pi,points); % set range of theta values

% Calculations
gauss_const = 1/sqrt(2*pi*det(Pxx)); % normal distribution constant
[R, T] = meshgrid(r_range, theta_range); % create (r,theta) pairs
pdf_a = zeros(size(R)); % preallocate pdf values

for i = 1:points
    for j = 1:points
        x = [R(i,j); T(i,j)]; % define operating point for pdf calculation
        l_sq = (x - mx)'/Pxx*(x - mx); % mahalanobis distance
        pdf_a(i,j) = gauss_const*exp(-0.5*l_sq); % calculate pdf values
    end
end

% plot pdf as a colored contour
figure
pcolor(R, T, pdf_a)
shading interp
hold on

ellipse = sigmabound(Pxx, mx, 3);

% plot ellipse over pdf distribution
plot(ellipse(1,:),ellipse(2,:),'r','LineWidth',2)
title('Polar Coordinates Probability Density Function With 3\sigma Bound - Sean Bohne','FontSize',8)
xlabel('r (distance units)')
ylabel('\theta (rad)')

%% Part C
% Initialization
x = linspace(-25,25,points); % x cartesian coordinates transformation
y = linspace(-25,25,points); % y cartesian coordinates transformation
[X, Y] = meshgrid(x, y); % create a grid for x and y values
pdf_cart_c = zeros(size(X)); % preallocate cartesian pdf values

% Calculations
for i = 1:points
    for j = 1:points
        if sqrt(X(i,j)^2 + Y(i,j)^2) <= 1e-2
            pdf_cart_c(i,j) = 0;
        else
            x_val = [sqrt(X(i,j)^2 + Y(i,j)^2); atan2(Y(i,j), X(i,j))]; % define operating point for pdf calculation
            l_sq = (x_val - mx)'/Pxx*(x_val - mx); % mahalanobis distance
            pdf_cart_c(i,j) = 1/sqrt((X(i,j)^2 + Y(i,j)^2))*gauss_const*exp(-0.5*l_sq); % calculate pdf values
        end
    end
end

% plot cartesian pdf as a colored contour
figure
pcolor(X, Y, pdf_cart_c)
shading interp
title('Cartesian Coordinates Probability Density Function With 3\sigma Bound - Sean Bohne','FontSize',8)
xlabel('x')
ylabel('y')
hold on

%% Part D
% EKF linearization
mz_EKF = [mx(1)*cos(mx(2)); mx(1)*sin(mx(2))];
H = [cos(mx(2)), -mx(1)*sin(mx(2)); sin(mx(2)), mx(1)*cos(mx(2))];
Pzz_EKF = H*Pxx*H';

% scaled unscented transform
alpha = 1*10^-3; % weighting factor parameter
beta = 2; % weighting factor parameter
kappa = 0; % weighting factor parameter
n = length(mx); % state dimension

[Pzz_sut, mz_sut] = SUT(Pxx, mx, alpha, beta, kappa, n);

% SUT 3 sigma ellipse
ellipse_sut = sigmabound(Pzz_sut, mz_sut, 3);

% EKF 3 sigma ellipse
ellipse_EKF = sigmabound(Pzz_EKF, mz_EKF, 3);

fprintf('---------------------------------------\n')
fprintf('EKF linearization mean approximation m_z = \n')
disp(mz_EKF)
fprintf('EKF linearization covariance approximation P_{zz} = \n')
disp(Pzz_EKF)
fprintf('---------------------------------------\n')
fprintf('SUT mean approximation m_z = \n')
disp(mz_sut)
fprintf('SUT covariance approximation P_{zz} = \n')
disp(Pzz_sut)
fprintf('---------------------------------------\n')

% plot ellipse over pdf distribution
figure
pcolor(X, Y, pdf_cart_c)
shading interp
title('Cartesian Coordinates Probability Density Function With 3\sigma Bound - Sean Bohne','FontSize',8)
xlabel('x')
ylabel('y')
hold on

plot(ellipse_sut(1,:),ellipse_sut(2,:),'r','LineWidth',2)
plot(ellipse_EKF(1,:),ellipse_EKF(2,:),'g','LineWidth',2)
legend('p(z)','SUT','EKF','Location','southwest')

%% Part E
% initialize gamma values
gamma = 0.9:-0.2:0.1;

for i = 1:length(gamma)
    pdf_e = zeros(points, points);
    pdf_cart_e = zeros(points, points);

    P = gamma(i)*Pxx; % scale covariance matrix by corresponding gamma

    % EKF linearization
    mz_EKF = [mx(1)*cos(mx(2)); mx(1)*sin(mx(2))];
    H = [cos(mx(2)), -mx(1)*sin(mx(2)); sin(mx(2)), mx(1)*cos(mx(2))];
    Pzz_EKF = H*P*H';
    
    % scaled unscented transform
    [Pzz_sut, mz_sut] = SUT(P, mx, alpha, beta, kappa, n);

    % original 3 sigma ellipse
    ellipse = sigmabound(P, mx, 3);

    % SUT 3 sigma ellipse
    ellipse_sut = sigmabound(Pzz_sut, mz_sut, 3);

    % EKF 3 sigma ellipse
    ellipse_EKF = sigmabound(Pzz_EKF, mz_EKF, 3);

    gauss_const = 1/sqrt(2*pi*det(P));

    for j = 1:points
        for k = 1:points
            x = [R(j,k); T(j,k)]; % define operating point for pdf calculation
            l_sq = (x - mx)'/P*(x - mx); % mahalanobis distance
            pdf_e(j,k) = gauss_const*exp(-0.5*l_sq); % calculate pdf values
        end
    end

    for j = 1:points
        for k = 1:points
            if sqrt(X(i,j)^2 + Y(i,j)^2) <= 1e-2
                pdf_cart_c(i,j) = 0;
            else
                x_val = [sqrt(X(j,k)^2 + Y(j,k)^2); atan2(Y(j,k), X(j,k))]; % define operating point for pdf calculation
                l_sq = (x_val - mx)'/P*(x_val - mx); % mahalanobis distance
                pdf_cart_e(j,k) = 1/sqrt((X(j,k)^2 + Y(j,k)^2))*gauss_const*exp(-0.5*l_sq); % calculate pdf values
            end
        end
    end
    
    % plot original pdf distribution from part A
    figure
    subplot(2,2,1)
    pcolor(R, T, pdf_e)
    shading interp
    hold on

    plot(ellipse(1,:),ellipse(2,:),'r','LineWidth',2)
    title(['Input p(x) - \gamma = ',num2str(gamma(i))])
    xlabel('r (distance units)')
    ylabel('\theta (rad)')

    % plot pdf distribution from part C
    subplot(2,2,2)
    pcolor(X, Y, pdf_cart_e)
    shading interp
    title(['Cartesian Transformation p(z) - \gamma = ',num2str(gamma(i))])
    xlabel('x')
    ylabel('y')
    
    % plot pdf distribution and 3 sigma ellipses from from part D
    subplot(2,2,3)
    pcolor(X, Y, pdf_cart_e)
    shading interp
    title(['SUT and EKF Approximations - \gamma = ',num2str(gamma(i))])
    xlabel('x')
    ylabel('y')
    hold on

    % plot 3 sigma ellipse over pdf distribution
    plot(ellipse_sut(1,:),ellipse_sut(2,:),'r','LineWidth',2)
    plot(ellipse_EKF(1,:),ellipse_EKF(2,:),'g','LineWidth',2)
    legend('p(z)','SUT','EKF','Location','southwest')
    sgtitle(['\gamma = ',num2str(gamma(i)),' - Sean Bohne'])
end

%% Part F
% vary SUT weighting parameters one by one starting with alpha to analyze
% its impact on pdf estimation performance
alpha = [0.001, 0.01, 0.1, 1];
beta = 2;
kappa = 0;

% recalculate pdf distribution for plotting
gauss_const = 1/sqrt(2*pi*det(Pxx));

for i = 1:points
    for j = 1:points
        if sqrt(X(i,j)^2 + Y(i,j)^2) <= 1e-2
            pdf_cart_c(i,j) = 0;
        else
            x_val = [sqrt(X(i,j)^2 + Y(i,j)^2); atan2(Y(i,j), X(i,j))]; % define operating point for pdf calculation
            l_sq = (x_val - mx)'/Pxx*(x_val - mx); % mahalanobis distance
            pdf_cart_f(i,j) = 1/sqrt((X(i,j)^2 + Y(i,j)^2))*gauss_const*exp(-0.5*l_sq); % calculate pdf values
        end
    end
end

% calculate SUT for each alpha and plot the distributions
figure
for i = 1:length(alpha)
    [Pzz_sut, mz_sut] = SUT(Pxx, mx, alpha(i), beta, kappa, n);
    ellipse_sut = sigmabound(Pzz_sut, mz_sut, 3);

    % plot pdf distribution and 3 sigma ellipses from from part D
    subplot(floor(sqrt(length(alpha))),ceil(sqrt(length(alpha))),i)
    pcolor(X, Y, pdf_cart_f)
    shading interp
    title(['\alpha = ',num2str(alpha(i)),', \beta = ',num2str(beta),', \kappa = ',num2str(kappa)])
    xlabel('x')
    ylabel('y')
    hold on

    % plot 3 sigma ellipse over pdf distribution
    plot(ellipse_sut(1,:),ellipse_sut(2,:),'r','LineWidth',2)
    legend('p(z)','SUT','Location','southwest')
    sgtitle('Varying \alpha - Sean Bohne')
end

% vary beta with other parameters fixed to analyze its impact
alpha = 1*10^-3;
beta = 0:3;
kappa = 0;

figure
for i = 1:length(beta)
    [Pzz_sut, mz_sut] = SUT(Pxx, mx, alpha, beta(i), kappa, n);
    ellipse_sut = sigmabound(Pzz_sut, mz_sut, 3);

    % plot pdf distribution and 3 sigma ellipses from from part D
    subplot(floor(sqrt(length(beta))),ceil(sqrt(length(beta))),i)
    pcolor(X, Y, pdf_cart_f)
    shading interp
    title(['\alpha = ',num2str(alpha),', \beta = ',num2str(beta(i)),', \kappa = ',num2str(kappa)])
    xlabel('x')
    ylabel('y')
    hold on

    % plot 3 sigma ellipse over pdf distribution
    plot(ellipse_sut(1,:),ellipse_sut(2,:),'r','LineWidth',2)
    legend('p(z)','SUT','Location','southwest')
    sgtitle('Varying \beta - Sean Bohne')
end

% vary kappa with other parameters fixed to analyze its impact
alpha = 1*10^-3;
beta = 2;
kappa = 0:3:9;

figure
for i = 1:length(kappa)
    [Pzz_sut, mz_sut] = SUT(Pxx, mx, alpha, beta, kappa(i), n);
    ellipse_sut = sigmabound(Pzz_sut, mz_sut, 3);

    % plot pdf distribution and 3 sigma ellipses from from part D
    subplot(floor(sqrt(length(kappa))),ceil(sqrt(length(kappa))),i)
    pcolor(X, Y, pdf_cart_f)
    shading interp
    title(['\alpha = ',num2str(alpha),', \beta = ',num2str(beta),', \kappa = ',num2str(kappa(i))])
    xlabel('x')
    ylabel('y')
    hold on

    % plot 3 sigma ellipse over pdf distribution
    plot(ellipse_sut(1,:),ellipse_sut(2,:),'r','LineWidth',2)
    legend('p(z)','SUT','Location','southwest')
    sgtitle('Varying \kappa - Sean Bohne')
end

% function to calculate any sigma bound ellipse
function ellipse = sigmabound(P, m, sigma)
    theta = linspace(0,2*pi,500); % define theta values for 3 sigma plotting
    
    % guarantee eigenvectors listed in descending order
    [V,D] = eig(P);
    [d, idx] = sort(diag(D), 'descend');
    D = diag(d);
    V = V(:,idx);
    
    lambda1 = D(1,1); % first eigenvalue
    lambda2 = D(2,2); % second eigenvalue
    
    a = sigma*sqrt(lambda1); % ellipse semimajor axis
    b = sigma*sqrt(lambda2); % ellipse semiminor axis
    
    ellipse_prime = [a*cos(theta); b*sin(theta)]; % define 3 sigma ellipse in principal frame
    ellipse = V*ellipse_prime + m; % transform ellipse to original frame
end

% function to calculate scaled unscented transform approximation
function [Pzz, mz] = SUT(P, m, alpha, beta, kappa, n)
    Z = zeros(n, 2*n + 1); % initialize transformed sigma points array
    
    % scaled unscented transform weighting parameters
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
    
    % transform sigma points through measurement transformation
    for i = 1:(2*n + 1)
        % extract r and theta components
        r = sig_points(1,i);
        theta = sig_points(2,i);
    
        % calculate x and y coordinate transformations
        x = r*cos(theta);
        y = r*sin(theta);
    
        Z(:,i) = [x; y]; % combine transformed sigma points
    end
    
    % calculate SUT estimated transformed mean
    mz = w_0m*Z(:,1);
    
    for i = 2:(2*n + 1)
        mz = mz + w_im*Z(:,i);
    end
    
    % calculate SUT estimated transformed covariance
    Pzz = w_0c*(Z(:,1) - mz)*(Z(:,1) - mz)';
    
    for i = 2:(2*n + 1)
        Pzz = Pzz + w_ic*(Z(:,i) - mz)*(Z(:,i) - mz)';
    end
end