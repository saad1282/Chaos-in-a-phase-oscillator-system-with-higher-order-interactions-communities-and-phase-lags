%% Lyapunov exponents and bifurcation diagram using QR factorization and polar equations
% We are using polar equations instead of complex questions so that we do
% not have to deal with real and imaginary parts separately. We can use
% polar equation becauselyapunov exponents are invariant under transformations 
% (refer to Transformation invariance of Lyapunov exponents by Ralf Eichhorn, 
% Stefan J. Linz, and Peter Hänggi)


% The first part of this codes calculates the local maxima and minima for
% bifurcation diagram whereas the lower half calculates lyapuno exponents.

tic
clear
clc

gamma_vec = 1.12:0.0005:1.42; %phase lag

time = 1000; % transient time
dt = 0.001;  % step size
num_iterations = time/dt; % number of transient timesteps
Delta1 = 1; % Delta from lorentzian distribution
time1 = 500; % steady state time
num_iterations1 = time1/dt; % number of iterations for steady state

r1_diff = 0;
r2_diff = 0;
tol = 2*1e-02;
K1 = 10;
K2 = 10;
alpha = -0.5;
lyapunov_exponent = zeros(length(gamma_vec), 3);

num_max = 33;
num_min = 33;

r1_max = transpose(zeros(length(gamma_vec), num_max));
r1_min = transpose(zeros(length(gamma_vec), num_min));
r2_max = transpose(zeros(length(gamma_vec), num_max));
r2_min = transpose(zeros(length(gamma_vec), num_min));

gamma = gamma_vec(1);

%% Jacobian Matrix used for QR factorization

DF11 = @(r1,r2,phi, gamma)(-1 + ((1-r1^2)/2) * (K1 * cos(gamma) + 3 * K2 * r1^2 * cos(gamma)...
    + alpha * K2 * r2^2 * cos(2 * phi - gamma)+ 2 * alpha * K2 * r1 * r2 * cos(phi+gamma))...
    - r1 * (K1 * r1 * cos(gamma) + K2 * r1^3 *cos(gamma) + alpha * K1 * r2 * cos(phi-gamma)...
    + alpha * K2 * r1 * r2^2 * cos(2*phi-gamma)...
    + alpha * K2 * r1^2 * r2 * cos(phi+gamma) + alpha^2 * K2 * r2^3 * cos(phi-gamma)));

DF12 = @(r1,r2,phi, gamma)(((1-r1^2)/2) * (alpha * K1 * cos(phi-gamma) + 3 * alpha^2 * K2 * r2^2 * cos(phi-gamma)...
    + 2 * alpha * K2 * r1* r2 * cos(2 * phi - gamma)+ alpha * K2 * r1^2 * cos(phi+gamma)));

DF13 = @(r1,r2,phi, gamma)(((1-r1^2)/2)*(-2 * alpha * K2 * r1 * r2^2 * sin(2*phi - gamma)...
    - alpha * K1 * r2* sin(phi - gamma)...
    - alpha^2 * K2 * r2^3 * sin(phi - gamma)...
    - alpha * K2 * r1^2 * r2 * sin(phi + gamma)));

DF21 = @(r1,r2,phi, gamma)(((1-r2^2)/2) * (alpha * K1 * cos(phi+gamma) + 3 * alpha^2 * K2 * r1^2 * cos(phi+gamma)...
    + 2 * alpha * K2 * r1 * r2 * cos(2 * phi + gamma)+ alpha * K2 * r2^2 * cos(phi-gamma)));

DF22 = @(r1,r2,phi, gamma)(-1 + ((1-r2^2)/2) * (K1 * cos(gamma) + 3 * K2 * r2^2 * cos(gamma)...
    + alpha * K2 * r1^2 * cos(2 * phi + gamma)+ 2 * alpha * K2 * r1 * r2 * cos(phi-gamma))...
    - r2 * (K1 * r2 * cos(gamma) + K2 * r2^3 *cos(gamma) + alpha * K1 * r1 * cos(phi+gamma)...
    + alpha * K2 * r1 * r2^2 * cos(phi-gamma)...
    + alpha * K2 * r1^2 * r2 * cos(2*phi+gamma) + cos(phi+gamma) * alpha^2 * K2 * r1^3));

DF23 = @(r1,r2,phi, gamma)(((1-r2^2)/2)*(-2 * alpha * K2 * r1^2 * r2 * sin(2*phi + gamma)...
    -alpha * K1 * r1* sin(phi + gamma)...
    -alpha^2 * K2 * r1^3 * sin(phi + gamma)...
    - alpha * K2 * r1 * r2^2 * sin(phi - gamma)));

DF31 = @(r1,r2,phi, gamma)(r1*sin(gamma)*(K1 + K2 + 2 * K2 * r1^2)...
    +alpha * K2 * r1 * r2^2 * sin(-2*phi + gamma)...
    + (alpha * r2 * sin(-phi + gamma))/(2*r1^2) *...
    (K1 * (r1^2 -1) - alpha * K2 * r2^2 + K2 * r1^2 * (-1 + r2^2 *(alpha-1)))...
    -(alpha/(2 * r2)) *...
    ((-K2 * (1 + 3 * r1^2) * r2^2 + K1 * (1+ r2^2) + 3 * alpha * K2 * r1^2 * (1+ r2^2))* sin(phi + gamma)...
    + 2 * K2 * r1 * r2 * (1 + r2^2) * sin(2 * phi + gamma)));

DF32 = @(r1,r2,phi, gamma)(-r2*sin(gamma)*(K1 + K2 + 2 * K2 * r2^2)...
    + alpha * K2 * r2 * (1+ r1^2) * sin(-2*phi + gamma)...
    + (alpha/(2*r1)) * sin(-phi + gamma) * (K1* (1 + r1^2) +...
    3 * alpha * K2 * r2^2 + K2 * r1^2 * (-1 + 3 * r2^2 * (alpha-1)))...
    + ((alpha * r1)/(2 * r2^2))* (sin(gamma + phi) * (K1 - K1 * r2^2 ...
    + K2 * (1 + r1^2) * r2^2 - alpha * K2 * r1^2 * (r2^2 - 1))...
    - 2 * K2 * r1 * r2^3 * sin(2*phi + gamma)));

DF33 = @(r1, r2, phi, gamma)((-alpha/(2 * r1 * r2)) * (2 * K2 * r1 * (1+ r1^2) * r2^3 * cos(2 * phi - gamma)...
    + r2^2 * cos(phi - gamma) * (K1 * (1 + r1^2) + alpha * K2 * r2^2 + K2 * r1^2 *(-1 + r2^2 * (alpha-1)))...
    + r1^2 * (cos(phi + gamma) * (-K2 * r2^2 * (1 + r1^2) + K1 * (1 + r2^2)...
    + alpha * K2 * r1^2 * (1 + r2^2))...
    + 2 * K2 * r1 * r2 * (1 + r2^2) * cos(2 * phi + gamma))));

DF = @(r1, r2, phi, gamma) ([DF11(r1, r2, phi, gamma) DF12(r1, r2, phi, gamma) DF13(r1, r2, phi, gamma);
    DF21(r1, r2, phi, gamma) DF22(r1, r2, phi, gamma) DF23(r1, r2, phi, gamma);
    DF31(r1, r2, phi, gamma) DF32(r1, r2, phi, gamma) DF33(r1, r2, phi, gamma)]);

r1_ini = rand();
r2_ini = rand();
phi_ini = pi*rand();


for s = 1: length(gamma_vec);
    
    s
    gamma = gamma_vec(s);
    r1_vec = zeros(num_iterations1,1);
    r2_vec = zeros(num_iterations1,1);
    phi_vec = zeros(num_iterations1,1);

    % transient solution
    for i = 1 : num_iterations
        [r1ans_bs, r2ans_bs, phi_ansbs] =  polar_bs(r1_ini, r2_ini, phi_ini, gamma);
        r1_ini = r1ans_bs;
        r2_ini = r2ans_bs;
        phi_ini = phi_ansbs;
    end

    % Extra Iteration

    r11 = r1_ini;
    r21 = r2_ini;
    phi1 = phi_ini;

    [r1new, r2new, phinew] =  polar_bs(r11, r21, phi1, gamma);

    cont = 1;
    index1 = 1;
    index2 = 1;
    index3 = 1;
    index4 = 1;

    % one new iteration to check if there is min/max

    while (cont)

        [r1new1, r2new1, phinew1] =  polar_bs(r1new, r2new, phinew, gamma);

        if ((r11 <= r1new) && (r1new >= r1new1))
            if (index1 <= num_max+1)
                r1_max(index1,s) = r1new;
                index1 = index1 + 1;
            end
        end

        if ((r11 >= r1new) && (r1new <= r1new1))
            if (index2 <= num_max+1)
                r1_min(index2,s) = r1new;
                index2 = index2 + 1;
            end
        end

        if ((r21 <= r2new) && (r2new >= r2new1))
            if (index3 <= num_max+1)
                r2_max(index3,s) = r2new;
                index3 = index3 + 1;
            end
        end

        if ((r21 >= r2new) && (r2new <= r2new1))
            if (index4 <= num_max+1)
                r2_min(index4,s) = r2new;
                index4 = index4 + 1;
            end
        end

        r11 = r1new;
        r1new = r1new1;
        r21 = r2new;
        r2new = r2new1;
        phi1 = phinew;
        phinew = phinew1;

        if (((index1) >= num_max) && ((index2) >= num_min) && ((index3) >= num_max) && ((index4) >= num_min))
            cont = 0;
        end
    end

    lyapunov_iter = zeros(num_iterations1,3);
    U = eye(3);

    %steady state & lyapunov exponents

    for t = 1 : num_iterations1

        [r1ans_bs, r2ans_bs, phi_ansbs] =  polar_bs(r1_ini, r2_ini, phi_ini, gamma);
        r1_vec(t) = r1ans_bs;
        r2_vec(t) = r2ans_bs;
        phi_vec(t) = phi_ansbs;
        r1_ini = r1ans_bs;
        r2_ini = r2ans_bs;
        phi_ini = phi_ansbs;

        U_n = (eye(3) + DF(r1_vec(t), r2_vec(t), phi_vec(t), gamma)*dt) * U; % This perturbs the unit sphere using Jacbobian matrix
        [Q,R] = qr(U_n); % QR factorization of perturbed spheres
        U = Q; % Directions of perturbed shape
        lyapunov_iter(t,:) = log(abs(diag(R))); % Growth rate of the object in each direction 

    end

    lyapunov_exponent(s, :) = sum(lyapunov_iter,1)/(dt*num_iterations1); % Averages the growrth rates for long time scales

end

figure(1);
plot(gamma_vec, lyapunov_exponent(:,1), 'b.-');
hold on;
plot(gamma_vec, lyapunov_exponent(:,2), 'g.-');
plot(gamma_vec, lyapunov_exponent(:,3), 'k.-');

figure(2);
plot(gamma_vec, r1_max(1:end-2,:), 'r.');
hold on;
plot(gamma_vec, r1_min(1:end-2,:), 'r.');
plot(gamma_vec, r2_max(1:end-2,:), 'b.');
plot(gamma_vec, r2_min(1:end-2,:), 'b.');
grid on;
grid minor;
ylabel('Min/ Max r_1, r_2');
xlabel('\gamma')


figure(2);
plot(gamma_vec, r1_max, 'r.');
hold on;
plot(gamma_vec, r1_min, 'b.');
xlabel('phase lags');
ylabel('min/max r_1');
title('local min/max r1 for community 1');
legend('local max', 'local min');

figure(3);
plot(gamma_vec, r2_max, 'r.');
hold on;
plot(gamma_vec, r2_min, 'b.');
xlabel('phase lags');
ylabel('min/max r_2');
title('local min/max r2 for community 2');
legend('local max', 'local min');

figure(4);
subplot(2,1,1);
plot(gamma_vec, lyapunov_exponent(:,1), 'b.-');
hold on;
plot(gamma_vec, lyapunov_exponent(:,2), 'g.-');
plot(gamma_vec, lyapunov_exponent(:,3), 'k.-');

subplot(2,1,2);
plot(gamma_vec, r1_max, 'r.');
xlabel('phase lags');
ylabel('min/max r_1');
title('local min/max r1 for community 1');
legend('local max', 'local min');

toc

%save('/directory/filename.mat', 'gamma_vec', 'lyapunov_exponent', 'K1', 'K2', 'alpha', 'eps', 'dt', 'r1_min', 'r1_max', 'r2_max', 'r2_min', 'num_max')

function[r1ans_bs, r2ans_bs, phi_ansbs] =  polar_bs(r1bs_ini, r2bs_ini, phibs_ini, gamma1)

K1 = 10;
K2 = 10;
alpha = -0.5;
Delta1 = 1;
dt = 0.001;

h1 = @(r1, r2, phi, gamma1) (-Delta1*r1 + ((1-r1^2)/2)*(K1*r1*cos(gamma1) + alpha*K1*r2*cos(phi-gamma1)...
    + r1^3 * K2 * cos(gamma1) + alpha * r1^2 * r2 * K2* cos(phi + gamma1) + alpha * r1* r2^2 * K2* cos(2*phi - gamma1)...
    + alpha^2 * K2 * r2^3 *cos(phi - gamma1)));

h2 = @(r1, r2, phi, gamma1) (-Delta1*r2 + ((1-r2^2)/2)*(K1*r2*cos(gamma1) + alpha*K1*r1*cos(phi+gamma1)...
    + r2^3 * K2 * cos(gamma1) + alpha *K2 *r1^2 * r2 * cos(2*phi + gamma1) + alpha *K2 * r1* r2^2 * cos(phi - gamma1)...
    + alpha^2 * K2 * r1^3 *cos(phi + gamma1)));

h3 = @(r1, r2, phi, gamma1) (((1+r2^2)/(2*r2))*(K1 * r2* sin(-gamma1) + alpha * K1 * r1* sin(-(phi+gamma1))...
    + K2* r2^3 * sin(-gamma1) + alpha* K2 * r1^2 *r2* sin(-(2*phi + gamma1))...
    + alpha * K2 * r1 * r2^2 * sin(phi -gamma1)+ alpha^2 * K2 * r1^3 * sin(-(phi+gamma1)))...
    -((1+r1^2)/(2*r1))*(K1* r1* sin(-gamma1) + alpha * K1* r2* sin(phi - gamma1)+ K2 * r1^3 * sin(-gamma1)...
    + alpha * K2 * r1^2 * r2* sin(-(phi+gamma1))...
    + alpha * K2 * r1* r2^2 * sin(2*phi - gamma1) + alpha^2 * K2 * r2^3 * sin(phi-gamma1)));

r1bs_k1 = h1(r1bs_ini,r2bs_ini, phibs_ini, gamma1);
r2bs_k1 = h2(r1bs_ini,r2bs_ini, phibs_ini, gamma1);
phibs_k1 = h3(r1bs_ini,r2bs_ini, phibs_ini, gamma1);

r1bs_int = r1bs_ini + (dt)*r1bs_k1;
r2bs_int = r2bs_ini + (dt)*r2bs_k1;
phibs_int = phibs_ini + (dt)*phibs_k1;

r1bs_k2 = h1(r1bs_int,r2bs_int, phibs_int, gamma1);
r2bs_k2 = h2(r1bs_int,r2bs_int, phibs_int, gamma1);
phibs_k2 = h3(r1bs_int,r2bs_int, phibs_int, gamma1);

r1ans_bs = r1bs_ini + (dt/2)*(r1bs_k1+r1bs_k2);
r2ans_bs = r2bs_ini + (dt/2)*(r2bs_k1+r2bs_k2);
phi_ansbs = phibs_ini + (dt/2)*(phibs_k1+phibs_k2);
end