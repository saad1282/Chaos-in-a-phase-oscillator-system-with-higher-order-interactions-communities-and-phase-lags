% This file generates timeseries plots of order parameters 
% using complex order parameters equations and polar equations. This file
% uses RK4 method. Heuns method is commented at the end, but included if
% one prefers to use Heuns method for better computation speed.

clear all
clc

gamma = 1.328; % phase lag value
gamma1 = gamma; % phase lag value for links
gamma2 = gamma; % phase lag value for triangles

time = 2000; % transient time
dt = 0.001; % timestep
num_iterations = time/dt; % transient iterations

time1 = 700;    % non-transient time
num_iterations1 = time1/dt; % steady state timesteps
r1_steady = zeros(num_iterations1,1);   % vector initialization for steady state r1 using complex equations
r2_steady = zeros(num_iterations1,1);   % vector initialization for steady state r2 using complex equations
phi_steady = zeros(num_iterations1,1);  % vector initialization for steady state phi using complex equations

r1_polar = zeros(num_iterations1,1);    % vector initialization for steady state r1 using polar equations
r2_polar = zeros(num_iterations1,1);    % vector initialization for steady state r2 using polar equations
phi_polar = zeros(num_iterations1,1);   % vector initialization for steady state phi using polar equations

num_ini_condition = 1; % number of initial conditions for r1, r2, phi

for i = 1 : num_ini_condition

    phi_ini = pi*rand();    % initial angle
    r1_ini = rand();    % initial r1
    r2_ini = rand();    % initial r2

    z1 = r1_ini;    % initial z1 (order parameter in complex form)
    z2 = r2_ini*exp(1i*phi_ini);    % initial z2 (order parameter in complex form)

    % r1bs_ini = r1_ini;
    % r2bs_ini = r2_ini;
    % phibs_ini = phi_ini;

    % transient simulation
    for t = 1 : num_iterations
        [z1ans, z2ans] = RK4(z1 , z2, gamma1, gamma2);  %complex equation simulation
        z1 = z1ans;
        z2 = z2ans;
        r1_complex_vec(t) = norm(z1);
        r2_complex_vec(t) = norm(z2);

        [r1ans, r2ans, phi_ans] = polar(r1_ini,r2_ini, phi_ini, gamma1, gamma2);    % polar equation simulation
        r1_ini = r1ans;
        r2_ini = r2ans;
        phi_ini = mod(phi_ans,2*pi);
        r1_vec(t) = r1ans;
        r2_vec(t) = r2ans;
        phi_vec(t) = phi_ans;
    end

    %steady state simulation
    for t1 = 1 : num_iterations1
        [z1ans, z2ans] = RK4(z1,z2, gamma1, gamma2);    % complex equation simulation
        r1_steady(t1) = norm(z1ans);
        r2_steady(t1) = norm(z2ans);
        phi_steady(t1) = mod(angle(z2ans) - angle(z1ans), 2*pi);
        z1 = z1ans;
        z2 = z2ans;

        [r1ans, r2ans, phi_ans] = polar(r1_ini,r2_ini, phi_ini, gamma1, gamma2);    % polar equation simulation

        r1_polar(t1) = r1ans;
        r2_polar(t1) = r2ans;
        phi_polar(t1) = mod(phi_ans, 2*pi);      
        r1_ini = r1ans;
        r2_ini = r2ans;
        phi_ini = phi_ans;
     end
end

% phase attractor plots
figure(1)
plot((r1_polar).*cos(phi_polar), (r2_polar).*sin(phi_polar));
xlabel('r1 \cdot cos{\Phi}')

% timeseries plot of r1 and r2 using polar equations
figure(2)
plot(dt*(num_iterations1-100000:num_iterations1), r1_polar(end-100000:end), 'b.', MarkerSize= 8);
hold on;
plot(dt*(num_iterations1-100000:num_iterations1), r2_polar(end-100000:end), 'r.', MarkerSize= 8);
xlabel('Time', 'FontSize', 22);
ylabel('r_1, r_2', 'FontSize', 22);
ylim([0 1])
grid on
grid minor
axes('Position', [0.5 0.7 0.2 0.2])
box on
plot(dt*(num_iterations1-30000:num_iterations1), (phi_polar(end-30000:end)), 'b.', MarkerSize=4)
xlabel('Time', 'FontSize',14)
ylabel('\Phi', 'FontSize',14)
hold off;

% timeseries plot of r1 and r2 using complex equations
figure(3)
plot(dt*(num_iterations1-100000:num_iterations1), r1_steady(end-100000:end), 'b.', MarkerSize= 8);
hold on;
plot(dt*(num_iterations1-100000:num_iterations1), r2_steady(end-100000:end), 'r.', MarkerSize= 8);
xlabel('Time', 'FontSize', 22);
ylabel('r_1, r_2', 'FontSize', 22);
grid on
grid minor

% simulation using complex equations
function[z1ans, z2ans] =  RK4(z1_ini, z2_ini, gamma1, gamma2)

K1 = 10;    % pairwise coupling strength
K2 = 10;    % triangular coupling strength
alpha = 0.5;    % parameter that determines positive and negative coupling
Delta1 = 1; % Delta in Lorentzian distribution of frequencies
dt = 0.001; % time step

% complex equations

f1 = @(z1,z2, gamma1, gamma2) (K1 * exp(-1i * gamma1) * (z1 + alpha * z2)...
    + K2 * exp(-1i * gamma2) * (z1^2 * conj(z1) + alpha * (z1)^2 * conj(z2) + alpha * (z2)^2 * conj(z1) + alpha^2 * (z2)^2 * conj(z2)));

f2 = @(z1,z2, gamma1, gamma2) (K1 * exp(-1i * gamma1) * (z2 + alpha * z1)...
    + K2 * exp(-1i * gamma2) * (z2^2 * conj(z2) + alpha * (z2)^2 * conj(z1) + alpha * (z1)^2 * conj(z2) + alpha^2 * (z1)^2 * conj(z1)));

H1 = f1(z1_ini,z2_ini,gamma1,gamma2);
H2 = f2(z1_ini,z2_ini,gamma1,gamma2);

z1_k1 = -Delta1 * z1_ini + (1/2) * (H1 - conj(H1)*z1_ini^2); % k1 RK4
z2_k1 = -Delta1 * z2_ini + (1/2) * (H2 - conj(H2)*z2_ini^2); % K1 RK4

z1_int1 = z1_ini + (dt/2)*z1_k1;
z2_int1 = z2_ini + (dt/2)*z2_k1;

H11 = f1(z1_int1,z2_int1,gamma1,gamma2);
H21 = f2(z1_int1,z2_int1,gamma1,gamma2);

z1_k2 = -Delta1 * z1_int1 + (1/2) * (H11 - conj(H11)*z1_int1^2);  %k2 RK4
z2_k2 = -Delta1 * z2_int1 + (1/2) * (H21 - conj(H21)*z2_int1^2);  %k2 RK4

z1_int2 = z1_ini + (dt/2)*z1_k2;
z2_int2 = z2_ini + (dt/2)*z2_k2;

H12 = f1(z1_int2,z2_int2,gamma1,gamma2);
H22 = f2(z1_int2,z2_int2,gamma1,gamma2);

z1_k3 = -Delta1 * z1_int2 + (1/2) * (H12 - conj(H12)*z1_int2^2);  %k3 RK4
z2_k3 = -Delta1 * z2_int2 + (1/2) * (H22 - conj(H22)*z2_int2^2);  %k3 RK4

z1_int3 = z1_ini + (dt)*z1_k3;
z2_int3 = z2_ini + (dt)*z2_k3;

H13 = f1(z1_int3,z2_int3,gamma1,gamma2);
H23 = f2(z1_int3,z2_int3,gamma1,gamma2);

z1_k4 = -Delta1 * z1_int3 + (1/2) * (H13 - conj(H13)*z1_int3^2);  %k4 RK4
z2_k4 = -Delta1 * z2_int3 + (1/2) * (H23 - conj(H23)*z2_int3^2);  %k4 RK4

z1ans = z1_ini + (dt/6)*(z1_k1 + 2*z1_k2 + 2*z1_k3 + z1_k4);
z2ans = z2_ini + (dt/6)*(z2_k1 + 2*z2_k2 + 2*z2_k3 + z2_k4);
end

% simulation using polar equations
function[r1ans, r2ans, phi_ans] =  polar(r1_ini, r2_ini, phi_ini, gamma1, gamma2)

K1 = 10; % pairwise coupling strength
K2 = 10;    % triangular coupling strength
alpha = 0.5;    % parameter that determines positive and negative coupling strength
Delta1 = 1; % Delta in Lorentzian distribution of frequencies 
dt = 0.001; % time step

% polar equations
g1 = @(r1,r2, phi, gamma1, gamma2)(-Delta1 * r1 + (1/2)* ...
    real(exp(-1i*gamma1)*(K1 * r1 + alpha * K1 * r2 * exp(1i*phi)...
    + K2 * r1^3 + alpha * K2 * r1^2 * r2 * exp(-1i*phi) + alpha * K2 *r1 *r2^2 * exp(2*1i*phi) + alpha^2 * K2 * r2^3 * exp(1i*phi))...
    - exp(1i*gamma1) * (K1 * r1^3 + alpha * K1 * r1^2 * r2 * exp(-1i*phi)+ K2 * r1^5 + alpha * K2 * r1^4 * r2 * exp(1i*phi)...
    + alpha * K2 * r1^3 * r2^2 * exp(-2*1i*phi) + alpha^2 * K2 * r1^2 * r2^3 * exp(-1i*phi)))); 

g2 = @(r1,r2, phi, gamma1, gamma2)(-Delta1 * r2 + (1/2) * real(exp(-1i*gamma1)*(K1 * r2 + alpha * K1 * r1 * exp(-1i*phi)...
    + K2 * r2^3 + alpha * K2 * r2^2 * r1 * exp(1i*phi) +alpha * K2* r2 * r1^2 * exp(-2*1i*phi) + alpha^2 * K2 * r1^3 * exp(-1i*phi))...
    -exp(1i*gamma1)*(K1 * r2^3 + alpha * K1 * r2^2 * r1 * exp(1i*phi)+ K2 * r2^5 + alpha * K2 * r2^4 * r1 * exp(-1i*phi)...
    + alpha * K2 *r2^3 *r1^2 * exp(2*1i*phi) + alpha^2 * K2 * r2^2 * r1^3 * exp(1i*phi))));

g3 = @(r1,r2, phi, gamma1, gamma2)((1/(2*r2)) * imag(exp(-1i*gamma1) * (K1 * r2 + alpha * K1 *r1 * exp(-1i*phi)...
    + K2 * r2^3 + alpha * K2 * r2^2 * r1 * exp(1i*phi) +alpha * K2 * r2 * r1^2 * exp(-2*1i*phi) + alpha^2 * K2 * r1^3 * exp(-1i*phi))...
    -exp(1i*gamma1) * (K1 * r2^3 + alpha * K1 * r2^2 * r1 * exp(1i*phi)+ K2 * r2^5 + alpha * K2 * r2^4 * r1*exp(-1i*phi)...
    + alpha * K2 * r2^3 *r1^2 * exp(2*1i*phi) + alpha^2 * K2 * r2^2 * r1^3 * exp(1i*phi)))...
    - (1/(2*r1))*imag(exp(-1i*gamma1) * (K1 * r1 + alpha * r2 * K1 * exp(1i*phi)...
    + K2 * r1^3 + alpha * K2 * r1^2 * r2 * exp(-1i*phi) +alpha * K2 * r1 * r2^2 * exp(2*1i*phi) +alpha^2 * K2 * r2^3 * exp(1i*phi))...
    -exp(1i*gamma1)*((K1 * r1^3 + alpha * K1 * r1^2 * r2 * exp(-1i*phi)+ K2 * r1^5 + alpha * K2 * r1^4 * r2*exp(1i*phi)...
    + alpha *K2 *r1^3 *r2^2 * exp(-2*1i*phi) + alpha^2 * K2 * r1^2 * r2^3 * exp(-1i*phi)))));

r1_k1 = g1(r1_ini,r2_ini, phi_ini, gamma1, gamma2); %k1 RK4
r2_k1 = g2(r1_ini,r2_ini, phi_ini, gamma1, gamma2); %k1 RK4
phi_k1 = g3(r1_ini,r2_ini, phi_ini, gamma1, gamma2);    %k1 RK4    

r1_int = r1_ini + (dt/2)*r1_k1;
r2_int = r2_ini + (dt/2)*r2_k1;
phi_int = phi_ini + (dt/2)*phi_k1;

r1_k2 = g1(r1_int,r2_int, phi_int, gamma1, gamma2); %k2 RK4
r2_k2 = g2(r1_int,r2_int, phi_int, gamma1, gamma2); %k2 RK4
phi_k2 = g3(r1_int,r2_int, phi_int, gamma1, gamma2);    %k2 RK4

r1_int1 = r1_ini + (dt/2)*r1_k2;
r2_int1 = r2_ini + (dt/2)*r2_k2;
phi_int1 = phi_ini + (dt/2)*phi_k2;

r1_k3 = g1(r1_int1,r2_int1, phi_int1, gamma1, gamma2);  % k3 RK4
r2_k3 = g2(r1_int1,r2_int1, phi_int1, gamma1, gamma2);  % k3 RK4
phi_k3 = g3(r1_int1,r2_int1, phi_int1, gamma1, gamma2); % k3 RK4

r1_int2 = r1_ini + (dt)*r1_k3;
r2_int2 = r2_ini + (dt)*r2_k3;
phi_int2 = phi_ini + (dt)*phi_k3;

r1_k4 = g1(r1_int2,r2_int2, phi_int2, gamma1, gamma2);  %k4 RK4
r2_k4 = g2(r1_int2,r2_int2, phi_int2, gamma1, gamma2);  %k4 RK4
phi_k4 = g3(r1_int2,r2_int2, phi_int2, gamma1, gamma2); %k4 RK4

r1ans = r1_ini + (dt/6)*(r1_k1 + 2* r1_k2 + 2* r1_k3 + r1_k4);  
r2ans = r2_ini + (dt/6)*(r2_k1 + 2* r2_k2 + 2* r2_k3 + r2_k4);
phi_ans = phi_ini + (dt/6)*(phi_k1 + 2* phi_k2 + 2* phi_k3 + phi_k4);
end

% %% In case we want to use Heuns Method
% 
% % simulation using complex equations
% function[z1ans, z2ans] =  Heuns_complex(z1_ini, z2_ini, gamma1, gamma2)
% 
% K1 = 10;    % pairwise coupling strength
% K2 = 10;    % triangular coupling strength
% alpha = 0.5;    % parameter that determines positive and negative coupling
% Delta1 = 1; % Delta in Lorentzian distribution of frequencies
% dt = 0.001; % time step
% 
% % complex equations
% 
% f1 = @(z1,z2, gamma1, gamma2) (K1 * exp(-1i * gamma1) * (z1 + alpha * z2)...
%     + K2 * exp(-1i * gamma2) * (z1^2 * conj(z1) + alpha * (z1)^2 * conj(z2) + alpha * (z2)^2 * conj(z1) + alpha^2 * (z2)^2 * conj(z2)));
% 
% f2 = @(z1,z2, gamma1, gamma2) (K1 * exp(-1i * gamma1) * (z2 + alpha * z1)...
%     + K2 * exp(-1i * gamma2) * (z2^2 * conj(z2) + alpha * (z2)^2 * conj(z1) + alpha * (z1)^2 * conj(z2) + alpha^2 * (z1)^2 * conj(z1)));
% 
% H1 = f1(z1_ini,z2_ini,gamma1,gamma2);
% H2 = f2(z1_ini,z2_ini,gamma1,gamma2);
% 
% z1_k1 = -Delta1 * z1_ini + (1/2) * (H1 - conj(H1)*z1_ini^2); 
% z2_k1 = -Delta1 * z2_ini + (1/2) * (H2 - conj(H2)*z2_ini^2); 
% 
% z1_int1 = z1_ini + (dt)*z1_k1;
% z2_int1 = z2_ini + (dt)*z2_k1;
% 
% H11 = f1(z1_int1,z2_int1,gamma1,gamma2);
% H21 = f2(z1_int1,z2_int1,gamma1,gamma2);
% 
% z1_k2 = -Delta1 * z1_int1 + (1/2) * (H11 - conj(H11)*z1_int1^2);  %k2 RK4
% z2_k2 = -Delta1 * z2_int1 + (1/2) * (H21 - conj(H21)*z2_int1^2);  %k2 RK4
% 
% z1ans = z1_ini + (dt/2)*(z1_k1 + z1_k2);
% z2ans = z2_ini + (dt/2)*(z2_k1 + z2_k2);
% end
% 
% % simulation using polar equations
% function[r1ans, r2ans, phi_ans] =  Heuns_polar(r1_ini, r2_ini, phi_ini, gamma1, gamma2)
% 
% K1 = 10; % pairwise coupling strength
% K2 = 10;    % triangular coupling strength
% alpha = 0.5;    % parameter that determines positive and negative coupling strength
% Delta1 = 1; % Delta in Lorentzian distribution of frequencies 
% dt = 0.001; % time step
% 
% % polar equations
% g1 = @(r1,r2, phi, gamma1, gamma2)(-Delta1 * r1 + (1/2)* ...
%     real(exp(-1i*gamma1)*(K1 * r1 + alpha * K1 * r2 * exp(1i*phi)...
%     + K2 * r1^3 + alpha * K2 * r1^2 * r2 * exp(-1i*phi) + alpha * K2 *r1 *r2^2 * exp(2*1i*phi) + alpha^2 * K2 * r2^3 * exp(1i*phi))...
%     - exp(1i*gamma1) * (K1 * r1^3 + alpha * K1 * r1^2 * r2 * exp(-1i*phi)+ K2 * r1^5 + alpha * K2 * r1^4 * r2 * exp(1i*phi)...
%     + alpha * K2 * r1^3 * r2^2 * exp(-2*1i*phi) + alpha^2 * K2 * r1^2 * r2^3 * exp(-1i*phi)))); 
% 
% g2 = @(r1,r2, phi, gamma1, gamma2)(-Delta1 * r2 + (1/2) * real(exp(-1i*gamma1)*(K1 * r2 + alpha * K1 * r1 * exp(-1i*phi)...
%     + K2 * r2^3 + alpha * K2 * r2^2 * r1 * exp(1i*phi) +alpha * K2* r2 * r1^2 * exp(-2*1i*phi) + alpha^2 * K2 * r1^3 * exp(-1i*phi))...
%     -exp(1i*gamma1)*(K1 * r2^3 + alpha * K1 * r2^2 * r1 * exp(1i*phi)+ K2 * r2^5 + alpha * K2 * r2^4 * r1 * exp(-1i*phi)...
%     + alpha * K2 *r2^3 *r1^2 * exp(2*1i*phi) + alpha^2 * K2 * r2^2 * r1^3 * exp(1i*phi))));
% 
% g3 = @(r1,r2, phi, gamma1, gamma2)((1/(2*r2)) * imag(exp(-1i*gamma1) * (K1 * r2 + alpha * K1 *r1 * exp(-1i*phi)...
%     + K2 * r2^3 + alpha * K2 * r2^2 * r1 * exp(1i*phi) +alpha * K2 * r2 * r1^2 * exp(-2*1i*phi) + alpha^2 * K2 * r1^3 * exp(-1i*phi))...
%     -exp(1i*gamma1) * (K1 * r2^3 + alpha * K1 * r2^2 * r1 * exp(1i*phi)+ K2 * r2^5 + alpha * K2 * r2^4 * r1*exp(-1i*phi)...
%     + alpha * K2 * r2^3 *r1^2 * exp(2*1i*phi) + alpha^2 * K2 * r2^2 * r1^3 * exp(1i*phi)))...
%     - (1/(2*r1))*imag(exp(-1i*gamma1) * (K1 * r1 + alpha * r2 * K1 * exp(1i*phi)...
%     + K2 * r1^3 + alpha * K2 * r1^2 * r2 * exp(-1i*phi) +alpha * K2 * r1 * r2^2 * exp(2*1i*phi) +alpha^2 * K2 * r2^3 * exp(1i*phi))...
%     -exp(1i*gamma1)*((K1 * r1^3 + alpha * K1 * r1^2 * r2 * exp(-1i*phi)+ K2 * r1^5 + alpha * K2 * r1^4 * r2*exp(1i*phi)...
%     + alpha *K2 *r1^3 *r2^2 * exp(-2*1i*phi) + alpha^2 * K2 * r1^2 * r2^3 * exp(-1i*phi)))));
% 
% r1_k1 = g1(r1_ini,r2_ini, phi_ini, gamma1, gamma2); 
% r2_k1 = g2(r1_ini,r2_ini, phi_ini, gamma1, gamma2); 
% phi_k1 = g3(r1_ini,r2_ini, phi_ini, gamma1, gamma2);       
% 
% r1_int = r1_ini + (dt)*r1_k1;
% r2_int = r2_ini + (dt)*r2_k1;
% phi_int = phi_ini + (dt)*phi_k1;
% 
% r1_k2 = g1(r1_int,r2_int, phi_int, gamma1, gamma2); 
% r2_k2 = g2(r1_int,r2_int, phi_int, gamma1, gamma2);
% phi_k2 = g3(r1_int,r2_int, phi_int, gamma1, gamma2);   
% 
% r1ans = r1_ini + (dt/2)*(r1_k1 + r1_k2);
% r2ans = r2_ini + (dt/2)*(r2_k1 + r2_k2);
% phi_ans = phi_ini + (dt/2)*(phi_k1 + phi_k2);
% end