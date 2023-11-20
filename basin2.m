% This file finds the basin of fixed point solution and periodic/chaotic
% solution. It takes the pairwise and triangular coupling strengths and
% phase lags value. Change the phase lags value to find the basin of
% attraction of periodic and chaotic solution. We run the simulation to
% find the steady state order parameters and if the difference between
% minimum and maximum values of the order parameters satisfy a threshold of
% 0.02, we conclude that it is a fixed point. If the difference is larger
% than the tolerance, we conclude that it is a periodic or chaotic solution
% depending on the phase lag values.

%% parameters and variables

K1 = 10;    % pairwise coupling strength
K2 = 10;    % triangular coupling strength
alpha = -0.5;   % parameter that determines positive/negative coupling strength
gamma = 1.406; % phase lag
gamma1 = gamma;
gamma2 = gamma;

time = 700; % transient time
dt = 0.001; % time step
num_iterations = time/dt;   % number of iterations for transient time
Delta1 = 1; % Delta from lorentzian distribution of frequencies

time1 = 700; % non transient time
num_iterations1 = time1/dt; % steady state time steps

%% ODES for order parameters for two communities z1 and z2

f1 = @(z1,z2, gamma1, gamma2) (K1 * exp(-1i * gamma1) * (z1 + alpha * z2)...
    + K2 * exp(-1i * gamma2) * (z1^2 * conj(z1) + alpha * (z1)^2 * conj(z2) + alpha * (z2)^2 * conj(z1) + alpha^2 * (z2)^2 * conj(z2)));

f2 = @(z1,z2, gamma1, gamma2) (K1 * exp(-1i * gamma1) * (z2 + alpha * z1)...
    + K2 * exp(-1i * gamma2) * (z2^2 * conj(z2) + alpha * (z2)^2 * conj(z1) + alpha * (z1)^2 * conj(z2) + alpha^2 * (z1)^2 * conj(z1)));

%% initial conditions for basin: Need to know that z1 = r1, z2 = r2* exp(1i * phi)

grid_num = 150; % total number of points on r1 and r2 direction

r1_ini = 0.01 : 1/grid_num : 1; % r1 initial values
r2_ini = 0.01 : 1/grid_num : 1; % r2 initial values

num_ini_condition = length(r1_ini)*length(r2_ini);  % total number of initial conditions
matrix = zeros(num_ini_condition,4); % initialization matrix with first column-r1_initial, second column-r2_initial, third column-phi_initial and fourth column- entries for basin 
angle_ini = pi*rand(); % initial condition for phi. we are using same phi for different r1 r2 values.
index_ini = 1; % initialization index for matrix
r1_value = zeros(num_iterations1/2,1);  % vector initialization for steady state solution
r2_value = zeros(num_iterations1/2,1);  % vector initialization for steady state solution

% fill in the entires of the matrix
for j = 1 : length(r1_ini)
    for m = 1 : length(r2_ini)
        matrix(index_ini,1) = r1_ini(j);
        matrix(index_ini,2) = r2_ini(m);
        matrix(index_ini,3) = angle_ini;
        index_ini = index_ini +1;
    end
end

% temporary variable for parfor loop
basin_ans = zeros(num_ini_condition,1);

% tolerance for fixed point solution
tol_variance = 2*1e-02;

%r1 = 0;
%r2 = 0;
%phi = 0;

% broadcast variable for parfor loop
r1_ini_vec = matrix(:,1);
r2_ini_vec = matrix(:,2);
phi_ini_vec = matrix(:,3);

% using 8 workers for par for loop
parfor (i = 1 : num_ini_condition, 8)

    r1 = r1_ini_vec(i);
    r2 = r2_ini_vec(i);
    phi = phi_ini_vec(i);

    z1 = r1;
    z2 = r2*exp(1i*phi);

    r1_value = zeros(num_iterations/2,1);
    r2_value = zeros(num_iterations/2,1);

    % transient solution using Heuns Method

    for t = 1 : num_iterations

        H1 = f1(z1,z2,gamma1,gamma2);
        H2 = f2(z1,z2,gamma1,gamma2);

        z1_k1 = -Delta1 * z1 + (1/2) * (H1 - conj(H1)*z1^2); 
        z2_k1 = -Delta1 * z2 + (1/2) * (H2 - conj(H2)*z2^2); 

        z1_int1 = z1 + (dt)*z1_k1;
        z2_int1 = z2 + (dt)*z2_k1;

        H11 = f1(z1_int1,z2_int1,gamma1,gamma2);
        H21 = f2(z1_int1,z2_int1,gamma1,gamma2);

        z1_k2 = -Delta1 * z1_int1 + (1/2) * (H11 - conj(H11)*z1_int1^2);  
        z2_k2 = -Delta1 * z2_int1 + (1/2) * (H21 - conj(H21)*z2_int1^2);  

        z1 = z1 + (dt/2)*(z1_k1 + z1_k2);
        z2 = z2 + (dt/2)*(z2_k1 + z2_k2);

    end
    
    % first we run the simulation for half of steady state timesteps to see
    % if we get fixed point, otherwise we will run the simulation for rest
    % of the steady state timesteps and check if we have fixed point.

    for t1 = 1 : num_iterations1/2

        H1 = f1(z1,z2,gamma1,gamma2);
        H2 = f2(z1,z2,gamma1,gamma2);

        z1_k1 = -Delta1 * z1 + (1/2) * (H1 - conj(H1)*z1^2); 
        z2_k1 = -Delta1 * z2 + (1/2) * (H2 - conj(H2)*z2^2); 

        z1_int1 = z1 + (dt)*z1_k1;
        z2_int1 = z2 + (dt)*z2_k1;

        H11 = f1(z1_int1,z2_int1,gamma1,gamma2);
        H21 = f2(z1_int1,z2_int1,gamma1,gamma2);

        z1_k2 = -Delta1 * z1_int1 + (1/2) * (H11 - conj(H11)*z1_int1^2);  
        z2_k2 = -Delta1 * z2_int1 + (1/2) * (H21 - conj(H21)*z2_int1^2);  

        z1 = z1 + (dt/2)*(z1_k1 + z1_k2);
        z2 = z2 + (dt/2)*(z2_k1 + z2_k2);

        r1_value(t1) = norm(z1);
        r2_value(t1) = norm(z2);

    end

    if ((max(r1_value)-min(r1_value)<tol_variance) && (max(r2_value)-min(r2_value)<tol_variance))
        basin = 0; % store fixed point as 0
    else % if not a fixed point, run rest of the simulation
        r1_value = zeros(num_iterations/2,1);
        r2_value = zeros(num_iterations/2,1);
        for t1 = 1 : num_iterations1/2

            r1_value = zeros(num_iterations1/2,1);
            r2_value = zeros(num_iterations1/2,1);

            H1 = f1(z1,z2,gamma1,gamma2);
            H2 = f2(z1,z2,gamma1,gamma2);

            z1_k1 = -Delta1 * z1 + (1/2) * (H1 - conj(H1)*z1^2); 
            z2_k1 = -Delta1 * z2 + (1/2) * (H2 - conj(H2)*z2^2); 

            z1_int1 = z1 + (dt)*z1_k1;
            z2_int1 = z2 + (dt)*z2_k1;

            H11 = f1(z1_int1,z2_int1,gamma1,gamma2);
            H21 = f2(z1_int1,z2_int1,gamma1,gamma2);

            z1_k2 = -Delta1 * z1_int1 + (1/2) * (H11 - conj(H11)*z1_int1^2); 
            z2_k2 = -Delta1 * z2_int1 + (1/2) * (H21 - conj(H21)*z2_int1^2);  

            z1 = z1 + (dt/2)*(z1_k1 + z1_k2);
            z2 = z2 + (dt/2)*(z2_k1 + z2_k2);

            r1_value(t1) = norm(z1);
            r2_value(t1) = norm(z2);
        end

        if ((max(r1_value)-min(r1_value)<tol_variance) && (max(r2_value)-min(r2_value)<tol_variance))
            basin = 0; % store fixed point as 0
        else
            basin  = 1; % store periodic/chaotic solution as 1
        end
    end

    %matrix(i,5) = basin;
    basin_ans(i) = basin;
end
matrix(:,4) = basin_ans;
%toc
%save('/rc_scratch/saad1282/basin_using_maxminfor1.406_Heuns_samephi1.mat', 'gamma', 'tol_variance', 'matrix', 'K1', 'K2')