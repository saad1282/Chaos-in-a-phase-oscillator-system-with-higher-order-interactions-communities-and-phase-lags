% %% This code generates bifurcation plot of the phase lagged oscilators

clear all
clc

%% Parameters

gamma1_vec = 1.43:0.0005:1.44;  % Phase lag values
gamma2_vec = gamma1_vec;    % Phase lag values
alpha = -0.5;   % alpha which determines coupling strengths
K1 = 10;    % pairwise coupling strengths
K2 = 10;    % triangular coupling strengths

Delta1 = 1;     % Delta for Lorentzian frequency distribution

transient_time = 3000;  % transient time
dt = 0.001;     % dt step for Heuns Method
transient_timesteps = transient_time/dt;    % Number of transient timesteps

num_max = 10;   % Number of local maxima
num_min = 10;   % Number of local minima

r1_max = transpose(zeros(length(gamma1_vec), num_max));     % Intialize the matrix for local maxima of community 1 order parameters
r1_min = transpose(zeros(length(gamma1_vec), num_min));     % Intialize the matrix for local minima of community 1 order parameters
r2_max = transpose(zeros(length(gamma1_vec), num_max));     % Intialize the matric for local maxima of community 2 order parameters
r2_min = transpose(zeros(length(gamma1_vec), num_min));     % Intialize the matric for local maxima of community 2 order parameters

num_initial = 1;    % Number of initial conditions (Repetition number)

%% complex equations
f1 = @(z1,z2, gamma1, gamma2) (K1 * exp(-1i * gamma1) * (z1 + alpha * z2))...
     + K2 * exp(-1i * gamma2) * (z1^2 * conj(z1) + alpha * (z1)^2 * conj(z2) + alpha * (z2)^2 * conj(z1) + alpha^2 * (z2)^2 * conj(z2));
 
f2 = @(z1,z2, gamma1, gamma2) (K1 * exp(-1i * gamma1) * (z2 + alpha * z1))...
     + K2 * exp(-1i * gamma2) * (z2^2 * conj(z2) + alpha * (z2)^2 * conj(z1) + alpha * (z1)^2 * conj(z2) + alpha^2 * (z1)^2 * conj(z1));
 
%% polar equations
 
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

%%
gamma1 = gamma1_vec(end);
gamma2 = gamma1_vec(end);

r1_diff = 0;
r2_diff = 0;
tol = 1e-2;

r1_ini = rand();
r2_ini = rand();
phi_ini = pi*rand();
 
z1 = r1_ini;
z2 = r2_ini*exp(1i*phi_ini); 
 
for j = 1 : num_initial
    for i = length(gamma1_vec) : -1: 1

        i
        gamma1 = gamma1_vec(i);
        gamma2 = gamma2_vec(i);

        % Run transient simulation
        for t = 1 : transient_timesteps
            H1 = f1(z1,z2,gamma1,gamma2);
            H2 = f2(z1,z2,gamma1,gamma2);

            z1_diff = -Delta1 * z1 + (1/2) * (H1 - conj(H1)*z1^2);
            z2_diff = -Delta1 * z2 + (1/2) * (H2 - conj(H2)*z2^2);

            z1_int = z1 + dt*z1_diff;
            z2_int = z2 + dt*z2_diff;

            H1_new = f1(z1_int,z2_int,gamma1,gamma2);
            H2_new = f2(z1_int,z2_int,gamma1,gamma2);

            z1_diff1 = -Delta1 * z1_int + (1/2) * (H1_new - conj(H1_new)*z1_int^2);
            z2_diff1 = -Delta1 * z2_int + (1/2) * (H2_new - conj(H2_new)*z2_int^2);

            z1 = z1 + (z1_diff + z1_diff1)*(dt/2);
            z2 = z2 + (z2_diff + z2_diff1)*(dt/2);
        end

        % Extra Iteration
        H1 = f1(z1,z2,gamma1,gamma2);
        H2 = f2(z1,z2,gamma1,gamma2);

        z1_diff = -Delta1 * z1 + (1/2) * (H1 - conj(H1)*z1^2);
        z2_diff = -Delta1 * z2 + (1/2) * (H2 - conj(H2)*z2^2);

        z1_int = z1 + dt*z1_diff;
        z2_int = z2 + dt*z2_diff;

        H1_new = f1(z1_int,z2_int,gamma1,gamma2);
        H2_new = f2(z1_int,z2_int,gamma1,gamma2);

        z1_diff1 = -Delta1 * z1_int + (1/2) * (H1_new - conj(H1_new)*z1_int^2);
        z2_diff1 = -Delta1 * z2_int + (1/2) * (H2_new - conj(H2_new)*z2_int^2);

        z1_new = z1 + (z1_diff + z1_diff1)*(dt/2);
        z2_new = z2 + (z2_diff + z2_diff1)*(dt/2);

        cont = 1;
        index1 = 1;
        index2 = 1;
        index3 = 1;
        index4 = 1;
        % One extra iteration and check if there are mins/max
        while (cont)

            H1 = f1(z1_new,z2_new,gamma1,gamma2);
            H2 = f2(z1_new,z2_new,gamma1,gamma2);

            z1_diff = -Delta1 * z1_new + (1/2) * (H1 - conj(H1)*z1_new^2);
            z2_diff = -Delta1 * z2_new + (1/2) * (H2 - conj(H2)*z2_new^2);

            z1_int = z1_new + dt*z1_diff;
            z2_int = z2_new + dt*z2_diff;

            H1_new = f1(z1_int,z2_int,gamma1,gamma2);
            H2_new = f2(z1_int,z2_int,gamma1,gamma2);

            z1_diff1 = -Delta1 * z1_int + (1/2) * (H1_new - conj(H1_new)*z1_int^2);
            z2_diff1 = -Delta1 * z2_int + (1/2) * (H2_new - conj(H2_new)*z2_int^2);

            z1_new1 = z1_new + (z1_diff + z1_diff1)*(dt/2);
            z2_new1 = z2_new + (z2_diff + z2_diff1)*(dt/2);

            if (norm(z1) <= norm(z1_new) && (norm(z1_new) >= norm(z1_new1)))
                if (index1 <= num_max+1)
                    r1_max(index1,i) = norm(z1_new);
                    index1 = index1 + 1;
                end
            end

            if ((norm(z1) >= norm(z1_new)) && (norm(z1_new) <= norm(z1_new1)))
                if (index2 <= num_min+1)
                    r1_min(index2,i) = norm(z1_new);
                    index2 = index2 + 1;
                end
            end

            if ((norm(z2) <= norm(z2_new)) && (norm(z2_new) >= norm(z2_new1)))
                if (index3 <= num_max+1)
                    r2_max(index3,i) = norm(z2_new);
                    index3 = index3 + 1;
                end
            end

            if ((norm(z2) >= norm(z2_new)) && (norm(z2_new) <= norm(z2_new1)))
                if (index4 <= num_min+1)
                    r2_min(index4,i) = norm(z2_new);
                    index4 = index4 + 1;
                end
            end

            z1 = z1_new;
            z1_new = z1_new1;
            z2 = z2_new;
            z2_new = z2_new1;

            if (((index1) >= num_max) && ((index2) >= num_min) && ((index3) >= num_max) && ((index4) >= num_min))
                cont = 0;
            end
        end
        z1 = z1 + rand*1e-02;
        z2 = z2 + rand*1e-02;
    end
end

%% If you want to use polar equations, run this section

% % for j = 1 : num_initial
% % 
% %     r1_ini = rand();
% %     r2_ini = rand();
% %     phi_ini = pi*rand();
% % 
% %     for i = 1 : length(gamma1_vec)
% %         gamma1 = gamma1_vec(i);
% %         gamma2 = gamma2_vec(i);
% % 
% %         for t = 1 : transient_timesteps
% % 
% %         r1_k1 = g1(r1_ini,r2_ini, phi_ini, gamma1, gamma2);
% %         r2_k1 = g2(r1_ini,r2_ini, phi_ini, gamma1, gamma2);
% %         phi_k1 = g3(r1_ini,r2_ini, phi_ini, gamma1, gamma2);
% % 
% %         r1_int = r1_ini + (dt)*r1_k1;
% %         r2_int = r2_ini + (dt)*r2_k1;
% %         phi_int = phi_ini + (dt)*phi_k1;
% % 
% %         r1_k2 = g1(r1_int,r2_int, phi_int, gamma1, gamma2);
% %         r2_k2 = g2(r1_int,r2_int, phi_int, gamma1, gamma2);
% %         phi_k2 = g3(r1_int,r2_int, phi_int, gamma1, gamma2);
% % 
% %         r1_ini = r1_ini + (dt/2)*(r1_k1 + r1_k2);
% %         r2_ini = r2_ini + (dt/2)*(r2_k1 + r2_k2);
% %         phi_ini = phi_ini + (dt/2)*(phi_k1 + phi_k2);
% %         end
% % 
% %        %% Extra Iteration
% % 
% %         r1_k1 = g1(r1_ini,r2_ini, phi_ini, gamma1, gamma2);
% %         r2_k1 = g2(r1_ini,r2_ini, phi_ini, gamma1, gamma2);
% %         phi_k1 = g3(r1_ini,r2_ini, phi_ini, gamma1, gamma2);
% % 
% %         r1_int = r1_ini + (dt)*r1_k1;
% %         r2_int = r2_ini + (dt)*r2_k1;
% %         phi_int = phi_ini + (dt)*phi_k1;
% % 
% %         r1_k2 = g1(r1_int,r2_int, phi_int, gamma1, gamma2);
% %         r2_k2 = g2(r1_int,r2_int, phi_int, gamma1, gamma2);
% %         phi_k2 = g3(r1_int,r2_int, phi_int, gamma1, gamma2);
% % 
% %         r1_new = r1_ini + (dt/2)*(r1_k1 + r1_k2);
% %         r2_new = r2_ini + (dt/2)*(r2_k1 + r2_k2);
% %         phi_new = phi_ini + (dt/2)*(phi_k1 + phi_k2);
% % 
% %         cont = 1;
% %         index1 = 1;
% %         index2 = 1;
% %         index3 = 1;
% %         index4 = 1;
% %         %% One extra iteration and check if there are mins/max
% %         while (cont)
% % 
% %         r1_k1 = g1(r1_new, r2_new, phi_new, gamma1, gamma2);
% %         r2_k1 = g2(r1_new, r2_new, phi_new, gamma1, gamma2);
% %         phi_k1 = g3(r1_new, r2_new, phi_new, gamma1, gamma2);
% % 
% %         r1_int = r1_new + (dt)*r1_k1;
% %         r2_int = r2_new + (dt)*r2_k1;
% %         phi_int = phi_new + (dt)*phi_k1;
% % 
% %         r1_k2 = g1(r1_int, r2_int, phi_int, gamma1, gamma2);
% %         r2_k2 = g2(r1_int, r2_int, phi_int, gamma1, gamma2);
% %         phi_k2 = g3(r1_int, r2_int, phi_int, gamma1, gamma2);
% % 
% %         r1_new1 = r1_new + (dt/2)*(r1_k1 + r1_k2);
% %         r2_new1 = r2_new + (dt/2)*(r2_k1 + r2_k2);
% %         phi_new1 = phi_new + (dt/2)*(phi_k1 + phi_k2);
% % 
% %             if (r1_ini <= r1_new && (r1_new >= r1_new1))
% %                 if (index1 <= num_max+1)
% %                 r1_max(index1,i) = r1_new;
% %                 index1 = index1 + 1;
% %                 end
% %             end
% % 
% %             if ((r1_ini >= r1_new) && (r1_new <= r1_new1))
% %                 if (index2 <= num_min+1)
% %                 r1_min(index2,i) = r1_new;
% %                 index2 = index2 + 1;
% %                 end
% %             end
% % 
% %             if ((r2_ini <= r2_new) && (r2_new >= r2_new1))
% %                 if (index3 <= num_max+1)
% %                 r2_max(index3,i) = r2_new;
% %                 index3 = index3 + 1;
% %                 end
% %             end
% % 
% %             if ((r2_ini >= r2_new) && (r2_new <= r2_new1))
% %                 if (index4 <= num_min+1)
% %                 r2_min(index4,i) = r2_new;
% %                 index4 = index4 + 1;
% %                 end
% %             end
% % 
% %             r1_ini = r1_new;
% %             r1_new = r1_new1;
% %             r2_ini = r2_new;
% %             r2_new = r2_new1;
% % 
% %             if (((index1) > num_max) && ((index2) > num_min) && ((index3) > num_max) && ((index4) > num_min))
% %                 cont = 0;
% %             end
% %         end
% %     end
% % end
 
%% Plot local maxima and minima
figure(1);
hold on;
plot(gamma1_vec,r1_max(1:10,:), 'r.', MarkerSize=1);
plot(gamma1_vec,r1_min(1:10,:), 'r.', MarkerSize=1);
plot(gamma1_vec,r2_max(1:10,:), 'b.', MarkerSize=1);
plot(gamma1_vec,r2_min(1:10,:), 'b.', MarkerSize=1);
hold off;

%% If you want to save the workspace

%save('yourdirectory/filename.mat', 'gamma1_vec', 'gamma2_vec', 'r1_max', 'r2_max', 'r1_min', 'r2_min', 'K1', 'K2', 'transient_time')