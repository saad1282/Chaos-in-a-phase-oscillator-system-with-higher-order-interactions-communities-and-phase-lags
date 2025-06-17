%% Simulation of the oscillator system before system reduction. This file uses the Kuramoto model 
% that describe the dynamics of phase oscillators

% $ \frac{d \theta_i^{\sigma}}{d t} = \omega_i^{\sigma} + \sum_{\sigma' = 1}^C \frac{K_{1}^{\sigma \sigma'}}{N_{\sigma'}} \sum_{j=1}^{N_{\sigma'}} \sin{(\theta_j^{\sigma'} - \theta_i^{\sigma} - \gamma_{\sigma \sigma'})} + \sum_{\sigma' = 1}^C \sum_{\sigma'' = 1}^C \frac{K_{2}^{\sigma \sigma' \sigma''}}{N_{\sigma'} N_{\sigma''}} \sum_{j=1}^{N_{\sigma'}} \sum_{l=1}^{N_{\sigma''}} \sin{(2\theta_j^{\sigma'} -\theta_l^{\sigma''} - \theta_i^{\sigma} - \gamma_{\sigma \sigma' \sigma''})}$ 

clear
clc

K1 = 10; % dyadic coupling strength
K2 = 10; % triadic coupling strength
alpha = -0.5; %alpha value that determines the coupling strengths
gamma_pair = 1.41; % phase lags that arises in pairwise interaction
gamma_tri = 1.41; % phase lags that arises in three-way interaction

time = 700; %transient time
dt = 0.001; % stepsize
num_iterations = time/dt; % number of time steps for transient state
time1 = 700; %steady state time
num_iterations1 = time1/dt; % number of time steps for steady state

N1 = 5000; % number of oscillators in community 1
N2 = N1; % number of oscillators in community 2

% Intrinsic frequency of oscillators
w1 = zeros(N1,1);

num_ini_conditions = 1;

for i = 1 : N1
    w1(i) = tan((i * pi)/N1 - (N1 + 1) * pi/(2 * N1));
end

w1 = w1(randperm(N1)); % natural frequency of oscillators in community 1
w2 = w1(randperm(N1)); % natural frequency of oscillators in community 2

r1_avg = zeros(num_ini_conditions,1);
r2_avg = zeros(num_ini_conditions,1);
angle_avg = zeros(num_ini_conditions,1);

for i = 1 : num_ini_conditions

    theta1_start = pi*ones(N1,1);
    theta2_start = pi*ones(N2,1);

    z1 = sum(exp(1i*theta1_start))/N1;
    z1_second = sum(exp(2i*theta1_start))/N1;
    z2 = sum(exp(1i*theta2_start))/N2;
    z2_second = sum(exp(2i*theta2_start))/N2;

    r1_timeseries = zeros(num_iterations1,1);
    r2_timeseries = zeros(num_iterations1,1);
    angle_z1 = zeros(num_iterations1,1);
    angle_z2 = zeros(num_iterations1,1);
    angle_diff = zeros(num_iterations1,1);

    for t = 1 : num_iterations

        theta1_diff1 = w1 + imag(exp(-1i*theta1_start)*(K1 * z1 * exp(-1i*gamma_pair)+...
            alpha * K1 * z2 * exp(-1i*gamma_pair)+...
            K2 * z1_second * conj(z1) * exp(-1i*gamma_tri)+...
            alpha * K2 * z1_second * conj(z2) * exp(-1i*gamma_tri)+...
            alpha * K2 * z2_second * conj(z1) * exp(-1i*gamma_tri)+...
            alpha^2 * K2 * z2_second * conj(z2) * exp(-1i*gamma_tri)));

        theta2_diff1 = w2 + imag(exp(-1i*theta2_start)*(K1 * z2 * exp(-1i*gamma_pair)+...
            alpha * K1 * z1 * exp(-1i*gamma_pair)+...
            K2 * z2_second * conj(z2) * exp(-1i*gamma_tri)+...
            alpha * K2 * z2_second * conj(z1) * exp(-1i*gamma_tri)+...
            alpha * K2 * z1_second * conj(z2) * exp(-1i*gamma_tri)+...
            alpha^2 * K2 * z1_second * conj(z1) * exp(-1i*gamma_tri)));

        theta1_start_midstep = theta1_start + dt*theta1_diff1;
        theta2_start_midstep = theta2_start + dt*theta2_diff1;

        z1_mid = sum(exp(1i*theta1_start_midstep))/N1;
        z2_mid = sum(exp(1i*theta2_start_midstep))/N2;
        z1_mid_second = sum(exp(2i*theta1_start_midstep))/N1;
        z2_mid_second = sum(exp(2i*theta2_start_midstep))/N2;

        theta1_diff2 = w1 + imag(exp(-1i*theta1_start_midstep)*(K1 * z1_mid * exp(-1i*gamma_pair)+...
            alpha * K1 * z2_mid * exp(-1i*gamma_pair)+...
            K2 * z1_mid_second * conj(z1_mid) * exp(-1i*gamma_tri)+...
            alpha * K2 * z1_mid_second * conj(z2_mid) * exp(-1i*gamma_tri)+...
            alpha * K2 * z2_mid_second * conj(z1_mid) * exp(-1i*gamma_tri)+...
            alpha^2 * K2 * z2_mid_second * conj(z2_mid) * exp(-1i*gamma_tri)));

        theta2_diff2 = w2 + imag(exp(-1i*theta2_start_midstep)*(K1 * z2_mid * exp(-1i*gamma_pair)+...
            alpha * K1 * z1_mid * exp(-1i*gamma_pair)+...
            K2 * z2_mid_second * conj(z2_mid) * exp(-1i*gamma_tri)+...
            alpha * K2 * z2_mid_second * conj(z1_mid) * exp(-1i*gamma_tri)+...
            alpha * K2 * z1_mid_second * conj(z2_mid) * exp(-1i*gamma_tri)+...
            alpha^2 * K2 * z1_mid_second * conj(z1_mid) * exp(-1i*gamma_tri)));

        theta1_start = theta1_start + (dt/2)*(theta1_diff1 + theta1_diff2);
        theta2_start = theta2_start + (dt/2)*(theta2_diff1 + theta2_diff2);

        z1 = sum(exp(1i*theta1_start))/N1;
        z1_second = sum(exp(2i*theta1_start))/N1;
        z2 = sum(exp(1i*theta2_start))/N2;
        z2_second = sum(exp(2i*theta2_start))/N2;
    end


    for t1 = 1 : num_iterations1

        theta1_diff1 = w1 + imag(exp(-1i*theta1_start)*(K1 * z1 * exp(-1i*gamma_pair)+...
            alpha * K1 * z2 * exp(-1i*gamma_pair)+...
            K2 * z1_second * conj(z1) * exp(-1i*gamma_tri)+...
            alpha * K2 * z1_second * conj(z2) * exp(-1i*gamma_tri)+...
            alpha * K2 * z2_second * conj(z1) * exp(-1i*gamma_tri)+...
            alpha^2 * K2 * z2_second * conj(z2) * exp(-1i*gamma_tri)));

        theta2_diff1 = w2 + imag(exp(-1i*theta2_start)*(K1 * z2 * exp(-1i*gamma_pair)+...
            alpha * K1 * z1 * exp(-1i*gamma_pair)+...
            K2 * z2_second * conj(z2) * exp(-1i*gamma_tri)+...
            alpha * K2 * z2_second * conj(z1) * exp(-1i*gamma_tri)+...
            alpha * K2 * z1_second * conj(z2) * exp(-1i*gamma_tri)+...
            alpha^2 * K2 * z1_second * conj(z1) * exp(-1i*gamma_tri)));

        theta1_start_midstep = theta1_start + dt*theta1_diff1;
        theta2_start_midstep = theta2_start + dt*theta2_diff1;

        z1_mid = sum(exp(1i*theta1_start_midstep))/N1;
        z2_mid = sum(exp(1i*theta2_start_midstep))/N2;
        z1_mid_second = sum(exp(2i*theta1_start_midstep))/N1;
        z2_mid_second = sum(exp(2i*theta2_start_midstep))/N2;

        theta1_diff2 = w1 + imag(exp(-1i*theta1_start_midstep)*(K1 * z1_mid * exp(-1i*gamma_pair)+...
            alpha * K1 * z2_mid * exp(-1i*gamma_pair)+...
            K2 * z1_mid_second * conj(z1_mid) * exp(-1i*gamma_tri)+...
            alpha * K2 * z1_mid_second * conj(z2_mid) * exp(-1i*gamma_tri)+...
            alpha * K2 * z2_mid_second * conj(z1_mid) * exp(-1i*gamma_tri)+...
            alpha^2 * K2 * z2_mid_second * conj(z2_mid) * exp(-1i*gamma_tri)));

        theta2_diff2 = w2 + imag(exp(-1i*theta2_start_midstep)*(K1 * z2_mid * exp(-1i*gamma_pair)+...
            alpha * K1 * z1_mid * exp(-1i*gamma_pair)+...
            K2 * z2_mid_second * conj(z2_mid) * exp(-1i*gamma_tri)+...
            alpha * K2 * z2_mid_second * conj(z1_mid) * exp(-1i*gamma_tri)+...
            alpha * K2 * z1_mid_second * conj(z2_mid) * exp(-1i*gamma_tri)+...
            alpha^2 * K2 * z1_mid_second * conj(z1_mid) * exp(-1i*gamma_tri)));

        theta1_start = theta1_start + (dt/2)*(theta1_diff1 + theta1_diff2);
        theta2_start = theta2_start + (dt/2)*(theta2_diff1 + theta2_diff2);

        z1 = sum(exp(1i*theta1_start))/N1;
        z1_second = sum(exp(2i*theta1_start))/N1;
        z2 = sum(exp(1i*theta2_start))/N2;
        z2_second = sum(exp(2i*theta2_start))/N2;

        r1_timeseries(t1) = norm(z1);
        r2_timeseries(t1) = norm(z2);

        angle_z1(t1) = angle(z1);
        angle_z2(t1) = angle(z2);

        angle_diff(t1) = angle(z1)-angle(z2);
    end

    r1_avg(i) = mean(r1_timeseries);
    r2_avg(i) = mean(r2_timeseries);
    angle_avg(i) = mean(angle_diff);
    angle_diff = angle(z1) - angle(z2);
end

%% plot r1 timeseries and r2 timeseries.
figure(1);
plot(dt*[num_iterations1-45000:10:num_iterations1], r1_timeseries(end-45000:10: end), 'b-.', LineWidth=5);
hold on;
plot(dt*[num_iterations1-45000:10:num_iterations1], r2_timeseries(end-45000:10:end), 'r.', Markersize=10);
grid on
grid minor
xlabel('Time', 'FontSize', 24);
ylabel('r_1, r_2', 'FontSize', 24);
ylim([0 1])
