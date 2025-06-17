% Lyapunov Dimension for our system. Plots Lyapunov Exponents and Lyapuov
% dimension if we had lyapunov exponents file

load('lyapunov_full_range_reduced_data.mat')
%load('lyapunov_full_range_newset.mat')
lyap_dim = zeros(length(gamma_vec),1);

for i = 1 : length(gamma_vec)

        if abs(lyapunov_exponent(i,1)) < 0.01
            lyapunov_exponent(i,1) = 0;
        end

    index = 1;
    lyap_sum = 0;

    while lyap_sum >= 0
        lyap_sum = lyap_sum + lyapunov_exponent(i,index);
        index = index + 1;
    end
    lyap_dim(i) = index - 1 + lyap_sum/(abs(lyapunov_exponent(i,index-1)));
end

%%
figure(1);
plot(gamma_vec(1:2:end), lyap_dim(1:2:end), 'b-', LineWidth=2);
ylabel('Lyapunov Dimension');
xlabel('\gamma');
xline(1.4078, 'g-', LineWidth=2)
xline(1.41084, 'g-', LineWidth=2)
xline(1.41128, 'g-', LineWidth=2)
xline(1.4164, 'g-', LineWidth=2)
xlim([1.4055 1.418])
ylim([0, 2.2]);
grid on;
grid minor;

%%
figure(2);
plot(gamma_vec(1,1:2:end), lyapunov_exponent(1:2:end,1), 'b-', LineWidth=2);
ylabel('Largest Lyapunov Exponent');
xlabel('\gamma');
xline(1.4078, 'g-', LineWidth=2)
xline(1.41084, 'g-', LineWidth=2)
xline(1.41128, 'g-', LineWidth=2)
xline(1.4164, 'g-', LineWidth=2)
xlim([1.4055 1.418])
grid on;
grid minor;
%%
figure(3);
plot(gamma_vec(1,1:2:end), lyapunov_exponent(1:2:end,1), 'b-', LineWidth=2);
hold on;
plot(gamma_vec(1,1:2:end), lyapunov_exponent(1:2:end,2), 'r-', LineWidth=2);
plot(gamma_vec(1,1:2:end), lyapunov_exponent(1:2:end,3), 'k-', LineWidth=2);
ylabel('Lyapunov Exponents');
xlabel('\gamma');
xline(1.4078, 'g-', LineWidth=2)
xline(1.41084, 'g-', LineWidth=2)
xline(1.41128, 'g-', LineWidth=2)
xline(1.4164, 'g-', LineWidth=2)
xlim([1.4055 1.418])
grid on;
grid minor;