%% Bifurcation for entire range. Fig. 3a

G = load('bifurcation_data_1_15_increasing_june5.mat');

figure(1);
plot(G.gamma1_vec, G.r1_max(1:2,:), 'r.', MarkerSize = 8);
hold on;
plot(G.gamma1_vec, G.r1_min(1:2,:), 'b.', MarkerSize = 8);
plot(G.gamma1_vec, G.r2_max(1:2,:), 'r.', MarkerSize = 8);
plot(G.gamma1_vec, G.r2_min(1:2,:), 'b.', MarkerSize = 8);
plot(G.gamma1_vec(810:838), G.r1_min(1:10,810:838), 'b.', MarkerSize = 8);
plot(G.gamma1_vec(810:838), G.r1_min(1:10,810:838), 'b.', MarkerSize = 8);
plot(G.gamma1_vec(810:838), G.r2_max(1:10,810:838), 'r.', MarkerSize = 8);
plot(G.gamma1_vec(810:838), G.r2_min(1:10,810:838), 'b.', MarkerSize = 8);
grid on
grid minor
xlabel('\gamma')
ylabel('Max/Min r_1, r_2')

%% Bifurcation diagram chaotic region. Fig 3b
a = -0.5;
K1 = 10;
K2 = 10;

% Anti-phase synchronized state function
fun1 = @(gamma)((1/2 - 1/(2*(1 + a)) * (K1/K2) + ((1/2 + 1/(2*(1 + a)) * K1/K2)^2 - 2/(1 - a^2) * sec(gamma)/K2).^(1/2)).^(1/2));
gamma1 = 1.4 : 0.0003: 1.44;
anti_phase = fun1(gamma1);

figure(4);
E = load('bifurcation_data_14_144.mat');
plot(gamma1, anti_phase, 'r-', 'LineWidth',2);
hold on;
plot(E.gamma1_vec(1:2:240), E.r1_max(1,1:2:240), 'k.', 'MarkerSize',5);
plot(E.gamma1_vec(1:2:240), E.r1_min(1,1:2:240), 'k.', 'MarkerSize',5);
plot(E.gamma1_vec(1:2:240), E.r2_max(1,1:2:240), 'k.', 'MarkerSize', 5);
plot(E.gamma1_vec(1:2:240), E.r2_min(1,1:2:240), 'k.', 'MarkerSize', 5);
plot(E.gamma1_vec(240:2:1250), E.r1_max(1:10,240:2:1250), 'k.', 'MarkerSize',5);
plot(E.gamma1_vec(240:2:1250), E.r1_min(1:10,240:2:1250), 'k.', 'MarkerSize',5);
plot(E.gamma1_vec(240:2:1250), E.r2_max(1:10,240:2:1250), 'k.', 'MarkerSize', 5);
plot(E.gamma1_vec(240:2:1250), E.r2_min(1:10,240:2:1250), 'k.', 'MarkerSize', 5);
plot(E.gamma1_vec(1251:2:end), E.r1_max(1,1251:2:end), 'k.', 'MarkerSize',5);
plot(E.gamma1_vec(1251:2:end), E.r1_min(1,1251:2:end), 'k.', 'MarkerSize',5);
plot(E.gamma1_vec(1251:2:end), E.r2_max(1,1251:2:end), 'k.', 'MarkerSize', 5);
plot(E.gamma1_vec(1251:2:end), E.r2_min(1,1251:2:end), 'k.', 'MarkerSize', 5);
grid on
grid minor
xlabel('\gamma');
ylabel('Max/Min r_1, r_2');
xlim([min(E.gamma1_vec), max(E.gamma1_vec)])

%% Bistability figure in the periodic regime Fig. 7b

figure(2);
C = load('bifurcation_data_1_12_decreasing.mat');
E = load('bifurcation_data_07_120_increasing.mat');
plot(C.gamma1_vec(1:3:end), C.r1_max(1:4,1:3:end), 'k^', 'MarkerFaceColor', [0.1,1,0.75]);
hold on;
plot(C.gamma1_vec(1:3:end), C.r1_min(1:4,1:3:end), 'k^', 'MarkerFaceColor', [0.1,1,0.75]);
plot(C.gamma1_vec(1:3:end), C.r2_max(1:4,1:3:end), 'k^', 'MarkerFaceColor', [0.1,1,0.75]);
plot(C.gamma1_vec(1:3:end), C.r2_min(1:4,1:3:end), 'k^', 'MarkerFaceColor', [0.1,1,0.75]);
plot(E.gamma1_vec(1:1:end), E.r1_max(1:4,1:1:end), 'k.', 'MarkerSize', 8);
plot(E.gamma1_vec(1:1:end), E.r1_min(1:4,1:1:end), 'k.', 'MarkerSize', 8);
plot(E.gamma1_vec(1:1:end), E.r2_max(1:4,1:1:end), 'k.', 'MarkerSize', 8);
plot(E.gamma1_vec(1:1:end), E.r2_min(1:4,1:1:end), 'k.', 'MarkerSize', 8);
xlim([1 1.2])
grid on
grid minor
xlabel('\gamma');
ylabel('Max/Min r_1, r_2');

%% Bistability figure in the periodic regime Fig. 7a
figure(3);
F = load('bifurcation_between_11_145_complex_equations_decreasing_perturbed.mat');
D = load('bifurcation_between11and143complex_equations_increasing_perturbed.mat');
plot(F.gamma1_vec(1:end), F.r1_max(1:10,1:end), 'k^', 'MarkerFaceColor', [0.1,1,0.75]);
hold on;
plot(F.gamma1_vec(1:end), F.r1_min(1:10,1:end), 'k^', 'MarkerFaceColor', [0.1,1,0.75]);
plot(F.gamma1_vec(1:end), F.r2_max(1:10,1:end), 'k^', 'MarkerFaceColor', [0.1,1,0.75]);
plot(F.gamma1_vec(1:end), F.r2_min(1:10,1:end), 'k^', 'MarkerFaceColor', [0.1,1,0.75]);
plot(D.gamma1_vec(1:end), D.r1_max(1:8,1:end), 'k.', 'MarkerSize',8);
plot(D.gamma1_vec(1:end), D.r1_min(1:8,1:end), 'k.', 'MarkerSize',8);
plot(D.gamma1_vec(1:end), D.r2_max(1:8,1:end), 'k.', 'MarkerSize',8);
plot(D.gamma1_vec(1:end), D.r2_min(1:8,1:end), 'k.', 'MarkerSize',8);
grid on
grid minor
xlabel('\gamma');
ylabel('Max/Min r_1, r_2');
xlim([1.412, 1.418])