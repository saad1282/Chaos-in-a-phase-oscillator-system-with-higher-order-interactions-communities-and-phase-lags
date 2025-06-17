%% Plots for basin of attraction for blanca files

% load the file
% matrix file gives grid points and their classification
% Row 1 is initial r1 value
% Row 2 is initial r2 value
% Row 3 is initial angle difference value
% Row 4 includes classification. 0 = Fixed point, 1 = Periodic/Chaotic
% solution
% For $\gamma = 1.406$, we should either get periodic or fixed point based
% on the timeseries plot
% For $\gamma = 1.416$, we should either get chaotic or fixed point based
% on the timeseries plot

%load('basin_using_maxminfor1406_ini_angle_pi6_150_grid_num_Heuns_samephi1_jun13.mat')
load('basin_using_maxminfor1.416_ini_angle_piby6_150_grid_num_Heuns_samephi1_check_symmetry.mat')
grid_num = 150;

% grid points
[X, Y] = meshgrid(0 :1/grid_num: 1, 0:1/grid_num:1);
Z = zeros(size(X));

x = 0.01 : 1/grid_num:1;
y = 0.01 : 1/grid_num:1;

for i = 1 : length(matrix)

    index = find(x == matrix(i,1));
    index1 = find(y == matrix(i,2));

    Z(index,index1) = matrix(i,4);
end
%% plot the basin
figure(1)
s = surf(X,Y,Z, 'Edgecolor', 'none');
colormap("sky")
xlabel('r_1');
ylabel('r_2');
xlim([0 0.99])
ylim([0 0.99])
%% 2D plot of the basin. We could just rotate figure 1 to get a nicer plot though

matrix_chaos = [];
matrix_fixed = [];
index1 = 1;
index2 = 1;

for k = 1 : length(matrix)
    if (matrix(k,4) == 1)
        matrix_chaos(index1,1) = matrix(k,1);
        matrix_chaos(index1,2) = matrix(k,2);
        index1 = index1 + 1;
    else
        matrix_fixed(index2,1) = matrix(k,1);
        matrix_fixed(index2,2) = matrix(k,2);
        index2 = index2 + 1;
    end
end

figure(3)
plot(matrix_chaos(:,1), matrix_chaos(:,2), 'r.');
hold on;
plot(matrix_fixed(:,1), matrix_fixed(:,2), 'b.');
xlabel('r_1');
ylabel('r_2');
title('Basin of attraction for \gamma=1.4155')