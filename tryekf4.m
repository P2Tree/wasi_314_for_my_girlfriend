% tryekf4
clear; clc;
% function [xf, p] = tryekf4(x, z, p0, q, r)
p0 = [1, 0; 0, 1];

T = 100;
Q = 1;
R = 1;
xf = [50; 30];
x = [50; 30];
p = p0;

x_out = zeros(2,T);
z_out = zeros(2,T);
x_est_out = zeros(2,T);

for k = 1:T
    x = 0.5 * x + (25 * x ./ (1 + x.^2)) + 8 * cos(1.2 * (k-1)) + sqrt(Q) * randn(2, 1);
    z = (x.^2 / 20) + sqrt(R) * randn(2, 1);
    
    xn = 0.5 * xf + 25 * xf ./ (1 + xf.^2) + 8 * cos(1.2 * (k-1));
    zn = xn.^2 / 20;
    
    % 相关矩阵
    F = diag(0.5 + 25 * (1 - xf.^2) ./ ((1 + xf.^2).^2));
    H = diag(xn / 10);    % 这里是xn的原因，与上边状态方程有关
    
    p = F * p * F' + Q;
    m = (H * p * H' + R);
    n = p * H';
    kg = n / m;
    
    xf = xn + kg * (z - zn);
    p = p - kg * H * p
    
    x_out(:, k) = x;
    z_out(:, k) = z;
    x_est_out(:, k) = xf;
end

figure
hold on;
plot(1:T, x_out, '-.k', 1:T, z_out, '*k', 1:T, x_est_out, '-or');