% 状态方程是：Ci=q/(2*pi*k)*1/(sqrt((xs-xi)^2+(ys-yi)^2)*exp(-V/2*k*(sqrt((xs-xi)^2+(ys-yi)^2))-deltax)
% 其中,q k 为常数，Ci xs ys xi yi V为变量
% deltax = (xs-xi)*cos(theta)+(ys-yi)*sin(theta)
% 其中,theta为变量，将二式带入一式，可得
% Ci = func(xs, ys, xi, yi, V, theta)
% 假设, (xs, ys) = (50, 30)
% 第一步，通过拟定一定量的(xi, yi, V, theta)来计算对应的Ci
% 第二步，认为(xs, ys)未知，通过已知的(xi, yi, V, theta, Ci)计算(xs, ys)
% 第三步，将得出的(xs, ys)带入滤波器，比较与真实值的均方误差

%{ 
程序中会调用一些文件，务必保证所有文件在一个文件夹，调用关系：
    wasi.m
        |--- pf3.m --> 内含绘图
        |--- upf.m --> 内含绘图
        |       |--- ukf3.m
        |               |--- ut.m
        |               |--- sigmas.m
        |--- ekpf.m --> 内含绘图
        |       |--- tryekf3.m : 重写的
        |--- mkpf.m --> 内含绘图
                |--- ukf3.m
                |       |--- 同上
                |--- tryekf3.m
        
%}
clear;clc;
% for cycle = 1:10
% clearvars -except cycle pf_rms_cyc ekpf_rms_cyc upf_rms_cyc mkpf_rms_cyc;
M = 100;    % 试验次数,M为偶数
m = 1;
N = 50;    % 采样次数
n = 1;
Xinoise = 0.01 * randn(1, N); % 位置观测噪声
Yinoise = 0.01 * randn(1, N);
Vnoise = 0.001 * randn(M, N);   % 速度观测噪声
Tnoise = 0.001 * randn(M, N);  % 角度观测噪声
Xs = 50 * ones(1,M); % 这两个参数的值是释放源的坐标
Ys = 30 * ones(1,M);
Rs = [Xs; Ys];
% 
wallx = 100;
wally = 100; % 假设墙壁尺寸为100*100
maxV = 100; % 假设最大风速为100，这是可以容忍的最小风速，如果小于100，pf粒子可能退化
maxtheta = 2*pi; % 最大风向角2pi
q = 500000; % 建模里边的常量，谨慎修改，可能导致第二步估计Ci的时候发散
k = 1000;
% 
xi = wallx * rand(1, N);
yi = wally * rand(1, N); % xi yi是某次N*M次试验过程中确定的
V = zeros(M, N);
theta = zeros(M, N);
Ci = zeros(M, N);
distance = zeros(1,N);
deltax = zeros(1, N);

%% 第一步
for m = 1:M
    % v theta需要M组试验每组重新测算
%     V(m,:) = maxV; % 如果取消屏蔽这条，而屏蔽下一条，意味着风速不变，，，（风速不变）
    V(m,:) = maxV * rand(1, N); % 反之，风速的变化服从高斯随机分布，0到maxV，，（风速变）
    theta(m,:) = maxtheta; % 如果取消屏蔽这条，而屏蔽下一条，意味着风向不变，，（风向不变）
%     theta(m,:) = maxtheta * rand(1, N); % 反之，风向的变化服从高斯随机分布，0到2pi，，（风向变）
    for n = 1:N
        distance(n) = sqrt((Xs(n)-xi(n))^2+(Ys(n)-yi(n))^2);
        deltax(n) = (Xs(n) - xi(n))*cos(theta(m, n)) + (Ys(n) - yi(n))*sin(theta(m, n));
        eee(n) = -V(m, n)/(2*k)*(distance(n))-deltax(n);
        Ci(m, n)=q/(2*pi*k)*1/(distance(n))*exp(eee(n));
    end
end

%% 第二步
Thold = 1000;
c = [Xs(1);Ys(1)];
options=optimset('fminsearch');
options.To1X=0.1;
options.Display='off';
    %rs = [Xs; Ys]; 为未知量，在该步中可以通过单纯形算法解多元非线性优化做参数估计
for m = 1:M
    ri = [xi + Xinoise; yi + Yinoise; V(m,:) + Vnoise(m); theta(m,:) + Tnoise(m)];
    [rs(m,:), ~, ~, ~] = fminsearch(@fun, c, options, ri, Ci(m,:)); %递归N次，估计一个最优参数rs
    % 下边是在做一个阈值滤波，限制可能导致的估计局部最小化问题
%     if rs(m,1) > (Xs(1,1)+Thold)
%         rs(m,1) = Xs(1,1)+Thold;
%     elseif rs(m,1) < (Xs(1,1)-Thold)
%         rs(m,1) = Xs(1,1)-Thold;
%     end
%     if rs(m,2) > (Ys(1,1)+Thold)
%         rs(m,2) = Ys(1,1) + Thold;
%     elseif rs(m,2) < (Ys(1,1) - Thold)
%         rs(m,2) = Ys(1,1)-Thold;
%     end
end
rs = rs';

%% 第三步
[pf_rs_est_out pf_error pf_rms] = pf3(Rs, rs, M, 100, 0.1, 1);
[ekpf_rs_est_out ekpf_error ekpf_rms] = ekpf(Rs, rs, M, 100, 0.1, 1);
[upf_rs_est_out upf_error upf_rms] = upf(Rs, rs, M, 100, 0.1, 1);
[mkpf_rs_est_out mkpf_error mkpf_rms] = mkpf(Rs, rs, M, 100, 0.1, 1);
% [ekf_rs_est_out, ekf_P0]=tryekf3(Rs, rs, M, 0.2, 0.3);

pf_rms
ekpf_rms
upf_rms
mkpf_rms

MEANxpf = mean(pf_rms)
VARxpf = var(pf_rms)
MEANxekpf = mean(ekpf_rms)
VARxekpf = var(ekpf_rms)
MEANxupf = mean(upf_rms)
VARxupf = var(upf_rms)
MEANxmkpf = mean(mkpf_rms)
VARxmkpf = var(mkpf_rms)

% pf_rms_cyc(:,cycle) = pf_rms
% ekpf_rms_cyc(:,cycle) = ekpf_rms
% upf_rms_cyc(:,cycle) = upf_rms
% mkpf_rms_cyc(:,cycle) = mkpf_rms
% end
% pf_rms_cyc_mean = mean(pf_rms_cyc, 2)
% ekpf_rms_cyc_mean = mean(ekpf_rms_cyc, 2)
% upf_rms_cyc_mean = mean(upf_rms_cyc, 2)
% mkpf_rms_cyc_mean = mean(mkpf_rms_cyc, 2)

%% 画图
%
m = 1:M;
figure(1);
subplot(211);
% plot(m, Rs(1,:), '.-k', m, rs(1,:), '*k', m, pf_rs_est_out(1,:) , '-mo', m, ukf_rs_est_out(1,:), '-r*',  m, ekpf_rs_est_out(1,:), '-go', m, upf_rs_est_out(1,:), '-b*', 'linewidth', 1);
plot(m, Rs(1,:), '.-k', m, rs(1,:), '*k', m, pf_rs_est_out(1,:) , '-mo', m, ekpf_rs_est_out(1,:), '-go', ...
    m, upf_rs_est_out(1,:), '-bo', m, mkpf_rs_est_out(1,:), '-ro', 'linewidth', 1);
title('x轴');
set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
xlabel('实验次数'); ylabel('状态估计值');
legend('true position',  'measurement position', 'PF estimate','EKPF estimate',  'UPF estimate', 'MKPF estimate');
% axis([0 M 40 60]);

subplot(212);
plot(m, Rs(2,:), '.-k', m, rs(2,:), '*k', m, pf_rs_est_out(2,:) , '-mo', m, ekpf_rs_est_out(2,:), '-go', ...
    m, upf_rs_est_out(2,:), '-bo', m, mkpf_rs_est_out(2,:), '-ro', 'linewidth', 1);
title('y轴');
set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
xlabel('实验次数'); ylabel('状态估计值');
% axis([0 M 20 40]);

figure(2);
hold on;
subplot(211);
% plot(m, pf_error(1,:), '-mo', m, ukf_error(1,:), '-r*',  m, ekpf_error(1,:), '-go', m, upf_error(1,:), '-b*', 'linewidth', 1);
plot(m, pf_error(1,:), '-mo', m, ekpf_error(1,:), '-go', m, upf_error(1,:), '-bo', m, mkpf_error(1,:), '-ro', 'linewidth', 1);
title('x轴估计误差');
set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
xlabel('实验次数'); ylabel('均方误差MSE');
legend('PF estimate', 'EKPF estimate',  'UPF estimate', 'MKPF estimate');
% axis([0 M 0 1]);

subplot(212);
% plot(m, pf_error(2,:), '-mo', m, ukf_error(2,:), '-r*',m, ekpf_error(2,:), '-go', m, upf_error(2,:), '-b*', 'linewidth', 1);
plot(m, pf_error(2,:), '-mo', m, ekpf_error(2,:), '-go', m, upf_error(2,:), '-bo', m, mkpf_error(2,:), '-ro', 'linewidth', 1);
title('y轴估计误差');
set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
xlabel('实验次数'); ylabel('均方误差MSE');
% axis([0 M 0 1]);
