% ekf3
% X 2*1
% Z 2*1
function [Xf, P0]=tryekf3(X, Z, P, Q, R)

% clc; clear;
% T = 100;
D = 2;  % 状态维度

% q = 0.05;
% r = 0.05;
% Q = [q,0; 0,q];    % 状态噪声方差
% R = [r,0; 0,r];      % 观测噪声方差
F = [1,0; 0,1]; % 状态转移矩阵，对于不同模型需要修改
H = [1,0; 0,1]; % 观测矩阵，对于不同模型需要修改

%% 状态与观测值生成
% X = zeros(D,T); % 真实状态初始化
% X(:,1) = [50; 30]; % (D*1)
% Z = zeros(D,T); % 观测初始化

% for t=2:T
%     X(:,t) = F*X(:,t-1); % + sqrt(Q)*randn(2,1);
% end
% for t=1:T
%     Z(:,t) = H*X(:,t) + sqrt(R)*randn(2,1);
% end

%% kalman filter
% Xf = zeros(2,T);
% Xf(:,1) = X(:,1);   % 初始化kalman估计值
Xf = X;
P0 = P;    % 初始化协方差矩阵
% for t = 2:T
%     Xn = F*Xf(:,t-1);   % 状态预测 X(k+1|k) = F * X(k|k)
    Xn = F*Xf;
    Zn = H*Xn;  % 观测预测
    
    % 如果是非线性状态方程或非线性观测方程，需要在这里计算雅克比矩阵
    % 更新F与H的雅克比矩阵形式
    % 更新结束
    
    P = F*P0*F' + Q;   % 协方差预测 P(k+1|k) = F * P(k|k) * F' + G * Q * G'
    
    Kg = P*H'/(H*P*H' + R);  % kalman gain计算
    
%     Xf(:,t) = Xn + Kg*(Z(:,t) - Zn);    % kalman状态更新
    Xf = Xn + Kg*(Z - Zn);
    P0 = (eye(D) - Kg*H) * P;   % 协方差更新
% end

%% 误差分析
% error = sqrt((X - Xf).^2);
% rms = sqrt( 1/T * sum((X - Xf).^2, 2) );

%% 画图
%{
figure
t = 1:T;
hold on; box on;
plot(t, X(1,:), '-k.', t, Xf(1,:), '-r+', t, Z(1,:), 'g*');
plot(t, X(2,:), '-k.', t, Xf(2,:), '-r+', t, Z(2,:), 'g*');
legend('真实状态', 'kalman 状态');

figure
hold on; box on;
plot(err_kalmanfilter, '-ks', 'MarkerFace', 'r');
xlabel('时间');
ylabel('误差');
%}
