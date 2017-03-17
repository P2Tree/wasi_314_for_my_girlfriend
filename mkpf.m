% 参考自：http://studentdavestutorials.weebly.com/particle-filter-with-matlab-code.html
% 粒子滤波器，指定状态矩阵与观测矩阵，如代码中描述
% 现在修改状态方程和观测方程
% 现在我要把它改成二维的，修改成功
% 现在将它变成函数
% r q 均为方差
function [x_est_out error rms]=pf3(Rs, rs, T, N, q, r)

%% 初始化变量

D = 2; % 状态维度
% q = 0.002; % 状态噪声（系统噪声），粒子更新时也需要噪声
% r = 0.003; % 观测噪声
% T = 100; % 采样周期
% N = 100; % 粒子数，粒子数增多会导致计算量增大

initx = Rs(:,1);
pukf = eye(D);     % 初始化状态协方差
pekf = eye(D);
f = @(x) [x(1); x(2)];     % 状态方程
h = @(x) [x(1); x(2)];      % 观测方程
Q = q*eye(D);
R = r*eye(D);

%% 初始化粒子， 先验粒子分布，高斯噪声下的状态初值
% 粒子状态初始估计协方差
P_w = zeros(2, N);

% 依赖于初始粒子的先验分布来构建第一组随机粒子
x_P = Rs;

% 初始化一些和输出有关的变量，添加第一项
z_out = zeros(2, T); % 真实的观测输出矩阵
x_out = zeros(2, T); % 真实的状态输出矩阵
x_est = zeros(2, 1); % 每一步的所有粒子状态均值矩阵
x_est_out = zeros(2, T); % 整个过程的粒子状态估计输出矩阵

%% 粒子滤波器过程
for t = 1:T
    % 更新真实状态，观测位置，这里可以使用非线性状态方程和观测方程
    x = Rs(:,t);% + q * randn(D,1);
    z = rs(:,t);% + r * randn(D,1);
    % 这里开始滤波器过程
    for i = 1:N
        %给定粒子先验，粒子更新过程会加入噪声，增加粒子获取的随机性
        x_P_update(:,i) = x_P(:,i) + q * randn(D,1);
%         x_P_update(:,i) = 0.5 * x_P(:,i) + 25 * x_P(:,i) ./ (1 + x_P(:,i).^2) + 8 * cos(1.2 * (t-1)) + q * randn(D, 1);
        % 用新更新的粒子更新观测量，由于我们已知观测方程，所以更新粒子的观测值不需要增加噪声
        z_update(:,i) = x_P_update(:,i);
%         z_update(:,i) = x_P_update(:,i);
        % 下边采用ukf -> ekf来更新
        [zukf, pukf] = ukf3(f, initx, pukf, h, z_update(:,i), Q, R);
        [zekf, pekf] = tryekf3(z, zukf, pekf, pukf, R);
        % 生成每个粒子的权重，权重更新时一种概率模型，但这里我们知道我们的观测误差是高斯分布，协方差为x_R
        P_w(:,i) = (1 / sqrt(2 * pi * r)) * exp(- (z - zekf).^2 / (2 * r));
%         P_w(:,i) = P_w(:,i) .* (mean(zekf,2)) + sqrt(pekf) * randn(2, 1);
    end
    
    % 判断粒子退化
    Psum = sum(P_w,2);
    if Psum(1) == 0 || Psum(2) == 0
        disp('mkpf算法，粒子退化，无法运算');
        break;
    end
    % 规范化粒子权重，使和为1
    for a = 1:D
        P_w(a,:) = P_w(a,:) ./ sum(P_w(a,:));
    end
    
    % 重采样过程：针对当前的新分布，随机采样来生成新的估计粒子
    % 根据粒子权重重采样，通过粒子权重生成概率分布的累积分布，概率性选择累计值并保留权重大的粒子估计
    for a = 1:D
        cum = cumsum(P_w(a,:));
        for i = 1:N
            x_P(a,i) = x_P_update(a, find(rand <= cum, 1));
        end
    end
    
    %最后的工作，比如计算期望和方差
    x_est = mean(x_P,2);
    
    % 保存状态和估计数据，用于输出
    x_out(:,t) = x;
    z_out(:,t) = z;
    x_est_out(:,t) = x_est ;
    
%     x_est_out(:,t) =x_est_out(:,t) + 0.45*randn(2,1);
end

% L = [0:1/100:0.99;0:1/100:0.99];
% x_est_out(:,t) = rand*(x_est_out(:,t))+(-0.5*randn(1,1));

error = sqrt((x_out - x_est_out).^2);
rms = sqrt( 1/T * sum((x_out - x_est_out).^2, 2) );

%% 绘图：最终结果
% figure(3);
% t = 1:T;
% for d = 1:D
%     subplot(D, 1, d);
%     plot(t, x_out(d,:), '.-b', t, x_est_out(d,:), '-.r', t, z_out(d,:), '*g', 'linewidth', 1);
% end
% set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
% xlabel('time stem'); ylabel('position');
% legend('true position', 'pf estimate', 'measurement position');
% 
% figure(4);
% hold on;
% plot(t, error(1,:), '.-b', t, error(2,:), '-.r', 'linewidth', 1);
% set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
% xlabel('time stem'); ylabel('error value');
% legend('rs_x', 'rs_y');
