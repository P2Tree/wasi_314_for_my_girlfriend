% ״̬�����ǣ�Ci=q/(2*pi*k)*1/(sqrt((xs-xi)^2+(ys-yi)^2)*exp(-V/2*k*(sqrt((xs-xi)^2+(ys-yi)^2))-deltax)
% ����,q k Ϊ������Ci xs ys xi yi VΪ����
% deltax = (xs-xi)*cos(theta)+(ys-yi)*sin(theta)
% ����,thetaΪ����������ʽ����һʽ���ɵ�
% Ci = func(xs, ys, xi, yi, V, theta)
% ����, (xs, ys) = (50, 30)
% ��һ����ͨ���ⶨһ������(xi, yi, V, theta)�������Ӧ��Ci
% �ڶ�������Ϊ(xs, ys)δ֪��ͨ����֪��(xi, yi, V, theta, Ci)����(xs, ys)
% �����������ó���(xs, ys)�����˲������Ƚ�����ʵֵ�ľ������

%{ 
�����л����һЩ�ļ�����ر�֤�����ļ���һ���ļ��У����ù�ϵ��
    wasi.m
        |--- pf3.m --> �ں���ͼ
        |--- upf.m --> �ں���ͼ
        |       |--- ukf3.m
        |               |--- ut.m
        |               |--- sigmas.m
        |--- ekpf.m --> �ں���ͼ
        |       |--- tryekf3.m : ��д��
        |--- mkpf.m --> �ں���ͼ
                |--- ukf3.m
                |       |--- ͬ��
                |--- tryekf3.m
        
%}
clear;clc;
% for cycle = 1:10
% clearvars -except cycle pf_rms_cyc ekpf_rms_cyc upf_rms_cyc mkpf_rms_cyc;
M = 100;    % �������,MΪż��
m = 1;
N = 50;    % ��������
n = 1;
Xinoise = 0.01 * randn(1, N); % λ�ù۲�����
Yinoise = 0.01 * randn(1, N);
Vnoise = 0.001 * randn(M, N);   % �ٶȹ۲�����
Tnoise = 0.001 * randn(M, N);  % �Ƕȹ۲�����
Xs = 50 * ones(1,M); % ������������ֵ���ͷ�Դ������
Ys = 30 * ones(1,M);
Rs = [Xs; Ys];
% 
wallx = 100;
wally = 100; % ����ǽ�ڳߴ�Ϊ100*100
maxV = 100; % ����������Ϊ100�����ǿ������̵���С���٣����С��100��pf���ӿ����˻�
maxtheta = 2*pi; % �������2pi
q = 500000; % ��ģ��ߵĳ����������޸ģ����ܵ��µڶ�������Ci��ʱ��ɢ
k = 1000;
% 
xi = wallx * rand(1, N);
yi = wally * rand(1, N); % xi yi��ĳ��N*M�����������ȷ����
V = zeros(M, N);
theta = zeros(M, N);
Ci = zeros(M, N);
distance = zeros(1,N);
deltax = zeros(1, N);

%% ��һ��
for m = 1:M
    % v theta��ҪM������ÿ�����²���
%     V(m,:) = maxV; % ���ȡ��������������������һ������ζ�ŷ��ٲ��䣬���������ٲ��䣩
    V(m,:) = maxV * rand(1, N); % ��֮�����ٵı仯���Ӹ�˹����ֲ���0��maxV���������ٱ䣩
    theta(m,:) = maxtheta; % ���ȡ��������������������һ������ζ�ŷ��򲻱䣬�������򲻱䣩
%     theta(m,:) = maxtheta * rand(1, N); % ��֮������ı仯���Ӹ�˹����ֲ���0��2pi����������䣩
    for n = 1:N
        distance(n) = sqrt((Xs(n)-xi(n))^2+(Ys(n)-yi(n))^2);
        deltax(n) = (Xs(n) - xi(n))*cos(theta(m, n)) + (Ys(n) - yi(n))*sin(theta(m, n));
        eee(n) = -V(m, n)/(2*k)*(distance(n))-deltax(n);
        Ci(m, n)=q/(2*pi*k)*1/(distance(n))*exp(eee(n));
    end
end

%% �ڶ���
Thold = 1000;
c = [Xs(1);Ys(1)];
options=optimset('fminsearch');
options.To1X=0.1;
options.Display='off';
    %rs = [Xs; Ys]; Ϊδ֪�����ڸò��п���ͨ���������㷨���Ԫ�������Ż�����������
for m = 1:M
    ri = [xi + Xinoise; yi + Yinoise; V(m,:) + Vnoise(m); theta(m,:) + Tnoise(m)];
    [rs(m,:), ~, ~, ~] = fminsearch(@fun, c, options, ri, Ci(m,:)); %�ݹ�N�Σ�����һ�����Ų���rs
    % �±�������һ����ֵ�˲������ƿ��ܵ��µĹ��ƾֲ���С������
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

%% ������
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

%% ��ͼ
%
m = 1:M;
figure(1);
subplot(211);
% plot(m, Rs(1,:), '.-k', m, rs(1,:), '*k', m, pf_rs_est_out(1,:) , '-mo', m, ukf_rs_est_out(1,:), '-r*',  m, ekpf_rs_est_out(1,:), '-go', m, upf_rs_est_out(1,:), '-b*', 'linewidth', 1);
plot(m, Rs(1,:), '.-k', m, rs(1,:), '*k', m, pf_rs_est_out(1,:) , '-mo', m, ekpf_rs_est_out(1,:), '-go', ...
    m, upf_rs_est_out(1,:), '-bo', m, mkpf_rs_est_out(1,:), '-ro', 'linewidth', 1);
title('x��');
set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
xlabel('ʵ�����'); ylabel('״̬����ֵ');
legend('true position',  'measurement position', 'PF estimate','EKPF estimate',  'UPF estimate', 'MKPF estimate');
% axis([0 M 40 60]);

subplot(212);
plot(m, Rs(2,:), '.-k', m, rs(2,:), '*k', m, pf_rs_est_out(2,:) , '-mo', m, ekpf_rs_est_out(2,:), '-go', ...
    m, upf_rs_est_out(2,:), '-bo', m, mkpf_rs_est_out(2,:), '-ro', 'linewidth', 1);
title('y��');
set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
xlabel('ʵ�����'); ylabel('״̬����ֵ');
% axis([0 M 20 40]);

figure(2);
hold on;
subplot(211);
% plot(m, pf_error(1,:), '-mo', m, ukf_error(1,:), '-r*',  m, ekpf_error(1,:), '-go', m, upf_error(1,:), '-b*', 'linewidth', 1);
plot(m, pf_error(1,:), '-mo', m, ekpf_error(1,:), '-go', m, upf_error(1,:), '-bo', m, mkpf_error(1,:), '-ro', 'linewidth', 1);
title('x��������');
set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
xlabel('ʵ�����'); ylabel('�������MSE');
legend('PF estimate', 'EKPF estimate',  'UPF estimate', 'MKPF estimate');
% axis([0 M 0 1]);

subplot(212);
% plot(m, pf_error(2,:), '-mo', m, ukf_error(2,:), '-r*',m, ekpf_error(2,:), '-go', m, upf_error(2,:), '-b*', 'linewidth', 1);
plot(m, pf_error(2,:), '-mo', m, ekpf_error(2,:), '-go', m, upf_error(2,:), '-bo', m, mkpf_error(2,:), '-ro', 'linewidth', 1);
title('y��������');
set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
xlabel('ʵ�����'); ylabel('�������MSE');
% axis([0 M 0 1]);
