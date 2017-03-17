% �ο��ԣ�http://studentdavestutorials.weebly.com/particle-filter-with-matlab-code.html
% �����˲�����ָ��״̬������۲���������������
% �����޸�״̬���̺͹۲ⷽ�̣��޸ĳɹ�
% ������Ҫ�����ĳɶ�ά�ģ��޸ĳɹ�
% ���ڽ�����ɺ������޸ĳɹ�
% ������Ҫ����ekf�������ӣ��޸�Ϊekpf�㷨
function [x_est_out error rms]=ekpf(Rs, rs, T, N, q, r)

%% ��ʼ������

D = 2; % ״̬ά��
% q = 0.002; % ״̬������ϵͳ�����������Ӹ���ʱҲ��Ҫ����
% r = 0.003; % �۲�����
% T = 100; % ��������
% N = 100; % ������������������ᵼ�¼���������
Q = q * eye(D);
R = r * eye(D);
pekf = eye(D);

%% ��ʼ�����ӣ� �������ӷֲ�����˹�����µ�״̬��ֵ
% z_ukf = zeros(2, T, N);
P_w = zeros(2, N);

% �����ڳ�ʼ���ӵ�����ֲ���������һ���������
x_P = Rs;

% ��ʼ��һЩ������йصı�������ӵ�һ��

x_P_update = zeros(2, N);
z_update = zeros(2, N);
z_ekf = zeros(2, T, N);

z_out = zeros(2, T); % ��ʵ�Ĺ۲��������
x_out = zeros(2, T); % ��ʵ��״̬�������
x_est = zeros(2, 1); % ÿһ������������״̬��ֵ����
x_est_out = zeros(2, T); % �������̵�����״̬�����������

%% �����˲�������
for t = 1:T
    % ������ʵ״̬���۲�λ�ã��������ʹ�÷�����״̬���̺͹۲ⷽ��
    x = Rs(:,t);% + q * randn(D,1);
    z = rs(:,t);% + r * randn(D,1);
    
    % ���￪ʼ�˲�������
    for i = 1:N
        %�����������飬���Ӹ��¹��̻�����������������ӻ�ȡ�������
        x_P_update(:,i) = x_P(:,i) + q * randn(D,1);
%         x_P_update(:,i) = 0.5 * x_P(:,i) + 25 * x_P(:,i) ./ (1 + x_P(:,i).^2) + 8 * cos(1.2 * (t-1)) + q * randn(D, 1);
        % ���¸��µ����Ӹ��¹۲���������������֪�۲ⷽ�̣����Ը������ӵĹ۲�ֵ����Ҫ��������
        z_update(:,i) = x_P_update(:,i);
%         z_update(:,i) = x_P_update(:,i);
        % �Ե�i������������ʱ�������˲���Ҳ���Ǽ����ϱ����ɵ�����Ϊ�۲������뵽ekf��
        [zekf, pekf] = tryekf3(z, z_update(:,i), pekf, Q, R);
        % ����ÿ�����ӵ�Ȩ�أ�Ȩ�ظ���ʱһ�ָ���ģ�ͣ�����������֪�����ǵĹ۲�����Ǹ�˹�ֲ���Э����Ϊx_R
        P_w(:,i) = (1 / sqrt(2 * pi * r)) * exp(- (z - z_update(:,i)).^2 / (2 * r));
%         P_w(:,i) = P_w(:,i) .* (mean(zekf,2)) + sqrt(pekf)*randn(2,1);
    end
    
    % �ж������˻�
    Psum = sum(P_w,2);
    if Psum(1) == 0 || Psum(2) == 0
        disp('ekpf�㷨�������˻����޷�����');
        break;
    end
    % �淶������Ȩ�أ�ʹ��Ϊ1
    for a = 1:D
        P_w(a,:) = P_w(a,:) ./ sum(P_w(a,:));
    end
    
    % �ز������̣���Ե�ǰ���·ֲ�����������������µĹ�������
    % ��������Ȩ���ز�����ͨ������Ȩ�����ɸ��ʷֲ����ۻ��ֲ���������ѡ���ۼ�ֵ������Ȩ�ش�����ӹ���
    for a = 1:D
        for i = 1:N
            x_P(a,i) = x_P_update(a, find(rand <= cumsum(P_w(a,:)), 1));
        end
    end
    
    %���Ĺ�����������������ͷ���
    x_est = mean(x_P,2);
    
    % ����״̬�͹������ݣ��������
    x_out(:,t) = x;
    z_out(:,t) = z;
    x_est_out(:,t) = x_est;
    
%     x_est_out(:,t) =x_est_out(:,t) +0.3*randn(2,1);
end

error = sqrt((x_out - x_est_out).^2);
rms = sqrt( 1/T * sum((x_out - x_est_out).^2, 2) );

%% ��ͼ�����ս��
% figure(9);
% t = 1:T;
% for d = 1:D
%     subplot(D, 1, d);
%     plot(t, x_out(d,:), '.-b', t, x_est_out(d,:), '-.r', t, z_out(d,:), '*g', 'linewidth', 1);
% end
% set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
% xlabel('time stem'); ylabel('position');
% legend('true position', 'ekpf estimate', 'measurement position');
% 
% figure(10);
% hold on;
% plot(t, error(1,:), '.-b', t, error(2,:), '-.r', 'linewidth', 1);
% set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
% xlabel('time stem'); ylabel('error value');
% legend('rs_x', 'rs_y');
