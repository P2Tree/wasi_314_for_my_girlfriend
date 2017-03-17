% ekf3
% X 2*1
% Z 2*1
function [Xf, P0]=tryekf3(X, Z, P, Q, R)

% clc; clear;
% T = 100;
D = 2;  % ״̬ά��

% q = 0.05;
% r = 0.05;
% Q = [q,0; 0,q];    % ״̬��������
% R = [r,0; 0,r];      % �۲���������
F = [1,0; 0,1]; % ״̬ת�ƾ��󣬶��ڲ�ͬģ����Ҫ�޸�
H = [1,0; 0,1]; % �۲���󣬶��ڲ�ͬģ����Ҫ�޸�

%% ״̬��۲�ֵ����
% X = zeros(D,T); % ��ʵ״̬��ʼ��
% X(:,1) = [50; 30]; % (D*1)
% Z = zeros(D,T); % �۲��ʼ��

% for t=2:T
%     X(:,t) = F*X(:,t-1); % + sqrt(Q)*randn(2,1);
% end
% for t=1:T
%     Z(:,t) = H*X(:,t) + sqrt(R)*randn(2,1);
% end

%% kalman filter
% Xf = zeros(2,T);
% Xf(:,1) = X(:,1);   % ��ʼ��kalman����ֵ
Xf = X;
P0 = P;    % ��ʼ��Э�������
% for t = 2:T
%     Xn = F*Xf(:,t-1);   % ״̬Ԥ�� X(k+1|k) = F * X(k|k)
    Xn = F*Xf;
    Zn = H*Xn;  % �۲�Ԥ��
    
    % ����Ƿ�����״̬���̻�����Թ۲ⷽ�̣���Ҫ����������ſ˱Ⱦ���
    % ����F��H���ſ˱Ⱦ�����ʽ
    % ���½���
    
    P = F*P0*F' + Q;   % Э����Ԥ�� P(k+1|k) = F * P(k|k) * F' + G * Q * G'
    
    Kg = P*H'/(H*P*H' + R);  % kalman gain����
    
%     Xf(:,t) = Xn + Kg*(Z(:,t) - Zn);    % kalman״̬����
    Xf = Xn + Kg*(Z - Zn);
    P0 = (eye(D) - Kg*H) * P;   % Э�������
% end

%% ������
% error = sqrt((X - Xf).^2);
% rms = sqrt( 1/T * sum((X - Xf).^2, 2) );

%% ��ͼ
%{
figure
t = 1:T;
hold on; box on;
plot(t, X(1,:), '-k.', t, Xf(1,:), '-r+', t, Z(1,:), 'g*');
plot(t, X(2,:), '-k.', t, Xf(2,:), '-r+', t, Z(2,:), 'g*');
legend('��ʵ״̬', 'kalman ״̬');

figure
hold on; box on;
plot(err_kalmanfilter, '-ks', 'MarkerFace', 'r');
xlabel('ʱ��');
ylabel('���');
%}
