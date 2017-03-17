function [x,P]=ukf3(fstate,x,P,hmeas,z,Q,R)
% UKF   Unscented Kalman Filter for nonlinear dynamic systems 
% [x,P]=ukf(f,x,P,h,z,Q,R) returns state estimate, x and state covariance, P  
% for nonlinear dynamic system (for simplicity, noises are assumed as additive): 
%           x_k+1 = f(x_k) + w_k 
%           z_k   = h(x_k) + v_k 
% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q 
%       v ~ N(0,R) meaning v is gaussian noise with covariance R 
% Inputs:   f: function handle for f(x) 
%           x: "a priori" state estimate 
%           P: "a priori" estimated state covariance 
%           h: fanction handle for h(x) 
%           z: current measurement 
%           Q: process noise covariance  
%           R: measurement noise covariance
% Output:   x: "a posteriori" state estimate 
%           P: "a posteriori" state covariance 
% 
% By Yi Cao at Cranfield University, 04/01/2008 % Modified by JD Liu 2010-4-20

L = numel(x);   % ״̬ά��
m = numel(z);  % �۲�ά��
alpha = 1e-3;
ki = 0;
beta = 2;
lambda = alpha^2 * (L+ki) - L;      % �߶�����
c = L+lambda;       % �߶�����
Wm = [lambda/c 0.5/c+zeros(1,2*L)];     % ��ֵȨ��
Wc = Wm;
Wc(1) = Wc(1) + (1-alpha^2 + beta);     % Э����Ȩ��
c = sqrt(c);
X = sigmas(x, P, c);        % Χ����x״̬��sigma��
[x1, X1, P1, X2] = ut(fstate, X, Wm, Wc, L, Q); % ��״̬��UT�任����
%  X1 = sigmas(x1, P1, c);  % Χ����X1״̬��sigma��
%  X2 = X1 - x1(:, ones(1, size(X1, 2)));  % X1ƫ��
[z1, Z1, P2, Z2] = ut(hmeas, X1, Wm, Wc, m, R); % �Բ�����UT�任����
P12 = X2 * diag(Wc) * Z2';  % ��Э����
K = P12 * inv(P2);
x = x1 + K * (z - z1);  % ״̬����
P = P1 - K * P12';    % Э�������
