clc;
clear all;
close all;

%%

% 시간 변수 정의
dt = 0.01; % 시간 간격 (초)
t_end = 30; % 시뮬레이션 끝 시간 (초)
t = 0:dt:t_end; % 시간 벡터

% 층 수와 자유도 정의
num_floors = 70; % 층 수
num_dofs = num_floors; % 각 층마다 1개의 자유도

% 시스템 매개변수 정의
% 질량 행렬 (kg)
mass_value = 10; % 각 층의 질량 (kg)
M = diag(ones(num_dofs,1) * mass_value); 

% 강성 행렬 (N/m)
stiffness_value = 1e5; % 각 층의 강성 (N/m)
K = diag(2*stiffness_value*ones(num_dofs,1)) - diag(stiffness_value*ones(num_dofs-1,1), 1) - diag(stiffness_value*ones(num_dofs-1,1), -1);

% 감쇠 행렬 (Ns/m)
damping_value = 1000; % 각 층의 감쇠 (Ns/m)
C = diag(2*damping_value*ones(num_dofs,1)) - diag(damping_value*ones(num_dofs-1,1), 1) - diag(damping_value*ones(num_dofs-1,1), -1);

% 층별 풍하중 정의
% 층별 하중을 다르게 설정하기 위해 각 층에 대해 서로 다른 크기의 하중을 설정
% 예를 들어, 하중을 각 층에 비례해서 다르게 설정
F_layers = zeros(num_floors, length(t)); % 층별 하중 초기화

% 각 층에 대해 하중 설정 (사인파 하중 예제)
frequency = 0.2; % 주파수 (Hz)
base_amplitude = 10000; % 기본 진폭 (N)
for layer = 1:num_floors
    % 층별로 하중의 크기를 다르게 설정 (층 수에 따라 진폭을 증가)
    amplitude = base_amplitude * (layer / num_floors); 
    F_layers(layer, :) = amplitude * sin(2 * pi * frequency * t); % 각 층의 사인파 하중
end

% 초기 조건
x = zeros(length(t), num_dofs); % 위치 초기화
v = zeros(length(t), num_dofs); % 속도 초기화
a = zeros(length(t), num_dofs); % 가속도 초기화

% Newmark-beta 방법 매개변수
beta = 0.25; % Newmark-beta 방법의 beta 값
gamma = 0.5; % Newmark-beta 방법의 gamma 값

% 초기 가속도 계산
F_current = F_layers(:, 1); % 초기 하중 벡터
a(1, :) = M \ (F_current - C * v(1, :)' - K * x(1, :)');

% 적분 수행
for i = 1:length(t) - 1
    % 시간 간격
    dt = t(i+1) - t(i);
    
    % 외력 벡터 (각 층별 하중 적용)
    F_current = F_layers(:, i+1);
    
    % 가속도, 속도, 위치 업데이트
    x(i+1, :) = x(i, :) + dt * v(i, :) + dt^2 * ((0.5 - beta) * a(i, :) + beta * a(i+1, :));
    v(i+1, :) = v(i, :) + dt * ((1 - gamma) * a(i, :) + gamma * a(i+1, :));
    a(i+1, :) = M \ (F_current - C * v(i+1, :)' - K * x(i+1, :)');
end

% 결과 시각화 (꼭대기층)
top_floor_index = num_dofs; % 가장 상단 층 인덱스 (70층)
figure;
plot(t, x(:, top_floor_index));
title(sprintf('제일 꼭대기층 (%d층)의 응답', top_floor_index));
xlabel('시간 (초)');
ylabel('위치 (m)');
grid on;