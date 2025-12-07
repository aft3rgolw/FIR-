%% 基于遗传算法(GA)的频率采样法 FIR 滤波器设计 - 最终增强版
% 功能：优化过渡点 T 值，对比优化前后的性能，输出详细指标，自动绘制过渡带宽。

clc; clear; close all;

%% 1. 参数设置
N = 65;         % 滤波器阶数
BW = 4;         % 通带采样点数 (0 ~ BW-1)
M = 1;          % 过渡点数量 (建议设为 1 或 2 进行观察)

% 遗传算法参数
popSize = 250;
maxGen = 500;       
crossoverProb = 0.9;
% 关闭GA绘图以专注于结果图，若想看收敛过程可去掉 'PlotFcn', []
options = optimoptions('ga', 'Display', 'off', 'PlotFcn', []); 

% 优化变量范围
lb = zeros(1, M);   
ub = ones(1, M);    

%% 2. 预计算：未优化的性能指标
fprintf('--------------------------------------------------\n');
fprintf('正在评估未优化(矩形窗)滤波器性能...\n');

% 构建未优化滤波器 (M=0, 无过渡点)
T_unopt = [];
[~, H_unopt_dB, w, ~] = get_filter_response(T_unopt, N, BW, 0);

% 计算未优化的阻带衰减
% 阻带起始位置估算：因为没有过渡点，阻带理论上从 BW 开始
% 但为了保险，我们在 BW 对应的频率后稍微偏移一点开始搜索旁瓣
w_stop_start_unopt = (BW + 0.5) * (2*pi/N); 
stop_indices_unopt = find(w >= w_stop_start_unopt & w <= pi);

if ~isempty(stop_indices_unopt)
    max_sidelobe_unopt = max(H_unopt_dB(stop_indices_unopt));
    As_unopt = -max_sidelobe_unopt; % 转为正数衰减值
else
    As_unopt = 0; % 异常情况
end

fprintf('未优化阻带衰减 (As_raw): %.2f dB\n', As_unopt);


%% 3. 运行遗传算法 (优化)
fprintf('--------------------------------------------------\n');
fprintf('正在运行遗传算法优化 %d 个过渡点...\n', M);

fitnessFunc = @(T) cost_SBMM(T, N, BW, M);

tic;
[T_opt, fval] = ga(fitnessFunc, M, [], [], [], [], lb, ub, [], options);
time_elapsed = toc;

As_opt = -fval; % 优化后的阻带衰减

fprintf('优化完成，耗时: %.2f 秒\n', time_elapsed);
fprintf('最优过渡点值 T: %s\n', mat2str(T_opt, 5));
fprintf('--------------------------------------------------\n');
fprintf('优化后阻带衰减 (As_opt): %.2f dB\n', As_opt);
fprintf('性能提升 (Gain):        %.2f dB\n', As_opt - As_unopt);

%% 4. 构建绘图数据与计算带宽

% 获取优化后的响应曲线
[~, H_opt_dB, ~, ~] = get_filter_response(T_opt, N, BW, M);

% 归一化频率 0~1
norm_w = w / pi; 

% --- 自动计算优化后的过渡带宽 ---
% 1. 寻找通带截止点 (-3dB)
idx_p = find(H_opt_dB <= -3, 1, 'first');
if isempty(idx_p), idx_p = BW; end 
wp = norm_w(idx_p);
val_p = H_opt_dB(idx_p);

% 2. 寻找阻带起始点 (首次达到 -As_opt)
% 从通带截止点往后找，找到第一个小于等于 -As_opt 的点
idx_s = find(H_opt_dB(idx_p:end) <= -As_opt, 1, 'first') + idx_p - 1;
if isempty(idx_s), idx_s = length(H_opt_dB); end
ws = norm_w(idx_s);
val_s = H_opt_dB(idx_s);

% 3. 计算带宽
trans_width = ws - wp;
fprintf('--------------------------------------------------\n');
fprintf('通带截止频率 (wp): %.4f PI\n', wp);
fprintf('阻带起始频率 (ws): %.4f PI\n', ws);
fprintf('过渡带宽 (Delta w): %.4f PI\n', trans_width);
fprintf('--------------------------------------------------\n');

%% 5. 绘图展示
figure('Color', 'w', 'Position', [200, 200, 1000, 600]);

hold on; box on; grid on;

% 1. 绘制未优化曲线 (灰色虚线)
plot(norm_w, H_unopt_dB, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, ...
    'DisplayName', ['未优化 (As=' num2str(As_unopt, '%.1f') 'dB)']);

% 2. 绘制优化后曲线 (蓝色实线)
plot(norm_w, H_opt_dB, 'b-', 'LineWidth', 1.5, ...
    'DisplayName', ['优化后 (M=' num2str(M) ')']);

% 3. 绘制 As 基准线 (红色虚线)
yline(-As_opt, 'r--', 'LineWidth', 1, ...
    'Label', ['As_{opt} = ' num2str(As_opt, '%.1f') 'dB'], ...
    'LabelHorizontalAlignment', 'right', 'DisplayName', '优化目标基准线');

% 4. 标注过渡带宽 (绿色区域 + 箭头)
ylim_vals = [-100 10]; 
fill([wp, ws, ws, wp], [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
    'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '过渡带区域');

% 绘制关键点
plot(wp, val_p, 'go', 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
plot(ws, val_s, 'ro', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');

% 绘制双向箭头
mid_y = -15; 
quiver(wp, mid_y, ws-wp, 0, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off'); 
quiver(ws, mid_y, wp-ws, 0, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off'); 
text((wp+ws)/2, mid_y+3, ['\Delta\omega = ' num2str(trans_width, '%.3f') '\pi'], ...
    'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold');

% 5. 绘制采样点 (Stem)
freq_indices_opt = 0:BW+M-1; 
vals_opt = [ones(1, BW), T_opt];
stem(freq_indices_opt * (2/N), 20*log10(vals_opt+eps), 'b', 'filled', 'MarkerSize', 4, ...
    'BaseValue', -100, 'DisplayName', '优化采样点');

% 图表修饰
title(['频率采样法 FIR 滤波器优化对比 (N=', num2str(N), ')']);
xlabel('归一化频率 (\times\pi rad/sample)');
ylabel('幅度 (dB)');
xlim([0 0.5]); 
ylim([-90 5]);
legend('Location', 'northeast');

hold off;

%% --- 辅助函数 ---

function cost = cost_SBMM(T, N, BW, M)
    [~, H_dB, ~, ~] = get_filter_response(T, N, BW, M);
    
    % 阻带起始位置：从最后一个过渡点(BW+M)之后开始算
    % 为了严谨，取 BW+M 对应的频率作为阻带硬边界
    stop_start_idx = BW + M; 
    w_stop_start = stop_start_idx * (2*pi/N);
    
    n_fft = 4096;
    w_axis = linspace(0, 2*pi, n_fft); % 0~2pi
    
    % 仅在阻带区域 (且在 0~pi 范围内) 找最大值
    stop_indices = find(w_axis >= w_stop_start & w_axis <= pi);
    
    if isempty(stop_indices)
        cost = 0; 
    else
        cost = max(H_dB(stop_indices));
    end
end

function [H_complex, H_dB, w, h_n] = get_filter_response(T, N, BW, M)
    % 构造 H(k)
    H_k = zeros(1, N);
    H_k(1 : BW) = 1;               
    if M > 0
        H_k(BW+1 : BW+M) = T;      
    end
    
    % 共轭对称
    for k = 1 : (N-1)/2
        idx = k + 1;
        sym_idx = N - k + 1;
        H_k(sym_idx) = H_k(idx);
    end
    
    h_n = real(ifft(H_k));
    h_n = fftshift(h_n); 
    
    n_fft = 4096;
    H_complex = fft(h_n, n_fft);
    H_mag = abs(H_complex);
    H_dB = 20 * log10(H_mag + eps);
    w = (0 : n_fft-1) * (2*pi / n_fft); % 0~2pi
    
    % 这里的 w 和 H_dB 是全频段的，用于内部计算
    % 外部绘图时会根据需要截取 0~pi
end