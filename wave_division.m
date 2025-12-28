clc;
clear all;

% % 读取csv文件
% data = readtable('024.14号线南约至打运区间-上行-小里程-K029+300-20230917-右轨-res.csv');
% % 把data存为矩阵
% data = table2array(data);
% save('024.14号线南约至打运区间-上行-小里程-K029+300-20230917-右轨-res.mat','data');

% 读取mat文件
data = load('024.14号线南约至打运区间-上行-小里程-K029+300-20230917-右轨-res.mat');

% 记录里程信息
dist = data.data(:,1);
dist = flipud(dist);

% 记录原始表面波形数据
wave = data.data(:,2);
wave = flipud(wave);

% 去噪函数
function wave_denoised = denoise_wave(wave, flag_median, medfilt1_size, ...
    flag_sgolay, sgolayfilt_size, ...
    flag_trend, mileage, elevation_smooth, ...
    flag_highpass, fc)
    % 去除异常值
    if flag_median == 1
        wave_denoised = wave;
    else
        window_size = medfilt1_size; % 窗口大小（奇数），可根据数据密度调整
        wave_denoised = medfilt1(wave, window_size, 'truncate');
    end

    % savitzky-golay 滤波
    if flag_sgolay == 1
        order = 2;        % 多项式阶数
        framelen = sgolayfilt_size;    % 窗口长度（奇数），值越大越平滑
        wave_denoised = sgolayfilt(wave_denoised, order, framelen);
    end
    % 去除趋势
    if flag_trend == 1
        p = polyfit(mileage, elevation_smooth, 1); % 1次线性拟合
        trend = polyval(p, mileage);
        wave_denoised = wave_denoised - trend;
    end
    % 高通滤波
    if flag_highpass == 1
        fc = fc; % 截止频率（归一化），值越小，保留的波长越长
        [b, a] = butter(4, fc, 'high');
        wave_denoised = filtfilt(b, a, wave_denoised);
    end
end

wave_denoised = denoise_wave(wave, 1, 21, 1, 21, 1, dist, wave, 1, 0.01);

%% 计算

% 适用于周期性较强的波磨
L = length(wave_denoised);        % 信号长度
Fs = 1 / (dist(2) - dist(1));   % 空间采样频率 (1/m)
Y = fft(wave_denoised);
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs * (0:(L/2)) / L; % 频率轴 (1/m)

% 找到主频（排除直流分量）
[pks_fft, locs_fft] = findpeaks(P1(2:end), 'SortStr', 'descend', 'NPeaks', 1);
if ~isempty(locs_fft)
    dominant_freq = f(locs_fft + 1); % +1 补偿之前排除的直流分量
    dominant_wavelength = 1 / dominant_freq; % 主波长 (m)
else
    dominant_wavelength = NaN;
end

% 显示频谱
figure;
plot(f, P1, 'b-');
title('Single-Sided Amplitude Spectrum of wave_denoised');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
grid on;

%% 使用带通滤波器
wave_bands = {
    [0.010, 0.030],
    [0.030, 0.100],
    [0.100, 0.300],
    [0.300, 1.000],
    [1.000, 3.000]
};
band_names = {'10-30mm', '30-100mm', '100-300mm', '300-1000mm', '1000-3000mm'};
num_bands = length(wave_bands);

% 关键修正1：确保采样间隔为正
dx = abs(mean(diff(dist)));     % 取绝对值保证采样间隔为正
Fs = 1/dx;                      % 空间采样频率 (1/m)
fprintf('实际平均采样间隔: %.4f m\n', dx);
fprintf('空间采样频率 Fs: %.1f cycles/m\n\n', Fs);

% 存储各波段滤波结果
elevation_bands = cell(num_bands, 1);
% FIR滤波器设计参数
filter_order = 200;  % 滤波器阶数，影响过渡带陡峭度

% figure('Position', [50, 50, 1400, 600]);

for i = 1:num_bands
    f_low = 1 / wave_bands{i}(2);
    f_high = 1 / wave_bands{i}(1);
    f_nyquist = Fs / 2;
    
    Wn = [f_low, f_high] / f_nyquist;
    b = fir1(filter_order, Wn, 'bandpass', hamming(filter_order+1));
    elevation_bands{i} = filtfilt(b, 1, wave_denoised);
end

% 第二步：绘制所有波段单独的子图
figure('Position', [50, 50, 1400, 600]);
for i = 1:num_bands
    subplot(2, 3, i);
    plot(dist, elevation_bands{i});
    xlabel('里程 (m)');
    ylabel('高程 (mm)');
    title(sprintf('波段 %d: %s', i, band_names{i}));
    grid on;
end

% 第六个子图：原始去噪信号
subplot(2, 3, 6);
plot(dist, wave_denoised);
xlabel('里程 (m)');
ylabel('高程 (mm)');
title('去噪后原始信号');
grid on;

% 第三步：绘制所有波段的叠加图
figure('Position', [100, 100, 1200, 400]);
hold on;
% 倒序画图，确保波段顺序正确
for i = num_bands:-1:1
    plot(dist, elevation_bands{i});
end
hold off;

xlabel('里程 (m)');
ylabel('高程 (mm)');
title('所有波段叠加图');
legend(band_names, 'Location', 'best');
grid on;
% 将elevation_bands转置后合并为矩阵
elevation_matrix = cell2mat(elevation_bands');

% 保存合并后的高程数据到mat文件
save('elevation_bands_merged.mat', 'elevation_matrix', '-v7.3');