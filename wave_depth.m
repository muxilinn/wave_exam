clc;
clear all;

data = load('elevation_bands_merged.mat');
elevation_matrix = data.elevation_matrix';
depth_matrix = zeros(size(elevation_matrix));
dx = 0.001; % 空间采样间隔 (m)

% 创建里程向量 - 这是你缺失的关键部分！
% 假设你的数据是等间距采样的，每个点间隔0.001米
num_points = size(elevation_matrix, 1);
dist = (0:num_points-1)' * dx; % 从0开始的等间距里程

wave_bands = {
    [0.010, 0.030],
    [0.030, 0.100],
    [0.100, 0.300],
    [0.300, 1.000],
    [1.000, 3.000]
};
band_names = {'10-30mm', '30-100mm', '100-300mm', '300-1000mm', '1000-3000mm'};
num_bands = length(wave_bands);

% 计算每个波段的深度
%% 对各波段计算波长和波深
results = table('Size', [num_bands, 5], 'VariableTypes', {'string', 'double', 'double', 'double', 'double'}, ...
                'VariableNames', {'Band', 'Avg_Wavelength_m', 'Std_Wavelength', 'Avg_Depth_mm', 'RMS_mm'});

for i = 1:num_bands
    signal = elevation_matrix(:,i); % 当前波段信号
    
    % --- 1. 参数调整 (关键修正) ---
    % 计算理论最小峰距（按点数）
    min_peak_distance_points = round((wave_bands{i}(1)*0.8) / dx);
    
    % 修正1：调整MinPeakDistance的逻辑，不再过度限制为2
    % 允许使用理论值，但设置一个合理上限（例如信号长度的1/10）和下限（至少5个点，避免太敏感）
    upper_limit = floor(length(signal) / 10);
    min_peak_distance = max(5, min(min_peak_distance_points, upper_limit));
    
    % 修正2：基于信号幅值（峰峰值）设置更合理的突显度阈值
    signal_peak2peak = peak2peak(signal);
    if signal_peak2peak == 0
        min_prominence = 1e-5; % 极小的默认值
    else
        % 使用峰峰值的1%~5%作为阈值，对微小信号更友好
        min_prominence = signal_peak2peak * 0.02; % 2%
    end
    fprintf('波段 %s: 峰峰值=%.6f, 使用突显度阈值=%.6f\n', band_names{i}, signal_peak2peak, min_prominence);
    
    % --- 2. 波长计算 (时域法) ---
    [pks, locs] = findpeaks(signal, 'MinPeakDistance', min_peak_distance, 'MinPeakProminence', min_prominence);
    
    if length(locs) >= 2
        wavelengths_m = diff(dist(locs)); % 相邻波峰间的里程差即为瞬时波长
        avg_wavelength = mean(wavelengths_m);
        std_wavelength = std(wavelengths_m);
    else
        avg_wavelength = NaN;
        std_wavelength = NaN;
        fprintf('波段 %s: 检测到的波峰数量不足（%d个），无法计算波长。尝试放宽参数...\n', band_names{i}, length(locs));
        % 备选方案：如果没有检测到，尝试不使用MinPeakProminence再试一次
        [pks, locs] = findpeaks(signal, 'MinPeakDistance', min_peak_distance);
        fprintf('波段 %s: 放宽后检测到波峰数=%d\n', band_names{i}, length(locs));
        if length(locs) >= 2
            wavelengths_m = diff(dist(locs));
            avg_wavelength = mean(wavelengths_m);
            std_wavelength = std(wavelengths_m);
        end
    end
    
    % --- 3. 波深计算 (稳健方法) ---
    [valleys, valley_locs] = findpeaks(-signal, 'MinPeakDistance', min_peak_distance, 'MinPeakProminence', min_prominence);
    valleys = -valleys;
    
    depths_mm = [];
    valid_peak_count = 0;
    
    for k = 1:length(locs)
        peak_loc = locs(k);
        left_valleys = valley_locs(valley_locs < peak_loc);
        right_valleys = valley_locs(valley_locs > peak_loc);
        
        if ~isempty(left_valleys) && ~isempty(right_valleys)
            left_idx = left_valleys(end);
            right_idx = right_valleys(1);
            depth = abs(signal(peak_loc) - ...
                       (signal(left_idx) + (signal(right_idx)-signal(left_idx)) * ...
                       (peak_loc-left_idx)/(right_idx-left_idx)));
            depths_mm = [depths_mm; depth];
            valid_peak_count = valid_peak_count + 1;
        end
    end
    
    if ~isempty(depths_mm)
        avg_depth = mean(depths_mm);
    else
        avg_depth = 0;
        fprintf('波段 %s: 无法计算有效的波深。\n', band_names{i});
    end
    
    band_rms = rms(signal);
    
    % --- 4. 存储结果 ---
    results.Band(i) = band_names{i};
    results.Avg_Wavelength_m(i) = avg_wavelength;
    results.Std_Wavelength(i) = std_wavelength;
    results.Avg_Depth_mm(i) = avg_depth;
    results.RMS_mm(i) = band_rms;
    
    % 输出调试信息
    fprintf('波段 %s: 检测到波峰数=%d, 有效波峰数=%d, 使用峰距=%d点\n', ...
            band_names{i}, length(locs), valid_peak_count, min_peak_distance);
end