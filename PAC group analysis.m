
clc; clear; close all;
fprintf('\n================ 正在计算全组平均时频耦合图 ================\n');


Dir_ON  = '';
Dir_OFF = '';
fileList_OFF = dir(fullfile(Dir_OFF, '*.mat'));
srate = 250; 
target_ch = 5; 

if isempty(fileList_OFF)
    error('错误：在 OFF 文件夹中未找到数据，请检查路径！');
end


fp_vec = 4:2:40;   
fa_vec = 40:5:110; 
n_surrogates = 500; 
min_shift = srate;  
num_files = length(fileList_OFF);


all_como_off = NaN(num_files, length(fa_vec), length(fp_vec));
all_como_on  = NaN(num_files, length(fa_vec), length(fp_vec));

fprintf('  -> 目标通道: Ch %d\n', target_ch);
fprintf('  -> 相位频率: %d-%d Hz (共 %d 个频点)\n', fp_vec(1), fp_vec(end), length(fp_vec));
fprintf('  -> 幅值频率: %d-%d Hz (共 %d 个频点)\n', fa_vec(1), fa_vec(end), length(fa_vec));
fprintf('  -> 提示: 分段滤波计算量较大，请耐心等待...\n\n');

for m = 1:num_files
    fprintf('     正在处理被试: %d / %d...\n', m, num_files);
    
   
    tmp_off = load(fullfile(Dir_OFF, fileList_OFF(m).name));
    [dat_off, ~] = helper_get_data(tmp_off);
    
  
    fileName_ON = regexprep(fileList_OFF(m).name, '(?i)off(?=\d*\.mat)', 'on'); 
    file_ON_path = fullfile(Dir_ON, fileName_ON);
    
    if ~exist(file_ON_path, 'file'), file_ON_path = fullfile(Dir_ON, fileList_OFF(m).name); end
    if ~exist(file_ON_path, 'file')
        fprintf('       - 未找到ON数据，跳过。\n');
        continue; 
    end 
    
    tmp_on = load(file_ON_path);
    [dat_on, ~] = helper_get_data(tmp_on);
    
    sig_off = double(dat_off{1,1}); if size(sig_off,1) > size(sig_off,2), sig_off = sig_off'; end
    sig_on  = double(dat_on{1,1});  if size(sig_on,1) > size(sig_on,2),   sig_on = sig_on'; end
    
    if size(sig_off, 1) < target_ch || size(sig_on, 1) < target_ch
        fprintf('       - 数据通道数不足，跳过。\n');
        continue;
    end
    
    dat_off_ch = sig_off(target_ch, :);
    dat_on_ch  = sig_on(target_ch, :);
    
    
    all_como_off(m, :, :) = fast_comodulogram(dat_off_ch, srate, fp_vec, fa_vec, n_surrogates, min_shift);
    all_como_on(m, :, :)  = fast_comodulogram(dat_on_ch, srate, fp_vec, fa_vec, n_surrogates, min_shift);
end


valid_idx = ~isnan(all_como_off(:, 1, 1)) & ~isnan(all_como_on(:, 1, 1));
avg_como_off = squeeze(mean(all_como_off(valid_idx, :, :), 1));
avg_como_on  = squeeze(mean(all_como_on(valid_idx, :, :), 1));
avg_como_diff = avg_como_off - avg_como_on; % OFF - ON

font_name = 'Arial';
axis_lw = 1.2;
roi_lw = 1.5;

fig_como = figure('Name', 'Group Average Comodulogram', 'Color', 'w', 'Position', [100 200 1300 400]);
custom_cmap = custom_redblue_cmap(256);


max_val = max(max(abs([avg_como_off, avg_como_on]))); 
clims = [-max_val, max_val];


ax1 = subplot(1, 3, 1);
contourf(fp_vec, fa_vec, avg_como_off, 100, 'linecolor', 'none'); 
caxis(clims); 
cb1 = colorbar; 
cb1.Label.String = 'PAC (Z score)'; % 更改为 Z score
cb1.Label.FontName = font_name;
colormap(ax1, custom_cmap); 
title('OFF State', 'FontName', font_name, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Phase Frequency (Hz)', 'FontName', font_name, 'FontSize', 11);
ylabel('Amplitude Frequency (Hz)', 'FontName', font_name, 'FontSize', 11);
set(gca, 'FontName', font_name, 'FontSize', 10, 'TickDir', 'out', 'LineWidth', axis_lw, 'Box', 'on', 'Layer', 'top');
rectangle('Position', [13, 50, 17, 30], 'EdgeColor', 'k', 'LineWidth', roi_lw, 'LineStyle', '--');
text(-0.2, 1.05, 'a', 'Units', 'normalized', 'FontName', font_name, 'FontSize', 14, 'FontWeight', 'bold');


ax2 = subplot(1, 3, 2);
contourf(fp_vec, fa_vec, avg_como_on, 100, 'linecolor', 'none');
caxis(clims); 
cb2 = colorbar; 
cb2.Label.String = 'PAC (Z score)';
cb2.Label.FontName = font_name;
colormap(ax2, custom_cmap);
title('ON State', 'FontName', font_name, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Phase Frequency (Hz)', 'FontName', font_name, 'FontSize', 11);
set(gca, 'FontName', font_name, 'FontSize', 10, 'TickDir', 'out', 'LineWidth', axis_lw, 'Box', 'on', 'Layer', 'top');
rectangle('Position', [13, 50, 17, 30], 'EdgeColor', 'k', 'LineWidth', roi_lw, 'LineStyle', '--');
text(-0.15, 1.05, 'b', 'Units', 'normalized', 'FontName', font_name, 'FontSize', 14, 'FontWeight', 'bold');


ax3 = subplot(1, 3, 3);
contourf(fp_vec, fa_vec, avg_como_diff, 100, 'linecolor', 'none');
cb3 = colorbar; 
cb3.Label.String = '\Delta PAC (Z score)';
cb3.Label.FontName = font_name;
max_diff = max(abs(avg_como_diff(:)));
caxis([-max_diff, max_diff]); 
colormap(ax3, custom_cmap);
title('Difference (OFF - ON)', 'FontName', font_name, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Phase Frequency (Hz)', 'FontName', font_name, 'FontSize', 11);
set(gca, 'FontName', font_name, 'FontSize', 10, 'TickDir', 'out', 'LineWidth', axis_lw, 'Box', 'on', 'Layer', 'top');
rectangle('Position', [13, 50, 17, 30], 'EdgeColor', 'k', 'LineWidth', roi_lw, 'LineStyle', '--');
text(-0.15, 1.05, 'c', 'Units', 'normalized', 'FontName', font_name, 'FontSize', 14, 'FontWeight', 'bold');

sgt = sgtitle(sprintf('Group Average PAC Comodulogram (Ch %d, N=%d)', target_ch, sum(valid_idx)));
sgt.FontName = font_name; sgt.FontSize = 14; sgt.FontWeight = 'bold';
fprintf('================ 时频耦合图绘制完成！ ================\n');


function pac_mat = fast_comodulogram(sig, srate, fp_vec, fa_vec, n_surr, min_shift)
    if all(isnan(sig)), pac_mat = NaN(length(fa_vec), length(fp_vec)); return; end
    
    
    is_valid = ~isnan(sig);
    d_mask = diff([0, is_valid, 0]);
    starts = find(d_mask == 1);
    ends = find(d_mask == -1) - 1;
    
    pac_mat = zeros(length(fa_vec), length(fp_vec));
    min_len = min_shift * 2; 
    for p = 1:length(fp_vec)
        f_p = fp_vec(p);
        [b_p, a_p] = butter(3, [max(1, f_p-2), f_p+2] / (srate/2), 'bandpass');
        
        for a = 1:length(fa_vec)
            f_a = fa_vec(a);
            bw_half = 25; 
            low_cut = max(1, f_a - bw_half);
            high_cut = min(srate/2 - 1, f_a + bw_half);
            
            if high_cut - low_cut < 10
                continue; 
            end
            [b_a, a_a] = butter(3, [low_cut, high_cut] / (srate/2), 'bandpass');
            
            all_phase = [];
            all_amp   = [];
            
            
            for k = 1:length(starts)
                seg = sig(starts(k):ends(k));
                if length(seg) >= min_len
                    ph_band = filtfilt(b_p, a_p, seg);
                    am_band = filtfilt(b_a, a_a, seg);
                    all_phase = [all_phase, angle(hilbert(ph_band))]; 
                    all_amp   = [all_amp, abs(hilbert(am_band))];     
            end
            
            if length(all_phase) < min_len
                pac_mat(a, p) = NaN; continue; 
            end
            
            exp_phase = exp(1i * all_phase);
            true_pac = abs(mean(all_amp .* exp_phase));
            
            surr_pac = zeros(1, n_surr);
            N_valid = length(all_amp);
            for k = 1:n_surr
                shift_val = randi([min_shift, N_valid - min_shift]);
                amp_shifted = circshift(all_amp, shift_val);
                surr_pac(k) = abs(mean(amp_shifted .* exp_phase));
            end
            
          
            m_surr = mean(surr_pac);
            std_surr = std(surr_pac);
            if std_surr == 0
                pac_mat(a, p) = true_pac - m_surr;  
            else
                pac_mat(a, p) = (true_pac - m_surr) / std_surr;
            end
        end
    end
end


function cmap = custom_redblue_cmap(m)
    if nargin < 1, m = 256; end 
    c1 = [0.1, 0.3, 0.7]; 
    c2 = [0.95, 0.95, 0.95]; 
    c3 = [0.8, 0.15, 0.2]; 
    
    half_m = ceil(m/2);
    r = [linspace(c1(1), c2(1), half_m)'; linspace(c2(1), c3(1), half_m)'];
    g = [linspace(c1(2), c2(2), half_m)'; linspace(c2(2), c3(2), half_m)'];
    b = [linspace(c1(3), c2(3), half_m)'; linspace(c2(3), c3(3), half_m)'];
    
    cmap = [r(1:m), g(1:m), b(1:m)];
end


function [d, name] = helper_get_data(S)
    if isfield(S, 'data_ica'), d = S.data_ica; name = 'data_ica';
    elseif isfield(S, 'data_new'), d = S.data_new; name = 'data_new';
    else, d = S.data; name = 'data'; end
end



























