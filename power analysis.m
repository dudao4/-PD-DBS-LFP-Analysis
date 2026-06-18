
clc; clear; close all;
fprintf('\n================ 正在计算并绘制 Z-scored PAC ================\n');


Dir_ON  = 'C:\Users\11560\Desktop\双极分析\on\Processed_Step1_Final\Manual_Cleaned';
Dir_OFF = 'C:\Users\11560\Desktop\双极分析\off\Processed_Step1_Final\Manual_Cleaned';
fileList_OFF = dir(fullfile(Dir_OFF, '*.mat'));
srate   = 250; 
n_target_chans = 8;

if isempty(fileList_OFF)
    error('错误：在 OFF 文件夹中未找到数据，请检查路径！');
end


f_phase = [13, 30]; 
f_amp   = [40, 80]; 
n_surrogates = 500; 
min_shift = srate;  
num_files = length(fileList_OFF);


all_pac_z_off = NaN(num_files, n_target_chans);
all_pac_z_on  = NaN(num_files, n_target_chans);


[b_p, a_p] = butter(4, f_phase / (srate/2), 'bandpass');
[b_a, a_a] = butter(4, f_amp / (srate/2), 'bandpass');

fprintf('  -> 正在提取并计算所有被试的 Z-scored PAC (13-30 Hz & 40-80 Hz)...\n');
for m = 1:num_files
 
    tmp_off = load(fullfile(Dir_OFF, fileList_OFF(m).name));
    [dat_off, ~] = helper_get_data(tmp_off);
    
    
  
    fileName_ON = strrep(lower(fileList_OFF(m).name), 'off', 'on'); 
    file_ON_path = fullfile(Dir_ON, fileName_ON);
    if ~exist(file_ON_path, 'file'), file_ON_path = fullfile(Dir_ON, fileList_OFF(m).name); end
    if ~exist(file_ON_path, 'file'), continue; end 
    
    tmp_on = load(file_ON_path);
    [dat_on, ~] = helper_get_data(tmp_on);
    
    sig_off = double(dat_off{1,1}); if size(sig_off,1) > size(sig_off,2), sig_off = sig_off'; end
    sig_on  = double(dat_on{1,1});  if size(sig_on,1) > size(sig_on,2),   sig_on = sig_on'; end
    
   
    for ch = 1:min([size(sig_off,1), size(sig_on,1), n_target_chans])
        all_pac_z_off(m, ch) = compute_zPAC(sig_off(ch, :), b_p, a_p, b_a, a_a, n_surrogates, min_shift, srate);
        all_pac_z_on(m, ch)  = compute_zPAC(sig_on(ch, :), b_p, a_p, b_a, a_a, n_surrogates, min_shift, srate);
    end
end


ch_labels = {'Right (1+3-)', 'Right (2+4-)', 'Right (3+2-)', 'Right (4+1-)', ...
             'Left (5+7-)',  'Left (6+8-)',  'Left (7+6-)',  'Left (8+5-)'};
subplot_pos = [3, 4, 7, 8, 1, 2, 5, 6]; 
ch_colors   = [0.8 0.2 0.2; 0.8 0.2 0.2; 0.8 0.2 0.2; 0.8 0.2 0.2; ... 
               0.2 0.6 0.2; 0.2 0.6 0.2; 0.2 0.6 0.2; 0.2 0.6 0.2];    


fig_indiv = figure('Name', 'Individual Channels PAC', 'Color', 'w', 'Position', [50 100 1400 700]);
sgtitle(sprintf('Hemispheric Mapping of Phase-Amplitude Coupling'), 'FontSize', 16, 'FontWeight', 'bold');
outlier_report = {}; 
for ch = 1:n_target_chans
    subplot(2, 4, subplot_pos(ch)); hold on;
    
    valid_off = all_pac_z_off(:, ch);
    valid_on  = all_pac_z_on(:, ch);
    valid_idx = ~isnan(valid_off) & ~isnan(valid_on);
    v_off_raw = valid_off(valid_idx); 
    v_on_raw  = valid_on(valid_idx);
    orig_indices = find(valid_idx); 
    
    
    if length(v_off_raw) >= 4 
        Q1_off = prctile(v_off_raw, 25); Q3_off = prctile(v_off_raw, 75); IQR_off = Q3_off - Q1_off;
        lb_off = Q1_off - 1.5 * IQR_off; ub_off = Q3_off + 1.5 * IQR_off;
        Q1_on  = prctile(v_on_raw, 25);  Q3_on  = prctile(v_on_raw, 75);  IQR_on  = Q3_on - Q1_on;
        lb_on  = Q1_on - 1.5 * IQR_on;   ub_on  = Q3_on + 1.5 * IQR_on;
        
        keep_idx = (v_off_raw >= lb_off & v_off_raw <= ub_off) & (v_on_raw >= lb_on & v_on_raw <= ub_on);
        outlier_idx = ~keep_idx;
        n_removed = sum(outlier_idx);
        
        if n_removed > 0
            bad_orig_inds = orig_indices(outlier_idx);
            for k = 1:length(bad_orig_inds)
                f_name = fileList_OFF(bad_orig_inds(k)).name;
                outlier_report{end+1} = sprintf('[%s] | 文件: %-25s | OFF: %6.2f | ON: %6.2f', ...
                    ch_labels{ch}, f_name, v_off_raw(outlier_idx(k)), v_on_raw(outlier_idx(k)));
            end
        end
        v_off = v_off_raw(keep_idx); v_on = v_on_raw(keep_idx);
    else
        v_off = v_off_raw; v_on = v_on_raw;
    end
    
    n_valid = length(v_off);
    if n_valid < 3, title(sprintf('%s (N不够)', ch_labels{ch})); continue; end
    
   
    p_val = signrank(v_off, v_on);
    mean_pac = [mean(v_off), mean(v_on)];
    sem_pac  = [std(v_off)/sqrt(n_valid), std(v_on)/sqrt(n_valid)];
    
   
    b = bar([1, 2], mean_pac, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1.2);
    b.CData(1,:) = [0.1 0.45 0.75]; b.CData(2,:) = [0.95 0.5 0.15];
    
    
    errorbar([1, 2], mean_pac, sem_pac, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);
    
   
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'OFF', 'ON'}); xlim([0.2 2.8]); grid on;
    if subplot_pos(ch) == 1 || subplot_pos(ch) == 5, ylabel('PAC (Z-score)', 'FontWeight', 'bold'); end
    
   
    y_max = max(max(v_off), max(v_on)); y_lim = ylim; 
    ylim([min(y_lim(1), -1), y_max + abs(y_max) * 0.3]); 
    
    title(sprintf('%s (n=%d)', ch_labels{ch}, n_valid), 'Color', ch_colors(ch,:));
    
  
    txt_y = y_max + abs(y_max) * 0.1;
    if p_val < 0.05
        sig_star = repmat('*', 1, sum(p_val < [0.05, 0.01, 0.001]));
        text(1.5, txt_y, sprintf('p = %.3f %s', p_val, sig_star), 'HorizontalAlignment', 'center', 'Color', 'r', 'FontWeight', 'bold');
        plot([1.1, 1.9], [txt_y*0.95, txt_y*0.95], '-k', 'LineWidth', 1.2);
    else
        text(1.5, txt_y, sprintf('p = %.3f', p_val), 'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
    end
end

fig_hemi = figure('Name', 'Hemispheric Summary', 'Color', 'w', 'Position', [300 200 800 500]);
sgtitle('Hemispheric PAC Summary', 'FontSize', 16, 'FontWeight', 'bold');

pac_Right_off = mean(all_pac_z_off(:, 1:4), 2, 'omitnan');
pac_Right_on  = mean(all_pac_z_on(:, 1:4),  2, 'omitnan');
pac_Left_off  = mean(all_pac_z_off(:, 5:8), 2, 'omitnan');
pac_Left_on   = mean(all_pac_z_on(:, 5:8),  2, 'omitnan');

hemi_data = {pac_Left_off, pac_Left_on; pac_Right_off, pac_Right_on};
hemi_names = {'Left Hemisphere (Avg Ch 5-8)', 'Right Hemisphere (Avg Ch 1-4)'};
hemi_colors = {[0.2 0.6 0.2], [0.8 0.2 0.2]};

for h = 1:2
    subplot(1, 2, h); hold on;
    
    v_off_raw = hemi_data{h, 1}; 
    v_on_raw  = hemi_data{h, 2};
    valid_idx = ~isnan(v_off_raw) & ~isnan(v_on_raw);
    v_off_raw = v_off_raw(valid_idx); v_on_raw  = v_on_raw(valid_idx);
    
    if length(v_off_raw) >= 4 
        Q1_off = prctile(v_off_raw, 25); Q3_off = prctile(v_off_raw, 75); IQR_off = Q3_off - Q1_off;
        lb_off = Q1_off - 1.5 * IQR_off; ub_off = Q3_off + 1.5 * IQR_off;
        Q1_on  = prctile(v_on_raw, 25);  Q3_on  = prctile(v_on_raw, 75);  IQR_on  = Q3_on - Q1_on;
        lb_on  = Q1_on - 1.5 * IQR_on;   ub_on  = Q3_on + 1.5 * IQR_on;
        keep_idx = (v_off_raw >= lb_off & v_off_raw <= ub_off) & (v_on_raw >= lb_on & v_on_raw <= ub_on);
        v_off = v_off_raw(keep_idx); v_on = v_on_raw(keep_idx);
    else
        v_off = v_off_raw; v_on = v_on_raw;
    end
    
    n_valid = length(v_off);
    if n_valid < 3, continue; end
    
    p_val = signrank(v_off, v_on);
    mean_pac = [mean(v_off), mean(v_on)];
    sem_pac  = [std(v_off)/sqrt(n_valid), std(v_on)/sqrt(n_valid)];
    
    b = bar([1, 2], mean_pac, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1.2);
    b.CData(1,:) = [0.1 0.45 0.75]; b.CData(2,:) = [0.95 0.5 0.15];
    errorbar([1, 2], mean_pac, sem_pac, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);
    plot([1, 2], [v_off, v_on]', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 1);
    scatter(ones(size(v_off)), v_off, 40, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(2*ones(size(v_on)), v_on, 40, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'OFF', 'ON'}); xlim([0.2 2.8]); grid on;
    ylabel('Mean PAC (Z-score)', 'FontWeight', 'bold');
    
  
    y_max = max(max(v_off), max(v_on)); y_lim = ylim; 
    ylim([min(y_lim(1), -1), y_max + abs(y_max) * 0.3]); 
    
    title(sprintf('%s\n(n=%d Patients)', hemi_names{h}, n_valid), 'Color', hemi_colors{h});
    
    txt_y = y_max + abs(y_max) * 0.1;
    if p_val < 0.05
        sig_star = repmat('*', 1, sum(p_val < [0.05, 0.01, 0.001]));
        text(1.5, txt_y, sprintf('p = %.3f %s', p_val, sig_star), 'HorizontalAlignment', 'center', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 12);
        plot([1.1, 1.9], [txt_y*0.95, txt_y*0.95], '-k', 'LineWidth', 1.5);
    else
        text(1.5, txt_y, sprintf('p = %.3f', p_val), 'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
    end
end


fprintf('\n================ 剔除的离群值汇总报告 ================\n');
if isempty(outlier_report)
    fprintf('  数据质量极佳，单个触点图未检测到离群值！\n');
else
    for i = 1:length(outlier_report), fprintf('  %s\n', outlier_report{i}); end
    fprintf('------------------------------------------------------\n');
    fprintf('  共计自动剔除了 %d 个通道异常数据对。\n', length(outlier_report));
end
fprintf('======================================================\n\n');

function z_pac = compute_zPAC(sig, b_p, a_p, b_a, a_a, n_surr, min_shift, fs)
    if all(isnan(sig)), z_pac = NaN; return; end
    
    
    is_valid = ~isnan(sig);
    d_mask = diff([0, is_valid, 0]);
    starts = find(d_mask == 1); 
    ends = find(d_mask == -1) - 1;
    
    all_phase = [];
    all_amp   = [];
    min_len   = fs; 
    
   
    for k = 1:length(starts)
        seg = sig(starts(k):ends(k));
        if length(seg) >= min_len
            ph_band = filtfilt(b_p, a_p, seg);
            am_band = filtfilt(b_a, a_a, seg);
            
            all_phase = [all_phase, angle(hilbert(ph_band))];
            all_amp   = [all_amp, abs(hilbert(am_band))];
        end
    end
    
  
    if length(all_phase) < min_shift * 2
        z_pac = NaN; return; 
    end
    
    N = length(all_amp);
    true_pac = abs(mean(all_amp .* exp(1i * all_phase)));
    surr_pac = zeros(1, n_surr);
    
   
    for k = 1:n_surr
        shift_val = randi([min_shift, N - min_shift]);
        amp_shifted = circshift(all_amp, shift_val);
        surr_pac(k) = abs(mean(amp_shifted .* exp(1i * all_phase)));
    end
    
    mu_surr = mean(surr_pac); std_surr = std(surr_pac);
    if std_surr == 0, z_pac = 0; else, z_pac = (true_pac - mu_surr) / std_surr; end
end

function [d, name] = helper_get_data(S)
    if isfield(S, 'data_ica'), d = S.data_ica; name = 'data_ica';
    elseif isfield(S, 'data_new'), d = S.data_new; name = 'data_new';
    else, d = S.data; name = 'data'; end
end










































clc; clear; close all;


Dir_ON  = 'C:\Users\11560\Desktop\双极分析\on\Processed_Step1_Final\Manual_Cleaned';
Dir_OFF = 'C:\Users\11560\Desktop\双极分析\off\Processed_Step1_Final\Manual_Cleaned';

srate   = 250;
win_len = 2 * srate; 
overlap = round(win_len * 0.5);
nfft    = 2^nextpow2(win_len);
beta_range  = [13, 30];   
total_range = [1, 100];    
n_target_chans = 8;      

ch_labels = {'Right (1+3-)', 'Right (2+4-)', 'Right (3+2-)', 'Right (4+1-)', ...
             'Left (5+7-)',  'Left (6+8-)',  'Left (7+6-)',  'Left (8+5-)'};


subplot_pos = [3, 4, 7, 8, 1, 2, 5, 6]; 

ch_colors   = [0.8 0.2 0.2; 0.8 0.2 0.2; 0.8 0.2 0.2; 0.8 0.2 0.2; ... 
               0.2 0.6 0.2; 0.2 0.6 0.2; 0.2 0.6 0.2; 0.2 0.6 0.2];    


States = {'OFF', 'ON'};
Dirs   = {Dir_OFF, Dir_ON};
Results = struct();

for st = 1:2
    current_state = States{st};
    current_dir   = Dirs{st};
    fileList = dir(fullfile(current_dir, '*.mat'));
    num_files = length(fileList);
    
    rel_power_mat = NaN(num_files, n_target_chans);
    cog_freq_mat  = NaN(num_files, n_target_chans);
    
    fprintf('>>> 正在提取 %s 状态下的数据指标...\n', current_state);
    
    for m = 1:num_files
        tmp = load(fullfile(current_dir, fileList(m).name));
        if isfield(tmp, 'data_ica'), data_struct = tmp.data_ica;
        elseif isfield(tmp, 'data_new'), data_struct = tmp.data_new;
        elseif isfield(tmp, 'data'), data_struct = tmp.data;
        else, error('文件 %s 中找不到变量！', fileList(m).name); end
        
        sig = double(data_struct{1,1});
        if size(sig, 1) > size(sig, 2), sig = sig'; end
        
        [~, f_axis] = pwelch(randn(win_len,1), win_len, overlap, nfft, srate);
        beta_idx  = find(f_axis >= beta_range(1)  & f_axis <= beta_range(2));
        total_idx = find(f_axis >= total_range(1) & f_axis <= total_range(2));
        f_beta    = f_axis(beta_idx);
        
        for ch = 1:min(size(sig,1), n_target_chans)
            ch_data = sig(ch, :);
            is_valid = ~isnan(ch_data);
            d_mask = diff([0, is_valid, 0]);
            starts = find(d_mask == 1); ends = find(d_mask == -1) - 1;
            
            psd_segments = [];
            for k = 1:length(starts)
                segment = ch_data(starts(k):ends(k));
                if length(segment) >= win_len
                    [p, ~] = pwelch(segment, win_len, overlap, nfft, srate);
                    psd_segments = [psd_segments, p];
                end
            end
            
            if ~isempty(psd_segments)
                ch_psd = mean(psd_segments, 2); 
                p_beta = ch_psd(beta_idx);
                p_total = sum(ch_psd(total_idx), 'all');
                rel_power_mat(m, ch) = (sum(p_beta, 'all') / p_total) * 100;
                cog_freq_mat(m, ch) = sum(f_beta .* p_beta(:), 'all') / sum(p_beta, 'all');
            end
        end
    end
    Results.(current_state).RelPower = rel_power_mat;
    Results.(current_state).CoG      = cog_freq_mat;
end

fprintf('\n================ 正在生成个体与通道可视化图表 ================\n');
metric = 'RelPower'; 
ylab = 'Relative Beta Power (%)';

figure('Name', 'Individual Changes', 'Color', 'w', 'Position', [50 100 1400 600]);
sgtitle('Relative Beta Power Distribution Across Contacts', 'FontSize', 16, 'FontWeight', 'bold');

mean_off = zeros(1, n_target_chans); mean_on = zeros(1, n_target_chans);
sem_off  = zeros(1, n_target_chans); sem_on  = zeros(1, n_target_chans);
p_vals   = zeros(1, n_target_chans);

for ch = 1:n_target_chans
    off_data = Results.OFF.(metric)(:, ch); on_data  = Results.ON.(metric)(:, ch);
    valid = ~isnan(off_data) & ~isnan(on_data);
    off_p = off_data(valid); on_p = on_data(valid);
    
    if length(off_p) < 3, continue; end
    
    mean_off(ch) = mean(off_p); mean_on(ch) = mean(on_p);
    sem_off(ch)  = std(off_p)/sqrt(length(off_p)); sem_on(ch) = std(on_p)/sqrt(length(on_p));
    [p_vals(ch), ~, ~] = signrank(off_p, on_p); 
    
    subplot(2, 4, subplot_pos(ch)); hold on;
    for i = 1:length(off_p)
        plot([1, 2], [off_p(i), on_p(i)], '-o', 'Color', [0.8 0.8 0.8 0.6], 'MarkerSize', 4);
    end
    plot([1, 2], [median(off_p), median(on_p)], '-k', 'LineWidth', 2.5);
    scatter([1, 2], [median(off_p), median(on_p)], 50, 'k', 'filled');
    
    set(gca, 'XTick', [1 2], 'XTickLabel', {'OFF', 'ON'}, 'FontSize', 10);
    title(sprintf('%s\n(n=%d, p=%.3f)', ch_labels{ch}, length(off_p), p_vals(ch)), 'Color', ch_colors(ch,:)); 
    grid on;
    if subplot_pos(ch) == 1 || subplot_pos(ch) == 5, ylabel(ylab); end
end

figure('Name', 'Bar Chart: Beta RelPower', 'Color', 'w', 'Position', [150 150 1200 500]);
hold on;
b = bar([mean_off', mean_on'], 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1);
b(1).CData = repmat([0.2 0.6 0.8], n_target_chans, 1); 
b(2).CData = repmat([0.9 0.6 0.2], n_target_chans, 1); 
x1 = b(1).XEndPoints; x2 = b(2).XEndPoints;
errorbar(x1, mean_off, sem_off, 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);
errorbar(x2, mean_on,  sem_on,  'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);

set(gca, 'XTick', 1:n_target_chans, 'XTickLabel', ch_labels, 'FontSize', 10);
xtickangle(45); 
ylabel('Mean Beta Relative Power (%) ± SEM', 'FontWeight', 'bold');
title('Comparison of Beta Relative Power Across All Group Channels (Mean ± SEM)', 'FontSize', 14, 'FontWeight', 'bold');
legend({'OFF State', 'ON State'}, 'Location', 'northeast'); grid on;

for ch = 1:n_target_chans
    if p_vals(ch) < 0.05
        max_h = max(mean_off(ch) + sem_off(ch), mean_on(ch) + sem_on(ch));
        y_line = max_h * 1.15;
        plot([x1(ch), x2(ch)], [y_line, y_line], '-k', 'LineWidth', 1.2);
        if p_vals(ch) < 0.001, star = '***'; elseif p_vals(ch) < 0.01, star = '**'; else, star = '*'; end
        text(mean([x1(ch), x2(ch)]), y_line * 1.05, star, 'FontSize', 14, 'HorizontalAlignment', 'center', 'Color', 'r');
    end
end
ylim([0, max(max(mean_off+sem_off), max(mean_on+sem_on)) * 1.3]); 


fprintf('\n================ 正在进行半球级别聚合分析 ================\n');

Right_off = mean(Results.OFF.(metric)(:, 1:4), 2, 'omitnan');
Right_on  = mean(Results.ON.(metric)(:, 1:4),  2, 'omitnan');
Left_off  = mean(Results.OFF.(metric)(:, 5:8), 2, 'omitnan');
Left_on   = mean(Results.ON.(metric)(:, 5:8),  2, 'omitnan');

hemi_data = {Left_off, Left_on; Right_off, Right_on};
hemi_names = {'Left Hemisphere (Avg Ch 5-8)', 'Right Hemisphere (Avg Ch 1-4)'};
hemi_colors = {[0.2 0.6 0.2], [0.8 0.2 0.2]};

fig_hemi = figure('Name', 'Hemispheric Summary', 'Color', 'w', 'Position', [300 200 800 500]);
sgtitle('Hemispheric Beta Power Summary', 'FontSize', 16, 'FontWeight', 'bold');

for h = 1:2
    subplot(1, 2, h); hold on;
    
    v_off_raw = hemi_data{h, 1}; v_on_raw  = hemi_data{h, 2};
    valid_idx = ~isnan(v_off_raw) & ~isnan(v_on_raw);
    v_off = v_off_raw(valid_idx); v_on  = v_on_raw(valid_idx);
    
    n_valid = length(v_off);
    if n_valid < 3, continue; end
    
    p_val = signrank(v_off, v_on);
    mean_pac = [mean(v_off), mean(v_on)];
    sem_pac  = [std(v_off)/sqrt(n_valid), std(v_on)/sqrt(n_valid)];
    
   
    b = bar([1, 2], mean_pac, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1.2);
    b.CData(1,:) = [0.2 0.6 0.8]; b.CData(2,:) = [0.9 0.6 0.2];
    errorbar([1, 2], mean_pac, sem_pac, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);
    
  
    plot([1, 2], [v_off, v_on]', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 1);
    scatter(ones(size(v_off)), v_off, 40, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(2*ones(size(v_on)), v_on, 40, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'OFF', 'ON'}); xlim([0.2 2.8]); grid on;
    if h == 1, ylabel('Average Relative Beta Power (%)', 'FontWeight', 'bold'); end
    
    y_max = max(max(v_off), max(v_on)); y_lim = ylim; ylim([0, y_max + (y_max * 0.3)]); 
    title(sprintf('%s\n(n=%d Patients)', hemi_names{h}, n_valid), 'Color', hemi_colors{h});
    
   
    txt_y = y_max + (y_max * 0.1);
    if p_val < 0.05
        if p_val < 0.001, sig_star = '***'; elseif p_val < 0.01, sig_star = '**'; else, sig_star = '*'; end
        text(1.5, txt_y, sprintf('p = %.3f %s', p_val, sig_star), 'HorizontalAlignment', 'center', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 12);
        plot([1.1, 1.9], [txt_y*0.95, txt_y*0.95], '-k', 'LineWidth', 1.5);
    else
        text(1.5, txt_y, sprintf('p = %.3f', p_val), 'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
    end
end
fprintf('================ 分析与绘图全部完成！ ================\n');


clc; 
fprintf('\n================ 正在进行基于半球平均的 Spearman 脑耦合分析 ================\n');


Left_off_mat  = Results.OFF.(metric)(:, 5:8);
Right_off_mat = Results.OFF.(metric)(:, 1:4);
Left_on_mat   = Results.ON.(metric)(:, 5:8);
Right_on_mat  = Results.ON.(metric)(:, 1:4);


Subj_L_off = mean(Left_off_mat, 2, 'omitnan');
Subj_R_off = mean(Right_off_mat, 2, 'omitnan');

Subj_L_on  = mean(Left_on_mat, 2, 'omitnan');
Subj_R_on  = mean(Right_on_mat, 2, 'omitnan');


Subj_L_del = Subj_L_on - Subj_L_off;
Subj_R_del = Subj_R_on - Subj_R_off;

fig_corr_avg = figure('Name', 'Averaged Hemispheric Coupling (Spearman)', 'Color', 'w', 'Position', [100 150 1500 500]);
sgtitle('Hemispheric Coupling Analysis (Spearman Rank Correlation on Patient Averages)', 'FontSize', 16, 'FontWeight', 'bold');

data_pairs = { {Subj_L_off, Subj_R_off}, {Subj_L_on, Subj_R_on}, {Subj_L_del, Subj_R_del} };
titles = {'OFF State', 'ON State', '\Delta Change (ON - OFF)'};
colors = {[0.2 0.6 0.8], [0.9 0.6 0.2], [0.5 0.5 0.5]};

for i = 1:3
    subplot(1, 3, i); hold on;
    
    x_raw = data_pairs{i}{1};
    y_raw = data_pairs{i}{2};
    
  
    valid_idx = ~isnan(x_raw) & ~isnan(y_raw);
    x_val = x_raw(valid_idx);
    y_val = y_raw(valid_idx);
    
    if length(x_val) > 2
      
        [rho_val, p_val] = corr(x_val, y_val, 'Type', 'Spearman');
        
       
        scatter(x_val, y_val, 60, colors{i}, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
        
       
        p_fit = polyfit(x_val, y_val, 1);
        x_fit = linspace(min(x_val)*0.9, max(x_val)*1.1, 100);
        y_fit = polyval(p_fit, x_fit);
        plot(x_fit, y_fit, '-k', 'LineWidth', 2);
        
     
        if p_val < 0.001, p_str = 'p < 0.001'; else, p_str = sprintf('p = %.3f', p_val); end
        txt = sprintf('\\rho = %.3f\n%s\n(N = %d Patients)', rho_val, p_str, length(x_val));
        
        if p_val < 0.05
            txt_color = 'r'; font_weight = 'bold';
        else
            txt_color = 'k'; font_weight = 'normal';
        end
        
        x_lims = xlim; y_lims = ylim;
        text(x_lims(1) + 0.05*diff(x_lims), y_lims(2) - 0.1*diff(y_lims), txt, ...
            'Color', txt_color, 'FontWeight', font_weight, 'FontSize', 12, ...
            'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'none');
    else
        text(0.5, 0.5, '有效样本量不足', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
    
   
    if i == 3
        xlabel('Average \Delta Left RelPower (%)', 'FontWeight', 'bold');
        ylabel('Average \Delta Right RelPower (%)', 'FontWeight', 'bold');
        xline(0, 'k--', 'LineWidth', 1, 'Alpha', 0.3);
        yline(0, 'k--', 'LineWidth', 1, 'Alpha', 0.3);
    else
        xlabel('Average Left RelPower (%)', 'FontWeight', 'bold');
        ylabel('Average Right RelPower (%)', 'FontWeight', 'bold');
    end
    title(titles{i}, 'Color', colors{i});
    grid on; axis square;
end

fprintf('================ Spearman 相关性分析与绘图完成！ ================\n');







































clc; clear; close all;
fprintf('');


Dir_ON  = 'C:\Users\11560\Desktop\双极分析\on\Processed_Step1_Final\Manual_Cleaned';
Dir_OFF = 'C:\Users\11560\Desktop\双极分析\off\Processed_Step1_Final\Manual_Cleaned';
fileList_OFF = dir(fullfile(Dir_OFF, '*.mat'));

srate   = 250;
win_len = 2 * srate; 
overlap = round(win_len * 0.5);
nfft    = 2^nextpow2(win_len);
beta_range  = [13, 30]; 
n_target_chans = 8;

if isempty(fileList_OFF), error('错误：在 OFF 文件夹中未找到数据，请检查路径！'); end


ch_labels = {'Right (1+3-)', 'Right (2+4-)', 'Right (3+2-)', 'Right (4+1-)', ...
             'Left (5+7-)',  'Left (6+8-)',  'Left (7+6-)',  'Left (8+5-)'};


subplot_pos = [3, 4, 7, 8, 1, 2, 5, 6]; 
ch_colors   = [0.8 0.2 0.2; 0.8 0.2 0.2; 0.8 0.2 0.2; 0.8 0.2 0.2; ... 
               0.2 0.6 0.2; 0.2 0.6 0.2; 0.2 0.6 0.2; 0.2 0.6 0.2];    


[~, f_axis] = pwelch(randn(win_len,1), win_len, overlap, nfft, srate);
f_len = length(f_axis);
num_files = length(fileList_OFF);


all_psd_off = NaN(num_files, n_target_chans, f_len);
all_psd_on  = NaN(num_files, n_target_chans, f_len);

fprintf('  -> 正在快速提取所有被试的完整频谱曲线...\n');
for m = 1:num_files
   
    tmp_off = load(fullfile(Dir_OFF, fileList_OFF(m).name));
    [dat_off, ~] = helper_get_data(tmp_off);
    
    fileName_ON = strrep(lower(fileList_OFF(m).name), 'off', 'on'); 
    file_ON_path = fullfile(Dir_ON, fileName_ON);
    if ~exist(file_ON_path, 'file'), file_ON_path = fullfile(Dir_ON, fileList_OFF(m).name); end
    if ~exist(file_ON_path, 'file'), continue; end 
    
    tmp_on = load(file_ON_path);
    [dat_on, ~] = helper_get_data(tmp_on);
    
    sig_off = double(dat_off{1,1}); if size(sig_off,1) > size(sig_off,2), sig_off = sig_off'; end
    sig_on  = double(dat_on{1,1});  if size(sig_on,1) > size(sig_on,2),   sig_on = sig_on'; end
    
    for ch = 1:min([size(sig_off,1), size(sig_on,1), n_target_chans])
      
        p_off = helper_calc_psd(sig_off(ch, :), win_len, overlap, nfft, srate);
        p_on  = helper_calc_psd(sig_on(ch, :), win_len, overlap, nfft, srate);
        
        all_psd_off(m, ch, :) = p_off;
        all_psd_on(m, ch, :)  = p_on;
    end
end


plot_idx = f_axis >= 2 & f_axis <= 90;
f_plot = f_axis(plot_idx);


fig_avg_psd = figure('Name', 'Group Average PSD', 'Color', 'w', 'Position', [50 100 1400 700]);
sgtitle(sprintf('Average power spectral density across all group contacts'), ...
    'FontSize', 16, 'FontWeight', 'bold');

for ch = 1:n_target_chans
    subplot(2, 4, subplot_pos(ch)); hold on;
    
   
    valid_mask = ~isnan(all_psd_off(:, ch, 10)) & ~isnan(all_psd_on(:, ch, 10)); 
    valid_off = squeeze(all_psd_off(valid_mask, ch, plot_idx));
    valid_on  = squeeze(all_psd_on(valid_mask, ch, plot_idx));
    
    if size(valid_off, 1) < 3
        title(sprintf('%s (有效样本过少)', ch_labels{ch})); continue; 
    end
    
    mean_off = mean(valid_off, 1); sem_off  = std(valid_off, 0, 1) / sqrt(size(valid_off, 1));
    mean_on  = mean(valid_on, 1);  sem_on   = std(valid_on, 0, 1) / sqrt(size(valid_on, 1));
    
    
    fill([beta_range(1) beta_range(2) beta_range(2) beta_range(1)], [-100 -100 100 100], ...
         [0.93 0.93 0.93], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
   
    color_off = [0.1 0.45 0.75]; color_on = [0.95 0.5 0.15];
    fill([f_plot', fliplr(f_plot')], [mean_off + sem_off, fliplr(mean_off - sem_off)], color_off, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    p_off_line = plot(f_plot, mean_off, 'Color', color_off, 'LineWidth', 2);
    
    fill([f_plot', fliplr(f_plot')], [mean_on + sem_on, fliplr(mean_on - sem_on)], color_on, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    p_on_line = plot(f_plot, mean_on, 'Color', color_on, 'LineWidth', 2);
    
   
    title(sprintf('%s (n=%d)', ch_labels{ch}, size(valid_off,1)), 'Color', ch_colors(ch,:));
    xlim([2, 90]); grid on; set(gca, 'Layer', 'top');
    if subplot_pos(ch) == 1 || subplot_pos(ch) == 5, ylabel('Power (dB/Hz)', 'FontWeight', 'bold'); end
    if subplot_pos(ch) == 1, legend([p_off_line, p_on_line], {'OFF', 'ON'}, 'Location', 'northeast'); end
    
    y_max = max(max(mean_off + sem_off), max(mean_on + sem_on));
    y_min = min(min(mean_off - sem_off), min(mean_on - sem_on));
    if y_max > y_min, ylim([y_min - 2, y_max + 5]); end
end


fprintf('  -> 正在计算半球聚合 PSD...\n');
fig_hemi_psd = figure('Name', 'Hemispheric Average PSD', 'Color', 'w', 'Position', [200 200 1000 500]);
sgtitle('Hemispheric Average PSD', 'FontSize', 16, 'FontWeight', 'bold');


psd_Right_off = squeeze(mean(all_psd_off(:, 1:4, plot_idx), 2, 'omitnan'));
psd_Right_on  = squeeze(mean(all_psd_on(:, 1:4, plot_idx),  2, 'omitnan'));
psd_Left_off  = squeeze(mean(all_psd_off(:, 5:8, plot_idx), 2, 'omitnan'));
psd_Left_on   = squeeze(mean(all_psd_on(:, 5:8, plot_idx),  2, 'omitnan'));

hemi_off_data = {psd_Left_off, psd_Right_off};
hemi_on_data  = {psd_Left_on, psd_Right_on};
hemi_names = {'Left Hemisphere (Avg Ch 5-8)', 'Right Hemisphere (Avg Ch 1-4)'};
hemi_colors_title = {[0.2 0.6 0.2], [0.8 0.2 0.2]};

for h = 1:2
    subplot(1, 2, h); hold on;
    
    valid_off_raw = hemi_off_data{h};
    valid_on_raw  = hemi_on_data{h};
    
   
    valid_mask = ~isnan(valid_off_raw(:, 10)) & ~isnan(valid_on_raw(:, 10)); 
    valid_off = valid_off_raw(valid_mask, :);
    valid_on  = valid_on_raw(valid_mask, :);
    
    mean_off = mean(valid_off, 1); sem_off  = std(valid_off, 0, 1) / sqrt(size(valid_off, 1));
    mean_on  = mean(valid_on, 1);  sem_on   = std(valid_on, 0, 1) / sqrt(size(valid_on, 1));
    
   
    fill([beta_range(1) beta_range(2) beta_range(2) beta_range(1)], [-100 -100 100 100], ...
         [0.93 0.93 0.93], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
   
    fill([f_plot', fliplr(f_plot')], [mean_off + sem_off, fliplr(mean_off - sem_off)], [0.1 0.45 0.75], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    p_off_line = plot(f_plot, mean_off, 'Color', [0.1 0.45 0.75], 'LineWidth', 2.5);
    
  
    fill([f_plot', fliplr(f_plot')], [mean_on + sem_on, fliplr(mean_on - sem_on)], [0.95 0.5 0.15], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    p_on_line = plot(f_plot, mean_on, 'Color', [0.95 0.5 0.15], 'LineWidth', 2.5);

    title(sprintf('%s\n(n=%d Patients)', hemi_names{h}, size(valid_off,1)), 'Color', hemi_colors_title{h});
    xlim([2, 90]); grid on; set(gca, 'Layer', 'top');
    if h == 1, ylabel('Power (dB/Hz)', 'FontWeight', 'bold'); end
    if h == 1, legend([p_off_line, p_on_line], {'OFF State', 'ON State'}, 'Location', 'northeast', 'FontSize', 11); end
    
    y_max = max(max(mean_off + sem_off), max(mean_on + sem_on));
    y_min = min(min(mean_off - sem_off), min(mean_on - sem_on));
    if y_max > y_min, ylim([y_min - 2, y_max + 3]); end
end

fprintf('================ 分析与绘图全部完成！ ================\n');

function [d, name] = helper_get_data(S)
    if isfield(S, 'data_ica'), d = S.data_ica; name = 'data_ica';
    elseif isfield(S, 'data_new'), d = S.data_new; name = 'data_new';
    elseif isfield(S, 'data'), d = S.data; name = 'data';
    else, error('数据结构不匹配'); end
end

function p_db = helper_calc_psd(sig, wl, ov, nf, fs)
    is_val = ~isnan(sig);
    starts = find(diff([0, is_val, 0]) == 1); 
    ends = find(diff([0, is_val, 0]) == -1) - 1;
    psds = [];
    for k = 1:length(starts)
        seg = sig(starts(k):ends(k));
        if length(seg) >= wl
            [p, ~] = pwelch(seg, wl, ov, nf, fs);
            psds = [psds, p];
        end
    end
    if isempty(psds), p_db = nan(nf/2+1, 1);
    else, p_db = 10 * log10(mean(psds, 2) + eps); end
end


































fprintf('\n================ 正在生成个体与通道可视化图表 ================\n');
metric = 'RelPower'; 
ylab = 'Relative Beta Power (%)';


mean_off = zeros(1, n_target_chans); mean_on = zeros(1, n_target_chans);
sem_off  = zeros(1, n_target_chans); sem_on  = zeros(1, n_target_chans);
p_vals   = zeros(1, n_target_chans);


figure('Name', 'Individual Changes', 'Color', 'w', 'Position', [50 100 1400 600]);
sgtitle('Relative Beta Power Distribution Across Contacts', 'FontSize', 16, 'FontWeight', 'bold');

fprintf('\n>>> 正在进行 IQR 离群值检测 (Factor = 1.5) <<<\n');

for ch = 1:n_target_chans
    off_data = Results.OFF.(metric)(:, ch); 
    on_data  = Results.ON.(metric)(:, ch);
    
   
    valid_mask = ~isnan(off_data) & ~isnan(on_data);
    idx_mapping = find(valid_mask); 
    off_p = off_data(valid_mask); 
    on_p = on_data(valid_mask);
    
    if length(off_p) < 3, continue; end
    
   
    q1_off = prctile(off_p, 25);
    q3_off = prctile(off_p, 75);
    iqr_off = q3_off - q1_off;
    lower_off = q1_off - 1.5 * iqr_off; upper_off = q3_off + 1.5 * iqr_off;

   
    q1_on = prctile(on_p, 25);
    q3_on = prctile(on_p, 75);
    iqr_on = q3_on - q1_on;
    lower_on = q1_on - 1.5 * iqr_on; upper_on = q3_on + 1.5 * iqr_on;
    
    non_outlier_mask = (off_p >= lower_off & off_p <= upper_off) & ...
                       (on_p >= lower_on & on_p <= upper_on);
                   
   
    outlier_idx = find(~non_outlier_mask);
    if ~isempty(outlier_idx)
        fprintf('通道 [%s] 检测到离群值并已剔除:\n', ch_labels{ch});
        for out_i = 1:length(outlier_idx)
            local_idx = outlier_idx(out_i);
            orig_idx = idx_mapping(local_idx); 
            fprintf('  -> 样本原始索引 (行号): %d | OFF值: %.2f%% (正常区间: [%.2f, %.2f]) | ON值: %.2f%% (正常区间: [%.2f, %.2f])\n', ...
                orig_idx, off_p(local_idx), lower_off, upper_off, on_p(local_idx), lower_on, upper_on);
        end
    end
    
    
    off_p = off_p(non_outlier_mask);
    on_p = on_p(non_outlier_mask);
    
    if length(off_p) < 3
        fprintf('  [Warning] 通道 %s 剔除离群值后有效样本不足 3 例，跳过计算。\n', ch_labels{ch});
        continue; 
    end
    
   
    mean_off(ch) = mean(off_p); mean_on(ch) = mean(on_p);
    sem_off(ch)  = std(off_p)/sqrt(length(off_p)); sem_on(ch) = std(on_p)/sqrt(length(on_p));
    [p_vals(ch), ~, ~] = signrank(off_p, on_p); 
    
  
    subplot(2, 4, subplot_pos(ch)); hold on;
    for i = 1:length(off_p)
        plot([1, 2], [off_p(i), on_p(i)], '-o', 'Color', [0.8 0.8 0.8 0.6], 'MarkerSize', 4);
    end
    plot([1, 2], [median(off_p), median(on_p)], '-k', 'LineWidth', 2.5);
    scatter([1, 2], [median(off_p), median(on_p)], 50, 'k', 'filled');
    
    set(gca, 'XTick', [1 2], 'XTickLabel', {'OFF', 'ON'}, 'FontSize', 10);
    title(sprintf('%s\n(n=%d, p=%.3f)', ch_labels{ch}, length(off_p), p_vals(ch)), 'Color', ch_colors(ch,:)); 
    grid on;
    if subplot_pos(ch) == 1 || subplot_pos(ch) == 5, ylabel(ylab); end
end
fprintf('>>> IQR 离群值处理完毕 <<<\n\n');


figure('Name', 'Bar Chart: Beta RelPower', 'Color', 'w', 'Position', [150 150 1200 500]);
hold on;
b = bar([mean_off', mean_on'], 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1);
b(1).CData = repmat([0.2 0.6 0.8], n_target_chans, 1); 
b(2).CData = repmat([0.9 0.6 0.2], n_target_chans, 1); 
x1 = b(1).XEndPoints; x2 = b(2).XEndPoints;
errorbar(x1, mean_off, sem_off, 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);
errorbar(x2, mean_on,  sem_on,  'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);
set(gca, 'XTick', 1:n_target_chans, 'XTickLabel', ch_labels, 'FontSize', 10);
xtickangle(45); 
ylabel('Mean Beta Relative Power (%) ± SEM', 'FontWeight', 'bold');
title('Comparison of Beta Relative Power Across All Group Channels (Mean ± SEM)', 'FontSize', 14, 'FontWeight', 'bold');
legend({'OFF State', 'ON State'}, 'Location', 'northeast'); grid on;

for ch = 1:n_target_chans
    if p_vals(ch) == 0, continue; end 
    if p_vals(ch) < 0.05
        max_h = max(mean_off(ch) + sem_off(ch), mean_on(ch) + sem_on(ch));
        y_line = max_h * 1.15;
        plot([x1(ch), x2(ch)], [y_line, y_line], '-k', 'LineWidth', 1.2);
        if p_vals(ch) < 0.001, star = '***'; elseif p_vals(ch) < 0.01, star = '**'; else, star = '*'; end
        text(mean([x1(ch), x2(ch)]), y_line * 1.05, star, 'FontSize', 14, 'HorizontalAlignment', 'center', 'Color', 'r');
    end
end
ylim([0, max(max(mean_off+sem_off), max(mean_on+sem_on)) * 1.3]);