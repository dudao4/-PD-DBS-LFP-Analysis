

clc; clear; close all;
fprintf('=========================================================================\n');
fprintf('  Hemi Delta (ON-OFF) vs Axial Improvement — Spearman Correlation\n');
fprintf('=========================================================================\n');


fprintf('\n--- Step 1: 加载中轴症状数据 ---\n');
clinical = readtable('_clinical_data.csv', 'TextType', 'string');
clinical_pinyin = clinical.pinyin;
clinical_improve_rate = clinical.improvement_rate;
fprintf('  共 %d 位患者有中轴评分\n', height(clinical));


Dir_OFF = fullfile('双极分析', 'off', 'Processed_Step1_Final', 'Manual_Cleaned');
Dir_ON  = fullfile('双极分析', 'on',  'Processed_Step1_Final', 'Manual_Cleaned');
fileList_OFF = dir(fullfile(Dir_OFF, '*.mat'));
srate = 250;
n_target_chans = 8;


f_phase = [13, 30];
f_amp   = [40, 80];
n_surrogates = 500;
min_shift = srate;


win_len = 2 * srate;
overlap = round(win_len * 0.5);
nfft    = 2^nextpow2(win_len);
beta_range  = [13, 30];
total_range = [1, 100];


[b_p, a_p] = butter(4, f_phase/(srate/2), 'bandpass');
[b_a, a_a] = butter(4, f_amp/(srate/2), 'bandpass');

num_files = length(fileList_OFF);
all_pac_z_off = NaN(num_files, n_target_chans);
all_pac_z_on  = NaN(num_files, n_target_chans);
all_rp_off    = NaN(num_files, n_target_chans);
all_rp_on     = NaN(num_files, n_target_chans);
file_stems    = repmat({''}, num_files, 1);


[~, f_axis] = pwelch(randn(win_len, 1), win_len, overlap, nfft, srate);
beta_idx  = find(f_axis >= beta_range(1)  & f_axis <= beta_range(2));
total_idx = find(f_axis >= total_range(1) & f_axis <= total_range(2));


fprintf('\n--- Step 2: 提取 PAC & Beta 相对功率 ---\n');
for m = 1:num_files
    tmp_off = load(fullfile(Dir_OFF, fileList_OFF(m).name));
    [dat_off, ~] = helper_get_data(tmp_off);


    fname_off = lower(fileList_OFF(m).name);
    fname_on  = strrep(fname_off, 'off', 'on');
    file_ON_path = fullfile(Dir_ON, fname_on);
    if ~exist(file_ON_path, 'file')
        file_ON_path = fullfile(Dir_ON, fileList_OFF(m).name);
    end
    if ~exist(file_ON_path, 'file'), continue; end


    if endsWith(fname_off, '2off2.mat')
        stem = extractBefore(fname_off, 'off2.mat');
    elseif endsWith(fname_off, '2off.mat')
        stem = extractBefore(fname_off, '2off.mat');
    else
        [~, stem] = fileparts(fname_off);
    end
    file_stems{m} = stem;

    tmp_on = load(file_ON_path);
    [dat_on, ~] = helper_get_data(tmp_on);

    sig_off = double(dat_off{1,1});
    if size(sig_off,1) > size(sig_off,2), sig_off = sig_off'; end
    sig_on  = double(dat_on{1,1});
    if size(sig_on,1) > size(sig_on,2),   sig_on = sig_on'; end

    n_ch = min([size(sig_off,1), size(sig_on,1), n_target_chans]);
    for ch = 1:n_ch
        all_pac_z_off(m, ch) = compute_zPAC(sig_off(ch,:), b_p, a_p, b_a, a_a, n_surrogates, min_shift, srate);
        all_pac_z_on(m, ch)  = compute_zPAC(sig_on(ch,:),  b_p, a_p, b_a, a_a, n_surrogates, min_shift, srate);
        all_rp_off(m, ch)    = compute_relpow(sig_off(ch,:), win_len, overlap, nfft, srate, beta_idx, total_idx);
        all_rp_on(m, ch)     = compute_relpow(sig_on(ch,:),  win_len, overlap, nfft, srate, beta_idx, total_idx);
    end

    if mod(m, 5) == 0
        fprintf('  已处理 %d/%d 文件...\n', m, num_files);
    end
end
fprintf('  LFP 指标提取完成。\n');


fprintf('\n--- Step 3: 患者匹配 ---\n');
[unique_stems, ~] = unique(file_stems);
unique_stems = unique_stems(~cellfun(@isempty, unique_stems));

matched_pinyin    = {};
matched_axial     = [];
LFP_pac_off_m = [];  LFP_pac_on_m = [];
LFP_rp_off_m  = [];  LFP_rp_on_m  = [];

for i = 1:length(unique_stems)
    stem = unique_stems{i};
    clin_idx = find(strcmpi(clinical_pinyin, stem), 1);
    if isempty(clin_idx), continue; end

    LFP_rows = find(strcmp(file_stems, stem));
    if isempty(LFP_rows), continue; end

    for r = 1:length(LFP_rows)
        row = LFP_rows(r);
        if all(~isnan(all_pac_z_off(row,:))) && all(~isnan(all_pac_z_on(row,:)))
            matched_pinyin{end+1} = stem;
            matched_axial(end+1)  = clinical_improve_rate(clin_idx);
            LFP_pac_off_m(end+1,:) = all_pac_z_off(row,:);
            LFP_pac_on_m(end+1,:)  = all_pac_z_on(row,:);
            LFP_rp_off_m(end+1,:)  = all_rp_off(row,:);
            LFP_rp_on_m(end+1,:)   = all_rp_on(row,:);
            fprintf('  匹配: %s -> axial improve = %.1f%%\n', stem, clinical_improve_rate(clin_idx));
            break;
        end
    end
end

n_subj = length(matched_axial);
fprintf('\n  共匹配 %d 位患者\n', n_subj);
if n_subj < 5
    error('匹配患者数不足 (<5)，无法进行相关性分析。');
end


fprintf('\n--- Step 4: 计算对称性指数及 Delta ---\n');


R_pac_off = mean(LFP_pac_off_m(:, 1:4), 2, 'omitnan');
L_pac_off = mean(LFP_pac_off_m(:, 5:8), 2, 'omitnan');
R_pac_on  = mean(LFP_pac_on_m(:, 1:4),  2, 'omitnan');
L_pac_on  = mean(LFP_pac_on_m(:, 5:8),  2, 'omitnan');

R_rp_off  = mean(LFP_rp_off_m(:, 1:4), 2, 'omitnan');
L_rp_off  = mean(LFP_rp_off_m(:, 5:8), 2, 'omitnan');
R_rp_on   = mean(LFP_rp_on_m(:, 1:4),  2, 'omitnan');
L_rp_on   = mean(LFP_rp_on_m(:, 5:8),  2, 'omitnan');


si_pac_off = (L_pac_off - R_pac_off) ./ (L_pac_off + R_pac_off);
si_pac_on  = (L_pac_on  - R_pac_on)  ./ (L_pac_on  + R_pac_on);
si_rp_off  = (L_rp_off  - R_rp_off)  ./ (L_rp_off  + R_rp_off);
si_rp_on   = (L_rp_on   - R_rp_on)   ./ (L_rp_on   + R_rp_on);


delta_pac = si_pac_off - si_pac_on;
delta_rp  = si_rp_off  - si_rp_on;

fprintf('  PAC SI OFF 范围: [%.3f, %.3f]\n', min(si_pac_off), max(si_pac_off));
fprintf('  PAC SI ON  范围: [%.3f, %.3f]\n', min(si_pac_on),  max(si_pac_on));
fprintf('  PAC Delta  范围: [%.3f, %.3f]\n', min(delta_pac),  max(delta_pac));
fprintf('  Beta SI OFF 范围: [%.3f, %.3f]\n', min(si_rp_off), max(si_rp_off));
fprintf('  Beta SI ON  范围: [%.3f, %.3f]\n', min(si_rp_on),  max(si_rp_on));
fprintf('  Beta Delta  范围: [%.3f, %.3f]\n', min(delta_rp),  max(delta_rp));


fprintf('\n--- Step 5: Spearman 相关分析 + Bootstrap 验证 ---\n');

n_boot = 2000;
rng(42);  


[rho_pac, p_pac] = corr(delta_pac, matched_axial(:), 'Type', 'Spearman');
rho_pac_boot = zeros(n_boot, 1);
for b = 1:n_boot
    idx = randi(n_subj, n_subj, 1);
    rho_pac_boot(b) = corr(delta_pac(idx), matched_axial(idx)', 'Type', 'Spearman');
end
ci_pac = prctile(rho_pac_boot, [2.5, 97.5]);


[rho_beta, p_beta] = corr(delta_rp, matched_axial(:), 'Type', 'Spearman');
rho_beta_boot = zeros(n_boot, 1);
for b = 1:n_boot
    idx = randi(n_subj, n_subj, 1);
    rho_beta_boot(b) = corr(delta_rp(idx), matched_axial(idx)', 'Type', 'Spearman');
end
ci_beta = prctile(rho_beta_boot, [2.5, 97.5]);


fprintf('\n  ===== 结果 =====\n');
fprintf('  PAC  Hemi Delta vs Axial:  rho = %+.4f,  p = %.5f,  95%%CI = [%+.4f, %+.4f]\n', ...
    rho_pac, p_pac, ci_pac(1), ci_pac(2));
fprintf('  Beta Hemi Delta vs Axial:  rho = %+.4f,  p = %.5f,  95%%CI = [%+.4f, %+.4f]\n', ...
    rho_beta, p_beta, ci_beta(1), ci_beta(2));

n_tests = 2; 
alpha_fdr = 0.05;


p_vals = [p_pac; p_beta];
metric_names = {'PAC'; 'Beta'};
[p_sorted, sort_idx] = sort(p_vals);
fdr_thresholds = (1:n_tests)' * alpha_fdr / n_tests;
is_sig_fdr = p_sorted <= fdr_thresholds;


max_sig = find(is_sig_fdr, 1, 'last');
if ~isempty(max_sig)
    is_sig_fdr(1:max_sig) = true;
end

fprintf('\n  --- 多重比较校正 (FDR, Benjamini-Hochberg, q = %.2f) ---\n', alpha_fdr);
fprintf('  检验次数: %d\n', n_tests);
fprintf('  %-6s  p_raw      rank   FDR_thr   FDR_sig\n', 'Metric');

for i = 1:n_tests
    orig_idx = sort_idx(i);
    fprintf('  %-6s  %.5f     %d      %.5f     %d\n', ...
        metric_names{orig_idx}, p_sorted(i), i, fdr_thresholds(i), is_sig_fdr(i));
end

% Report
fh_pac = is_sig_fdr(find(strcmp(metric_names, 'PAC')));
fh_beta = is_sig_fdr(find(strcmp(metric_names, 'Beta')));

if fh_pac
    fprintf('\n  *** PAC Hemi Delta FDR 显著 (q < %.2f) ***\n', alpha_fdr);
else
    fprintf('\n  PAC Hemi Delta FDR 不显著\n');
end
if fh_beta
    fprintf('  *** Beta Hemi Delta FDR 显著 (q < %.2f) ***\n', alpha_fdr);
else
    fprintf('  Beta Hemi Delta FDR 不显著\n');
end


fprintf('\n--- Step 6: 绘制散点图 ---\n');


delta_pac_col     = delta_pac(:);
delta_rp_col      = delta_rp(:);
matched_axial_col = matched_axial(:);

data_pairs = { ...
    {delta_pac_col, matched_axial_col}, ... 
    {delta_rp_col,  matched_axial_col}       
};


rho_vals = [rho_pac, rho_beta];
p_vals   = [p_pac, p_beta];


colors_corr = {[0.2 0.6 0.8], [0.9 0.5 0.2]}; 
titles_corr = {'PAC Hemi Delta vs Axial Improvement', 'Beta Hemi Delta vs Axial Improvement'};
x_labels    = {'PAC Hemi Delta (ON-OFF)  [L-R difference change]', 'Beta Power Hemi Delta (ON-OFF)  [L-R difference change]'};


fig_scatter = figure('Name', 'Hemi Delta vs Axial Analysis', 'Color', 'w', 'Position', [100 200 1200 500]);

for i = 1:2
    subplot(1, 2, i); hold on;
    
    x_raw = data_pairs{i}{1};
    y_raw = data_pairs{i}{2};
    
    
    min_len = min(length(x_raw), length(y_raw));
    x_raw = x_raw(1:min_len);
    y_raw = y_raw(1:min_len);
   
    valid_idx = ~isnan(x_raw) & ~isnan(y_raw);
    x_val = x_raw(valid_idx);
    y_val = y_raw(valid_idx);
    
    if length(x_val) > 2
        
        xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
        
        
        scatter(x_val, y_val, 70, colors_corr{i}, 'filled', ...
            'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.75);
        
        
        p_fit = polyfit(x_val, y_val, 1);
        x_fit = linspace(min(x_val)*1.1, max(x_val)*1.1, 100);
        y_fit = polyval(p_fit, x_fit);
        plot(x_fit, y_fit, '-k', 'LineWidth', 1.5);
        
    
        if p_vals(i) < 0.001
            p_str = 'p < 0.001'; 
        else
            p_str = sprintf('p = %.3f', p_vals(i)); 
        end
        
        
        txt = sprintf('\\rho = %.3f\n%s\n(N = %d)', rho_vals(i), p_str, length(x_val));
        
    
        txt_color = [0.3 0.3 0.3]; 
        if p_vals(i) < 0.05, txt_color = 'r'; end 
        
        x_lims = xlim; y_lims = ylim;
        text(x_lims(1) + 0.05*diff(x_lims), y_lims(2) - 0.12*diff(y_lims), txt, ...
            'Color', txt_color, 'FontWeight', 'bold', 'FontSize', 12, ...
            'BackgroundColor', [1 1 1 0.75], 'EdgeColor', 'none');
    else
        text(0.5, 0.5, '有效样本量不足', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
    
   
    xlabel(x_labels{i}, 'FontWeight', 'bold', 'FontSize', 11);
    ylabel('Axial Improvement Rate (%)', 'FontWeight', 'bold', 'FontSize', 11);
    title(titles_corr{i}, 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold');
    grid on; 
    axis square; 
end


saveas(fig_scatter, 'Fig_HemiDelta_Combined_vs_Axial.png');
fprintf('  已保存: Fig_HemiDelta_Combined_vs_Axial.png\n');



fig3 = figure('Name', 'Bootstrap Distribution', 'Color', 'w', 'Position', [350 50 1000 450]);
sgtitle('Bootstrap Distribution of Spearman \\rho (Hemi Delta vs Axial)', ...
    'FontSize', 13, 'FontWeight', 'bold');

boot_data = {rho_pac_boot, rho_beta_boot};

for i = 1:2
    subplot(1, 2, i); hold on;
    histogram(boot_data{i}, 35, 'FaceColor', colors_corr{i}, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    
    
    xline(rho_vals(i), 'r-', 'LineWidth', 2.5);
    xline(0, 'k-', 'LineWidth', 1);
    
  
    title(sprintf('%s (\\rho = %+.3f)', titles_corr{i}(1:4), rho_vals(i)), 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('Spearman \\rho', 'FontWeight', 'bold'); 
    ylabel('Frequency', 'FontWeight', 'bold'); 
    grid on;
    axis square;
end

saveas(fig3, 'Fig_HemiDelta_Bootstrap.png');
fprintf('  已保存: Fig_HemiDelta_Bootstrap.png\n');
fprintf('================ 左右脑数据相关性及 Bootstrap 分析完成 ================\n');

fprintf('\n--- Step 7: 导出结果 ---\n');
T = table({'PAC'; 'Beta'}, ...
    [rho_pac; rho_beta], ...
    [p_pac; p_beta], ...
    [ci_pac(1); ci_beta(1)], ...
    [ci_pac(2); ci_beta(2)], ...
    [n_subj; n_subj], ...
    'VariableNames', {'Metric', 'Spearman_rho', 'p_value', 'CI95_lo', 'CI95_hi', 'N'});
writetable(T, 'Hemi_Delta_Axial_Results.csv');
fprintf('  已保存: Hemi_Delta_Axial_Results.csv\n');


patient_table = table(matched_pinyin', matched_axial', delta_pac, delta_rp, ...
    'VariableNames', {'pinyin', 'Axial_ImproveRate', 'PAC_HemiDelta', 'Beta_HemiDelta'});
writetable(patient_table, 'Hemi_Delta_PerPatient.csv');
fprintf('  已保存: Hemi_Delta_PerPatient.csv\n');

fprintf('\n==================== Hemi Delta 分析全部完成！ ====================\n');


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

    mu_surr = mean(surr_pac);
    std_surr = std(surr_pac);
    if std_surr == 0, z_pac = 0; else, z_pac = (true_pac - mu_surr) / std_surr; end
end

function rel_pow = compute_relpow(sig, win_len, overlap, nfft, fs, beta_idx, total_idx)
    if all(isnan(sig)), rel_pow = NaN; return; end

    is_valid = ~isnan(sig);
    d_mask = diff([0, is_valid, 0]);
    starts = find(d_mask == 1);
    ends = find(d_mask == -1) - 1;

    psd_segments = [];
    for k = 1:length(starts)
        seg = sig(starts(k):ends(k));
        if length(seg) >= win_len
            [p, ~] = pwelch(seg, win_len, overlap, nfft, fs);
            psd_segments = [psd_segments, p];
        end
    end

    if isempty(psd_segments)
        rel_pow = NaN; return;
    end

    ch_psd   = mean(psd_segments, 2);
    p_beta   = ch_psd(beta_idx);
    p_total  = sum(ch_psd(total_idx), 'all');

    if p_total == 0
        rel_pow = NaN;
    else
        rel_pow = (sum(p_beta, 'all') / p_total) * 100;
    end
end

function [d, name] = helper_get_data(S)
    if isfield(S, 'data_ica'), d = S.data_ica; name = 'data_ica';
    elseif isfield(S, 'data_new'), d = S.data_new; name = 'data_new';
    elseif isfield(S, 'data'), d = S.data; name = 'data';
    else, d = S.data; name = 'data'; end
end