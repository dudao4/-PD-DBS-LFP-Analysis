
clc; clear; close all;
fprintf('=======================================================================\n');
fprintf('  Brady_Hemi: Rigidity+Bradykinesia Rate vs Hemi LFP Delta\n');
fprintf('=======================================================================\n');


fprintf('\n--- Step 1: 加载临床数据 ---\n');
brady_clinical = readtable('_bradykinesia_clinical_data.csv', 'TextType', 'string');
clinical_pinyin = brady_clinical.pinyin;
R_rate = brady_clinical.R_brady_improve_rate;
L_rate = brady_clinical.L_brady_improve_rate;
fprintf('  共 %d 位患者\n', height(brady_clinical));


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
    if mod(m, 5) == 0, fprintf('  已处理 %d/%d 文件...\n', m, num_files); end
end
fprintf('  LFP 指标提取完成。\n');

fprintf('\n--- Step 3: 患者匹配 ---\n');
[unique_stems, ~] = unique(file_stems);
unique_stems = unique_stems(~cellfun(@isempty, unique_stems));

matched_pinyin  = {};
matched_R_rate  = [];
matched_L_rate  = [];
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
            matched_R_rate(end+1)  = R_rate(clin_idx);
            matched_L_rate(end+1)  = L_rate(clin_idx);
            LFP_pac_off_m(end+1,:) = all_pac_z_off(row,:);
            LFP_pac_on_m(end+1,:)  = all_pac_z_on(row,:);
            LFP_rp_off_m(end+1,:)  = all_rp_off(row,:);
            LFP_rp_on_m(end+1,:)   = all_rp_on(row,:);
            fprintf('  %s: R_rate=%.1f%%, L_rate=%.1f%%\n', stem, R_rate(clin_idx), L_rate(clin_idx));
            break;
        end
    end
end

n_subj = length(matched_R_rate);
fprintf('\n  共匹配 %d 位患者\n', n_subj);
if n_subj < 5, error('匹配患者数不足 (<5)。'); end


fprintf('\n--- Step 4: 计算半球 Δ (ON-OFF) ---\n');
R_pac_off = mean(LFP_pac_off_m(:, 1:4), 2, 'omitnan');
L_pac_off = mean(LFP_pac_off_m(:, 5:8), 2, 'omitnan');
R_pac_on  = mean(LFP_pac_on_m(:, 1:4),  2, 'omitnan');
L_pac_on  = mean(LFP_pac_on_m(:, 5:8),  2, 'omitnan');
R_rp_off  = mean(LFP_rp_off_m(:, 1:4), 2, 'omitnan');
L_rp_off  = mean(LFP_rp_off_m(:, 5:8), 2, 'omitnan');
R_rp_on   = mean(LFP_rp_on_m(:, 1:4),  2, 'omitnan');
L_rp_on   = mean(LFP_rp_on_m(:, 5:8),  2, 'omitnan');

delta_R_pac = R_pac_on - R_pac_off;
delta_L_pac = L_pac_on - L_pac_off;
delta_R_rp  = R_rp_off  - R_rp_on;
delta_L_rp  = L_rp_off  - L_rp_on;

fprintf('\n--- Step 5: Spearman 相关 + Bootstrap ---\n');
n_boot = 2000;
rng(42);


test_configs = {
    delta_R_pac, matched_L_rate,  'R_PAC_vs_L_Rate';
    delta_L_pac, matched_R_rate,  'L_PAC_vs_R_Rate';
    delta_R_rp,  matched_L_rate,  'R_Beta_vs_L_Rate';
    delta_L_rp,  matched_R_rate,  'L_Beta_vs_R_Rate';
};

n_tests = 4;
rho_vals = zeros(n_tests, 1);
p_vals   = zeros(n_tests, 1);
ci_lo    = zeros(n_tests, 1);
ci_hi    = zeros(n_tests, 1);
boot_cells = cell(n_tests, 1);
metric_names = test_configs(:, 3);

for t = 1:n_tests
    LFP_delta = test_configs{t, 1};
    clin_out  = test_configs{t, 2};
    [rho_vals(t), p_vals(t)] = corr(LFP_delta, clin_out(:), 'Type', 'Spearman');
    rho_boot = zeros(n_boot, 1);
    for b = 1:n_boot
        idx = randi(n_subj, n_subj, 1);
        rho_boot(b) = corr(LFP_delta(idx), clin_out(idx)', 'Type', 'Spearman');
    end
    ci_lo(t) = prctile(rho_boot, 2.5);
    ci_hi(t) = prctile(rho_boot, 97.5);
    boot_cells{t} = rho_boot;
    fprintf('  %-22s: rho = %+.4f,  p = %.5f,  95%%CI = [%+.4f, %+.4f]\n', ...
        metric_names{t}, rho_vals(t), p_vals(t), ci_lo(t), ci_hi(t));
end

alpha_fdr = 0.05;
[p_sorted, sort_idx] = sort(p_vals);
fdr_thresholds = (1:n_tests)' * alpha_fdr / n_tests;
is_sig_fdr = p_sorted <= fdr_thresholds;
max_sig = find(is_sig_fdr, 1, 'last');
if ~isempty(max_sig), is_sig_fdr(1:max_sig) = true;
else, is_sig_fdr(:) = false; end

fprintf('\n  --- FDR 校正 (BH, n=%d, q=%.2f) ---\n', n_tests, alpha_fdr);
fprintf('  %-22s  p_raw      rank   FDR_thr   FDR_sig\n', 'Metric');
for i = 1:n_tests
    orig_idx = sort_idx(i);
    fprintf('  %-22s  %.5f     %d      %.5f     %d\n', ...
        metric_names{orig_idx}, p_sorted(i), i, fdr_thresholds(i), is_sig_fdr(i));
end
fprintf('\n  FDR 结果:\n');
for t = 1:n_tests
    if is_sig_fdr(find(sort_idx == t))
        fprintf('  *** %-20s FDR 显著 ***\n', metric_names{t});
    else
        fprintf('  %-20s FDR 不显著\n', metric_names{t});
    end
end


fprintf('\n--- Step 6: 绘制散点图与统计图 ---\n');


delta_R_pac_col = delta_R_pac(:);
delta_L_pac_col = delta_L_pac(:);
delta_R_rp_col  = delta_R_rp(:);
delta_L_rp_col  = delta_L_rp(:);
matched_L_rate_col = matched_L_rate(:);
matched_R_rate_col = matched_R_rate(:);


data_pairs = { ...
    {delta_R_pac_col, matched_L_rate_col}, ... 
    {delta_L_pac_col, matched_R_rate_col}, ... 
    {delta_R_rp_col,  matched_L_rate_col}, ... 
    {delta_L_rp_col,  matched_R_rate_col}  ... 
};

colors_boot = {[0.2 0.6 0.8], [0.2 0.8 0.6], [0.9 0.5 0.2], [1.0 0.6 0.2]};
x_labels    = {'R-Hemi PAC \Delta', 'L-Hemi PAC \Delta', 'R-Hemi Beta \Delta', 'L-Hemi Beta \Delta'};
y_labels    = {'L-Brady+Rig Improvement Rate (%)', 'R-Brady+Rig Improvement Rate (%)', ...
               'L-Brady+Rig Improvement Rate (%)', 'R-Brady+Rig Improvement Rate (%)'};
titles_corr = {'R-Hemi PAC vs L-Brady+Rig Rate', 'L-Hemi PAC vs R-Brady+Rig Rate', ...
               'R-Hemi Beta vs L-Brady+Rig Rate', 'L-Hemi Beta vs R-Brady+Rig Rate'};



fig1 = figure('Name', 'PAC', 'Color', 'w', 'Position', [50 200 1100 480]);
sgtitle(sprintf('PAC \\Delta vs Bradykinesia+Rigidity Improve Rate (N=%d)', n_subj), ...
    'FontSize', 13, 'FontWeight', 'bold');

for t = 1:2 
    subplot(1, 2, t); hold on;
    
    x_raw = data_pairs{t}{1}; y_raw = data_pairs{t}{2};
    min_len = min(length(x_raw), length(y_raw));
    x_raw = x_raw(1:min_len); y_raw = y_raw(1:min_len);
    
    valid_idx = ~isnan(x_raw) & ~isnan(y_raw);
    x_val = x_raw(valid_idx); y_val = y_raw(valid_idx);
    
    if length(x_val) > 2
        xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
        scatter(x_val, y_val, 65, colors_boot{t}, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.75);
        
        p_fit = polyfit(x_val, y_val, 1);
        x_fit = linspace(min(x_val)*1.1, max(x_val)*1.1, 100);
        y_fit = polyval(p_fit, x_fit);
        plot(x_fit, y_fit, '-k', 'LineWidth', 1.5);
        
        if p_vals(t) < 0.001, p_str = 'p < 0.001'; else, p_str = sprintf('p = %.3f', p_vals(t)); end
        txt = sprintf('\\rho = %.3f\n%s\n(N = %d)', rho_vals(t), p_str, length(x_val));
        
        txt_color = [0.3 0.3 0.3]; if p_vals(t) < 0.05, txt_color = 'r'; end
        x_lims = xlim; y_lims = ylim;
        text(x_lims(1) + 0.05*diff(x_lims), y_lims(2) - 0.12*diff(y_lims), txt, ...
            'Color', txt_color, 'FontWeight', 'bold', 'FontSize', 11, ...
            'BackgroundColor', [1 1 1 0.75], 'EdgeColor', 'none');
    else
        text(0.5, 0.5, '有效样本量不足', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
    xlabel(x_labels{t}, 'FontWeight', 'bold', 'FontSize', 10);
    ylabel(y_labels{t}, 'FontWeight', 'bold', 'FontSize', 10);
    title(titles_corr{t}, 'Color', 'k', 'FontSize', 11, 'FontWeight', 'bold');
    grid on; axis square;
end
saveas(fig1, 'Fig_BradyPAC_Total.png');
fprintf('  已保存: Fig_BradyPAC_Total.png\n');



fig2 = figure('Name', 'Beta', 'Color', 'w', 'Position', [100 250 1100 480]);
sgtitle(sprintf('Beta \\Delta vs Bradykinesia+Rigidity Improve Rate (N=%d)', n_subj), ...
    'FontSize', 13, 'FontWeight', 'bold');

for t = 3:4 
    subplot(1, 2, t-2); hold on; 
    
    x_raw = data_pairs{t}{1}; y_raw = data_pairs{t}{2};
    min_len = min(length(x_raw), length(y_raw));
    x_raw = x_raw(1:min_len); y_raw = y_raw(1:min_len);
    
    valid_idx = ~isnan(x_raw) & ~isnan(y_raw);
    x_val = x_raw(valid_idx); y_val = y_raw(valid_idx);
    
    if length(x_val) > 2
        xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
        scatter(x_val, y_val, 65, colors_boot{t}, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.75);
        
        p_fit = polyfit(x_val, y_val, 1);
        x_fit = linspace(min(x_val)*1.1, max(x_val)*1.1, 100);
        y_fit = polyval(p_fit, x_fit);
        plot(x_fit, y_fit, '-k', 'LineWidth', 1.5);
        
        if p_vals(t) < 0.001, p_str = 'p < 0.001'; else, p_str = sprintf('p = %.3f', p_vals(t)); end
        txt = sprintf('\\rho = %.3f\n%s\n(N = %d)', rho_vals(t), p_str, length(x_val));
        
        txt_color = [0.3 0.3 0.3]; if p_vals(t) < 0.05, txt_color = 'r'; end
        x_lims = xlim; y_lims = ylim;
        text(x_lims(1) + 0.05*diff(x_lims), y_lims(2) - 0.12*diff(y_lims), txt, ...
            'Color', txt_color, 'FontWeight', 'bold', 'FontSize', 11, ...
            'BackgroundColor', [1 1 1 0.75], 'EdgeColor', 'none');
    else
        text(0.5, 0.5, '有效样本量不足', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
    xlabel(x_labels{t}, 'FontWeight', 'bold', 'FontSize', 10);
    ylabel(y_labels{t}, 'FontWeight', 'bold', 'FontSize', 10);
    title(titles_corr{t}, 'Color', 'k', 'FontSize', 11, 'FontWeight', 'bold');
    grid on; axis square;
end
saveas(fig2, 'Fig_BradyBeta_Total.png');
fprintf('  已保存: Fig_BradyBeta_Total.png\n');



fig3 = figure('Name', 'Bootstrap', 'Color', 'w', 'Position', [50 50 1000 750]);
sgtitle('Bootstrap Distributions', 'FontSize', 13, 'FontWeight', 'bold');

for t = 1:4
    subplot(2, 2, t); hold on;
    histogram(boot_cells{t}, 30, 'FaceColor', colors_boot{t}, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    
    xline(rho_vals(t), 'r-', 'LineWidth', 2.5);
    xline(0, 'k-', 'LineWidth', 1);
    
    clean_name = strrep(metric_names{t}, '_', '\_');
    title(sprintf('%s (\\rho = %+.3f)', clean_name, rho_vals(t)), 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('Spearman \\rho', 'FontWeight', 'bold'); ylabel('Freq', 'FontWeight', 'bold'); 
    grid on; axis square;
end
saveas(fig3, 'Fig_BradyTotal_Bootstrap.png');
fprintf('  已保存: Fig_BradyTotal_Bootstrap.png\n');



fig4 = figure('Name', 'Combined', 'Color', 'w', 'Position', [150 150 650 450]);
hold on;

for t = 1:4
    bar(t, rho_vals(t), 'FaceColor', colors_boot{t}, 'EdgeColor', 'k', 'LineWidth', 1.2);
    
    if p_vals(t) < 0.001, p_str = 'p < 0.001'; else, p_str = sprintf('p = %.3f', p_vals(t)); end
    sig_str = sprintf('\\rho = %.3f\n%s', rho_vals(t), p_str);
    
    y_offset = 0.05 * sign(rho_vals(t));
    if rho_vals(t) >= 0, v_align = 'bottom'; else, v_align = 'top'; end
    
    text(t, rho_vals(t) + y_offset, sig_str, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', v_align, ...
        'FontWeight', 'bold', 'FontSize', 9);
end

set(gca, 'XTick', 1:4, 'XTickLabel', strrep(metric_names, '_', '\_'), 'FontSize', 9, 'FontWeight', 'bold');
xtickangle(25);
yline(0, 'k-', 'LineWidth', 1.2);
ylabel('Spearman \\rho', 'FontWeight', 'bold', 'FontSize', 11);
title(sprintf('Hemi \\Delta vs Brady+Rig Rate (N=%d)', n_subj), 'FontWeight', 'bold', 'FontSize', 12);
grid on; axis square;

saveas(fig4, 'Fig_BradyTotal_Combined.png');
fprintf('  已保存: Fig_BradyTotal_Combined.png\n');
fprintf('================ 左右半脑症状相关性分析全部完成 ================\n');


fprintf('\n--- Step 7: 导出 ---\n');
T = table(metric_names, rho_vals, p_vals, ci_lo, ci_hi, repmat(n_subj, n_tests, 1), ...
    'VariableNames', {'Metric', 'Spearman_rho', 'p_value', 'CI95_lo', 'CI95_hi', 'N'});
writetable(T, 'Brady_Hemi_Results.csv');
fprintf('  Brady_Hemi_Results.csv\n');

patient_table = table(matched_pinyin', matched_R_rate', matched_L_rate', ...
    delta_R_pac, delta_L_pac, delta_R_rp, delta_L_rp, ...
    'VariableNames', {'pinyin', 'R_BradyRig_Rate', 'L_BradyRig_Rate', ...
    'R_PAC_Delta', 'L_PAC_Delta', 'R_Beta_Delta', 'L_Beta_Delta'});
writetable(patient_table, 'Brady_Hemi_PerPatient.csv');
fprintf('  Brady_Hemi_PerPatient.csv\n');

fprintf('\n==================== 分析完成！ ====================\n');


function plot_scatter_sub(ax, x_data, y_data, rho_val, p_val, ci_lo, ci_hi, xlbl, tt, mc)
    axes(ax); hold on;
    scatter(x_data, y_data, 60, mc, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    p_fit = polyfit(x_data, y_data, 1);
    xr = linspace(min(x_data)*1.1, max(x_data)*1.1, 100);
    plot(xr, polyval(p_fit, xr), '-k', 'LineWidth', 1.5);
    xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
    xlabel(xlbl, 'FontWeight', 'bold', 'FontSize', 10);
    ylabel('Improve Rate (%)', 'FontWeight', 'bold', 'FontSize', 10);
    title(tt, 'FontSize', 11, 'FontWeight', 'bold');
    if p_val < 0.001, ps = 'p < 0.001';
    elseif p_val < 0.01, ps = sprintf('p = %.4f', p_val);
    else, ps = sprintf('p = %.3f', p_val); end
    tc = ternary(p_val < 0.05, 'r', [0.3 0.3 0.3]);
    xl2 = xlim; yl2 = ylim;
    text(xl2(1)+0.05*diff(xl2), yl2(2)-0.08*diff(yl2), ...
        sprintf('\\rho=%.3f %s\nCI[%.2f,%.2f]', rho_val, ps, ci_lo, ci_hi), ...
        'Color', tc, 'FontWeight', 'bold', 'FontSize', 9);
    grid on;
end

function s = ternary(cond, t, f)
    if cond, s = t; else, s = f; end
end

function z_pac = compute_zPAC(sig, b_p, a_p, b_a, a_a, n_surr, min_shift, fs)
    if all(isnan(sig)), z_pac = NaN; return; end
    iv = ~isnan(sig); dm = diff([0, iv, 0]);
    starts = find(dm == 1); ends = find(dm == -1) - 1;
    ap = []; aa = []; ml = fs;
    for k = 1:length(starts)
        seg = sig(starts(k):ends(k));
        if length(seg) >= ml
            ap = [ap, angle(hilbert(filtfilt(b_p, a_p, seg)))];
            aa = [aa, abs(hilbert(filtfilt(b_a, a_a, seg)))];
        end
    end
    if length(ap) < min_shift*2, z_pac = NaN; return; end
    N = length(aa);
    tp = abs(mean(aa .* exp(1i*ap)));
    sp = zeros(1, n_surr);
    for k = 1:n_surr
        sp(k) = abs(mean(circshift(aa, randi([min_shift,N-min_shift])) .* exp(1i*ap)));
    end
    mu_s = mean(sp); sd_s = std(sp);
    if sd_s == 0, z_pac = 0; else, z_pac = (tp - mu_s)/sd_s; end
end

function rp = compute_relpow(sig, wl, ol, nf, fs, bi, ti)
    if all(isnan(sig)), rp = NaN; return; end
    iv = ~isnan(sig); dm = diff([0, iv, 0]);
    starts = find(dm == 1); ends = find(dm == -1) - 1;
    psds = [];
    for k = 1:length(starts)
        seg = sig(starts(k):ends(k));
        if length(seg) >= wl
            [p, ~] = pwelch(seg, wl, ol, nf, fs);
            psds = [psds, p];
        end
    end
    if isempty(psds), rp = NaN; return; end
    cpsd = mean(psds, 2); pt = sum(cpsd(ti), 'all');
    if pt == 0, rp = NaN; else, rp = sum(cpsd(bi), 'all')/pt*100; end
end

function [d, name] = helper_get_data(S)
    if isfield(S, 'data_ica'), d = S.data_ica; name = 'data_ica';
    elseif isfield(S, 'data_new'), d = S.data_new; name = 'data_new';
    elseif isfield(S, 'data'), d = S.data; name = 'data';
    else, d = S.data; name = 'data'; end
end