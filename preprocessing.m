clc; clear; close all;

%% 1. 设置路径
RawDir = 'C:\Users\11560\Desktop\双极分析\on\'; 
SaveDir = fullfile(RawDir, 'Processed_Step1_Final'); 
if ~exist(SaveDir, 'dir'), mkdir(SaveDir); end
DirList = dir(fullfile(RawDir, '*.mat'));
FileNames = {DirList.name};

%% 2. 参数设置
target_srate   = 250;      
winsor_uv      = 1000; 
fallback_srate = 1000; 

%% 3. 数据处理
for m = 1:length(FileNames)
    fullPath = fullfile(RawDir, FileNames{m});
    fprintf('正在处理 %d/%d: %s ... \n', m, length(FileNames), FileNames{m});
    
  
    vars = whos('-file', fullPath);
    if ismember('data_ica', {vars.name})
        load(fullPath, 'data_ica');
        current_data = data_ica;
    elseif ismember('data', {vars.name})
        load(fullPath, 'data');
        current_data = data;
    else
        warning('未找到 data 或 data_ica，跳过：%s', FileNames{m});
        continue;
    end
  
    if size(current_data, 2) < 2
        warning('文件格式不符，未发现 {1,2} 索引，跳过：%s', FileNames{m});
        continue;
    end
    
    rawdata = double(current_data{1,2}); 
    if size(rawdata, 1) > size(rawdata, 2), rawdata = rawdata'; end
  
    srate_orig = [];
    if length(current_data) >= 5
        srate_val = current_data{1,5};
        if ischar(srate_val), srate_orig = str2double(regexp(srate_val, '[\d\.]+', 'match', 'once'));
        elseif isnumeric(srate_val), srate_orig = double(srate_val); end
    end
    if isempty(srate_orig) || isnan(srate_orig) || srate_orig <= 0
        srate_orig = fallback_srate; 
    end
    
  
    rawdata(rawdata > winsor_uv) = winsor_uv;
    rawdata(rawdata < -winsor_uv) = -winsor_uv;
    
    
    if abs(srate_orig - target_srate) > 1
        [P, Q] = rat(target_srate / srate_orig);
        data_resampled = resample(rawdata', P, Q)';
        srate = target_srate; 
    else
        data_resampled = rawdata;
        srate = srate_orig;
    end
    
    
    fprintf('  [工频去噪] 执行 Cleanline...\n');
    try
        EEG = eeg_emptyset();
        EEG.data = data_resampled;
        EEG.srate = srate;
        EEG.pnts = size(data_resampled, 2);
        EEG.nbchan = size(data_resampled, 1);
        EEG.trials = 1;
        
        EEG = pop_cleanline(EEG, 'bandwidth', 2, 'chanlist', 1:EEG.nbchan, ...
            'computepower', 0, 'linefreqs', [50 100], 'normSpectrum', 0, ...
            'p', 0.01, 'pad', 2, 'plotfigures', 0, 'scanforlines', 0, ...
            'sigtype', 'Channels', 'tau', 100, 'verb', 0, 'winsize', 4, 'winstep', 1);
            
        data_final = EEG.data;
    catch
        data_final = data_resampled;
        fprintf('  [Warning] Cleanline 失败，保留原信号。\n');
    end
    
    %% 4. 保存结果
    
    data_to_save = cell(1, 5);
    data_to_save{1,1} = data_final;
    data_to_save{1,5} = srate;
    
   
    if length(current_data) >= 3, data_to_save{1,3} = current_data{1,3}; end
    
    data = data_to_save;
    save(fullfile(SaveDir, FileNames{m}), 'data');
    clear data current_data data_final rawdata; 
end

fprintf('------------------------------------------------\n');
fprintf('处理完成！已从 data{1,2} 提取并完成后续清洗。\n');




























clc; clear; close all;

%% 1. 设置路径
Group1Dir = 'C:\Users\11560\Desktop\双极分析\on\Processed_Step1_Final';
SaveDir = fullfile(Group1Dir, 'Manual_Cleaned'); 
if ~exist(SaveDir, 'dir'), mkdir(SaveDir); end
DirGroup1 = dir(fullfile(Group1Dir, '*.mat'));
FileNamesGroup1 = {DirGroup1.name};

fprintf('======================================================\n');
fprintf('  启动手动剔除模式 (全兼容版：支持 data_ica)\n');
fprintf('  1. 【左键】点击起点 -> 再点终点 -> 标红剔除 (相应波形变淡)\n');
fprintf('  2. 【↑ 上箭头】放大波形 (Zoom In)\n');
fprintf('  3. 【↓ 下箭头】缩小波形 (Zoom Out)\n');
fprintf('  4. 【r 键】重置视图 (Reset View)\n');
fprintf('  5. 【Enter】或【Esc】保存并处理下一个文件\n');
fprintf('======================================================\n');

%% 2. 循环处理
for m = 1:length(FileNamesGroup1)
    fullPath = fullfile(Group1Dir, FileNamesGroup1{m});
    fprintf('正在加载文件 %d/%d: %s ...\n', m, length(FileNamesGroup1), FileNamesGroup1{m});
    
    tmp = load(fullPath);
    
    % ========================================================
    % 【修复 1】：智能安全读取 data，优先识别 data_ica
    % ========================================================
    if isfield(tmp, 'data_ica')
        data = tmp.data_ica;
        var_name = 'data_ica';
    elseif isfield(tmp, 'data_new')
        data = tmp.data_new;
        var_name = 'data_new';
    elseif isfield(tmp, 'data')
        data = tmp.data;
        var_name = 'data';
    else
        error('文件中找不到 data_ica, data_new 或 data 变量，请检查！');
    end
    
    rawdata = double(data{1,1}); 
    if size(rawdata, 1) > size(rawdata, 2), rawdata = rawdata'; end
    [nChans, nPnts] = size(rawdata);
    
    srate_raw = data{1,5};
    if ischar(srate_raw)
        srate = str2double(regexp(srate_raw, '[\d\.]+', 'match', 'once'));
    else
        srate = double(srate_raw); 
    end
    
    time_axis = (0:nPnts-1)/srate;
    
    % ========================================================
    % 作为独立变量读取掩码
    % ========================================================
    if isfield(tmp, 'artifact_mask')
        artifact_mask = tmp.artifact_mask; 
    else
        artifact_mask = false(1, nPnts); 
    end
    
    %% --- 交互式剔除界面 ---
    f = figure('Name', ['Manual Rejection: ' FileNamesGroup1{m}], ...
               'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8], 'Color', 'w');
    
    max_amp = max(abs(rawdata(:)), [], 'omitnan');
    if isempty(max_amp) || max_amp == 0, max_amp = 10; end
    base_ylim = [-max_amp, max_amp] * 1.1; 
    current_ylim = base_ylim; 
    
    finished = false;
    while ~finished
        clf; 
        axes('Position', [0.05 0.05 0.9 0.9]); 
        
        % --- 图层 1 (底层背景)：绘制所有原始波形 ---
        plot(time_axis, rawdata, 'Color', [0.85 0.85 0.85], 'LineWidth', 0.4); hold on;
        
        % --- 图层 2 (顶层前景)：绘制只有有效信号的波形 ---
        data_good_only = rawdata;
        data_good_only(:, artifact_mask) = NaN;
        plot(time_axis, data_good_only, 'Color', [0.2 0.5 0.8 0.7], 'LineWidth', 0.6); 
        
        % --- 画红框 (显示被标记为伪迹的区域) ---
        ylim(current_ylim); 
        y_lims = get(gca, 'YLim');
        patch_y = [y_lims(1), y_lims(2), y_lims(2), y_lims(1)];
        
        d_mask = diff([0, artifact_mask, 0]);
        starts = find(d_mask == 1);
        ends   = find(d_mask == -1) - 1;
        
        for k = 1:length(starts)
            x_start = time_axis(starts(k));
            x_end   = time_axis(ends(k));
            fill([x_start, x_end, x_end, x_start], patch_y, 'r', ...
                 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
        end
        
        title({['文件: ' FileNamesGroup1{m} ' (变量: ' var_name ')'], ...
               '【左键】标红剔除 (剔除段波形将变淡) | 【↑】放大 | 【↓】缩小 | 【r】重置 | 【Enter】保存'}, ...
               'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
        xlabel('Time (s)', 'FontSize', 11); 
        ylabel('Voltage (\muV)', 'FontSize', 11);
        axis tight; 
        ylim(current_ylim); 
        grid on;
        
        % --- 2. 交互部分 ---
        try
            [x_in, ~, button] = ginput(1); 
        catch
            button = []; 
        end
        
        y_center = mean(current_ylim);
        y_range = diff(current_ylim);
        
        if isempty(button) || button == 13 || button == 27 
            finished = true;
        elseif button == 30 % ↑ 放大
            current_ylim = y_center + [-y_range/2, y_range/2] * 0.7; 
        elseif button == 31 % ↓ 缩小
            current_ylim = y_center + [-y_range/2, y_range/2] * 1.4; 
        elseif button == 114 || button == 82 % r / R 重置
            current_ylim = base_ylim;
        elseif button == 1 % 左键
            xline(x_in, 'r-', 'LineWidth', 1.5); 
            [x_out, ~] = ginput(1); 
            
            t_start = min(x_in, x_out);
            t_end   = max(x_in, x_out);
            idx_start = max(1, round(t_start * srate));
            idx_end   = min(nPnts, round(t_end * srate));
            
            artifact_mask(idx_start:idx_end) = true;
            fprintf('  -> 标红伪迹区间: %.2fs - %.2fs\n', t_start, t_end);
        end
    end
    
    if ishandle(f), close(f); end 
    
    %% 3. 保存
    % ========================================================
    % 【修复 2】：动态变量名保存，绝不篡改原数据结构
    % ========================================================
    if strcmp(var_name, 'data_ica')
        data_ica = data;
        save(fullfile(SaveDir, FileNamesGroup1{m}), 'data_ica', 'artifact_mask');
    elseif strcmp(var_name, 'data_new')
        data_new = data;
        save(fullfile(SaveDir, FileNamesGroup1{m}), 'data_new', 'artifact_mask');
    else
        save(fullfile(SaveDir, FileNamesGroup1{m}), 'data', 'artifact_mask');
    end
    
    fprintf('  文件 %s 已保存 (包含独立的 artifact_mask 变量)。\n\n', FileNamesGroup1{m});
end
fprintf('全部手动剔除完成！\n');