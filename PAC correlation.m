
fprintf('\n================ 正在进行左右半球 PAC Spearman 相关性分析 (基于半球平均) ================\n');

left_keep_idx  = [5, 6, 7, 8]; 
right_keep_idx = [1, 2, 3, 4]; 

Left_off_mat  = all_pac_z_off(:, left_keep_idx);
Right_off_mat = all_pac_z_off(:, right_keep_idx);
Left_on_mat   = all_pac_z_on(:, left_keep_idx);
Right_on_mat  = all_pac_z_on(:, right_keep_idx);

Subj_L_off = mean(Left_off_mat, 2, 'omitnan');
Subj_R_off = mean(Right_off_mat, 2, 'omitnan');
Subj_L_on  = mean(Left_on_mat, 2, 'omitnan');
Subj_R_on  = mean(Right_on_mat, 2, 'omitnan');


fig_pac_corr = figure('Name', 'Averaged PAC Spearman Correlation', 'Color', 'w', 'Position', [100 150 1000 500]);
sgtitle('Hemispheric Spearman Correlation of Average PAC (OFF vs. ON)', 'FontSize', 16, 'FontWeight', 'bold');


states_corr = {'OFF State', 'ON State'};
data_pairs  = { {Subj_L_off, Subj_R_off}, {Subj_L_on, Subj_R_on} };
colors_corr = {[0.2 0.6 0.8], [0.9 0.6 0.2]}; 

for i = 1:2
    subplot(1, 2, i); hold on;
    
    x_raw = data_pairs{i}{1};
    y_raw = data_pairs{i}{2};
    

    valid_idx = ~isnan(x_raw) & ~isnan(y_raw);
    x_val = x_raw(valid_idx);
    y_val = y_raw(valid_idx);
    
    if length(x_val) > 2
       
        [rho, p_val] = corr(x_val, y_val, 'Type', 'Spearman');
     
        scatter(x_val, y_val, 60, colors_corr{i}, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
        
       
        p_fit = polyfit(x_val, y_val, 1);
        x_fit = linspace(min(x_val), max(x_val), 100);
        y_fit = polyval(p_fit, x_fit);
        plot(x_fit, y_fit, '-k', 'LineWidth', 2);
        

        if p_val < 0.001, p_str = 'p < 0.001'; else, p_str = sprintf('p = %.3f', p_val); end
        txt = sprintf('\\rho = %.3f\n%s\n(N = %d)', rho, p_str, length(x_val));
        
        txt_color = 'k'; 
        if p_val < 0.05, txt_color = 'r'; end
        
        x_lims = xlim; y_lims = ylim;
        text(x_lims(1) + 0.05*diff(x_lims), y_lims(2) - 0.1*diff(y_lims), txt, ...
            'Color', txt_color, 'FontWeight', 'bold', 'FontSize', 12, ...
            'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'none');
    else
        text(0.5, 0.5, '有效样本量不足', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
    

    xlabel('Average Left PAC', 'FontWeight', 'bold');
    ylabel('Average Right PAC', 'FontWeight', 'bold');
    title(states_corr{i}, 'Color', colors_corr{i}, 'FontSize', 14);
    grid on; axis square;
end
fprintf('================ 左右脑 PAC 相关性分析完成 (仅展示 OFF/ON) ================\n');