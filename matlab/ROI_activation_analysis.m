% Episodic memory task for adolescents
% :: ROI activation - - analysis (figures and statistics)
% :: code written by Kahyun Choi


%% define parameters

% load subject list
subject_path = 'C:\subject_info_directory';
load(fullfile(subject_path, 'subject_list_for_GLM.mat')); % 'bhv_list', 'sbj_list', 'sbj_nums', 'num_sbj'

% load ROI activation files
dir_metric = 'C:\ROI_metric_directory';
dir_seg = 'C:\native_ROI_directory';
dir_beta = 'C:\native_beta_directory';
load(fullfile(dir_metric, sprintf('metric_neural_N%d.mat', num_sbj))); % metric_list_enc, flag_name_mri, roi_name_all
disp('Metric data loaded');

% behavior
load(fullfile(subject_path, 'behavioral_data.mat')); % 'em_acc', 'sbj_age', 'sbj_sex'

load(fullfile(subject_path, 'survey_data.mat'), 'sbj_svy', 'svy_label');

% ETI
dir_survey = fullfile(subject_path, 'ETI_survey.xlsx');
sbj_eti = readtable(dir_survey);
sbj_eti = sbj_eti(sbj_nums, :);
sbj_eti = table2struct(sbj_eti); % ALL, GT, PA, EA, Q1~21
eti_ea = [sbj_eti.EA];


%% hippocampal activation: group difference

roi_idx_L = [ 1 4 7 ]; % Hippocampus, ant, post
roi_idx_R = [ 18 21 24 ];
roi_idx_B = [ 35 38 41 ];
roi_idx_list = {roi_idx_L, roi_idx_R, roi_idx_B};

con_i = 2; % 2:control
flag_i = 1; % 1:all

group_idx = (eti_ea > 0);
tmp_group = group_idx;
activation_ylim = [-2.5 3];
plot_color = [0.3 0.3 0.3; 223/255 153/255 153/255];

%
fig_draw = true;
if fig_draw
    figure;
    set(gcf, 'Position', [100 100 420 150]);

    tmp_category = categorical(group_idx, logical([0 1]), {'non-abused', 'abused'});

    % - subregions -
    s_names = {'Ant', 'Post'};
    for s_i = 1:2
        subplot(1,6,[2*(s_i-1)+1 2*(s_i)]);
    
        data_all = metric_list_enc{con_i}{flag_i}{roi_idx_B(s_i+1)};

        fig_box = boxchart(data_all, 'GroupByColor', tmp_category); hold on
        for g_i = 1:length(fig_box)
            fig_box(g_i).BoxEdgeColor = 'k';
            fig_box(g_i).BoxFaceColor = plot_color(g_i, :);
            fig_box(g_i).BoxFaceAlpha = 0.7;
            fig_box(g_i).MarkerStyle = 'none';
            fig_box(g_i).LineWidth = 0.8;
        end
        fig_dot = []; dot_density = 3; % 0, 1, 2, 3
        bar_width_half = 0.3 - dot_density * 0.07;
    
        for graph_i = 1:size(fig_box, 1)
            curr_data = data_all(tmp_group == graph_i-1);
    
            y_draw = sort(curr_data); 
            y_draw = y_draw(~isnan(y_draw));
            y_unique = unique(y_draw);
            num_repeated = arrayfun(@(x) sum(x==y_draw),y_unique);
            x_draw = arrayfun(@(x) linspace(-bar_width_half, bar_width_half, x - mod(x,2)),num_repeated,'UniformOutput',false);
            x_draw(mod(num_repeated,2)==1) = cellfun(@(x) [x 0], x_draw(mod(num_repeated,2)==1),'UniformOutput',false);
            x_draw = cell2mat(reshape(x_draw,1,[])) + graph_i/2 + 0.05;
        end
        yline(0, 'k:')
        ylim(activation_ylim); ylabel('beta coefficient');
        xticklabels(s_names{s_i});
        [~, p] = ttest2(data_all(tmp_group == 0), data_all(tmp_group > 0));
    
        c_stats = meanEffectSize(data_all(tmp_group == 0), data_all(tmp_group > 0), 'Effect', 'Cohen');
        fprintf('\n%s: Cohens d = %.3f\n', s_names{s_i}, c_stats.Effect);
        
        hold off
        box off
        set(gca,'LineWidth',0.6);
        set(gca,'FontName','Arial','FontSize',6, 'FontWeight','normal')

        if s_i > 1
            ax = gca;
            set(ax.YAxis, 'Visible', 'off');
        end
    end
end

% ---------- statistics ---------- %
stats_on = true;
if stats_on
    fprintf('<< hippocampal activation >>\n');
    fprintf(' - no emotional abuse -\n');
    for i = 1:3
        [~,p,~,stats] = ttest(metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(1)}(group_idx == 0));
        fprintf('%10s: T(%d)=%.3f, p=%.3f\n', roi_name_all{roi_idx_list{i}(1)}, stats.df, stats.tstat, p);
    end
    disp(' ');
    
    fprintf(' - emotional abuse -\n');
    for i = 1:3
        [~,p,~,stats] = ttest(metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(1)}(group_idx == 1));
        fprintf('%10s: T(%d)=%.3f, p=%.3f\n', roi_name_all{roi_idx_list{i}(1)}, stats.df, stats.tstat, p);
    end
    disp(' ');
    
    fprintf(' - group diff -\n');
    for i = 1:3
        [~,p,~,stats] = ttest2(metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(1)}(group_idx == 0), metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(1)}(group_idx == 1));
        fprintf('%10s: T(%d)=%.3f, p=%.3f\n', roi_name_all{roi_idx_list{i}(1)}, stats.df, stats.tstat, p);
    end
    disp(' ');
    
    fprintf('<< hp subregion activation >>\n');
    fprintf(' - no emotional abuse -\n');
    for i = 1:3
        for j = 2:3
            [~,p,~,stats] = ttest(metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(j)}(group_idx == 0));
            fprintf('%10s: T(%d)=%.3f, p=%.3f\n', roi_name_all{roi_idx_list{i}(j)}, stats.df, stats.tstat, p);
        end
    end
    disp(' ');
    
    fprintf(' - emotional abuse -\n');
    for i = 1:3
        for j = 2:3
            [~,p,~,stats] = ttest(metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(j)}(group_idx == 1));
            fprintf('%10s: T(%d)=%.3f, p=%.3f\n', roi_name_all{roi_idx_list{i}(j)}, stats.df, stats.tstat, p);
        end
    end
    disp(' ');
    
    fprintf(' - group diff -\n');
    for i = 1:3
        for j = 2:3
            [~,p,~,stats] = ttest2(metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(j)}(group_idx == 0), metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(j)}(group_idx == 1));
            fprintf('%10s: T(%d)=%.3f, p=%.3f\n', roi_name_all{roi_idx_list{i}(j)}, stats.df, stats.tstat, p);
        end
    end
    disp(' ');
end

% ---------- paired comparison ---------- %
for g_i = 0:1
    for hemis_i = 1:3
        ant_beta = metric_list_enc{con_i}{flag_i}{roi_idx_list{hemis_i}(2)};
        post_beta = metric_list_enc{con_i}{flag_i}{roi_idx_list{hemis_i}(3)};
        [~, p, ~, stats] = ttest(ant_beta(group_idx == g_i), post_beta(group_idx == g_i));
        fprintf('group %d, hemisphere %d, Ant vs Post: t(%d)=%.3f, p=%.4f\n', ...
            g_i, hemis_i, stats.df, stats.tstat, p);
    end
end

% ---------- Linear mixed-effects model ---------- %
lme_analysis = true;
if lme_analysis
    tmp_beta = [metric_list_enc{2}{1}{roi_idx_list{1}(1)} metric_list_enc{2}{1}{roi_idx_list{2}(1)}];
    tmp_type = [ones(1, length(group_idx)) ones(1, length(group_idx))*2];
    tmp_age = sbj_age; tmp_age = [tmp_age tmp_age];
    tmp_sex = double(sbj_sex); tmp_sex = [tmp_sex tmp_sex];
    tmp_acc = em_acc{6}; tmp_acc = [tmp_acc tmp_acc]; % full EM
    tmp_eti = eti_ea; tmp_eti = [tmp_eti tmp_eti]; % continuous
    tmp_sbj = 1:length(group_idx); tmp_sbj = [tmp_sbj tmp_sbj];
    tmp_names = {'beta', 'hemis', 'age', 'sex', 'acc', 'eti', 'sbj'};
    tmp_cell = {tmp_beta', tmp_type', tmp_age', tmp_sex', tmp_acc', tmp_eti', tmp_sbj'};
    tmp_array = cell2mat(tmp_cell);
    tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
    tmp_table.hemis = categorical(tmp_table.hemis);

    %
    lme = fitlme(tmp_table, 'beta~hemis*eti+(1|sbj)')
    anova(lme)

    % - print results
    fprintf('<< LME results >>\n');
    for i = 2:length(lme.CoefficientNames)
        tmp = ' ';
        if lme.Coefficients{i,6} < 0.05
            tmp = '* ';
        end
        fprintf('%2s%18s: β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', tmp, lme.Coefficients{i,1}, ...
            lme.Coefficients{i,2}, lme.Coefficients{i,3}, double(lme.Coefficients{i,5}), lme.Coefficients{i,4}, lme.Coefficients{i,6})
    end
end



%% hippocampal activation (subregion): non-abused

roi_idx_L = [ 1 4 7 ]; % Hippocampus, ant, post
roi_idx_R = [ 18 21 24 ];
roi_idx_B = [ 35 38 41 ];
roi_idx_list = {roi_idx_L, roi_idx_R, roi_idx_B};

con_i = 2; % 2:control
flag_i = 1; % 1:all
hemis_i = 3; % 3:bilateral

group_idx = (eti_ea > 0);

% ---------- statistics ---------- %
stats_on = true;
if stats_on
    fprintf('<< hippocampal correlation w/ age >>\n');
    fprintf(' - no emotional abuse -\n');
    for i = 1:3
        for j = 1:3
            tmp_data = metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(j)};
            [r,p] = corr(sbj_age(group_idx == 0)', tmp_data(group_idx == 0)');
            fprintf('%10s: r=%.3f, p=%.3f\n', roi_name_all{roi_idx_list{i}(j)}, r, p);
        end
    end
    disp(' ');
end


% ---------- Linear mixed-effects model ---------- %
lme_analysis = true;
if lme_analysis            
    tmp_beta = arrayfun(@(x) metric_list_enc{2}{1}{x}, roi_idx_list{hemis_i}(2:3), 'uni', 0); % enc-ctrl, all
    tmp_beta = [tmp_beta{1}(group_idx == 0) tmp_beta{2}(group_idx == 0)];
    tmp_type = [ones(1, sum(group_idx == 0)) ones(1, sum(group_idx == 0))*2];
    tmp_age = sbj_age(group_idx == 0); tmp_age = [tmp_age tmp_age];    
    tmp_sex = double(sbj_sex(group_idx == 0)); tmp_sex = [tmp_sex tmp_sex];
    tmp_acc = em_acc{6}(group_idx == 0); tmp_acc = [tmp_acc tmp_acc]; % full EM
    tmp_eti = eti_ea(group_idx == 0); tmp_eti = [tmp_eti tmp_eti]; % continuous
    tmp_sbj = 1:sum(group_idx == 0); tmp_sbj = [tmp_sbj tmp_sbj];
    tmp_names = {'beta', 'subregion', 'age', 'sex', 'acc', 'eti', 'sbj'};
    tmp_cell = {tmp_beta', tmp_type', tmp_age', tmp_sex', tmp_acc', tmp_eti', tmp_sbj'};
    tmp_array = cell2mat(tmp_cell);
    tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
    tmp_table.subregion = categorical(tmp_table.subregion);

    %
    lme = fitlme(tmp_table, 'beta~subregion*age+(1|sbj)')
    anova(lme)

    % - print results
    fprintf('<< LME results >>\n');
    for i = 2:length(lme.CoefficientNames)
        tmp = ' ';
        if lme.Coefficients{i,6} < 0.05
            tmp = '* ';
        end
        fprintf('%2s%18s: β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', tmp, lme.Coefficients{i,1}, ...
            lme.Coefficients{i,2}, lme.Coefficients{i,3}, double(lme.Coefficients{i,5}), lme.Coefficients{i,4}, lme.Coefficients{i,6})
    end
end


%% hippocampal activation (subregion): abused

roi_idx_L = [ 1 4 7 ]; % Hippocampus, ant, post
roi_idx_R = [ 18 21 24 ];
roi_idx_B = [ 35 38 41 ];
roi_idx_list = {roi_idx_L, roi_idx_R, roi_idx_B};

con_i = 2; % 2:control
flag_i = 1; % 1:all
hemis_i = 3; % 3:bilateral

group_idx = (eti_ea > 0);

% ---------- stats ---------- %
stats_on = true;
if stats_on
    fprintf('<< hippocampal correlation w/ age >>\n');
    fprintf(' - emotional abuse -\n');
    for i = 1:3
        for j = 1:3
            tmp_data = metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(j)};
            [r,p] = corr(sbj_age(group_idx == 1)', tmp_data(group_idx == 1)');
            fprintf('%10s: r=%.3f, p=%.3f\n', roi_name_all{roi_idx_list{i}(j)}, r, p);
        end
    end
    disp(' ');
end

% ---------- Linear mixed-effects model ---------- %
lme_analysis = true;
if lme_analysis            
    tmp_beta = arrayfun(@(x) metric_list_enc{2}{1}{x}, roi_idx_list{hemis_i}(2:3), 'uni', 0); % enc-ctrl, all
    tmp_beta = [tmp_beta{1}(group_idx == 1) tmp_beta{2}(group_idx == 1)];
    tmp_type = [ones(1, sum(group_idx == 1)) ones(1, sum(group_idx == 1))*2];
    tmp_age = sbj_age(group_idx == 1); tmp_age = [tmp_age tmp_age];
    tmp_sex = double(sbj_sex(group_idx == 1)); tmp_sex = [tmp_sex tmp_sex];
    tmp_acc = em_acc{6}(group_idx == 1); tmp_acc = [tmp_acc tmp_acc]; % full EM
    tmp_eti = eti_ea(group_idx == 1); tmp_eti = [tmp_eti tmp_eti]; % continuous
    tmp_sbj = 1:sum(group_idx == 1); tmp_sbj = [tmp_sbj tmp_sbj];
    tmp_names = {'beta', 'subregion', 'age', 'sex', 'acc', 'eti', 'sbj'};
    tmp_cell = {tmp_beta', tmp_type', tmp_age', tmp_sex', tmp_acc', tmp_eti', tmp_sbj'};
    tmp_array = cell2mat(tmp_cell);
    tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
    tmp_table.subregion = categorical(tmp_table.subregion);
    
    %
    lme = fitlme(tmp_table, 'beta~subregion*age+(1|sbj)')
    anova(lme)

    % - print results
    fprintf('<< LME results >>\n');
    for i = 2:length(lme.CoefficientNames)
        tmp = ' ';
        if lme.Coefficients{i,6} < 0.05
            tmp = '* ';
        end
        fprintf('%2s%18s: β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', tmp, lme.Coefficients{i,1}, ...
            lme.Coefficients{i,2}, lme.Coefficients{i,3}, double(lme.Coefficients{i,5}), lme.Coefficients{i,4}, lme.Coefficients{i,6})
    end
end


%% hippocampal activation (subregion): abuse interaction

roi_idx_L = [ 1 4 7 ]; % Hippocampus, ant, post
roi_idx_R = [ 18 21 24 ];
roi_idx_B = [ 35 38 41 ];
roi_idx_list = {roi_idx_L, roi_idx_R, roi_idx_B};

con_i = 2; % 2:control
flag_i = 1; % 1:all
hemis_i = 3; % 3:bilateral

group_idx = (eti_ea > 0);

% ---------- Linear mixed-effects model ---------- %
lme_analysis = true;
if lme_analysis            
    tmp_beta = arrayfun(@(x) metric_list_enc{2}{1}{x}, roi_idx_list{hemis_i}(2:3), 'uni', 0); % enc-ctrl, all
    tmp_beta = [tmp_beta{1} tmp_beta{2}];
    tmp_type = [ones(1, length(group_idx)) ones(1, length(group_idx))*2];
    tmp_age = sbj_age; tmp_age = [tmp_age tmp_age];
    tmp_sex = double(sbj_sex); tmp_sex = [tmp_sex tmp_sex];
    tmp_acc = em_acc{6}; tmp_acc = [tmp_acc tmp_acc]; % full EM
    tmp_eti = eti_ea; tmp_eti = [tmp_eti tmp_eti]; % continuous
    tmp_sbj = 1:length(group_idx); tmp_sbj = [tmp_sbj tmp_sbj];
    tmp_names = {'beta', 'subregion', 'age', 'sex', 'acc', 'eti', 'sbj'};
    tmp_cell = {tmp_beta', tmp_type', tmp_age', tmp_sex', tmp_acc', tmp_eti', tmp_sbj'};
    tmp_array = cell2mat(tmp_cell);
    tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
    tmp_table.subregion = categorical(tmp_table.subregion);
    
    %
    lme = fitlme(tmp_table, 'beta~subregion*age*eti+(1|sbj)')
    tbl = anova(lme)

    % - print results
    fprintf('<< LME results >>\n');
    for i = 2:length(lme.CoefficientNames)
        tmp = ' ';
        if lme.Coefficients{i,6} < 0.05
            tmp = '* ';
        end
        fprintf('%2s%18s: β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', tmp, lme.Coefficients{i,1}, ...
            lme.Coefficients{i,2}, lme.Coefficients{i,3}, double(lme.Coefficients{i,5}), lme.Coefficients{i,4}, lme.Coefficients{i,6})
    end
end


%% relative activation of subregions

hp_index = [ 1 18 35 ];
sub_index = [ 4 7 ; 21 24 ; 38 41 ];
roi_idx = [ 1 2 ];

roi_idx_L = [ 1 4 7 ]; % Hippocampus, ant, post
roi_idx_R = [ 18 21 24 ];
roi_idx_B = [ 35 38 41 ];
roi_idx_list = {roi_idx_L, roi_idx_R, roi_idx_B};

con_i = 2; % 2:control
flag_i = 1; % 1:all
hemis_i = 3; % 3:bilateral

group_idx = (eti_ea > 0);
plot_color = [0.3 0.3 0.3; 223/255 153/255 153/255];

relative_ylim = [-0.1 4.8];

% rank ratio
sub_rank = cell(1, 3); hp_rank = cell(1, 3); sub_relative = cell(1, 3);
for hemis_i = 1:3
    hp_rank{hemis_i} = tiedrank(metric_list_enc{con_i}{flag_i}{hp_index(hemis_i)});

    for roi_i = 1:length(sub_index(hemis_i, :))
        sub_rank{hemis_i}{roi_i} = tiedrank(metric_list_enc{con_i}{flag_i}{sub_index(hemis_i, roi_i)});
    end

    sub_relative{hemis_i} = cellfun(@(x) x ./ hp_rank{hemis_i}, sub_rank{hemis_i}, 'uni', 0);
end


% ---------- group comparison ---------- %
hemis_i = 3;
for roi_i = 1:2
    data_all = sub_relative{hemis_i}{roi_i};
    [~, p, ~, stats] = ttest2(data_all(group_idx == 0), data_all(group_idx == 1));
    fprintf('relative ROI %d: t(%d)=%.3f, p=%.4f\n', roi_i, stats.df, stats.tstat, p);
end

% ---------- figure ---------- %
fig_draw = true;
if fig_draw
    for g_i = 1:2
        figure;
        set(gcf, 'Position', [100 100 500 100]);

        tmp_position = [1 4.1 7.2 10.3];
        for r_i = 1:length(roi_idx)
            subplot(1,12,[tmp_position(r_i) tmp_position(r_i)+1.2]);

            data_all = sub_relative{hemis_i}{roi_idx(r_i)};

            data = data_all(group_idx == g_i - 1);
            data_bhv = sbj_age(group_idx == g_i - 1);

            x_reg = reshape(data_bhv,[],1); y_reg = reshape(data,[],1);
            nan_reg = isnan(x_reg) | isnan(y_reg); x_reg(nan_reg) = []; y_reg(nan_reg) = [];
            [r, p] = corr(x_reg, y_reg);
            fig_scatter = scatter(x_reg, y_reg, 6, 'filled'); hold on
            fig_line = lsline;
            mdl = fitlm(x_reg, y_reg);
            xx = linspace(min(x_reg), max(x_reg), 1000); [~, ci] = predict(mdl, xx');
            fill([xx, flip(xx)]', [ci(:,1)', flip(ci(:,2))']', [.5 .5 .5], 'FaceAlpha', 0.2, 'linestyle', 'none');
            fig_scatter.MarkerFaceColor = plot_color(g_i, :);
            fig_scatter.SizeData = 6;
            fig_line.LineWidth = 1.2;
            if p < 0.05; fig_line.Color = [1 0.5 0.5]; end

            xlim([9 19]); xticks(10:2:18); xlabel('age');

            ylabel('relative activation');
            ylim(relative_ylim); yticks(0:1:4);

            set(gca,'XTickLabelRotation',0);
            set(gca,'LineWidth',0.6);
            set(gca,'FontName','Arial','FontSize',6, 'FontWeight','normal')
        end
    end
end


% ---------- Linear mixed-effect model ---------- %
lme_analysis = true;
if lme_analysis            
    tmp_beta = arrayfun(@(x) sub_relative{hemis_i}{x}, 1:2, 'uni', 0); % enc-ctrl, all
    tmp_beta = [tmp_beta{1} tmp_beta{2}];
    tmp_type = [ones(1, length(group_idx)) ones(1, length(group_idx))*2];
    tmp_age = sbj_age; tmp_age = [tmp_age tmp_age];
    tmp_sex = double(sbj_sex); tmp_sex = [tmp_sex tmp_sex];
    tmp_acc = em_acc{6}; tmp_acc = [tmp_acc tmp_acc]; % full EM
    tmp_eti = eti_ea; tmp_eti = [tmp_eti tmp_eti]; % continuous
    tmp_sbj = 1:length(group_idx); tmp_sbj = [tmp_sbj tmp_sbj];
    tmp_names = {'beta', 'subregion', 'age', 'sex', 'acc', 'eti', 'sbj'};
    tmp_cell = {tmp_beta', tmp_type', tmp_age', tmp_sex', tmp_acc', tmp_eti', tmp_sbj'};
    tmp_array = cell2mat(tmp_cell);
    tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
    tmp_table.subregion = categorical(tmp_table.subregion);
    
    %
    lme = fitlme(tmp_table, 'beta~subregion*age*eti+(1|sbj)')
    tbl = anova(lme)

    % - print results
    fprintf('<< LME results >>\n');
    for i = 2:length(lme.CoefficientNames)
        tmp = ' ';
        if lme.Coefficients{i,6} < 0.05
            tmp = '* ';
        end
        fprintf('%2s%18s: β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', tmp, lme.Coefficients{i,1}, ...
            lme.Coefficients{i,2}, lme.Coefficients{i,3}, double(lme.Coefficients{i,5}), lme.Coefficients{i,4}, lme.Coefficients{i,6})
    end
end


%% peak location: hippocampal activation peak along longitudinal axis (fig. 3e)

% ROI indices: big / hbt
hp_idx = [1 4];
hbt_idx = [1 2 3; 6 7 8]; % head, body, tail (L/R)

hemis_label = {'L', 'R', 'B'};
num_hemis = 2;

peak_bw_mm = 8; % smoothing bandwidth (mm)

trial_count_enc = nan(1, num_sbj);
peak_pos = nan(num_sbj, num_hemis);
peak_pos_comb = nan(num_sbj, 1);
hpc_len_mm = nan(num_sbj, num_hemis);
hpc_nvox = nan(num_sbj, num_hemis);
peak_beta = nan(num_sbj, num_hemis);

for sbj_i = 1:num_sbj
    fprintf('[%d/%d] %s\n', sbj_i, num_sbj, sbj_list{sbj_i});

    % load seg
    seg = load(fullfile(dir_seg, sprintf('%s.mat', sbj_list{sbj_i})));
    hdr = seg.sbj.seg_fit.header;
    vol_size = hdr.ImageSize;
    T = hdr.Transform.T;

    % load beta
    bdata = load(fullfile(dir_beta, sprintf('%s.mat', sbj_list{sbj_i})));
    n_trials = length(bdata.glm.beta);
    trial_count_enc(sbj_i) = sum(arrayfun(@(x) ~isempty(x.enc), bdata.glm.beta));

    beta_all = zeros(vol_size);
    trial_count = 0;
    for trial_i = 1:n_trials
        tmp_beta = bdata.glm.beta(trial_i).enc_ctrl;
        if ~isempty(tmp_beta)
            beta_all = beta_all + double(tmp_beta);
            trial_count = trial_count + 1;
        end
    end
    if trial_count == 0
        beta_all = [];
    else
        beta_all = beta_all / trial_count;
    end
    if isempty(beta_all); continue; end

    % per hemisphere
    proj_norm_tmp = {[], []};
    beta_tmp = {[], []};
    bw_norm_tmp = [nan, nan];
    for hemis_i = 1:num_hemis
        idx_hpc = seg.sbj.seg_fit.hpc.big.roi_list_nn{hp_idx(hemis_i)};
        idx_head = seg.sbj.seg_fit.hpc.hbt.roi_list_nn{hbt_idx(hemis_i, 1)};
        idx_body = seg.sbj.seg_fit.hpc.hbt.roi_list_nn{hbt_idx(hemis_i, 2)};
        idx_tail = seg.sbj.seg_fit.hpc.hbt.roi_list_nn{hbt_idx(hemis_i, 3)};

        if isempty(idx_hpc) || isempty(idx_head) || isempty(idx_body) || isempty(idx_tail); continue; end

        % voxel coordinates (mm)
        mm_hpc = vox2mm(idx_hpc, vol_size, T);
        mm_head = vox2mm(idx_head, vol_size, T);
        mm_tail = vox2mm(idx_tail, vol_size, T);

        % long axis (PCA)
        [coeff, ~, ~] = pca(mm_hpc);
        ax = coeff(:, 1);
        if mean(mm_head * ax) > mean(mm_tail * ax)
            ax = -ax;
        end

        all_proj = mm_hpc * ax;
        head_proj = min(all_proj);
        tail_proj = max(all_proj);
        ax_range = tail_proj - head_proj;
        hpc_len_mm(sbj_i, hemis_i) = ax_range;
        hpc_nvox(sbj_i, hemis_i) = length(idx_hpc);
        proj_norm_tmp{hemis_i} = (mm_hpc * ax - head_proj) / ax_range;
        beta_tmp{hemis_i} = double(beta_all(idx_hpc));
        bw_norm_tmp(hemis_i) = peak_bw_mm / ax_range;
        beta_vals = double(beta_all(idx_hpc));
        proj_vals = mm_hpc * ax;
        x_eval = linspace(head_proj, head_proj + ax_range, 200);
        kernel_weight = exp(-0.5 * ((x_eval' - proj_vals') / peak_bw_mm).^2);
        profile = (kernel_weight * beta_vals) ./ sum(kernel_weight, 2);
        [peak_beta(sbj_i, hemis_i), peak_idx] = max(profile);
        peak_pos(sbj_i, hemis_i) = (x_eval(peak_idx) - head_proj) / ax_range;
    end

    % bilateral peak
    if ~isempty(proj_norm_tmp{1}) && ~isempty(proj_norm_tmp{2})
        proj_bil = [proj_norm_tmp{1}(:); proj_norm_tmp{2}(:)];
        beta_bil = [beta_tmp{1}(:); beta_tmp{2}(:)];
        bw_bil = mean(bw_norm_tmp, 'omitnan');
        x_eval = linspace(0, 1, 200);
        kernel_weight = exp(-0.5 * ((x_eval' - proj_bil(:)') / bw_bil).^2);
        profile = (kernel_weight * beta_bil(:)) ./ sum(kernel_weight, 2);
        [~, peak_idx] = max(profile);
        peak_pos_comb(sbj_i) = x_eval(peak_idx);
    end
end

% bilateral
peak_pos(:, 3) = peak_pos_comb;
hpc_len_mm = [hpc_len_mm, mean(hpc_len_mm, 2, 'omitnan')];
hpc_nvox = [hpc_nvox, mean(hpc_nvox, 2, 'omitnan')];
peak_beta = [peak_beta, mean(peak_beta, 2, 'omitnan')];


%% peak location: figures and stats


% ---------- parameters ---------- %
ap_data = peak_pos;  ap_label = 'peak';
h_plot = 3;  % 1=L, 2=R, 3=bilateral

group_idx = (eti_ea > 0);
group_flag = (eti_ea == 0)';
group_names = {'non-abused', 'abused'};
group_colors = {[0.35 0.35 0.35], [0.78 0.25 0.25]};
color1 = group_colors{1};
color2 = group_colors{2};

num_perm = 5000;

% ---------- density ---------- %
ap = ap_data(:, h_plot);
valid = ~isnan(ap);
x_eval = linspace(0, 1, 300);

g1_valid = valid &  group_flag(:);
g2_valid = valid & ~group_flag(:);
ap_g1 = ap(g1_valid);
ap_g2 = ap(g2_valid);

[~, ~, kde_bw] = ksdensity(ap(valid), x_eval, 'BoundaryCorrection', 'reflection');
f_g1 = ksdensity(ap_g1, x_eval, 'BoundaryCorrection', 'reflection', 'Bandwidth', kde_bw);
f_g2 = ksdensity(ap_g2, x_eval, 'BoundaryCorrection', 'reflection', 'Bandwidth', kde_bw);
obs_diff = f_g1 - f_g2;

% permutation
ap_pooled = ap(g1_valid | g2_valid);
n1 = sum(g1_valid);
perm_diff = zeros(num_perm, length(x_eval));
rng(42);
for perm_i = 1:num_perm
    idx_shuf = randperm(length(ap_pooled));
    fp1 = ksdensity(ap_pooled(idx_shuf(1:n1)),     x_eval, 'BoundaryCorrection', 'reflection', 'Bandwidth', kde_bw);
    fp2 = ksdensity(ap_pooled(idx_shuf(n1+1:end)), x_eval, 'BoundaryCorrection', 'reflection', 'Bandwidth', kde_bw);
    perm_diff(perm_i, :) = fp1 - fp2;
end

% pointwise threshold (two-sided p < .05)
thr_pos_pw = prctile(perm_diff, 97.5, 1);
thr_neg_pw = prctile(perm_diff,  2.5, 1);

% cluster mass
max_mass_null = zeros(num_perm, 1);
for perm_i = 1:num_perm
    d = perm_diff(perm_i, :);
    for sign_val = [1 -1]
        if sign_val == 1; mask = d > thr_pos_pw;
        else;             mask = d < thr_neg_pw; end
        s = find(diff([false, mask]) == 1);
        e = find(diff([mask, false]) == -1);
        for k = 1:length(s)
            max_mass_null(perm_i) = max(max_mass_null(perm_i), sum(abs(d(s(k):e(k)))));
        end
    end
end
mass_thr = prctile(max_mass_null, 95);

% significant clusters
sig_pos = false(size(obs_diff));
sig_neg = false(size(obs_diff));
s = find(diff([false, obs_diff > thr_pos_pw]) == 1);
e = find(diff([obs_diff > thr_pos_pw, false]) == -1);
for k = 1:length(s)
    if sum(abs(obs_diff(s(k):e(k)))) >= mass_thr
        sig_pos(s(k):e(k)) = true;
    end
end
s = find(diff([false, obs_diff < thr_neg_pw]) == 1);
e = find(diff([obs_diff < thr_neg_pw, false]) == -1);
for k = 1:length(s)
    if sum(abs(obs_diff(s(k):e(k)))) >= mass_thr
        sig_neg(s(k):e(k)) = true;
    end
end

sig_pos_sel = sig_pos;
sig_neg_sel = sig_neg;

% ---------- statistics ---------- %
[~, p_tt, ~, stats_tt] = ttest2(ap_g1, ap_g2);
fprintf('\n========== AP density group statistics [%s, %s HPC] ==========\n', ap_label, hemis_label{h_plot});
fprintf('%-6s  mean=%.3f  std=%.3f  N=%d\n', group_names{1}, mean(ap_g1), std(ap_g1), length(ap_g1));
fprintf('%-6s  mean=%.3f  std=%.3f  N=%d\n', group_names{2}, mean(ap_g2), std(ap_g2), length(ap_g2));
fprintf('t-test: t(%.0f)=%.3f  p=%.4f\n', stats_tt.df, stats_tt.tstat, p_tt);

report_range = @(s, x) sprintf('%.3f-%.3f', min(x(s)), max(x(s)));
fprintf('[correction method: cluster]\n');
for g_i = 1:2
    if g_i == 1; sig_sel = sig_pos_sel; group1 = group_names{1}; group2 = group_names{2};
    else;        sig_sel = sig_neg_sel; group1 = group_names{2}; group2 = group_names{1}; end
    if any(sig_sel)
        fprintf('%s > %s  (cluster): %s\n', group1, group2, report_range(sig_sel, x_eval));
        starts = find(diff([false, sig_sel]) == 1);
        ends = find(diff([sig_sel, false]) == -1);
        for cluster_i = 1:length(starts)
            mass = sum(abs(obs_diff(starts(cluster_i):ends(cluster_i))));
            p_cluster = (sum(max_mass_null >= mass) + 1)/(num_perm + 1);
            fprintf('cluster %d: corrected p=%.4f\n', cluster_i, p_cluster);
        end
    end
end
if ~any(sig_pos_sel) && ~any(sig_neg_sel)
    fprintf('No significant region (cluster)\n');
end


% ---------- figure ---------- %
figure;
set(gcf, 'Position', [100 100 200 200]);

% - density -
ax1 = subplot('Position', [0.16 0.30 0.80 0.62]);
hold(ax1, 'on');

y_max = max([f_g1, f_g2]) * 1.05;

for tmp_k = 1:2
    if tmp_k == 1; sig_mask = sig_pos_sel(:)';
    else;           sig_mask = sig_neg_sel(:)'; end
    sig_color = [1 0.75 0];

    if ~any(sig_mask); continue; end
    s_idx = find(diff([false, sig_mask]) == 1);
    e_idx = find(diff([sig_mask, false]) == -1);
    for k = 1:length(s_idx)
        patch([x_eval(s_idx(k)) x_eval(e_idx(k)) x_eval(e_idx(k)) x_eval(s_idx(k))], ...
            [0 0 y_max y_max], sig_color, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', 'Parent', ax1);
    end
end

fill([x_eval, fliplr(x_eval)], [f_g1, zeros(1,length(x_eval))], color1, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
fill([x_eval, fliplr(x_eval)], [f_g2, zeros(1,length(x_eval))], color2, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
plot(x_eval, f_g1, 'Color', color1, 'LineWidth', 1.2);
plot(x_eval, f_g2, 'Color', color2, 'LineWidth', 1.2);

ylabel('Density');
xlim([0 1]); ylim([0 y_max]);
set(ax1, 'XTickLabel', {}, 'LineWidth', 0.6, 'FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'normal');
box off;

% - strip -
ax2 = subplot('Position', [0.16 0.08 0.80 0.17]);
hold(ax2, 'on');

dot_jitter = 0.06;
group_y1 = 1.5;
group_y2 = 0.5;

rng(1);
jitter1 = dot_jitter * randn(length(ap_g1), 1);
jitter2 = dot_jitter * randn(length(ap_g2), 1);
scatter(ap_g1, group_y1 + jitter1, 6, color1, 'filled', 'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor', 'none');
scatter(ap_g2, group_y2 + jitter2, 6, color2, 'filled', 'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor', 'none');

for g_i = 1:2
    if g_i == 1; grp_data = ap_g1; color = color1; ym = group_y1;
    else;        grp_data = ap_g2; color = color2; ym = group_y2; end
    data_mean = mean(grp_data);
    data_sem = std(grp_data) / sqrt(length(grp_data));
    plot([data_mean-data_sem, data_mean+data_sem], [ym, ym], '-', 'Color', color, 'LineWidth', 1.2);
    plot([data_mean-data_sem, data_mean-data_sem], [ym-0.05, ym+0.05], '-', 'Color', color, 'LineWidth', 1.2);
    plot([data_mean+data_sem, data_mean+data_sem], [ym-0.05, ym+0.05], '-', 'Color', color, 'LineWidth', 1.2);
    scatter(data_mean, ym, 30, color, 'd', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1);
end

xlim([0 1]);
ylim([0.0 2.1]);
yticks([group_y2, group_y1]);
yticklabels({'', ''});
set(ax2, 'LineWidth', 0.6, 'FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'normal', 'YDir', 'normal');
box off;

linkaxes([ax1, ax2], 'x');

% ---------- statistics ---------- %
ok = ~isnan(ap);
idx0 = (group_idx == 0) & ok';
idx1 = (group_idx == 1) & ok';

fprintf('\n<< AP position group diff >>\n');
[~, p, ~, stats] = ttest2(ap(idx0), ap(idx1));
fprintf(' AP peak: non-abused = %.3f ± %.3f, abused = %.3f ± %.3f\n', mean(ap(idx0)), std(ap(idx0)), mean(ap(idx1)), std(ap(idx1)));
fprintf('          - T(%d)=%.3f, p=%.3f\n', stats.df, stats.tstat, p);


%% peak location: age correlation and interaction
% Requires peak_pos from the peak calculation section.

h_plot = 3;
num_boot = 1000;
plot_color = [0.3 0.3 0.3; 223/255 153/255 153/255];
peak_ylim = [-0.05 1.05];

ap = peak_pos(:, h_plot);
group_idx = (eti_ea > 0);
xvar = sbj_age(:);
figure;
set(gcf, 'Position', [100 100 500 130]);
rng(42);
for g_i = 0:1
    ax = subplot(1, 2, g_i+1);
    ok = (group_idx(:) == g_i) & ~isnan(ap) & ~isnan(xvar);
    x_reg = xvar(ok); y_reg = ap(ok);
    [tau, p] = corr(x_reg, y_reg, 'Type', 'Kendall', 'Tail', 'both');
    xx = linspace(min(x_reg), max(x_reg), 100)';
    boot_fit = nan(num_boot, length(xx));
    for boot_i = 1:num_boot
        idx = randi(length(x_reg), length(x_reg), 1);
        coef = polyfit(x_reg(idx), y_reg(idx), 1);
        boot_fit(boot_i, :) = polyval(coef, xx)';
    end
    ci = prctile(boot_fit, [2.5 97.5], 1);
    fill(ax, [xx; flipud(xx)], [ci(1,:)'; flipud(ci(2,:)')], plot_color(g_i+1,:), ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold(ax, 'on');
    scatter(ax, x_reg, y_reg, 6, plot_color(g_i+1,:), 'filled');
    coef = polyfit(x_reg, y_reg, 1);
    plot(ax, xx, polyval(coef, xx), 'Color', plot_color(g_i+1,:), 'LineWidth', 1.2);
    xlim(ax, [9 19]); ylim(ax, peak_ylim);
    xlabel(ax, 'age'); ylabel(ax, 'AP peak position');
    set(ax, 'FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'normal', 'LineWidth', 0.6);
    box(ax, 'off');
    fprintf('group %d: tau=%.3f, p=%.4f\n', g_i, tau, p);
end

% ---------- linear model ---------- %
tmp_names = {'ap', 'age', 'abuse', 'sex'};
tmp_cell = {ap, sbj_age(:), double(group_idx(:)), double(sbj_sex(:))};
tmp_array = cell2mat(tmp_cell);
tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
lm = fitlm(tmp_table, 'ap~age*abuse')
lm_sex = fitlm(tmp_table, 'ap~age*abuse+sex')


%% hippocampal activation: demographic and psychological factors

roi_idx = 35; % bilateral whole HP
group_idx = (eti_ea > 0);
cesd = sbj_svy{ismember(svy_label, 'CES_D')}(:);
stais = sbj_svy{ismember(svy_label, 'STAI_S')}(:);
data_all = metric_list_enc{2}{1}{roi_idx}(:);
factor_list = {cesd >= median(cesd), stais >= median(stais), ...
    sbj_age(:) >= mean(sbj_age), double(sbj_sex(:))};
for factor_i = 1:length(factor_list)
    [p, tbl] = anovan(data_all, {group_idx(:), factor_list{factor_i}}, ...
        'model', 'interaction', 'display', 'off');
    fprintf('factor %d: abuse p=%.4f, factor p=%.4f, interaction p=%.4f\n', factor_i, p);
    disp(tbl);
end


%% hippocampal activation: memory correlation

roi_idx = [35 38 41]; % bilateral HP, anterior, posterior
acc_idx = [6 5 4]; % full EM, where-when, what-when
group_idx = (eti_ea > 0);
for acc_i = acc_idx
    for roi_i = roi_idx
        data_all = metric_list_enc{2}{1}{roi_i};
        for g_i = 0:1
            [r, p] = corr(em_acc{acc_i}(group_idx == g_i)', data_all(group_idx == g_i)');
            fprintf('acc %d, %s, group %d: r=%.3f, p=%.4f\n', ...
                acc_i, roi_name_all{roi_i}, g_i, r, p);
        end
    end
end


%% retrieval activation and encoding-retrieval similarity

roi_idx = [35 38 41]; % bilateral HP, anterior, posterior
load(fullfile(dir_metric, sprintf('metric_neural_N%d.mat', num_sbj)), 'metric_list_ret', 'metric_list_act');
load(fullfile(dir_metric, sprintf('pattern_sim_ERS_N%d.mat', num_sbj)), 'metric_ERS');
load(fullfile(dir_metric, sprintf('pattern_sim_ERS_enact_N%d.mat', num_sbj)), 'metric_ERS_enact');

data_list = {metric_list_ret{2}{1}, metric_list_act{3}{1}, metric_ERS{7}{1}, metric_ERS_enact{7}{1}};
data_names = {'imagery activation', 'reenactment activation', 'imagery ERS', 'reenactment ERS'};
group_idx = (eti_ea > 0);
bf_rscale = sqrt(2)/2; % JZS; ERS metric 7 = same-trial minus run-matched other-trial
for data_i = 1:length(data_list)
    for roi_i = roi_idx
        data_all = data_list{data_i}{roi_i};
        non_data = data_all(group_idx == 0);
        abused_data = data_all(group_idx == 1);
        [~, p, ~, stats] = ttest2(non_data, abused_data);
        n_eff = numel(non_data)*numel(abused_data)/(numel(non_data)+numel(abused_data));
        bf10 = jzs_bf10(stats.tstat, stats.df, n_eff, bf_rscale);
        fprintf('%s, %s: t(%d)=%.3f, p=%.4f, BF10=%.3f, BF01=%.3f\n', ...
            data_names{data_i}, roi_name_all{roi_i}, stats.df, stats.tstat, p, bf10, 1/bf10);
    end
end


%% hippocampal volume: ICV correction, group difference and age

volume_index = [1 2 3]; % whole HP, anterior, posterior
load(fullfile(dir_metric, sprintf('volume_N%d.mat', num_sbj)), 'volume_b', 'volume_icv', 'volume_label');
group_idx = (eti_ea > 0);
volume_residual = nan(num_sbj, length(volume_index));
for roi_i = 1:length(volume_index)
    mdl = fitlm(volume_icv(:), volume_b(volume_index(roi_i), :)');
    volume_residual(:, roi_i) = mdl.Residuals.Raw;
    data_all = volume_residual(:, roi_i);
    [~, p, ~, stats] = ttest2(data_all(group_idx == 0), data_all(group_idx == 1));
    fprintf('%s: t(%d)=%.3f, p=%.4f\n', volume_label{volume_index(roi_i)}, stats.df, stats.tstat, p);
    for g_i = 0:2
        idx = (group_idx(:) == g_i);
        if g_i == 2; idx = true(num_sbj, 1); end
        [r, p] = corr(sbj_age(idx)', data_all(idx));
        fprintf('group %d: r=%.3f, p=%.4f\n', g_i, r, p);
    end
end
% ---------- Linear mixed-effects model ---------- %
tmp_names = {'vol', 'subregion', 'age', 'eti', 'sbj'};
tmp_cell = {reshape(volume_residual(:, 2:3), [], 1), repelem((1:2)', num_sbj), ...
    repmat(sbj_age(:), 2, 1), repmat(eti_ea(:), 2, 1), repmat((1:num_sbj)', 2, 1)};
tmp_array = cell2mat(tmp_cell);
tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
tmp_table.subregion = categorical(tmp_table.subregion);
lme = fitlme(tmp_table, 'vol~subregion*age*eti+(1|sbj)')
anova(lme)


%% hippocampal activation: sex and valid-trial sensitivity
% Requires sub_relative and trial_count_enc from the preceding calculations.

group_idx = (eti_ea > 0);
for data_i = 1:3
    if data_i == 1
        tmp_beta = [metric_list_enc{2}{1}{1}(:); metric_list_enc{2}{1}{18}(:)];
        model_formula = 'beta~region*eti';
    elseif data_i == 2
        tmp_beta = [metric_list_enc{2}{1}{38}(:); metric_list_enc{2}{1}{41}(:)];
        model_formula = 'beta~region*age*eti';
    else
        tmp_beta = [sub_relative{3}{1}(:); sub_relative{3}{2}(:)];
        model_formula = 'beta~region*age*eti';
    end
    tmp_names = {'beta', 'region', 'age', 'eti', 'sex', 'tnum', 'sbj'};
    tmp_cell = {tmp_beta, repelem((1:2)', num_sbj), repmat(sbj_age(:), 2, 1), ...
        repmat(eti_ea(:), 2, 1), repmat(double(sbj_sex(:)), 2, 1), ...
        repmat(trial_count_enc(:), 2, 1), repmat((1:num_sbj)', 2, 1)};
    tmp_array = cell2mat(tmp_cell);
    tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
    tmp_table.region = categorical(tmp_table.region);
    lme_sex = fitlme(tmp_table, [model_formula '+sex+(1|sbj)'])
    anova(lme_sex)
    lme_trials = fitlme(tmp_table, [model_formula '+tnum+(1|sbj)'])
    anova(lme_trials)
end

tmp_beta = [metric_list_enc{2}{1}{38}(:); metric_list_enc{2}{1}{41}(:)];
tmp_names = {'beta', 'region', 'age', 'abuse', 'sex', 'sbj'};
tmp_cell = {tmp_beta, repelem((1:2)', num_sbj), repmat(sbj_age(:), 2, 1), ...
    repmat(double(group_idx(:)), 2, 1), repmat(double(sbj_sex(:)), 2, 1), repmat((1:num_sbj)', 2, 1)};
tmp_array = cell2mat(tmp_cell);
tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
tmp_table.region = categorical(tmp_table.region);
for g_i = 0:1
    lme_sex = fitlme(tmp_table(tmp_table.abuse == g_i, :), 'beta~region*age+sex+(1|sbj)')
    anova(lme_sex)
end

%% local functions

function mm = vox2mm(idx, vol_size, T)
    [vx, vy, vz] = ind2sub(vol_size, idx(:));
    vox_h = [vx, vy, vz, ones(numel(idx), 1)];
    mm = (vox_h * T);
    mm = mm(:, 1:3);
end

function bf10 = jzs_bf10(tstat, df, n_eff, rscale)
    a = n_eff * rscale^2;
    
    integrand = @(z) sqrt(z.^2 ./ (z.^2 + a)) .* ...
        (1 + (tstat^2 / df) .* z.^2 ./ (z.^2 + a)).^(-(df + 1) / 2) .* ...
        sqrt(2 / pi) .* exp(-z.^2 / 2);
    
    m1 = integral(integrand, 0, Inf);
    m0 = (1 + tstat^2 / df)^(-(df + 1) / 2);
    
    bf10 = m1 / m0;
end
