% Episodic memory task for adolescents
% :: ROI activation - - analysis (figures and statistics)
% :: code written by Kahyun Choi


%% define parameters

% load subject list
subject_path = 'C:\subject_info_directory';
bhv_nums = load(fullfile(subject_path, 'subject_list_for_GLM.mat')); % 'bhv_list', 'sbj_list', 'sbj_nums', 'num_sbj'

% load ROI activation files
dir_seg = 'C:\ROI_data_directory';
load(fullfile(dir_seg, sprintf('metric_neural_N%d.mat', num_sbj))); % metric_list_enc, flag_name_mri, roi_name_all
disp('Metric data loaded');

% behavior
load(fullfile(subject_path, 'behavioral_data.mat')); % 'em_acc', 'sbj_age', 'sbj_sex'

% ETI
dir_survey = fullfile(subject_path, 'ETI_survey.xlsx');
sbj_eti = readtable(dir_survey);
sbj_eti = sbj_eti(sbj_nums, :);
sbj_eti = table2struct(sbj_eti); % ALL, GT, PA, EA, Q1~21
eti_ea = [sbj_eti.EA];


%% Hippocampal activation: group difference

%%%%%%%%%%%%%%%%
roi_idx_L = [ 1 8 9 10 ]; % Hippocampus, CA3, CA1, subiculum
roi_idx_R = [ 17 24 25 26 ];
roi_idx_B = [ 33 40 41 42 ];
roi_idx_list = {roi_idx_L, roi_idx_R, roi_idx_B};

con_i = 2; % 2:control
flag_i = 1; % 1:all

group_idx = (eti_ea > 0);
plot_color = [0.3 0.3 0.3; 223/255 153/255 153/255];
%%%%%%%%%%%%%%%%

%
fig_draw = true;
if fig_draw
    
    figure;
    set(gcf, 'Position', [100 100 250 320]);

    % --- Hippocampus ---
    subplot(2,3,[1 2]);
    data_all = metric_list_enc{con_i}{flag_i}{roi_idx_B(1)};
    data = {data_all(group_idx == 0), data_all(group_idx == 1)};

    avg = cellfun(@nanmean,data);
    err = cellfun(@(x) nanstd(x)/sqrt(length(x)), data);

    %
    hold on
    fig_bar = bar(avg);
    set(fig_bar, 'FaceColor', 'flat');
    fig_bar.CData = plot_color;
    fig_bar.EdgeColor = 'none';

    %
    fig_err = errorbar(avg, err);
    fig_err.Color = [0 0 0];
    fig_err.LineStyle = 'none';
    fig_err.LineWidth = 1.2;
    fig_err.CapSize = 0;

    %
    y_text = abs(max(avg)) + 0.15;
    for data_i = 1:length(data)
        [~, p, ~, stats] = ttest(data{data_i});
        mark = '';
        if p < .1; mark = '†'; end
        if p < .05; mark = '*'; end
        if p < .01; mark = '**'; end
        if p < .001; mark = '***'; end
        if p < .05 
            text(data_i, y_text, mark, 'FontSize',7, 'HorizontalAlignment','center','VerticalAlignment','middle')
        end
    end

    [~, p1] = ttest(data{1});
    [~, p2] = ttest(data{2});
    [~, p12, ~, stats] = ttest2(data{1}, data{2});

    xticks(1.5); xticklabels('HP');
    ylim([-0.25 0.55]);
    ylabel('beta coefficient', 'FontSize', 7);
    set(gca,'LineWidth',0.8);
    set(gca,'FontName','Arial','FontSize',7, 'FontWeight','bold')

    % --- subfields ---
    subplot(2,3,[4 6]);

    data_all = arrayfun(@(x) metric_list_enc{con_i}{flag_i}{x}, roi_idx_B(2:4), 'uni', 0);
    data = {data_all{1}(group_idx == 0), data_all{1}(group_idx == 1), ...
            data_all{2}(group_idx == 0), data_all{2}(group_idx == 1), ...
            data_all{3}(group_idx == 0), data_all{3}(group_idx == 1)};

    avg = cellfun(@nanmean,data);
    err = cellfun(@(x) nanstd(x)/sqrt(length(x)), data);

    %
    hold on
    b1 = bar(1, avg(1), 'FaceColor', plot_color(1, :), 'EdgeColor', 'none');
    b2 = bar(2, avg(2), 'FaceColor', plot_color(2, :), 'EdgeColor', 'none');
    b3 = bar(3.5, avg(3), 'FaceColor', plot_color(1, :), 'EdgeColor', 'none');
    b4 = bar(4.5, avg(4), 'FaceColor', plot_color(2, :), 'EdgeColor', 'none');
    b5 = bar(6, avg(5), 'FaceColor', plot_color(1, :), 'EdgeColor', 'none');
    b6 = bar(7, avg(6), 'FaceColor', plot_color(2, :), 'EdgeColor', 'none');

    %
    err1 = errorbar(1, avg(1), err(1));
    err1.Color = [0 0 0]; err1.LineStyle = 'none'; err1.LineWidth = 1; err1.CapSize = 0;
    err2 = errorbar(2, avg(2), err(2));
    err2.Color = [0 0 0]; err2.LineStyle = 'none'; err2.LineWidth = 1; err2.CapSize = 0;
    err3 = errorbar(3.5, avg(3), err(3));
    err3.Color = [0 0 0]; err3.LineStyle = 'none'; err3.LineWidth = 1; err3.CapSize = 0;
    err4 = errorbar(4.5, avg(4), err(4));
    err4.Color = [0 0 0]; err4.LineStyle = 'none'; err4.LineWidth = 1; err4.CapSize = 0;
    err5 = errorbar(6, avg(5), err(5));
    err5.Color = [0 0 0]; err5.LineStyle = 'none'; err5.LineWidth = 1; err5.CapSize = 0;
    err6 = errorbar(7, avg(6), err(6));
    err6.Color = [0 0 0]; err6.LineStyle = 'none'; err6.LineWidth = 1; err6.CapSize = 0;

    %
    y_text = 0.5; 
    x_coords = [1 2 3.5 4.5 6 7];
    for data_i = 1:length(data)
        [~, p, ~, stats] = ttest(data{data_i});
        mark = '';
        if p < .1; mark = '†'; end
        if p < .05; mark = '*'; end
        if p < .01; mark = '**'; end
        if p < .001; mark = '***'; end
        if p < .1 
            text(x_coords(data_i), y_text, mark, 'FontSize',7, 'HorizontalAlignment','center','VerticalAlignment','middle')
        end
    end
    label_list = {'CA23DG', 'CA1', 'Subiculum'};
    xticks([1.5 4 6.5]); xticklabels(label_list);    
    ylim([-0.25 0.55]);
    ylabel('beta coefficient');

    set(gca,'LineWidth',0.8);
    set(gca,'FontName','Arial','FontSize',7, 'FontWeight','bold')
end


% ---------- statistics ---------- %
stats_on = true;
if stats_on
    fprintf('<< hippocampal activation >>\n');
    fprintf(' - non-abused -\n');
    for i = 1:3
        [~,p,~,stats] = ttest(metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(1)}(group_idx == 0));
        fprintf('%10s: T(%d)=%.3f, p=%.3f\n', roi_name_all{roi_idx_list{i}(1)}, stats.df, stats.tstat, p);
    end
    disp(' ');
    
    fprintf(' - abused -\n');
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
    
    fprintf('<< hp subfield activation >>\n');
    fprintf(' - non-abused -\n');
    for i = 1:3
        for j = 2:4
            [~,p,~,stats] = ttest(metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(j)}(group_idx == 0));
            fprintf('%10s: T(%d)=%.3f, p=%.3f\n', roi_name_all{roi_idx_list{i}(j)}, stats.df, stats.tstat, p);
        end
    end
    disp(' ');
    
    fprintf(' - abused -\n');
    for i = 1:3
        for j = 2:4
            [~,p,~,stats] = ttest(metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(j)}(group_idx == 1));
            fprintf('%10s: T(%d)=%.3f, p=%.3f\n', roi_name_all{roi_idx_list{i}(j)}, stats.df, stats.tstat, p);
        end
    end
    disp(' ');
    
    fprintf(' - group diff -\n');
    for i = 1:3
        for j = 2:4
            [~,p,~,stats] = ttest2(metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(j)}(group_idx == 0), metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(j)}(group_idx == 1));
            fprintf('%10s: T(%d)=%.3f, p=%.3f\n', roi_name_all{roi_idx_list{i}(j)}, stats.df, stats.tstat, p);
        end
    end
    disp(' ');
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

    %
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
            lme.Coefficients(i,2), lme.Coefficients(i,3), double(lme.Coefficients(i,5)), lme.Coefficients(i,4), lme.Coefficients(i,6))
    end
end



%% Hippocampal activation (subfield): non-abused

%%%%%%%%%%%%%%%%
roi_idx_L = [ 1 8 9 10 ]; % Hippocampus, CA3, CA1, subiculum
roi_idx_R = [ 17 24 25 26 ];
roi_idx_B = [ 33 40 41 42 ];
roi_idx_list = {roi_idx_L, roi_idx_R, roi_idx_B};

con_i = 2; % 2:control
flag_i = 1; % 1:all
hemis_i = 3; % 3:bilateral

group_idx = (eti_ea > 0);
%%%%%%%%%%%%%%%%

% ---------- statistics ---------- %
stats_on = true;
if stats_on
    fprintf('<< hippocampal correlation w/ age >>\n');
    fprintf(' - non-abused -\n');
    for i = 1:3
        for j = 1:4
            tmp_data = metric_list_enc{con_i}{flag_i}{roi_idx_list{i}(j)};
            [r,p] = corr(sbj_age(tmp_group == 0)', tmp_data(tmp_group == 0)');
            fprintf('%10s: r=%.3f, p=%.3f\n', roi_name_all{roi_idx_list{i}(j)}, r, p);
        end
    end
    disp(' ');
end


% ---------- Linear mixed-effects model ---------- %
lme_analysis = true;
if lme_analysis            
    tmp_beta = arrayfun(@(x) metric_list_enc{2}{1}{x}, roi_idx_list{hemis_i}(2:4), 'uni', 0); % enc-ctrl, all
    tmp_beta = [tmp_beta{1}(group_idx == 0) tmp_beta{2}(group_idx == 0) tmp_beta{3}(group_idx == 0)];

    tmp_type = [ones(1, sum(group_idx == 0)) ones(1, sum(group_idx == 0))*2 ones(1, sum(group_idx == 0))*3];

    tmp_age = sbj_age(group_idx == 0); tmp_age = [tmp_age tmp_age tmp_age];    
    tmp_sex = double(sbj_sex(group_idx == 0)); tmp_sex = [tmp_sex tmp_sex tmp_sex];
    tmp_acc = em_acc{6}(group_idx == 0); tmp_acc = [tmp_acc tmp_acc tmp_acc]; % full EM
    tmp_eti = eti_ea(group_idx == 0); tmp_eti = [tmp_eti tmp_eti tmp_eti]; % continuous

    tmp_sbj = 1:sum(group_idx == 0); tmp_sbj = [tmp_sbj tmp_sbj tmp_sbj];

    %
    tmp_names = {'beta', 'subfield', 'age', 'sex', 'acc', 'eti', 'sbj'};
    tmp_cell = {tmp_beta', tmp_type', tmp_age', tmp_sex', tmp_acc', tmp_eti', tmp_sbj'};
    tmp_array = cell2mat(tmp_cell);
    tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
    tmp_table.subfield = categorical(tmp_table.subfield);

    %
    lme = fitlme(tmp_table, 'beta~subfield*age+(1|sbj)')
    anova(lme)

    % - print results
    fprintf('<< LME results >>\n');
    for i = 2:length(lme.CoefficientNames)
        tmp = ' ';
        if lme.Coefficients{i,6} < 0.05
            tmp = '* ';
        end
        fprintf('%2s%18s: β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', tmp, lme.Coefficients{i,1}, ...
            lme.Coefficients(i,2), lme.Coefficients(i,3), double(lme.Coefficients(i,5)), lme.Coefficients(i,4), lme.Coefficients(i,6))
    end
end


%% Hippocampal activation (subfield): abused

%%%%%%%%%%%%%%%%
roi_idx_L = [ 1 8 9 10 ]; % Hippocampus, CA3, CA1, subiculum
roi_idx_R = [ 17 24 25 26 ];
roi_idx_B = [ 33 40 41 42 ];
roi_idx_list = {roi_idx_L, roi_idx_R, roi_idx_B};

con_i = 2; % 2:control
flag_i = 1; % 1:all
hemis_i = 3; % 3:bilateral

group_idx = (eti_ea > 0);
%%%%%%%%%%%%%%%%

% ---------- stats ---------- %
stats_on = true;
if stats_on
    fprintf('<< hippocampal correlation w/ age >>\n');
    fprintf(' - abused -\n');
    for i = 1:3
        for j = 1:4
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
    tmp_beta = arrayfun(@(x) metric_list_enc{2}{1}{x}, roi_idx_list{hemis_i}(2:4), 'uni', 0); % enc-ctrl, all
    tmp_beta = [tmp_beta{1}(group_idx == 1) tmp_beta{2}(group_idx == 1) tmp_beta{3}(group_idx == 1)];

    tmp_type = [ones(1, sum(group_idx == 1)) ones(1, sum(group_idx == 1))*2 ones(1, sum(group_idx == 1))*3];

    tmp_age = sbj_age(group_idx == 1); tmp_age = [tmp_age tmp_age tmp_age];
    tmp_sex = double(sbj_sex(group_idx == 1)); tmp_sex = [tmp_sex tmp_sex tmp_sex];
    tmp_acc = em_acc{acc_i}(group_idx == 1); tmp_acc = [tmp_acc tmp_acc tmp_acc]; % full EM
    tmp_eti = eti_ea(group_idx == 1); tmp_eti = [tmp_eti tmp_eti tmp_eti]; % continuous

    tmp_sbj = 1:sum(group_idx == 1); tmp_sbj = [tmp_sbj tmp_sbj tmp_sbj];

    %
    tmp_names = {'beta', 'subfield', 'age', 'sex', 'acc', 'eti', 'sbj'};
    tmp_cell = {tmp_beta', tmp_type', tmp_age', tmp_sex', tmp_acc', tmp_eti', tmp_sbj'};
    tmp_array = cell2mat(tmp_cell);
    tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
    tmp_table.subfield = categorical(tmp_table.subfield);
    
    %
    lme = fitlme(tmp_table, 'beta~subfield*age+(1|sbj)')
    anova(lme)

    % - print results
    fprintf('<< LME results >>\n');
    for i = 2:length(lme.CoefficientNames)
        tmp = ' ';
        if lme.Coefficients{i,6} < 0.05
            tmp = '* ';
        end
        fprintf('%2s%18s: β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', tmp, lme.Coefficients{i,1}, ...
            lme.Coefficients(i,2), lme.Coefficients(i,3), double(lme.Coefficients(i,5)), lme.Coefficients(i,4), lme.Coefficients(i,6))
    end
end


%% Hippocampal activation (subfield): abuse interaction

%%%%%%%%%%%%%%%%
roi_idx_L = [ 1 8 9 10 ]; % Hippocampus, CA3, CA1, subiculum
roi_idx_R = [ 17 24 25 26 ];
roi_idx_B = [ 33 40 41 42 ];
roi_idx_list = {roi_idx_L, roi_idx_R, roi_idx_B};

con_i = 2; % 2:control
flag_i = 1; % 1:all
hemis_i = 3; % 3:bilateral

group_idx = (eti_ea > 0);
%%%%%%%%%%%%%%%%

% ---------- Linear mixed-effects model ---------- %
lme_analysis = true;
if lme_analysis            
    tmp_beta = arrayfun(@(x) metric_list_enc{2}{1}{x}, roi_idx_list{hemis_i}(2:4), 'uni', 0); % enc-ctrl, all
    tmp_beta = [tmp_beta{1} tmp_beta{2} tmp_beta{3}];

    tmp_type = [ones(1, length(group_idx)) ones(1, length(group_idx))*2 ones(1, length(group_idx))*3];

    tmp_age = sbj_age; tmp_age = [tmp_age tmp_age tmp_age];
    tmp_sex = double(sbj_sex); tmp_sex = [tmp_sex tmp_sex tmp_sex];
    tmp_acc = em_acc{6}; tmp_acc = [tmp_acc tmp_acc tmp_acc]; % full EM
    tmp_eti = eti_ea; tmp_eti = [tmp_eti tmp_eti tmp_eti]; % continuous

    tmp_sbj = 1:length(group_idx); tmp_sbj = [tmp_sbj tmp_sbj tmp_sbj];

    %
    tmp_names = {'beta', 'subfield', 'age', 'sex', 'acc', 'eti', 'sbj'};
    tmp_cell = {tmp_beta', tmp_type', tmp_age', tmp_sex', tmp_acc', tmp_eti', tmp_sbj'};
    tmp_array = cell2mat(tmp_cell);
    tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
    tmp_table.subfield = categorical(tmp_table.subfield);
    
    %
    lme = fitlme(tmp_table, 'beta~subfield*age*eti+(1|sbj)')
    tbl = anova(lme)

    % - print results
    fprintf('<< LME results >>\n');
    for i = 2:length(lme.CoefficientNames)
        tmp = ' ';
        if lme.Coefficients{i,6} < 0.05
            tmp = '* ';
        end
        fprintf('%2s%18s: β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', tmp, lme.Coefficients{i,1}, ...
            lme.Coefficients(i,2), lme.Coefficients(i,3), double(lme.Coefficients(i,5)), lme.Coefficients(i,4), lme.Coefficients(i,6))
    end
end


%% Hippocampal activation (relative activation of subfields)

%%%%%%%%%%%%%%%%
hp_index = [ 1 17 33 ];
sub_index = [ 4 7 8 9 10 ; 20 23 24 25 26 ; 36 39 40 41 42 ];
roi_idx = [3 4 5];

roi_idx_L = [ 1 8 9 10 ]; % Hippocampus, CA3, CA1, subiculum
roi_idx_R = [ 17 24 25 26 ];
roi_idx_B = [ 33 40 41 42 ];
roi_idx_list = {roi_idx_L, roi_idx_R, roi_idx_B};

con_i = 2; % 2:control
flag_i = 1; % 1:all
hemis_i = 3; % 3:bilateral

group_idx = (eti_ea > 0);
plot_color = [0.3 0.3 0.3; 223/255 153/255 153/255];
%%%%%%%%%%%

% rank ratio
sub_sort = cell(1, 3); sub_rank = cell(1, 3); hp_sort = cell(1, 3); hp_rank = cell(1, 3);
for hemis_i = 1:3
    hp_rank{hemis_i} = tiedrank(metric_list_enc{con_i}{flag_i}{hp_index(hemis_i)});

    for roi_i = 1:length(sub_index(hemis_i, :))
        sub_rank{hemis_i}{roi_i} = tiedrank(metric_list_enc{con_i}{flag_i}{sub_index(hemis_i, roi_i)});
    end

    sub_relative{hemis_i} = cellfun(@(x) x ./ hp_rank{hemis_i}, sub_rank{hemis_i}, 'uni', 0);
end


% ---------- figure ---------- %
fig_draw = true;
if fig_draw
    for g_i = 1:2
        figure;
        set(gcf, 'Position', [100 100 500 100]);

        tmp_position = [1 4.5 7.6 10.7];
        for r_i = 1:4
            subplot(1,12,[tmp_position(r_i) tmp_position(r_i)+1.2]);

            if r_i == 1
                data_all = var_list_all{var_i}{flag_i}{hp_index(hemis_i)};
            else
                data_all = sub_relative{hemis_i}{roi_idx(r_i - 1)};
            end
            data = data_all(group_idx == g_i - 1);
            data_bhv = sbj_age(group_idx == g_i - 1);

            [r, p] = corr(reshape(data_bhv,[],1), reshape(data,[],1));
            fig_scatter = scatter(data_bhv, data, 6, 'filled');
            fig_scatter.MarkerFaceColor = plot_color(g_i, :);
            fig_line = lsline;
            fig_line.LineWidth = 1.2;

            mdl = fitlm(reshape(data_bhv,[],1), reshape(data,[],1));
            x = linspace(min(reshape(data_bhv,[],1)), max(reshape(data_bhv,[],1)), 1000);
            [y,ci] = predict(mdl,x');
            hold on
            fill([x,flip(x)]', [[ci(:,1)]', flip(ci(:,2))']', [.5 .5 .5],'FaceAlpha',0.2,'linestyle','none'); 

            xlim([9 19]); xticks(10:2:18); xlabel('age');

            if r_i == 1
                ylabel('beta coefficient');
                if g_i == 1
                    ylim([-2 3.2]); yticks(-3:1:3); 
                elseif g_i ==2
                    ylim([-3.2 3.2]); yticks(-3:1:3); 
                end

            elseif r_i == 2
                ylabel('relative activation');
                if g_i == 1
                    ylim([-0.1 3.1]); yticks(0:1:3); 
                elseif g_i == 2
                    ylim([-1 15]); yticks(0:5:15); 
                end

            elseif r_i == 3
                ylabel('relative activation');
                if g_i == 1
                    ylim([-0.1 4.1]); yticks(0:1:4); 
                elseif g_i == 2
                    ylim([-0.5 8.5]); yticks(0:2:8); 
                end

            elseif r_i == 4
                ylabel('relative activation');
                if g_i == 1
                    ylim([-1 21]); yticks(0:5:20); 
                elseif g_i == 2
                    ylim([-2 30]); yticks(0:10:30); 
                end
            end
            box off
            set(gca,'XTickLabelRotation',0);
            set(gca,'LineWidth',0.8);
            set(gca,'FontName','Arial','FontSize',6, 'FontWeight','bold')
        end
    end
end


% ---------- Linear mixed-effect model ---------- %
lme_analysis = true;
if lme_analysis            
    tmp_beta = arrayfun(@(x) sub_relative{hemis_i}{x}, 1:3, 'uni', 0); % enc-ctrl, all
    tmp_beta = [tmp_beta{1} tmp_beta{2} tmp_beta{3}];

    tmp_type = [ones(1, length(group_idx)) ones(1, length(group_idx))*2 ones(1, length(group_idx))*3];

    tmp_age = sbj_age; tmp_age = [tmp_age tmp_age tmp_age];
    tmp_sex = double(sbj_sex); tmp_sex = [tmp_sex tmp_sex tmp_sex];
    tmp_acc = em_acc{6}; tmp_acc = [tmp_acc tmp_acc tmp_acc]; % full EM
    tmp_eti = eti_ea; tmp_eti = [tmp_eti tmp_eti tmp_eti]; % continuous

    tmp_sbj = 1:length(group_idx); tmp_sbj = [tmp_sbj tmp_sbj tmp_sbj];

    %
    tmp_names = {'beta', 'subfield', 'age', 'sex', 'acc', 'eti', 'sbj'};
    tmp_cell = {tmp_beta', tmp_type', tmp_age', tmp_sex', tmp_acc', tmp_eti', tmp_sbj'};
    tmp_array = cell2mat(tmp_cell);
    tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
    tmp_table.subfield = categorical(tmp_table.subfield);
    
    %
    lme = fitlme(tmp_table, 'beta~subfield*age*eti+(1|sbj)')
    tbl = anova(lme)

    % - print results
    fprintf('<< LME results >>\n');
    for i = 2:length(lme.CoefficientNames)
        tmp = ' ';
        if lme.Coefficients{i,6} < 0.05
            tmp = '* ';
        end
        fprintf('%2s%18s: β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', tmp, lme.Coefficients{i,1}, ...
            lme.Coefficients(i,2), lme.Coefficients(i,3), double(lme.Coefficients(i,5)), lme.Coefficients(i,4), lme.Coefficients(i,6))
    end
end

