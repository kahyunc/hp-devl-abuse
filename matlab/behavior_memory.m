% Episodic memory task for adolescents
% :: behavior - analysis (memory)
% :: code written by Kahyun Choi


%% define parameters

% load subject list
subject_path = 'C:\subject_info_directory';
load(fullfile(subject_path, 'subject_list_for_GLM.mat')); % 'bhv_list', 'sbj_list', 'sbj_nums', 'num_sbj'

% behavior
load(fullfile(subject_path, 'behavioral_data.mat')); % 'em_acc', 'sbj_age', 'sbj_sex'
load(fullfile(subject_path, 'survey_data.mat')); % 'sbj_svy', 'svy_label'
load(fullfile(subject_path, 'em_plot.mat')); % 'em_color', 'em_label'

% ETI
dir_survey = fullfile(subject_path, 'ETI_survey.xlsx');
sbj_eti = readtable(dir_survey);
sbj_eti = sbj_eti(sbj_nums, :);
sbj_eti = table2struct(sbj_eti); % ALL, GT, PA, EA, Q1~21
eti_ea = [sbj_eti.EA];
eti_pa = [sbj_eti.PA];


%% average memory performance

group_idx = 1:num_sbj; group_name = '(all)';
% group_idx = (eti_ea == 0); group_name = '(no emotional abuse)';
% group_idx = (eti_ea > 0); group_name = '(emotional abuse)';

label_order = [1 2 4 5 3 6];

fprintf('\n----- %s -----\n', group_name);
for acc_i = 1:length(label_order)
    curr_idx = label_order(acc_i);
    fprintf('%11s: M=%.3f, SD=%.3f\n', em_label{curr_idx}, mean(em_acc{curr_idx}(group_idx)), std(em_acc{curr_idx}(group_idx)));
end
fprintf('\t\t   (range: %.3f-%.3f)\n\n', min(em_acc{6}(group_idx)), max(em_acc{6}(group_idx)));

% ---------- figure ---------- %
fig_draw = true;
if fig_draw                 
    figure;
    set(gcf, 'Position', [100 100 300 180]);
    
    tmpData = {em_acc{1}, em_acc{2}, em_acc{4}, em_acc{5}, em_acc{3}, em_acc{6}};
    tmpData = cellfun(@(x) x(group_idx), tmpData, 'uni', 0);
    
    avg = cellfun(@(x) mean(x, 'omitnan'),tmpData);
    err = cellfun(@(x) std(x, 0, 'omitnan')/sqrt(length(x)), tmpData);

    %
    hold on
    b1 = bar(1, avg(1), 'FaceColor', [160/255 185/255 212/255], 'EdgeColor', 'none');
    b3 = bar(2, avg(3), 'FaceColor', [160/255 185/255 212/255], 'EdgeColor', 'none');
    b2 = bar(3.5, avg(2), 'FaceColor', [171/255 201/255 191/255], 'EdgeColor', 'none');
    b4 = bar(4.5, avg(4), 'FaceColor', [171/255 201/255 191/255], 'EdgeColor', 'none');
    b6 = bar(6, avg(6), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');

    %
    err1 = errorbar(1, avg(1), err(1));
    err1.Color = [0 0 0]; err1.LineStyle = 'none'; err1.LineWidth = 1; err1.CapSize = 0;
    err3 = errorbar(2, avg(3), err(3));
    err3.Color = [0 0 0]; err3.LineStyle = 'none'; err3.LineWidth = 1; err3.CapSize = 0;
    err2 = errorbar(3.5, avg(2), err(2));
    err2.Color = [0 0 0]; err2.LineStyle = 'none'; err2.LineWidth = 1; err2.CapSize = 0;
    err4 = errorbar(4.5, avg(4), err(4));
    err4.Color = [0 0 0]; err4.LineStyle = 'none'; err4.LineWidth = 1; err4.CapSize = 0;
    err6 = errorbar(6, avg(6), err(6));
    err6.Color = [0 0 0]; err6.LineStyle = 'none'; err6.LineWidth = 1; err6.CapSize = 0;

    %
    % requires hatchfill2 (MATLAB File Exchange; not redistributed here)
    hatchfill2(b3, 'single', 'HatchAngle', 60, 'HatchDensity', 15, 'HatchColor', 'k');
    hatchfill2(b4, 'single', 'HatchAngle', 60, 'HatchDensity', 15, 'HatchColor', 'k');
    hatchfill2(b6, 'single', 'HatchAngle', 60, 'HatchDensity', 15, 'HatchColor', 'k');
 
    tmpLabel = {'what', 'what-when', 'where', 'where-when', 'full EM'};
    xticks([1 2 3.5 4.5 6]); xticklabels(tmpLabel); 
    
    ylim([0 1]); ylabel('accuracy');
    ylim([0.5 1]); ylabel('accuracy');

    set(gca,'LineWidth',0.6);
    set(gca,'FontName','Arial','FontSize',7, 'FontWeight','normal')
    disp(' ');
end


%% correlation with age (No emotional abuse, Emotional abuse)

tmp_group = (eti_ea > 0);
tmp_category = categorical(tmp_group, logical([0 1]), {'no-EA', 'EA'});

% ---------- figure ---------- %
fig_draw = true;
if fig_draw                 
    %% full EM
    tmp_data = em_acc{6};
    group_color = [0.3 0.3 0.3; 0.7 0.7 0.7];

    figure;
    set(gcf, 'Position', [100 100 600 130]);

    % - box plot -
    subplot(1,4,1);

    fig_box = boxchart(tmp_data, 'GroupByColor', tmp_category); hold on
    for g_i = 1:length(fig_box)
        fig_box(g_i).BoxEdgeColor = 'k';
        fig_box(g_i).BoxFaceColor = group_color(g_i, :);
        fig_box(g_i).BoxFaceAlpha = 1;
        fig_box(g_i).MarkerStyle = 'none';
        fig_box(g_i).LineWidth = 0.8;
    end
    fig_dot = []; dot_density = 3; % 0, 1, 2, 3
    bar_width_half = 0.3 - dot_density * 0.07;

    for graph_i = 1:size(tmp_data, 2)
        curr_data = tmp_data(tmp_group == graph_i-1);

        y_draw = sort(curr_data); 
        y_draw = y_draw(~isnan(y_draw));
        y_unique = unique(y_draw);
        num_repeated = arrayfun(@(x) sum(x == y_draw), y_unique);
        
        x_draw = arrayfun(@(x) linspace(-bar_width_half, bar_width_half, x - mod(x,2)),num_repeated, 'uni', 0);
        x_draw(mod(num_repeated, 2) == 1) = cellfun(@(x) [x 0], x_draw(mod(num_repeated,2) == 1), 'uni', 0);
        x_draw = cell2mat(reshape(x_draw, 1, [])) + graph_i/2 + 0.05;
            
        fig_dot_temp = scatter(x_draw, y_draw, 4,' k','filled', 'HandleVisibility', 'off');
        fig_dot_temp.MarkerFaceAlpha = 0.2;
        fig_dot = [fig_dot fig_dot_temp];
    end
    ylim([0 1.05]); yticks(0:0.2:1); ylabel('accuracy');
    xticklabels(em_label{6});
    [~, p] = ttest2(tmp_data(tmp_group == 0), tmp_data(tmp_group > 0));
    hold off   
    box off
    set(gca, 'LineWidth', 0.6);
    set(gca, 'FontName', 'Arial', 'FontSize', 6, 'FontWeight', 'normal');

    % - correlation - 
    subplot(1,4,2.3); % no emotional abuse
    [r, p] = corr(reshape(sbj_age(tmp_group == 0),[],1), reshape(tmp_data(tmp_group == 0),[],1));
    fig_scatter = scatter(sbj_age(tmp_group == 0), tmp_data(tmp_group == 0), 6, 'filled');
    fig_scatter.MarkerFaceColor = group_color(1,:);
    fig_line = lsline;
    fig_line.LineWidth = 1.2;
    mdl = fitlm(reshape(sbj_age(tmp_group == 0),[],1), reshape(tmp_data(tmp_group == 0),[],1));
    x = linspace(min(reshape(sbj_age(tmp_group == 0),[],1)), max(reshape(sbj_age(tmp_group == 0),[],1)), 1000);
    [y,ci] = predict(mdl,x');
    hold on
    fill([x,flip(x)]', [[ci(:,1)]', flip(ci(:,2))']', [.5 .5 .5],'FaceAlpha',0.2,'linestyle','none'); 
    if p < 0.05; fig_line.Color = [1 0.5 0.5]; end
    xlabel('age', 'FontSize', 7); ylabel('accuracy', 'FontSize', 7); xlim([9 19]); xticks(10:2:18); 
    ylim([0 1.05]);
    set(gca, 'FontSize', 6);
    set(gca,'LineWidth',0.8);

    subplot(1,4,3.4); % emotional abuse
    [r, p] = corr(reshape(sbj_age(tmp_group == 1),[],1), reshape(tmp_data(tmp_group == 1),[],1));
    fig_scatter = scatter(sbj_age(tmp_group == 1), tmp_data(tmp_group == 1), 6, 'filled');
    fig_scatter.MarkerFaceColor = group_color(2,:);
    fig_line = lsline;
    fig_line.LineWidth = 1.2;
    mdl = fitlm(reshape(sbj_age(tmp_group == 1),[],1), reshape(tmp_data(tmp_group == 1),[],1));
    x = linspace(min(reshape(sbj_age(tmp_group == 1),[],1)), max(reshape(sbj_age(tmp_group == 1),[],1)), 1000);
    [y,ci] = predict(mdl,x');
    hold on
    fill([x,flip(x)]', [[ci(:,1)]', flip(ci(:,2))']', [.5 .5 .5],'FaceAlpha',0.2,'linestyle','none'); 
    if p < 0.05; fig_line.Color = [1 0.5 0.5]; end
    xlabel('age', 'FontSize', 7); ylabel('accuracy', 'FontSize', 7); xlim([9 19]); xticks(10:2:18); 
    ylim([0 1.05]);
    set(gca, 'LineWidth', 0.6);
    set(gca, 'FontName', 'Arial', 'FontSize', 6, 'FontWeight', 'normal');


    %% where-when
    tmp_data = em_acc{5};
    group_color = [76 120 105; 171 201 191]./255;

    figure;
    set(gcf, 'Position', [100 100 600 130]);

    % - box plot -
    subplot(1,4,1);

    fig_box = boxchart(tmp_data, 'GroupByColor', tmp_category); hold on
    for g_i = 1:length(fig_box)
        fig_box(g_i).BoxEdgeColor = 'k';
        fig_box(g_i).BoxFaceColor = group_color(g_i, :);
        fig_box(g_i).BoxFaceAlpha = 1;
        fig_box(g_i).MarkerStyle = 'none';
        fig_box(g_i).LineWidth = 0.8;
    end
    fig_dot = []; dot_density = 3; % 0, 1, 2, 3
    bar_width_half = 0.3 - dot_density * 0.07;

    for graph_i = 1:size(tmp_data, 2)
        curr_data = tmp_data(tmp_group == graph_i-1);

        y_draw = sort(curr_data); 
        y_draw = y_draw(~isnan(y_draw));
        y_unique = unique(y_draw);
        num_repeated = arrayfun(@(x) sum(x == y_draw),y_unique);
        
        x_draw = arrayfun(@(x) linspace(-bar_width_half, bar_width_half, x - mod(x,2)),num_repeated, 'uni', 0);
        x_draw(mod(num_repeated, 2) == 1) = cellfun(@(x) [x 0], x_draw(mod(num_repeated,2) == 1), 'uni', 0);
        x_draw = cell2mat(reshape(x_draw, 1, [])) + graph_i/2 + 0.05;
            
        fig_dot_temp = scatter(x_draw, y_draw, 4,' k','filled', 'HandleVisibility', 'off');
        fig_dot_temp.MarkerFaceAlpha = 0.2;
        fig_dot = [fig_dot fig_dot_temp];
    end
    ylim([0 1.05]); yticks(0:0.2:1); ylabel('accuracy');
    xticklabels(em_label{5});
    [~, p] = ttest2(tmp_data(tmp_group == 0), tmp_data(tmp_group > 0));
    hold off   
    box off
    set(gca, 'LineWidth', 0.6);
    set(gca, 'FontName', 'Arial', 'FontSize', 6, 'FontWeight', 'normal');

    % - correlation - 
    subplot(1,4,2.3); % no emotional abuse
    [r, p] = corr(reshape(sbj_age(tmp_group == 0),[],1), reshape(tmp_data(tmp_group == 0),[],1));
    fig_scatter = scatter(sbj_age(tmp_group == 0), tmp_data(tmp_group == 0), 6, 'filled');
    fig_scatter.MarkerFaceColor = group_color(1,:);
    fig_line = lsline;
    fig_line.LineWidth = 1.2;
    mdl = fitlm(reshape(sbj_age(tmp_group == 0),[],1), reshape(tmp_data(tmp_group == 0),[],1));
    x = linspace(min(reshape(sbj_age(tmp_group == 0),[],1)), max(reshape(sbj_age(tmp_group == 0),[],1)), 1000);
    [y,ci] = predict(mdl,x');
    hold on
    fill([x,flip(x)]', [[ci(:,1)]', flip(ci(:,2))']', [.5 .5 .5],'FaceAlpha',0.2,'linestyle','none'); 
    if p < 0.05; fig_line.Color = [1 0.5 0.5]; end
    xlabel('age'); ylabel('accuracy'); xlim([9 19]); xticks(10:2:18); 
    ylim([0 1.05]);
    set(gca, 'LineWidth', 0.6);
    set(gca, 'FontName', 'Arial', 'FontSize', 6, 'FontWeight', 'normal');
    
    subplot(1,4,3.4); % emotional abuse
    [r, p] = corr(reshape(sbj_age(tmp_group == 1),[],1), reshape(tmp_data(tmp_group == 1),[],1));
    fig_scatter = scatter(sbj_age(tmp_group == 1), tmp_data(tmp_group == 1), 6, 'filled');
    fig_scatter.MarkerFaceColor = group_color(2,:);
    fig_line = lsline;
    fig_line.LineWidth = 1.2;
    mdl = fitlm(reshape(sbj_age(tmp_group == 1),[],1), reshape(tmp_data(tmp_group == 1),[],1));
    x = linspace(min(reshape(sbj_age(tmp_group == 1),[],1)), max(reshape(sbj_age(tmp_group == 1),[],1)), 1000);
    [y,ci] = predict(mdl,x');
    hold on
    fill([x,flip(x)]', [[ci(:,1)]', flip(ci(:,2))']', [.5 .5 .5],'FaceAlpha',0.2,'linestyle','none'); 
    if p < 0.05; fig_line.Color = [1 0.5 0.5]; end
    xlabel('age', 'FontSize', 7); ylabel('accuracy', 'FontSize', 7); xlim([9 19]); xticks(10:2:18); 
    ylim([0 1.05]);
    set(gca, 'LineWidth', 0.6);
    set(gca, 'FontName', 'Arial', 'FontSize', 6, 'FontWeight', 'normal');
    


    %% what-when
    tmp_data = em_acc{4};
    group_color = [70 111 156; 160 185 212]./255;

    figure;
    set(gcf, 'Position', [100 100 600 130]);

    % - box plot -
    subplot(1,4,1);
    fig_box = boxchart(tmp_data, 'GroupByColor', tmp_category); hold on
    for g_i = 1:length(fig_box)
        fig_box(g_i).BoxEdgeColor = 'k';
        fig_box(g_i).BoxFaceColor = group_color(g_i, :);
        fig_box(g_i).BoxFaceAlpha = 1;
        fig_box(g_i).MarkerStyle = 'none';
        fig_box(g_i).LineWidth = 0.8;
    end
    fig_dot = []; dot_density = 3; % 0, 1, 2, 3
    bar_width_half = 0.3 - dot_density * 0.07;

    for graph_i = 1:size(tmp_data, 2)
        curr_data = tmp_data(tmp_group == graph_i-1);

        y_draw = sort(curr_data); 
        y_draw = y_draw(~isnan(y_draw));
        y_unique = unique(y_draw);
        num_repeated = arrayfun(@(x) sum(x == y_draw),y_unique);
        
        x_draw = arrayfun(@(x) linspace(-bar_width_half, bar_width_half, x - mod(x,2)),num_repeated, 'uni', 0);
        x_draw(mod(num_repeated, 2) == 1) = cellfun(@(x) [x 0], x_draw(mod(num_repeated,2) == 1), 'uni', 0);
        x_draw = cell2mat(reshape(x_draw, 1, [])) + graph_i/2 + 0.05;
            
        fig_dot_temp = scatter(x_draw, y_draw, 4,' k','filled', 'HandleVisibility', 'off');
        fig_dot_temp.MarkerFaceAlpha = 0.2;
        fig_dot = [fig_dot fig_dot_temp];
    end
    ylim([0 1.05]); yticks(0:0.2:1); ylabel('accuracy');
    xticklabels(em_label{4});
    [~, p] = ttest2(tmp_data(tmp_group == 0), tmp_data(tmp_group > 0));
    hold off   
    box off
    set(gca, 'LineWidth', 0.6);
    set(gca, 'FontName', 'Arial', 'FontSize', 6, 'FontWeight', 'normal');

    % - correlation - 
    subplot(1,4,2.3); % no emotional abuse
    [r, p] = corr(reshape(sbj_age(tmp_group == 0),[],1), reshape(tmp_data(tmp_group == 0),[],1));
    fig_scatter = scatter(sbj_age(tmp_group == 0), tmp_data(tmp_group == 0), 6, 'filled');
    fig_scatter.MarkerFaceColor = group_color(1,:);
    fig_line = lsline;
    fig_line.LineWidth = 1.2;

    mdl = fitlm(reshape(sbj_age(tmp_group == 0),[],1), reshape(tmp_data(tmp_group == 0),[],1));
    x = linspace(min(reshape(sbj_age(tmp_group == 0),[],1)), max(reshape(sbj_age(tmp_group == 0),[],1)), 1000);
    [y,ci] = predict(mdl,x');
    hold on
    fill([x,flip(x)]', [[ci(:,1)]', flip(ci(:,2))']', [.5 .5 .5],'FaceAlpha',0.2,'linestyle','none'); 
    if p < 0.05; fig_line.Color = [1 0.5 0.5]; end
    xlabel('age'); ylabel('accuracy'); xlim([9 19]); xticks(10:2:18); 
    ylim([0 1.05]);
    set(gca, 'LineWidth', 0.6);
    set(gca, 'FontName', 'Arial', 'FontSize', 6, 'FontWeight', 'normal');
    
    subplot(1,4,3.4); % emotional abuse
    [r, p] = corr(reshape(sbj_age(tmp_group == 1),[],1), reshape(tmp_data(tmp_group == 1),[],1));
    fig_scatter = scatter(sbj_age(tmp_group == 1), tmp_data(tmp_group == 1), 6, 'filled');
    fig_scatter.MarkerFaceColor = group_color(2,:);
    fig_line = lsline;
    fig_line.LineWidth = 1.2;
    mdl = fitlm(reshape(sbj_age(tmp_group == 1),[],1), reshape(tmp_data(tmp_group == 1),[],1));
    x = linspace(min(reshape(sbj_age(tmp_group == 1),[],1)), max(reshape(sbj_age(tmp_group == 1),[],1)), 1000);
    [y,ci] = predict(mdl,x');
    hold on
    fill([x,flip(x)]', [[ci(:,1)]', flip(ci(:,2))']', [.5 .5 .5],'FaceAlpha',0.2,'linestyle','none'); 
    if p < 0.05; fig_line.Color = [1 0.5 0.5]; end
    xlabel('age'); ylabel('accuracy'); xlim([9 19]); xticks(10:2:18); 
    ylim([0 1.05]);
    set(gca, 'LineWidth', 0.6);
    set(gca, 'FontName', 'Arial', 'FontSize', 6, 'FontWeight', 'normal');

end

% ---------- stats ---------- %
fprintf('\n----- no emotional abuse -----\n');
label_order = [6 1 2 4 5 3];
for acc_i = 1:length(label_order)
    curr_idx = label_order(acc_i);
    [r,p] = corr(sbj_age(tmp_group == 0)', em_acc{curr_idx}(tmp_group == 0)');
    fprintf('%11s: r=%.3f, p=%.4f\n', em_label{curr_idx}, r, p);
end
disp(' ');

fprintf('\n----- emotional abuse -----\n');
label_order = [6 1 2 4 5 3];
for acc_i = 1:length(label_order)
    curr_idx = label_order(acc_i);
    [r,p] = corr(sbj_age(tmp_group == 1)', em_acc{curr_idx}(tmp_group == 1)');
    fprintf('%11s: r=%.3f, p=%.4f\n', em_label{curr_idx}, r, p);
end
disp(' ');


%% abuse-by-age interaction: trial-wise (full EM, where-when, what-when)

% ---------- Linear mixed-effects model ---------- %
lme_analysis = true;
if lme_analysis                 
    
    acc_i = [ 6 ]; % full EM
    % acc_i = [ 5 ]; % where-when
    % acc_i = [ 4 ]; % what-when
    % acc_i = [ 3 ]; % what-where
    
    % - trial-wise -
    tmp = cellfun(@(x) x.organized.trial_acc(acc_i, :), bhv_list, 'uni', 0);
    tmp_acc = cell2mat(tmp);
    tmp_trial = repmat(1:8, 1, length(bhv_list));
    tmp_age = sbj_age; tmp_age = reshape(repmat(tmp_age, 8, 1), 1, length(bhv_list)*8);
    tmp_eti = eti_ea; tmp_eti = reshape(repmat(tmp_eti, 8, 1), 1, length(bhv_list)*8); % continuous
    tmp_sex = double(sbj_sex); tmp_sex = reshape(repmat(tmp_sex, 8, 1), 1, length(bhv_list)*8);
    tmp_sbj = 1:length(bhv_list); tmp_sbj = reshape(repmat(tmp_sbj, 8, 1), 1, length(bhv_list)*8);
    tmp_names = {'acc', 'trial', 'age', 'eti', 'sex', 'sbj'};
    tmp_cell = {tmp_acc', tmp_trial', tmp_age', tmp_eti', tmp_sex', tmp_sbj'};
    tmp_array = cell2mat(tmp_cell);
    tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
   
    lme = fitlme(tmp_table, 'acc~age*eti+(trial|sbj)')
    anova(lme)

    % - print results
    fprintf('<< LME results >>\n');
    for i = 2:4
        fprintf('%8s: β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', lme.Coefficients{i,1}, ...
            lme.Coefficients{i,2}, lme.Coefficients{i,3}, double(lme.Coefficients{i,5}), lme.Coefficients{i,4}, lme.Coefficients{i,6})
    end
end


%% abuse-by-age-by-component interaction: where-when, what-when in one model

% ---------- Linear mixed-effects model ---------- %
lme_analysis = true;
if lme_analysis                 
    acc_idx = [ 4 5 ];
    
    tmp_acc = [em_acc{acc_idx(1)} em_acc{acc_idx(2)}];
    tmp_type = [ones(1, length(em_acc{1})) ones(1, length(em_acc{1}))*2];
    tmp_age = sbj_age; tmp_age = [tmp_age tmp_age];
    tmp_eti = eti_ea; tmp_eti = [tmp_eti tmp_eti];    
    tmp_sex = double(sbj_sex); tmp_sex = [tmp_sex tmp_sex];
    tmp_sbj = 1:length(em_acc{1}); tmp_sbj = [tmp_sbj tmp_sbj];
    tmp_names = {'acc', 'type', 'age', 'eti', 'sex', 'sbj'};
    tmp_cell = {tmp_acc', tmp_type', tmp_age', tmp_eti', tmp_sex', tmp_sbj'};
    tmp_array = cell2mat(tmp_cell);
    tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
    tmp_table.type = categorical(tmp_table.type);

    lme = fitlme(tmp_table, 'acc~type*age*eti+(type|sbj)')
    anova(lme)

    % - print results
    fprintf('<< LME results >>\n');
    for i = 2:length(lme.CoefficientNames)
        tmp = ' ';
        if lme.Coefficients{i,6} < 0.05
            tmp = '* ';
        end
        fprintf('%2s%12s: β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', tmp, lme.Coefficients{i,1}, ...
            lme.Coefficients{i,2}, lme.Coefficients{i,3}, double(lme.Coefficients{i,5}), lme.Coefficients{i,4}, lme.Coefficients{i,6})
    end
end


%% error patterns (binding preserved ratio)


tmp_group = (eti_ea > 0);
tmp_category = categorical(tmp_group, logical([0 1]), {'no-EA', 'EA'});

% set-up: error bias
whatwhere_ans = cellfun(@(x) x.expInfo.animal.*10 + x.expInfo.location, bhv_list, 'uni', 0);
whatwhere_resp = cellfun(@(x) x.reenact.animal.*10 + x.reenact.location, bhv_list, 'uni', 0);

whatwhere_list = cell(1, length(whatwhere_ans));
for sbj_i = 1:length(whatwhere_ans)
    for x_i = 1:size(whatwhere_ans{sbj_i}, 1)
        for y_i = 1:size(whatwhere_ans{sbj_i}, 2)
            tmp = sum(whatwhere_ans{sbj_i}(x_i, :) == whatwhere_resp{sbj_i}(x_i, y_i));
            whatwhere_list{sbj_i}(x_i, y_i) = tmp;
        end
    end
end
whatwhen_list = cellfun(@(x) x.reenact.animalAcc, bhv_list, 'uni', 0);
wherewhen_list = cellfun(@(x) x.reenact.locationAcc, bhv_list, 'uni', 0);
fullem_list = cellfun(@(x) x.reenact.animalAcc + x.reenact.locationAcc, bhv_list, 'uni', 0);

%
x_binding = []; x_fullem = []; o_whatwhen = []; o_wherewhen = []; o_whatwhere = []; 
for sbj_i = 1:length(fullem_list)
    wrong_flag = (fullem_list{sbj_i} < 2);
    x_fullem(sbj_i) = sum((wrong_flag==1), 'all');
    
    o_whatwhen(sbj_i)  = sum(whatwhen_list{sbj_i}(wrong_flag) == 1);
    o_wherewhen(sbj_i) = sum(wherewhen_list{sbj_i}(wrong_flag) == 1);
    o_whatwhere(sbj_i) = sum(whatwhere_list{sbj_i}(wrong_flag) == 1);
    
    tmp = sum(whatwhen_list{sbj_i}(wrong_flag) == 0 & wherewhen_list{sbj_i}(wrong_flag) == 0 & whatwhere_list{sbj_i}(wrong_flag) == 0);
    x_binding(sbj_i) = tmp;
end

%
r_whatwhen  = arrayfun(@(num) o_whatwhen(num)/x_fullem(num), 1:length(fullem_list));
r_wherewhen = arrayfun(@(num) o_wherewhen(num)/x_fullem(num), 1:length(fullem_list));
r_whatwhere = arrayfun(@(num) o_whatwhere(num)/x_fullem(num), 1:length(fullem_list));
r_binding_x = arrayfun(@(num) x_binding(num)/x_fullem(num), 1:length(fullem_list));

error_r_list = {r_whatwhen, r_wherewhen, r_whatwhere, r_binding_x};
tmpLabel = {'what-when', 'where-when', 'what-where', 'failure'};

group_list = {(eti_ea == 0), (eti_ea > 0)}; group_name = 'eti';
group_label = {'no emotional abuse', 'emotional abuse'};
group_color = [0.3 0.3 0.3; 223/255 153/255 153/255];

% --- Figure 1: Group difference (box plots) ---
draw_fig1 = true;
if draw_fig1
    type_colors = {[70 111 156; 160 185 212]./255, ...    % what-when (blue)
                   [76 120 105; 171 201 191]./255, ...    % where-when (green)
                   [115 75 140; 185 160 210]./255, ...    % what-where (purple)
                   [0.3 0.3 0.3; 0.7 0.7 0.7]};           % failure  (gray)

    figure;
    set(gcf, 'Position', [100 100 650 150]);

    fprintf('\n== Group difference (t-test) ==\n');
    for acc_i = 1:length(error_r_list)
        subplot(1, 4, acc_i);

        curr_group_color = type_colors{acc_i};

        fig_box = boxchart(error_r_list{acc_i}(:), 'GroupByColor', tmp_category); hold on
        for g_i = 1:length(fig_box)
            fig_box(g_i).BoxEdgeColor = 'k';
            fig_box(g_i).BoxFaceColor = curr_group_color(g_i, :);
            fig_box(g_i).BoxFaceAlpha = 1;
            fig_box(g_i).MarkerStyle = 'none';
            fig_box(g_i).LineWidth = 0.8;
        end

        fig_dot = []; dot_density = 3;
        bar_width_half = 0.3 - dot_density * 0.07;
        for graph_i = 1:2
            curr_data = error_r_list{acc_i}(tmp_group == graph_i-1);
            y_draw = sort(curr_data(:));
            y_draw = y_draw(~isnan(y_draw));
            y_unique = unique(y_draw);
            num_repeated = arrayfun(@(x) sum(x==y_draw), y_unique);
            x_draw = arrayfun(@(x) linspace(-bar_width_half, bar_width_half, x - mod(x,2)), num_repeated, 'UniformOutput', false);
            x_draw(mod(num_repeated,2)==1) = cellfun(@(x) [x 0], x_draw(mod(num_repeated,2)==1), 'UniformOutput', false);
            x_draw = cell2mat(reshape(x_draw,1,[])) + graph_i/2 + 0.05;
            fig_dot_temp = scatter(x_draw, y_draw', 4, ' k', 'filled', 'HandleVisibility', 'off');
            fig_dot_temp.MarkerFaceAlpha = 0.2;
            fig_dot = [fig_dot fig_dot_temp];
        end
        d1 = error_r_list{acc_i}(group_list{1});
        d2 = error_r_list{acc_i}(group_list{2});
        [~, p, ~, stats] = ttest2(d1(:), d2(:));
        fprintf('%s: t(%d) = %.3f, p = %.4f\n', tmpLabel{acc_i}, stats.df, stats.tstat, p);

        xticklabels(tmpLabel(acc_i));
        ylabel('ratio'); ylim([-0.05 1.05]);
        hold off;
        box off;
        set(gca, 'LineWidth', 0.6);
        set(gca, 'FontName', 'Arial', 'FontSize', 6, 'FontWeight', 'normal');
    end
end


% --- Figure 2: Group correlation (age) ---
draw_fig2 = true;
if draw_fig2
    all_correct_idx = (x_fullem == 0);

    group_idx = (eti_ea == 0); group_names = {'no emotional abuse', 'emotional abuse'};

    tmpData = {r_whatwhen, r_wherewhen, r_whatwhere, r_binding_x};

    fprintf('\n== Age correlation per group ==\n');
    for fig_i = 1:length(tmpData)
        figure;
        set(gcf, 'Position', [100 100 280 140]);

        currData = tmpData{fig_i};

        subplot(1, 2, 1);
        tmp_idx = (~all_correct_idx & group_idx);
        x_reg = reshape(sbj_age(tmp_idx),[],1); y_reg = reshape(currData(tmp_idx),[],1);
        nan_reg = isnan(x_reg) | isnan(y_reg); x_reg(nan_reg) = []; y_reg(nan_reg) = [];
        [r, p] = corr(x_reg, y_reg);
        fig_scatter = scatter(x_reg, y_reg, 6, 'filled'); hold on
        fig_line = lsline;
        mdl = fitlm(x_reg, y_reg);
        xx = linspace(min(x_reg), max(x_reg), 1000); [~, ci] = predict(mdl, xx');
        fill([xx, flip(xx)]', [ci(:,1)', flip(ci(:,2))']', [.5 .5 .5], 'FaceAlpha', 0.2, 'linestyle', 'none');
        fig_scatter.MarkerFaceColor = group_color(1,:);
        fig_scatter.SizeData = 6;
        fig_line.LineWidth = 1.2;
        if p < 0.05; fig_line.Color = [1 0.5 0.5]; end
        ylabel(tmpLabel{fig_i});  xlabel('age');
        xlim([9 19]);  xticks(10:2:18); ylim([-0.05 1.05]);
        box off;
        set(gca, 'LineWidth', 0.6);
        set(gca, 'FontName', 'Arial', 'FontSize', 6, 'FontWeight', 'normal');
        fprintf('%s - %s: r = %.3f, p = %.4f\n', tmpLabel{fig_i}, group_names{1}, r, p);

        subplot(1, 2, 2);
        tmp_idx = (~all_correct_idx & ~group_idx);
        x_reg = reshape(sbj_age(tmp_idx),[],1); y_reg = reshape(currData(tmp_idx),[],1);
        nan_reg = isnan(x_reg) | isnan(y_reg); x_reg(nan_reg) = []; y_reg(nan_reg) = [];
        [r, p] = corr(x_reg, y_reg);
        fig_scatter = scatter(x_reg, y_reg, 6, 'filled'); hold on
        fig_line = lsline;
        mdl = fitlm(x_reg, y_reg);
        xx = linspace(min(x_reg), max(x_reg), 1000); [~, ci] = predict(mdl, xx');
        fill([xx, flip(xx)]', [ci(:,1)', flip(ci(:,2))']', [.5 .5 .5], 'FaceAlpha', 0.2, 'linestyle', 'none');
        fig_scatter.MarkerFaceColor = group_color(2,:);
        fig_scatter.SizeData = 6;
        fig_line.LineWidth = 1.2;
        if p < 0.05; fig_line.Color = [1 0.5 0.5]; end
        ylabel('');  xlabel('age');
        xlim([9 19]);  xticks(10:2:18); ylim([-0.05 1.05]);
        box off;
        set(gca, 'LineWidth', 0.6);
        set(gca, 'FontName', 'Arial', 'FontSize', 6, 'FontWeight', 'normal');
        fprintf('%s - %s: r = %.3f, p = %.4f\n', tmpLabel{fig_i}, group_names{2}, r, p);
    end
end


%% behavioral sensitivity: anxiety, depression, physical abuse and sex

acc_idx = [6 5 4]; % full EM, where-when, what-when
cesd = sbj_svy{ismember(svy_label, 'CES_D')}(:);
stais = sbj_svy{ismember(svy_label, 'STAI_S')}(:);
model_formulas = {'acc~age*eti+cesd+stais+(trial|sbj)', ...
    'acc~age*eti+pa+(trial|sbj)', 'acc~age*eti+age*pa+(trial|sbj)', ...
    'acc~age*eti+sex+(trial|sbj)', 'acc~age*pa+(trial|sbj)'};
for acc_i = acc_idx
    tmp_acc = cell2mat(cellfun(@(x) x.organized.trial_acc(acc_i, :), bhv_list, 'uni', 0));
    tmp_table = table(tmp_acc(:), repmat((1:8)', num_sbj, 1), ...
        repelem(sbj_age(:), 8), repelem(eti_ea(:), 8), repelem(eti_pa(:), 8), ...
        repelem(double(sbj_sex(:)), 8), repelem(cesd, 8), repelem(stais, 8), ...
        repelem((1:num_sbj)', 8), ...
        'VariableNames', {'acc', 'trial', 'age', 'eti', 'pa', 'sex', 'cesd', 'stais', 'sbj'});
    for model_i = 1:length(model_formulas)
        fprintf('\naccuracy %d: %s\n', acc_i, model_formulas{model_i});
        lme = fitlme(tmp_table, model_formulas{model_i})
        anova(lme)
    end
end


%% component model: covariate sensitivity

tmp_table = table([em_acc{4}(:); em_acc{5}(:)], categorical(repelem((1:2)', num_sbj)), ...
    repmat(sbj_age(:), 2, 1), repmat(eti_ea(:), 2, 1), repmat(eti_pa(:), 2, 1), ...
    repmat(double(sbj_sex(:)), 2, 1), repmat(cesd, 2, 1), repmat(stais, 2, 1), ...
    repmat((1:num_sbj)', 2, 1), ...
    'VariableNames', {'acc', 'type', 'age', 'eti', 'pa', 'sex', 'cesd', 'stais', 'sbj'});
model_formulas = {'acc~type*age*eti+cesd+stais+(type|sbj)', ...
    'acc~type*age*eti+pa+(type|sbj)', 'acc~type*age*eti+type*age*pa+(type|sbj)', ...
    'acc~type*age*eti+sex+(type|sbj)'};
for model_i = 1:length(model_formulas)
    lme = fitlme(tmp_table, model_formulas{model_i})
    anova(lme)
end


%% physical abuse: model comparison

acc_idx = [6 5 4];
age_c = sbj_age - mean(sbj_age);
for acc_i = acc_idx
    tmp_acc = cell2mat(cellfun(@(x) x.organized.trial_acc(acc_i, :), bhv_list, 'uni', 0));
    tmp_table = table(tmp_acc(:), repmat((1:8)', num_sbj, 1), ...
        repelem(age_c(:), 8), repelem(eti_ea(:), 8), repelem(eti_pa(:), 8), ...
        repelem((1:num_sbj)', 8), 'VariableNames', {'acc', 'trial', 'age', 'eti', 'pa', 'sbj'});
    lme_base = fitlme(tmp_table, 'acc~age*eti+(trial|sbj)', 'FitMethod', 'ML');
    lme_pa = fitlme(tmp_table, 'acc~age*eti+pa+(trial|sbj)', 'FitMethod', 'ML');
    lme_pa_interaction = fitlme(tmp_table, 'acc~age*eti+age*pa+(trial|sbj)', 'FitMethod', 'ML');
    fprintf('\naccuracy %d: PA model comparison\n', acc_i);
    disp(compare(lme_base, lme_pa));
    disp(compare(lme_base, lme_pa_interaction));
    disp([lme_base.ModelCriterion; lme_pa.ModelCriterion; lme_pa_interaction.ModelCriterion]);
end


%% exposure groups: none, physical only, emotional with or without physical abuse

acc_idx = [6 5 4];
abuse_group = ones(num_sbj, 1);
abuse_group(eti_ea(:) == 0 & eti_pa(:) > 0) = 2;
abuse_group(eti_ea(:) > 0) = 3;
for acc_i = acc_idx
    data_all = em_acc{acc_i}(:);
    [p, tbl] = anova1(data_all, abuse_group, 'off');
    disp(tbl);
    fprintf('accuracy %d, group ANOVA: p=%.4f\n', acc_i, p);
    for g_i = 1:3
        idx = (abuse_group == g_i);
        [r, p] = corr(sbj_age(idx)', data_all(idx));
        fprintf('group %d: r=%.3f, p=%.4f\n', g_i, r, p);
    end
    tmp_acc = cell2mat(cellfun(@(x) x.organized.trial_acc(acc_i, :), bhv_list, 'uni', 0));
    tmp_table = table(tmp_acc(:), repmat((1:8)', num_sbj, 1), ...
        repelem(age_c(:), 8), categorical(repelem(abuse_group, 8)), ...
        repelem((1:num_sbj)', 8), 'VariableNames', {'acc', 'trial', 'age', 'group', 'sbj'});
    lme = fitlme(tmp_table, 'acc~age*group+(trial|sbj)')
    anova(lme)
end


%% age subgroups and physical abuse group comparison

acc_idx = [6 5 4];
age_group = (sbj_age >= mean(sbj_age));
for acc_i = acc_idx
    for age_i = 0:1
        idx0 = (eti_ea == 0) & (age_group == age_i);
        idx1 = (eti_ea > 0) & (age_group == age_i);
        [~, p, ~, stats] = ttest2(em_acc{acc_i}(idx0), em_acc{acc_i}(idx1));
        fprintf('accuracy %d, age group %d: t(%d)=%.3f, p=%.4f\n', acc_i, age_i, stats.df, stats.tstat, p);
    end
end
[~, p, ~, stats] = ttest2(em_acc{6}(eti_pa == 0), em_acc{6}(eti_pa > 0));
fprintf('physical abuse, full EM: t(%d)=%.3f, p=%.4f\n', stats.df, stats.tstat, p);
for g_i = 0:1
    idx = ((eti_pa > 0) == g_i);
    [r, p] = corr(sbj_age(idx)', em_acc{6}(idx)');
    fprintf('physical abuse group %d: r=%.3f, p=%.4f\n', g_i, r, p);
end

