% Episodic memory task for adolescents
% :: seed-based functional connectivity - analysis (statistics)
% :: code written by Kahyun Choi


%% define parameters

% load subject list
subject_path = 'C:\subject_info_directory';
bhv_nums = load(fullfile(subject_path, 'subject_list_for_GLM.mat')); % 'bhv_list', 'sbj_list', 'sbj_nums', 'num_sbj'

% behavior
load(fullfile(subject_path, 'behavioral_data.mat')); % 'em_acc', 'sbj_age', 'sbj_sex'

% ETI
dir_survey = fullfile(subject_path, 'ETI_survey.xlsx');
sbj_eti = readtable(dir_survey);
sbj_eti = sbj_eti(sbj_nums, :);
sbj_eti = table2struct(sbj_eti); % ALL, GT, PA, EA, Q1~21
eti_ea = [sbj_eti.EA];


% connectivity
dir_conn_first = 'C:\connectivity_first_level_directory'; % directory of CONN first-level analysis results

Z_surf = cell(1,3);
for con_i = 1:3
    roi_file_name = sprintf('resultsROI_Condition%.3d.mat', con_i);

    %%% roi_file: DOF(n2), SE(sbj*n2), Z(n1*n2*sbj), names(n1), names2(n2), regressors(n1), xyz(n2)
    results = load(fullfile(dir_conn_first, roi_file_name));

    %%% ROI
    % seed ROI (names): hippocampus
    target_name = {'LaHPC', 'LpHPC', 'LCA1', 'LCA23DG', 'LSub', ...
               'RaHPC', 'RpHPC', 'RCA1', 'RCA23DG', 'RSub', ...
               'BaHPC', 'BpHPC', 'BCA1', 'BCA23DG', 'BSub'};

    seed_hp_idx1 = find(contains(results.names, 'Hippocampus'));
    seed_hp_idx2 = cellfun(@(x) find(contains(results.names, x), 1), target_name);
    roi_seed_idx = [seed_hp_idx1, seed_hp_idx2];

    roi_seed_names = {'LHP', 'RHP', 'LaHP', 'LpHP', 'LCA1', 'LCA23DG', 'LSub', ...
                      'RaHP', 'RpHP', 'RCA1', 'RCA23DG', 'RSub', 'RCA1-A', ...
                      'BaHP', 'BpHP', 'BCA1', 'BCA23DG', 'BSub', 'BCA1-A'};

    % target ROI (names2)
    exclude_name = {'Cerebral', 'Vent', 'Matter', 'Brain-Stem', 'CSF', 'choroid', 'hypointensities', 'Optic-Chiasm', 'CC_', 'QC', 'HPC'};
    exclude_idx = cellfun(@(x) find(contains(results.names2, x)), exclude_name, 'uni', 0);
    exclude_idx = cell2mat(exclude_idx);
    exclude_idx = [exclude_idx, roi_seed_idx];

    roi_target_idx = setdiff(1:length(results.names2), exclude_idx);
    roi_target_names = results.names2(roi_target_idx);

    %%% Z value (transformed from Pearson's r correlation)
    Z_surf{con_i} = results.Z(roi_seed_idx, roi_target_idx, :);
end
roi_target_names = cellfun(@(x) x(9:end), roi_target_names, 'uni', 0);

tmp_label = {'Left-', '-lh-'}; roi_target_idx_left = cell2mat(cellfun(@(x) find(contains(roi_target_names, x)), tmp_label, 'uni', 0));
tmp_label = {'Right-', '-rh-'}; roi_target_idx_right = cell2mat(cellfun(@(x) find(contains(roi_target_names, x)), tmp_label, 'uni', 0));

% connectivity - bilateral seed
avg_roi_name = ['B', roi_seed_names{1}(2:end)];
roi_seed_names = {roi_seed_names{:}, avg_roi_name};
for con_i = 1:length(Z_surf)
    Z_surf{con_i}(end+1, :, :) = mean(Z_surf{con_i}(1:2, :, :), 1);
end


%% CA23DG connectivity: correlation (non-abused)

%%%%%%%%%%%%%%%%
con_i = 3; % 3: encoding
seed_i = 28; % 28: B.CA23DG, 27: B.CA1, 29: B.Subiculum, 36: B.Hippocampus

group_idx = (eti_ea > 0);
%%%%%%%%%%%%%%%%

% ---------- stats (correlation) ---------- %
stats_on = true;
if stats_on
    fprintf('<< CA23DG connectivity - correlation >>\n');
    fprintf(' - non-abused -\n');

    sig_idx = [];
    for roi_i = 1:length(roi_target_names)
        data = squeeze(Z_surf{con_i}( seed_i, roi_i, : ));

        % - age -
        [r1,p1] = corr(sbj_age(group_idx == 0)', data(group_idx == 0));
    
        % - full em -
        [r2,p2] = corr(em_acc{6}(group_idx == 0)', data(group_idx == 0));

        % - print results -
        if p1 < 0.05 || p2 < 0.05
            sig_idx = [sig_idx roi_i];

            fprintf('%30s: age) r=%7.3f, p=%7.3f', roi_target_names{roi_i}, r1, p1);
            fprintf('\tEM) r=%7.3f, p=%7.3f\n', r2, p2);
        end
    end
    disp(' ');

    % - check hemsipheric interactions -
    if true
        sig_idx_sort = [];
        for roi_i = sig_idx
            tmp_l = find(roi_target_idx_left == roi_i);
            tmp_r = find(roi_target_idx_right == roi_i);
            tmp_b = find(roi_target_idx_avg == roi_i);
    
            sig_idx_sort = [sig_idx_sort tmp_l tmp_r tmp_b];
        end
        sig_idx_sort = unique(sig_idx_sort);
    
        idx_age = []; idx_h_age = []; idx_acc = []; idx_h_acc = []; 
        tbl_h_age = []; tbl_h_acc = []; count = 0;
        for roi_i = sig_idx_sort
            count = count+1;
            tmp_conn_l = squeeze(Z_surf{con_i}( seed_i, roi_target_idx_left(roi_i), : ));
            tmp_conn_r = squeeze(Z_surf{con_i}( seed_i, roi_target_idx_right(roi_i), : ));
            tmp_conn = [tmp_conn_l(group_idx == 0)' tmp_conn_r(group_idx == 0)'];

            tmp_hemis = [ones(1, sum(group_idx == 0)) ones(1, sum(group_idx == 0))*2];

            tmp_eti = eti_ea(group_idx == 0); tmp_eti = [tmp_eti tmp_eti];
            tmp_age = sbj_age(group_idx == 0); tmp_age = [tmp_age tmp_age];
            tmp_acc = em_acc{6}(group_idx == 0); tmp_acc = [tmp_acc tmp_acc]; % full EM

            tmp_sbj = 1:sum(group_idx==0); tmp_sbj = [tmp_sbj tmp_sbj];

            %
            tmp_names = {'conn', 'hemis', 'age', 'acc', 'eti', 'sbj'};
            tmp_cell = {tmp_conn', tmp_hemis', tmp_age', tmp_acc', tmp_eti', tmp_sbj'};
            tmp_array = cell2mat(tmp_cell);
            tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
            tmp_table.hemis = categorical(tmp_table.hemis);

            % - hemisphere by age interaction -
            lme = fitlme(tmp_table, 'conn~hemis*age+(hemis|sbj)');
            if lme.Coefficients{3,6} < 0.05
                idx_age = [idx_age roi_i];
                disp(roi_target_names{roi_target_idx_avg(roi_i)}(5:end));
                anova(lme)
            end
            if lme.Coefficients{4,6} < 0.05
                idx_h_age = [idx_h_age roi_i];
            end
            tbl_h_age{count} = lme.Coefficients;
            
            % - hemisphere by acc interaction -
            lme = fitlme(tmp_table, 'conn~hemis*acc+(hemis|sbj)');
            if lme.Coefficients{3,6} < 0.05
                idx_acc = [idx_acc roi_i];
                disp(roi_target_names{roi_target_idx_avg(roi_i)}(5:end));
                anova(lme)
            end
            if lme.Coefficients{4,6} < 0.05
                idx_h_acc = [idx_h_acc roi_i];
            end
            tbl_h_acc{count} = lme.Coefficients;
        end

        % - print results -
        fprintf(' -> hemispheric interaction?\n');
        fprintf('  # main effects of age\n');
        for i = 1:length(idx_age)
            curr_idx = find(sig_idx_sort == idx_age(i));
            curr_tbl = tbl_h_age{curr_idx};
            fprintf('%30s: age) β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', roi_target_names{roi_target_idx_avg(idx_age(i))}(5:end), ...
                                            curr_tbl{3,2}, curr_tbl{3,3}, curr_tbl{3,5}, curr_tbl{3,4}, curr_tbl{3,6})
        end
        disp(' ');

        fprintf('  # interaction effects between age and hemisphere\n');
        for i = 1:length(idx_h_age)
            curr_idx = find(sig_idx_sort == idx_h_age(i));
            curr_tbl = tbl_h_age{curr_idx};
            fprintf('%30s: age) β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', roi_target_names{roi_target_idx_avg(idx_h_age(i))}(5:end), ...
                                            curr_tbl{3,2}, curr_tbl{3,3}, curr_tbl{3,5}, curr_tbl{3,4}, curr_tbl{3,6})
        end
        disp(' ');

        fprintf('  # main effects of acc\n');
        for i = 1:length(idx_acc)
            curr_idx = find(sig_idx_sort == idx_acc(i));
            curr_tbl = tbl_h_acc{curr_idx};
            fprintf('%30s: acc) β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', roi_target_names{roi_target_idx_avg(idx_acc(i))}(5:end), ...
                                            curr_tbl{3,2}, curr_tbl{3,3}, curr_tbl{3,5}, curr_tbl{3,4}, curr_tbl{3,6})
        end
        disp(' ');

        fprintf('  # interaction effects between acc and hemisphere\n');
        for i = 1:length(idx_h_acc)
            curr_idx = find(sig_idx_sort == idx_h_acc(i));
            curr_tbl = tbl_h_acc{curr_idx};
            fprintf('%30s: acc) β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', roi_target_names{roi_target_idx_avg(idx_h_acc(i))}(5:end), ...
                                            curr_tbl{3,2}, curr_tbl{3,3}, curr_tbl{3,5}, curr_tbl{3,4}, curr_tbl{3,6})
        end
        disp(' ');
    end
end


%% CA23DG connectivity: correlation (abused)

%%%%%%%%%%%%%%%%
con_i = 3; % 3: encoding
seed_i = 28; % 28: B.CA23DG, 27: B.CA1, 29: B.Subiculum, 36: B.Hippocampus

group_idx = (eti_ea > 0);
%%%%%%%%%%%%%%%%

% ---------- stats (correlation) ---------- %
stats_on = true;
% stats_on = false;
if stats_on
    fprintf('<< CA23DG connectivity - correlation >>\n');
    fprintf(' - abused -\n');

    sig_idx = [];
    for roi_i = 1:length(roi_target_names)
        data = squeeze(Z_surf{con_i}( seed_i, roi_i, : ));

        % - age -
        [r1,p1] = corr(sbj_age(group_idx == 1)', data(group_idx == 1));
    
        % - full em -
        [r2,p2] = corr(em_acc{6}(group_idx == 1)', data(group_idx == 1));

        % - print results -
        if p1 < 0.05 || p2 < 0.05
            sig_idx = [sig_idx roi_i];

            fprintf('%30s: age) r=%7.3f, p=%7.3f', roi_target_names{roi_i}, r1, p1);
            fprintf('\tEM) r=%7.3f, p=%7.3f\n', r2, p2);
        end
    end
    disp(' ');

    % - check hemsipheric interactions -
    if true
        sig_idx_sort = [];
        for roi_i = sig_idx
            tmp_l = find(roi_target_idx_left == roi_i);
            tmp_r = find(roi_target_idx_right == roi_i);
            tmp_b = find(roi_target_idx_avg == roi_i);
    
            sig_idx_sort = [sig_idx_sort tmp_l tmp_r tmp_b];
        end
        sig_idx_sort = unique(sig_idx_sort);
    
        idx_age = []; idx_h_age = []; idx_acc = []; idx_h_acc = []; 
        tbl_h_age = []; tbl_h_acc = []; count = 0;
        for roi_i = sig_idx_sort
            count = count+1;
            tmp_conn_l = squeeze(Z_surf{con_i}( seed_i, roi_target_idx_left(roi_i), : ));
            tmp_conn_r = squeeze(Z_surf{con_i}( seed_i, roi_target_idx_right(roi_i), : ));
            tmp_conn = [tmp_conn_l(group_idx == 1)' tmp_conn_r(group_idx == 1)'];

            tmp_hemis = [ones(1, sum(group_idx == 1)) ones(1, sum(group_idx == 1))*2];

            tmp_eti = eti_ea(group_idx == 1); tmp_eti = [tmp_eti tmp_eti];
            tmp_age = sbj_age(group_idx == 1); tmp_age = [tmp_age tmp_age];
            tmp_acc = em_acc{6}(group_idx == 1); tmp_acc = [tmp_acc tmp_acc]; % full EM

            tmp_sbj = 1:sum(group_idx == 1); tmp_sbj = [tmp_sbj tmp_sbj];

            %
            tmp_names = {'conn', 'hemis', 'age', 'acc', 'eti', 'sbj'};
            tmp_cell = {tmp_conn', tmp_hemis', tmp_age', tmp_acc', tmp_eti', tmp_sbj'};
            tmp_array = cell2mat(tmp_cell);
            tmp_table = array2table(tmp_array, 'VariableNames', tmp_names);
            tmp_table.hemis = categorical(tmp_table.hemis);

            % - hemisphere by age interaction -
            lme = fitlme(tmp_table, 'conn~hemis*age+(hemis|sbj)');
            if lme.Coefficients{3,6} < 0.05
                idx_age = [idx_age roi_i];
            end
            if lme.Coefficients{4,6} < 0.05
                idx_h_age = [idx_h_age roi_i];
            end
            tbl_h_age{count} = lme.Coefficients;
            
            % - hemisphere by acc interaction -
            lme = fitlme(tmp_table, 'conn~hemis*acc+(hemis|sbj)');
            if lme.Coefficients{3,6} < 0.05
                idx_acc = [idx_acc roi_i];
            end
            if lme.Coefficients{4,6} < 0.05
                idx_h_acc = [idx_h_acc roi_i];
            end
            tbl_h_acc{count} = lme.Coefficients;
        end

        % - print results -
        fprintf(' -> hemispheric interaction?\n');
        fprintf('  # main effects of age\n');
        for i = 1:length(idx_age)
            curr_idx = find(sig_idx_sort == idx_age(i));
            curr_tbl = tbl_h_age{curr_idx};
            fprintf('%30s: age) β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', roi_target_names{roi_target_idx_avg(idx_age(i))}(5:end), ...
                                            curr_tbl{3,2}, curr_tbl{3,3}, curr_tbl{3,5}, curr_tbl{3,4}, curr_tbl{3,6})
        end
        disp(' ');

        fprintf('  # interaction effects between age and hemisphere\n');
        for i = 1:length(idx_h_age)
            curr_idx = find(sig_idx_sort == idx_h_age(i));
            curr_tbl = tbl_h_age{curr_idx};
            fprintf('%30s: age) β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', roi_target_names{roi_target_idx_avg(idx_h_age(i))}(5:end), ...
                                            curr_tbl{4,2}, curr_tbl{4,3}, curr_tbl{4,5}, curr_tbl{4,4}, curr_tbl{4,6})
        end
        disp(' ');

        fprintf('  # main effects of acc\n');
        for i = 1:length(idx_acc)
            curr_idx = find(sig_idx_sort == idx_acc(i));
            curr_tbl = tbl_h_acc{curr_idx};
            fprintf('%30s: acc) β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', roi_target_names{roi_target_idx_avg(idx_acc(i))}(5:end), ...
                                            curr_tbl{3,2}, curr_tbl{3,3}, curr_tbl{3,5}, curr_tbl{3,4}, curr_tbl{3,6})
        end
        disp(' ');

        fprintf('  # interaction effects between acc and hemisphere\n');
        for i = 1:length(idx_h_acc)
            curr_idx = find(sig_idx_sort == idx_h_acc(i));
            curr_tbl = tbl_h_acc{curr_idx};
            fprintf('%30s: acc) β=%.3f, SE=%.3f, t(%d)=%.3f, p=%.3f\n', roi_target_names{roi_target_idx_avg(idx_h_acc(i))}(5:end), ...
                                            curr_tbl{4,2}, curr_tbl{4,3}, curr_tbl{4,5}, curr_tbl{4,4}, curr_tbl{4,6})
        end
        disp(' ');
    end

    % - interaction with abuse -
    fprintf(' - interaction with abuse -\n');
    for roi_i = sig_idx
        data = squeeze(Z_surf{con_i}( seed_i, roi_i, : ));

        % -- full EM --
        [~, tbl, stats, ~] = anovan(data, {eti_ea, em_acc{6}}, 'continuous', [1 2], 'varnames', {'eti', 'acc'}, 'model', 'full', 'display', 'off');

        % -- display results --
        fprintf('%30s: full EM ) F=%7.3f, p=%7.3f\n', roi_target_names{roi_i}, tbl{4,6}, tbl{4,7});
    end
end

