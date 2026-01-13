% Episodic memory task for adolescents
% :: behavior - analysis (demographics)
% :: code written by Kahyun Choi


%% define parameters

% load subject list
subject_path = 'C:\subject_info_directory';
bhv_nums = load(fullfile(subject_path, 'subject_list_for_GLM.mat')); % 'bhv_list', 'sbj_list', 'sbj_nums', 'num_sbj'

% behavior
load(fullfile(subject_path, 'behavioral_data.mat')); % 'em_acc', 'sbj_age', 'sbj_sex'
load(fullfile(subject_path, 'survey_data.mat')); % 'sbj_svy', 'svy_label'

% ETI
dir_survey = fullfile(subject_path, 'ETI_survey.xlsx');
sbj_eti = readtable(dir_survey);
sbj_eti = sbj_eti(sbj_nums, :);
sbj_eti = table2struct(sbj_eti); % ALL, GT, PA, EA, Q1~21
eti_ea = [sbj_eti.EA];


%% Participant demographics (Non-abused, Abused)

tmp_group = (eti_ea > 0);

sbj_ses = sbj_svy{ismember(svy_label, 'SES')};
sbj_iq = sbj_svy{ismember(svy_label, 'IQ')};
sbj_stais = sbj_svy{ismember(svy_label, 'STAI_S')};
sbj_cesd = sbj_svy{ismember(svy_label, 'CES_D')};

tmp_list = {sbj_age, sbj_sex, sbj_ses, sbj_iq, sbj_stais, sbj_cesd};
tmp_name = {'age', 'sex', 'SES', 'IQ', 'STAI-S', 'CES-D'};

disp('<< Participant demographics >>');

% age
[~,p,~,stats] = ttest2(tmp_list{1}(tmp_group==0), tmp_list{1}(tmp_group==1));
fprintf('%8s: non - %.2f(%.2f), abused - %.2f(%.2f), T=%.3f, p=%.3f\n', tmp_name{1}, ...
    nanmean(tmp_list{1}(tmp_group == 0)), nanstd(tmp_list{1}(tmp_group == 0)), nanmean(tmp_list{1}(tmp_group == 1)), nanstd(tmp_list{1}(tmp_group == 1)), ...
    stats.tstat, p);

% sex
[~,p,~,stats] = ttest2(tmp_list{2}(tmp_group==0), tmp_list{2}(tmp_group==1));
fprintf('%8s: non - %.2f(%.2f), abused - %.2f(%.2f), T=%.3f, p=%.3f\n', tmp_name{2}, ...
    sum(tmp_list{2}(tmp_group==0)==1), sum(tmp_list{2}(tmp_group==0)==1)/sum(tmp_group==0), sum(tmp_list{2}(tmp_group==1)==1), sum(tmp_list{2}(tmp_group==1)==1)/sum(tmp_group==1), ...
    stats.tstat, p);

% svy
for i = 3:6
    [~,p,~,stats] = ttest2(tmp_list{i}(tmp_group==0), tmp_list{i}(tmp_group==1));
    fprintf('%8s: non - %.2f(%.2f), abused - %.2f(%.2f), T=%.3f, p=%.3f\n', tmp_name{i}, ...
        nanmean(tmp_list{i}(tmp_group == 0)), nanstd(tmp_list{i}(tmp_group == 0)), nanmean(tmp_list{i}(tmp_group == 1)), nanstd(tmp_list{i}(tmp_group == 1)), ...
        stats.tstat, p);
end
disp(' ');


%% Participant demographics - family environment

tmp_group = (eti_ea > 0);

disp('<< Participant demographics: family environment >>');
% svy
for i = 34:42
    [~,p,~,stats] = ttest2(sbj_svy{i}(tmp_group==0), sbj_svy{i}(tmp_group==1));
    fprintf('%12s: non - %.2f(%.2f), abused - %.2f(%.2f), T=%.3f, p=%.3f\n', svy_label{i}, ...
        nanmean(sbj_svy{i}(tmp_group == 0)), nanstd(sbj_svy{i}(tmp_group == 0)), nanmean(sbj_svy{i}(tmp_group == 1)), nanstd(sbj_svy{i}(tmp_group == 1)), ...
        stats.tstat, p);
end
disp(' ');


%% Correlation table

sbj_iq = sbj_svy{ismember(svy_label, 'IQ')};
sbj_stais = sbj_svy{ismember(svy_label, 'STAI_S')};
sbj_cesd = sbj_svy{ismember(svy_label, 'CES_D')};

sbj_data = {sbj_age', sbj_stais, sbj_cesd, sbj_iq, em_acc{6}', em_acc{5}', em_acc{4}'};
sbj_data = cell2mat(sbj_data);
tbl_label = {'age', 'stais', 'cesd', 'iq', 'fullem', 'wherew', 'whatw'};
sbj_table = array2table(sbj_data, 'VariableNames', tbl_label);

%
[r, p] = corrcoef(sbj_data, 'rows', 'pairwise');
idx_tril = tril(r); r(idx_tril == 0) = NaN; p(idx_tril == 0) = NaN;
idx_ones = (r==1); r(idx_ones) = NaN; p(idx_ones) = NaN;

% - table -
disp('<< correlation table >>');
disp(' - r values - ');
disp(round(r, 3))
disp(' - p values - ');
disp(p)

% - STAIS & CESD - 
fprintf('%8s & %8s: r=%.3f, p=%.3f\n', tbl_label{2}, tbl_label{5}, r(5,2), p(5,2));
fprintf('%8s & %8s: r=%.3f, p=%.3f\n', tbl_label{3}, tbl_label{5}, r(5,3), p(5,3));
fprintf('%8s & %8s: r=%.3f, p=%.3f\n', tbl_label{2}, tbl_label{6}, r(6,2), p(6,2));
fprintf('%8s & %8s: r=%.3f, p=%.3f\n', tbl_label{3}, tbl_label{6}, r(6,3), p(6,3));
fprintf('%8s & %8s: r=%.3f, p=%.3f\n', tbl_label{2}, tbl_label{7}, r(7,2), p(7,2));
fprintf('%8s & %8s: r=%.3f, p=%.3f\n', tbl_label{3}, tbl_label{7}, r(7,3), p(7,3));


