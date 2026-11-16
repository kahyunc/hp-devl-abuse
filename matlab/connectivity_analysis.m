% Episodic memory task for adolescents
% :: seed-based functional connectivity - second level (SPM12)
% :: code written by Kahyun Choi


%% define parameters

WORKING_DIRECTORY = pwd;

IN_PATH_FIRST = 'C:\connectivity_output_directory'; % directory of connectivity output 

% load subject list
subject_path = 'C:\subject_info_directory';
load(fullfile(subject_path, 'subject_list_for_GLM.mat')); % 'bhv_list', 'sbj_list', 'sbj_nums', 'num_sbj'

% setup
OUT_PATH = fullfile(IN_PATH_FIRST, 'connectivity_2nd_level');
if ~exist(OUT_PATH, 'dir'); mkdir(OUT_PATH); end

% contrasts
condition_i = 2; % 1: rest, 2: encoding
cond_name = {'rest', 'encoding'};

% ROIs
roi_file_name = sprintf('resultsROI_Condition%.3d.mat', condition_i);
load(fullfile(IN_PATH_FIRST, roi_file_name), 'names');

seed_names = {'B_HPC', 'B_aHPC', 'B_pHPC'};
roi_seed_idx = cellfun(@(x) find(strcmp(names, x), 1), seed_names);
roi_seed_names = names(roi_seed_idx);
roi_seed_names = cellfun(@(x) strrep(x, '_', '-'), roi_seed_names, 'uni', 0);

% behavior
load(fullfile(subject_path, 'behavioral_data.mat')); % 'em_acc', 'sbj_age', 'sbj_sex'

% ETI
dir_survey = fullfile(subject_path, 'ETI_survey.xlsx');
sbj_eti = readtable(dir_survey);
sbj_eti = sbj_eti(sbj_nums, :);
sbj_eti = table2struct(sbj_eti); % ALL, GT, PA, EA, Q1~21
eti_ea = [sbj_eti.EA];


%% run 2nd level: condition contrast - group difference 

roi_idx = 1:length(roi_seed_idx);

target_group = {find(eti_ea == 0), find(eti_ea > 0)}; 
target_group_name = 'eti';

group_idx = (1:length(bhv_list)); group_name = '';

group_2nd_name = ['group_N', num2str(length(bhv_list))];
GROUP_OUT_PATH = ['../results/connectivity/', group_2nd_name];

for roi_i = roi_idx
    cd(WORKING_DIRECTORY)

    out_dir = fullfile(GROUP_OUT_PATH, cond_name{condition_i}, sprintf('%s%s', target_group_name, group_name), roi_seed_names{roi_i});
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    file_names1 = arrayfun(@(x) sprintf('BETA_Subject%.3d_Condition%.3d_Source%.3d.nii', x, condition_i, roi_seed_idx(roi_i)), target_group{1}, 'uni', 0);
    in_dir1 = cellfun(@(x) fullfile(IN_PATH_FIRST, x), file_names1, 'uni', 0);
    file_names2 = arrayfun(@(x) sprintf('BETA_Subject%.3d_Condition%.3d_Source%.3d.nii', x, condition_i, roi_seed_idx(roi_i)), target_group{2}, 'uni', 0);
    in_dir2 = cellfun(@(x) fullfile(IN_PATH_FIRST, x), file_names2, 'uni', 0);


    clear matlabbatch

    matlabbatch{1}.spm.stats.factorial_design.dir = {out_dir};

    matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = in_dir1(:);
    matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = in_dir2(:);

    matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {'C:\mask_directory\mask_ICV.nii'};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = '-';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'abused > non-abused';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 0;

    matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    result_contrasts = [1 2];
    for result_i = 1:length(result_contrasts)
        matlabbatch{4}.spm.stats.results.conspec(result_i).titlestr = '';
        matlabbatch{4}.spm.stats.results.conspec(result_i).contrasts = result_contrasts(result_i);
        matlabbatch{4}.spm.stats.results.conspec(result_i).threshdesc = 'none';
        matlabbatch{4}.spm.stats.results.conspec(result_i).thresh = 0.001;
        matlabbatch{4}.spm.stats.results.conspec(result_i).extent = 15;
        matlabbatch{4}.spm.stats.results.conspec(result_i).conjunction = 1;
        matlabbatch{4}.spm.stats.results.conspec(result_i).mask.none = 1;
    end
    matlabbatch{4}.spm.stats.results.units = 1;
    matlabbatch{4}.spm.stats.results.export{1}.ps = true;
    matlabbatch{4}.spm.stats.results.export{2}.tspm.basename = 'unc001_k15';


    %% run batch
    batch = matlabbatch;

    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_jobman('run',batch);

    cd(WORKING_DIRECTORY)
end


%% run 2nd level: multiple regression (correlation)

%%%%%%%%%%%%%%%%
roi_idx = 1:length(roi_seed_idx);

target_bhv = sbj_age; target_names = {'age'};
% target_bhv = em_acc{6}; target_names = {'fullem'};
% target_bhv = em_acc{5}; target_names = {'wherewhen'};
% target_bhv = em_acc{4}; target_names = {'whatwhen'};

% group_idx = 1:length(bhv_list); group_name = '';
group_idx = find(eti_ea == 0); group_name = 'no-EA_';
% group_idx = find(eti_ea > 0); group_name = 'EA_';
%%%%%%%%%%%%%%%%

correlation_name = sprintf('corr_N%d', length(bhv_list));
CORR_OUT_PATH = ['../results/connectivity/', correlation_name];
if ~exist(CORR_OUT_PATH, 'dir'); mkdir(CORR_OUT_PATH); end

for roi_i = roi_idx
    cd(WORKING_DIRECTORY)

    file_names = arrayfun(@(x) sprintf('BETA_Subject%.3d_Condition%.3d_Source%.3d.nii', x, condition_i, roi_seed_idx(roi_i)), group_idx, 'uni', 0);
    in_dir = cellfun(@(x) fullfile(IN_PATH_FIRST, x), file_names, 'uni', 0);
    cov = target_bhv(:, group_idx);
    clear matlabbatch

    corr_sign = 'both';
    out_dir = fullfile(CORR_OUT_PATH, cond_name{condition_i}, sprintf('%s%s', group_name, target_names{1}), roi_seed_names{roi_i});
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end



    n_cov = size(cov, 1);

    matlabbatch{1}.spm.stats.factorial_design.dir = {out_dir};
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = in_dir(:);
    for i = 1:n_cov
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).c = cov(i, :);
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).cname = target_names{i};
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).iCC = 1;
    end
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {'C:\mask_directory\mask_ICV.nii'};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    if strcmp(corr_sign, 'all')
        contrast_idx = 1:n_cov;
    else
        contrast_idx = 1;
    end

    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

    con_idx = 0;
    for cov_i = contrast_idx
        col_idx = cov_i + 1;
        con_weights = zeros(1, n_cov + 1);
        if strcmp(corr_sign, 'neg')
            con_idx = con_idx + 1;
            con_weights(col_idx) = -1;
            matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.name = [target_names{cov_i}, '_neg'];
            matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.weights = con_weights;
            matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';
        else
            con_idx = con_idx + 1;
            con_weights(col_idx) = 1;
            matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.name = [target_names{cov_i}, '_pos'];
            matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.weights = con_weights;
            matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';
            if strcmp(corr_sign, 'both') || strcmp(corr_sign, 'all')
                con_idx = con_idx + 1;
                con_weights(col_idx) = -1;
                matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.name = [target_names{cov_i}, '_neg'];
                matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.weights = con_weights;
                matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';
            end
        end
    end
    matlabbatch{3}.spm.stats.con.delete = 0;

    matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    result_contrasts = 1:con_idx;
    for result_i = 1:length(result_contrasts)
        matlabbatch{4}.spm.stats.results.conspec(result_i).titlestr = '';
        matlabbatch{4}.spm.stats.results.conspec(result_i).contrasts = result_contrasts(result_i);
        matlabbatch{4}.spm.stats.results.conspec(result_i).threshdesc = 'none';
        matlabbatch{4}.spm.stats.results.conspec(result_i).thresh = 0.001;
        matlabbatch{4}.spm.stats.results.conspec(result_i).extent = 15;
        matlabbatch{4}.spm.stats.results.conspec(result_i).conjunction = 1;
        matlabbatch{4}.spm.stats.results.conspec(result_i).mask.none = 1;
    end
    matlabbatch{4}.spm.stats.results.units = 1;
    matlabbatch{4}.spm.stats.results.export{1}.ps = true;
    matlabbatch{4}.spm.stats.results.export{2}.tspm.basename = 'unc001_k15';


    %% run batch
    batch = matlabbatch;

    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_jobman('run',batch);

    cd(WORKING_DIRECTORY)
end


%% run 2nd level: condition contrast - ANOVA (2-way interaction)

roi_idx = 1:length(roi_seed_idx);

var1 = double(eti_ea == 0);  var1_name = 'eti';

% var2 = sbj_age;            var2_name = 'age';
var2 = em_acc{6};            var2_name = 'fullem';
% var2 = em_acc{5};            var2_name = 'wherew';
% var2 = em_acc{4};            var2_name = 'whatw';

group_idx = (1:length(bhv_list)); group_name = '';

anova_name = ['anova_N', num2str(length(bhv_list))];
ANOVA_OUT_PATH = ['../results/connectivity/', anova_name];
if ~exist(ANOVA_OUT_PATH, 'dir'); mkdir(ANOVA_OUT_PATH); end

for roi_i = roi_idx
    cd(WORKING_DIRECTORY)

    file_names = arrayfun(@(x) sprintf('BETA_Subject%.3d_Condition%.3d_Source%.3d.nii', x, condition_i, roi_seed_idx(roi_i)), group_idx, 'uni', 0);
    in_dir = cellfun(@(x) fullfile(IN_PATH_FIRST, x), file_names, 'uni', 0);

    grp_v1 = var1(group_idx);
    grp_v2 = var2(group_idx);

    main_cov  = [grp_v1; grp_v2];
    cov_names = {var1_name, var2_name, [var1_name '#' var2_name]};
    interact_cov = grp_v1 .* grp_v2;

    out_dir = fullfile(ANOVA_OUT_PATH, cond_name{condition_i}, sprintf('interaction_%s%s#%s', group_name, var1_name, var2_name), roi_seed_names{roi_i});
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end


    clear matlabbatch
    n_main = size(main_cov, 1);
    n_interact = size(interact_cov, 1);
    n_cov = n_main + n_interact;

    all_cov = [main_cov; interact_cov];
    if ~iscell(cov_names)
        cov_names = {cov_names};
    end

    matlabbatch{1}.spm.stats.factorial_design.dir = {out_dir};
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = in_dir(:);
    for i = 1:n_cov
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).c = all_cov(i, :);
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).cname = cov_names{i};
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).iCC = 0;
    end
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {'C:\mask_directory\mask_ICV.nii'};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val','{}',{2},'.','val','{}',{1},'.','val','{}',{1}), substruct('.','spmmat'));

    con_idx = 0;

    f_weights = [zeros(n_interact, n_main+1), eye(n_interact)];
    interact_name = strjoin(cov_names(n_main+1:end), '+');

    con_idx = con_idx + 1;
    matlabbatch{3}.spm.stats.con.consess{con_idx}.fcon.name = ['F_' interact_name];
    matlabbatch{3}.spm.stats.con.consess{con_idx}.fcon.weights = f_weights;
    matlabbatch{3}.spm.stats.con.consess{con_idx}.fcon.sessrep = 'none';

    for interact_i = 1:n_interact
        col_idx = n_main + 1 + interact_i;
        con_weights = zeros(1, n_cov+1);

        con_idx = con_idx + 1;
        con_weights(col_idx) = 1;
        matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.name = [cov_names{n_main+interact_i}, '_pos'];
        matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.weights = con_weights;
        matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';

        con_idx = con_idx + 1;
        con_weights(col_idx) = -1;
        matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.name = [cov_names{n_main+interact_i}, '_neg'];
        matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.weights = con_weights;
        matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';
    end
    matlabbatch{3}.spm.stats.con.delete = 0;

    matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val','{}',{3},'.','val','{}',{1},'.','val','{}',{1}), substruct('.','spmmat'));
    result_contrasts = 2:con_idx;
    for result_i = 1:length(result_contrasts)
        matlabbatch{4}.spm.stats.results.conspec(result_i).titlestr = '';
        matlabbatch{4}.spm.stats.results.conspec(result_i).contrasts = result_contrasts(result_i);
        matlabbatch{4}.spm.stats.results.conspec(result_i).threshdesc = 'none';
        matlabbatch{4}.spm.stats.results.conspec(result_i).thresh = 0.001;
        matlabbatch{4}.spm.stats.results.conspec(result_i).extent = 15;
        matlabbatch{4}.spm.stats.results.conspec(result_i).conjunction = 1;
        matlabbatch{4}.spm.stats.results.conspec(result_i).mask.none = 1;
    end
    matlabbatch{4}.spm.stats.results.units = 1;
    matlabbatch{4}.spm.stats.results.export{1}.ps = true;
    matlabbatch{4}.spm.stats.results.export{2}.tspm.basename = 'unc001_k15';


    %% run batch
    batch = matlabbatch;

    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_jobman('run',batch);

    cd(WORKING_DIRECTORY)
end

%% connectivity window sensitivity: original 25-s clusters

window_paths = {'C:\connectivity_25s_output', 'C:\connectivity_18s_output', 'C:\connectivity_19s_output'};
condition_idx = [3 1 1]; % original CONN projects: 25s / 18s / 19s
window_subdir = {'HPC', '', ''};
seed_names = {'B_aHPC', 'B_pHPC'};
DIR_CLUSTER = 'C:\original_25s_cluster_masks'; % separate seed folders
window_seed_idx = cell(1, 3);
for window_i = 1:3
    tmp = load(fullfile(window_paths{window_i}, ...
        sprintf('resultsROI_Condition%.3d.mat', condition_idx(window_i))), 'names');
    window_seed_idx{window_i} = cellfun(@(x) find(strcmp(tmp.names, x), 1), seed_names);
end

for seed_i = 1:length(seed_names)
    mask_list = dir(fullfile(DIR_CLUSTER, seed_names{seed_i}, '*.nii'));
    for mask_i = 1:length(mask_list)
        mask = niftiread(fullfile(mask_list(mask_i).folder, mask_list(mask_i).name)) > 0;
        for window_i = 1:3
            fc = nan(num_sbj, 1);
            for sbj_i = 1:num_sbj
                beta_file = fullfile(window_paths{window_i}, window_subdir{window_i}, ...
                    sprintf('BETA_Subject%.3d_Condition%.3d_Source%.3d.nii', ...
                    sbj_i, condition_idx(window_i), window_seed_idx{window_i}(seed_i)));
                beta = niftiread(beta_file);
                fc(sbj_i) = mean(beta(mask));
            end
            tmp_table = table(fc, double(eti_ea(:) == 0), em_acc{6}(:), ...
                'VariableNames', {'fc', 'group', 'acc'});
            mdl = fitlm(tmp_table, 'fc~group*acc');
            fprintf('\n%s, %s, window %d\n', seed_names{seed_i}, mask_list(mask_i).name, window_i);
            disp(mdl.Coefficients('group:acc', :));
        end
    end
end
