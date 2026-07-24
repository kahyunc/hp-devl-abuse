% Episodic memory task for adolescents
% :: seed-based functional connectivity - second level (SPM12)
% :: code written by Kahyun Choi


%% define parameters

WORKING_DIRECTORY = pwd;

IN_PATH_FIRST = 'C:\connectivity_output_directory'; % directory of connectivity output 

d = dir(IN_PATH_FIRST); d_names = {d.name};
sbj_list = d_names(contains(d_names, 'A'));
sbj_nums = cellfun(@(x) str2double(x(2:end))-100, sbj_list);


% load subject list
subject_path = 'C:\subject_info_directory';
bhv_nums = load(fullfile(subject_path, 'subject_list_for_GLM.mat')); % 'bhv_list', 'sbj_list', 'sbj_nums', 'num_sbj'


% setup
OUT_PATH = fullfile(IN_PATH_FIRST, 'connectivity_2nd_level');
if ~exist(OUT_PATH, 'dir'); mkdir(OUT_PATH); end


% contrasts
condition_i = 2; % 1: rest, 2: encoding
cond_name = {'rest', 'encoding'};


% ROIs
roi_file_name = sprintf('resultsROI_Condition%.3d.mat', condition_i);
load(fullfile(IN_PATH_FIRST, roi_file_name), 'names');

roi_seed_idx = [165 166 167 168 169 170 171 172 173]; % hippocampus
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

%%%%%%%%%%%%%%%%
roi_idx = 1:length(roi_seed_idx);

target_group = {find(eti_ea == 0), find(eti_ea > 0)}; 
target_group_name = 'eti';

group_idx = (1:length(bhv_list)); group_name = '';
%%%%%%%%%%%%%%%%

group_2nd_name = ['group_N', num2str(length(bhv_list))];
GROUP_OUT_PATH = ['../results/glm_norm_smooth/', group_2nd_name];

for roi_i = roi_idx
    cd(WORKING_DIRECTORY)

    out_dir = fullfile(GROUP_OUT_PATH, cond_name{condition_i}, sprintf('%s%s', target_group_name, group_name), roi_seed_names{roi_i});
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    file_names1 = arrayfun(@(x) sprintf('BETA_Subject%.3d_Condition%.3d_Source%.3d.nii', x, condition_i, roi_seed_idx(roi_i)), target_group{1}, 'uni', 0);
    in_dir1 = cellfun(@(x) fullfile(IN_PATH_FIRST, 'HPC', x), file_names1, 'uni', 0);
    file_names2 = arrayfun(@(x) sprintf('BETA_Subject%.3d_Condition%.3d_Source%.3d.nii', x, condition_i, roi_seed_idx(roi_i)), target_group{2}, 'uni', 0);
    in_dir2 = cellfun(@(x) fullfile(IN_PATH_FIRST, 'HPC', x), file_names2, 'uni', 0);

    %% batch
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
    matlabbatch{3}.spm.stats.con.delete = 0;
    
    matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
    matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
    matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
    matlabbatch{4}.spm.stats.results.conspec.extent = 0;
    matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{4}.spm.stats.results.units = 1;
    matlabbatch{4}.spm.stats.results.export{1}.ps = true;
    %
    matlabbatch{4}.spm.stats.results.export{2}.tspm.basename = 'FWE';
    
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
group_idx = find(eti_ea == 0); group_name = 'ETIX_';
% group_idx = find(eti_ea > 0); group_name = 'ETIO_';
%%%%%%%%%%%%%%%%

correlation_name = sprintf('corr_N%d', length(bhv_list));
CORR_OUT_PATH = ['../results/glm_norm_smooth/', correlation_name];
if ~exist(CORR_OUT_PATH, 'dir'); mkdir(CORR_OUT_PATH); end

for roi_i = roi_idx
    cd(WORKING_DIRECTORY)

    file_names = arrayfun(@(x) sprintf('BETA_Subject%.3d_Condition%.3d_Source%.3d.nii', x, condition_i, roi_seed_idx(roi_i)), group_idx, 'uni', 0);
    in_dir = cellfun(@(x) fullfile(IN_PATH_FIRST, 'HPC', x), file_names, 'uni', 0);
    cov = target_bhv(:, group_idx);

    corr_sign = 'both';
    out_dir = fullfile(CORR_OUT_PATH, cond_name{condition_i}, sprintf('%s%s', group_name, target_names), roi_seed_names{roi_i});
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end


    %% batch
    n_cov = size(cov, 1);

    matlabbatch{1}.spm.stats.factorial_design.dir = {out_dir};
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = in_dir(:);
    for i = 1:n_cov
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).c     = cov(i, :);
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).cname = target_names{i};
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).iCC   = 1; % overall mean centering
    end
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint   = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov       = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none  = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im          = 0;
    matlabbatch{1}.spm.stats.factorial_design.masking.em          = {'C:\mask_directory\mask_ICV.nii'};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit      = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm     = 1;

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', ...
        substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    if strcmp(corr_sign, 'all')
        contrast_idx = 1:n_cov;  % generate contrasts for all regressors
    else
        contrast_idx = 1;         % only variable of interest (first covariate)
    end
    
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', ...
        substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    
    con_idx = 0;
    for ci = contrast_idx
        col = ci + 1;  % +1 for intercept
        w = zeros(1, n_cov + 1);
        if strcmp(corr_sign, 'neg')
            con_idx = con_idx + 1;
            w(col) = -1;
            matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.name    = [target_names{ci}, '_neg'];
            matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.weights = w;
            matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';
        else  % 'pos', 'both', 'all'
            con_idx = con_idx + 1;
            w(col) = 1;
            matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.name    = [target_names{ci}, '_pos'];
            matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.weights = w;
            matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';
            if strcmp(corr_sign, 'both') || strcmp(corr_sign, 'all')
                con_idx = con_idx + 1;
                w(col) = -1;
                matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.name    = [target_names{ci}, '_neg'];
                matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.weights = w;
                matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';
            end
        end
    end
    matlabbatch{3}.spm.stats.con.delete = 0;

    matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', ...
        substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{4}.spm.stats.results.conspec.titlestr    = '';
    matlabbatch{4}.spm.stats.results.conspec.contrasts   = 1;
    matlabbatch{4}.spm.stats.results.conspec.threshdesc  = 'FWE';
    matlabbatch{4}.spm.stats.results.conspec.thresh      = 0.05;
    matlabbatch{4}.spm.stats.results.conspec.extent      = 0;
    matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{4}.spm.stats.results.conspec.mask.none   = 1;
    matlabbatch{4}.spm.stats.results.units               = 1;
    matlabbatch{4}.spm.stats.results.export{1}.ps        = true;
    matlabbatch{4}.spm.stats.results.export{2}.tspm.basename = 'FWE';
    
    batch = matlabbatch;
        
    %% run batch
    batch = matlabbatch;

    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_jobman('run',batch);

    cd(WORKING_DIRECTORY)
end




%% run 2nd level: condition contrast - ANOVA (2-way interaction)

%%%%%%%%%%%%%%%%
roi_idx = 1:length(roi_seed_idx);

var1 = double(eti_ea == 0);  var1_name = 'eti';

% var2 = sbj_age;            var2_name = 'age';
var2 = em_acc{6};            var2_name = 'fullem';
% var2 = em_acc{5};            var2_name = 'wherew';
% var2 = em_acc{4};            var2_name = 'whatw';

group_idx = (1:length(bhv_list)); group_name = '';
%%%%%%%%%%%%%%%%

anova_name = ['anova_N', num2str(length(bhv_list))];
ANOVA_OUT_PATH = ['../results/glm_norm_smooth/', anova_name];
if ~exist(ANOVA_OUT_PATH, 'dir'); mkdir(ANOVA_OUT_PATH); end

for roi_i = roi_idx
    cd(WORKING_DIRECTORY)

    file_names = arrayfun(@(x) sprintf('BETA_Subject%.3d_Condition%.3d_Source%.3d.nii', x, condition_i, roi_seed_idx(roi_i)), group_idx, 'uni', 0);
    in_dir = cellfun(@(x) fullfile(IN_PATH_FIRST, 'HPC', x), file_names, 'uni', 0);

    grp_v1 = var1(group_idx);
    grp_v2 = var2(group_idx);

    main_cov  = [grp_v1; grp_v2];
    cov_names = {var1_name, var2_name, [var1_name '#' var2_name]};
    interact_cov = grp_v1 .* grp_v2;

    out_dir = fullfile(ANOVA_OUT_PATH, cond_name{condition_i}, sprintf('interaction_%s%s#%s', group_name, var1_name, var2_name), roi_seed_names{roi_i});
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end


    %% batch
    n_main     = size(main_cov,     1);
    n_interact = size(interact_cov, 1);
    n_cov      = n_main + n_interact;
    
    all_cov = [main_cov; interact_cov];    
    if ~iscell(cov_names)
        cov_names = {cov_names};
    end
    
    matlabbatch{1}.spm.stats.factorial_design.dir = {out_dir};
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = in_dir(:);
    for i = 1:n_cov
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).c     = all_cov(i, :);
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).cname = cov_names{i};
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).iCC   = 0;
    end
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint   = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov       = struct('c',{}, 'cname',{}, 'iCFI',{}, 'iCC',{});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files',{}, 'iCFI',{}, 'iCC',{});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none  = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im          = 0;
    matlabbatch{1}.spm.stats.factorial_design.masking.em          = {'C:\mask_directory\mask_ICV.nii'};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit      = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm     = 1;

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', ...
        substruct('.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    % Design columns: intercept(1), main(2..n_main+1), interact(n_main+2..n_cov+1)
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', ...
        substruct('.','val','{}',{2},'.','val','{}',{1},'.','val','{}',{1}), substruct('.','spmmat'));
    
    con_idx = 0;
    
    % F-contrast: jointly tests all interaction terms (n_interact df)
    % rows = n_interact, cols = n_cov+1 (intercept + all covariates)
    f_weights = [zeros(n_interact, n_main+1), eye(n_interact)];
    interact_name = strjoin(cov_names(n_main+1:end), '+');
    
    con_idx = con_idx + 1;
    matlabbatch{3}.spm.stats.con.consess{con_idx}.fcon.name    = ['F_' interact_name];
    matlabbatch{3}.spm.stats.con.consess{con_idx}.fcon.weights = f_weights;
    matlabbatch{3}.spm.stats.con.consess{con_idx}.fcon.sessrep = 'none';
    
    % t-contrasts: pos and neg for each interaction term individually
    for ii = 1:n_interact
        col = n_main + 1 + ii;   % column in design matrix (+1 for intercept)
        w   = zeros(1, n_cov+1);
    
        con_idx  = con_idx + 1;
        w(col)   = 1;
        matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.name    = [cov_names{n_main+ii}, '_pos'];
        matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.weights = w;
        matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';
    
        con_idx  = con_idx + 1;
        w(col)   = -1;
        matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.name    = [cov_names{n_main+ii}, '_neg'];
        matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.weights = w;
        matlabbatch{3}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';
    end
    matlabbatch{3}.spm.stats.con.delete = 0;

    matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', ...
        substruct('.','val','{}',{3},'.','val','{}',{1},'.','val','{}',{1}), substruct('.','spmmat'));
    matlabbatch{4}.spm.stats.results.conspec.titlestr    = '';
    matlabbatch{4}.spm.stats.results.conspec.contrasts   = 1;
    matlabbatch{4}.spm.stats.results.conspec.threshdesc  = 'FWE';
    matlabbatch{4}.spm.stats.results.conspec.thresh      = 0.05;
    matlabbatch{4}.spm.stats.results.conspec.extent      = 0;
    matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{4}.spm.stats.results.conspec.mask.none   = 1;
    matlabbatch{4}.spm.stats.results.units               = 1;
    matlabbatch{4}.spm.stats.results.export{1}.ps        = true;
    matlabbatch{4}.spm.stats.results.export{2}.tspm.basename = 'FWE';

        
    %% run batch
    batch = matlabbatch;

    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_jobman('run',batch);

    cd(WORKING_DIRECTORY)
end
