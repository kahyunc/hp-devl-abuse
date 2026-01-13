% Episodic memory task for adolescents
% :: GLM second level (SPM12)
% :: code written by Kahyun Choi


%% define parameters

WORKING_DIRECTORY = pwd;

IN_PATH_FIRST = 'C:\GLM_first_level_directory'; % directory of GLM first level output data

d = dir(IN_PATH_FIRST); d_names = {d.name};
sbj_list = d_names(contains(d_names, 'A'));
sbj_nums = cellfun(@(x) str2double(x(2:end))-100, sbj_list);


% load subject list
subject_path = 'C:\subject_info_directory';
bhv_nums = load(fullfile(subject_path, 'subject_list_for_GLM.mat')); % 'bhv_list', 'sbj_list', 'sbj_nums', 'num_sbj'

excluded_sbj = setdiff(sbj_nums, bhv_nums.sbj_nums);
included_sbj = intersect(sbj_nums, bhv_nums.sbj_nums);

excluded_idx = cell2mat(arrayfun(@(num) find(sbj_nums == excluded_sbj(num)), 1:length(excluded_sbj), 'uni', 0));
sbj_nums(excluded_idx) = []; sbj_list(excluded_idx) = [];


% setup
glm_2nd_name = sprintf('2nd_level_N%d', length(included_sbj));
OUT_PATH = ['../results/glm_norm_smooth/', glm_2nd_name];
if ~exist(OUT_PATH, 'dir'); mkdir(OUT_PATH); end


% contrasts
con_names = {'enc','ret', 'enc-fix','ret-fix', 'enc-ctrl','ret-ctrl', 'enc-ctrl-fix','ret-ctrl-fix', 'reenact','interf'};
con_list = 1:length(con_names);
run_list = {'M1', 'M2', 'M3', 'M4'};


% behavior
load(fullfile(subject_path, 'subject_list_for_GLM.mat'), 'bhv_list'); 
load(fullfile(subject_path, 'behavioral_data.mat')); % 'em_acc', 'sbj_age', 'sbj_sex'


% ETI
dir_survey = fullfile(subject_path, 'ETI_survey.xlsx');
sbj_eti = readtable(dir_survey);
sbj_eti = sbj_eti(sbj_nums, :);
sbj_eti = table2struct(sbj_eti); % ALL, GT, PA, EA, Q1~21
eti_ea = [sbj_eti.EA];


%% run 2nd level: condition contrast

%%%%%%%%%%%%%%%%
target_con = {'enc-ctrl'};

% group_idx = (1:length(included_sbj)); group_name = '';
group_idx = (eti_ea == 0); group_name = '_ETIX';
% group_idx = (eti_ea > 0); group_name = '_ETIO';
%%%%%%%%%%%%%%%%

excluded_idx = cell2mat(arrayfun(@(num) find(sbj_nums == excluded_sbj(num)), 1:length(excluded_sbj), 'uni', 0));
sbj_nums(excluded_idx) = [];
sbj_list(excluded_idx) = [];

for con_i = 1:length(target_con)
    cd(WORKING_DIRECTORY) 
    
    tmp_sbj = sbj_list(group_idx);
    
    out_dir = fullfile(OUT_PATH, group_name, target_con{con_i});
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    in_path = fullfile(IN_PATH, '_organized');
    file_list = cellfun(@(x) fullfile(in_path, x, sprintf('%s.nii', target_con{con_i})), tmp_sbj, 'uni', 0);

    % check if empty
    empty_idx = cellfun(@(x) exist(x, 'file'), file_list);
    file_list(empty_idx == 0) = [];
    in_dir = file_list;

    
    %% batch
    matlabbatch{1}.spm.stats.factorial_design.dir = {out_dir};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = in_dir(:);
    if isempty(covars)
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    else
        for c_i = 1:length(covars)
            matlabbatch{1}.spm.stats.factorial_design.cov(c_i).c = covars{c_i}{2};
            matlabbatch{1}.spm.stats.factorial_design.cov(c_i).cname = covars{c_i}{1};
            matlabbatch{1}.spm.stats.factorial_design.cov(c_i).iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(c_i).iCC = 1;
        end
    end
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
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
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
    matlabbatch{4}.spm.stats.results.export{2}.tspm.basename = 'FWE';

    %% run batch
    batch = matlabbatch;

    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_jobman('run',batch);

    cd(WORKING_DIRECTORY)
end


%% run 2nd level: condition contrast - group difference 

%%%%%%%%%%%%%%%%
target_con = {'enc-ctrl'};
target_group = {(eti_ea == 0), (eti_ea > 0)}; 
target_group_name = 'eti';

group_idx = (1:length(included_sbj)); group_name = '';
%%%%%%%%%%%%%%%%

group_2nd_name = ['group_N', num2str(length(included_sbj))];
GROUP_OUT_PATH = ['../results/glm_norm_smooth/', group_2nd_name];

for con_i = 1:length(target_con)
    cd(WORKING_DIRECTORY)

    in_path = fullfile(IN_PATH, '_organized');

    out_dir = fullfile(GROUP_OUT_PATH, sprintf('%s%s', target_group_name, group_name), target_con{con_i});
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    in_dir1 = cellfun(@(x) fullfile(in_path, x, sprintf('%s.nii', target_con{con_i})), sbj_list(target_group{1} & group_idx), 'uni', 0);
    in_dir2 = cellfun(@(x) fullfile(in_path, x, sprintf('%s.nii', target_con{con_i})), sbj_list(target_group{2} & group_idx), 'uni', 0);

    %% batch
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


%% run 2nd level: condition contrast - ANOVA (2-way interaction)

%%%%%%%%%%%%%%%%
target_con = {'enc-ctrl'};

var1 = double(eti_ea == 0);  var1_name = 'eti';
var2 = em_acc{4};            var2_name = 'whatw';

group_idx = (1:length(included_sbj)); group_name = '';
%%%%%%%%%%%%%%%%

var_list_all = {var1, var2, var1.*var2};

cov_name = {var1_name, var2_name, [var1_name '#' var2_name]};
cov = cell2mat(var_list_all');

correlation_name = ['anova_N', num2str(length(included_sbj)), '_case', CASE_NUM];
CORR_OUT_PATH = ['../results/glm_norm_smooth/', correlation_name];
if ~exist(CORR_OUT_PATH, 'dir'); mkdir(CORR_OUT_PATH); end

for con_i = 1:length(target_con)
    cd(WORKING_DIRECTORY)

    tmp_sbj = sbj_list(group_idx);
    cov = cov(:, group_idx);

    in_dir = cellfun(@(x) fullfile(IN_PATH, '_organized', x, sprintf('%s.nii', target_con{con_i})), tmp_sbj, 'uni', 0);
    out_dir = fullfile(CORR_OUT_PATH, target_con{con_i}, sprintf('interaction_%s%s#%s', group_name, var1_name, var2_name));
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    corr_sign = 'full';

    %% batch
    weight_matrix = zeros(1, size(cov,1)+1);
    
    if ~iscell(cov_name)
        cov_name = {cov_name};
    end
    
    matlabbatch{1}.spm.stats.factorial_design.dir = {out_dir};
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = in_dir(:);
    if size(cov,1) == 1
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = cov;
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = cov_name{1};
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1;
    else
        for i = 1:size(cov,1)
            matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).c = cov(i, :);
            matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).cname = cov_name{i};
            matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).iCC = 1;
        end
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
    
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = '-';
    if strcmp(c_sign, 'pos')
        tmp_weights = weight_matrix; tmp_weights(2) = 1;
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = tmp_weights;
    elseif strcmp(c_sign, 'neg')
        tmp_weights = weight_matrix; tmp_weights(2) = -1;
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = tmp_weights;
    elseif strcmp(c_sign, 'full')
        tmp_weights = weight_matrix; tmp_weights(end) = 1;
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = tmp_weights;
    end
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
    matlabbatch{4}.spm.stats.results.export{2}.tspm.basename = 'FWE';
        
    %% run batch
    batch = matlabbatch;

    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_jobman('run',batch);

    cd(WORKING_DIRECTORY)
end



%% Post-hoc: significant clusters from ANOVA

%%%%%%%%%%%%%%%%
plot_color = [70 111 156; 160 185 212]./255;  % non-abused vs. abused - whatwhen

group_idx = (eti_ea > 0); group_names = {'non', 'abused'}; group_type = 'ETI EA';

beta_name = 'enc-ctrl';

data_behav = em_acc{4};
%%%%%%%%%%%%%%%%

DIR_BETA = 'C:\GLM_first_level_directory';

mask_path = 'C:\GLM_second_level_ANOVA_output_directory';
mask_names = {'1_L_SMA', '2_R_SMA', '3_R_ITG', '4_L_IN', '5_L_CAU', '6_R_CAU', '7_L_MFG'};

fprintf('<< cluster correlation >>\n');
for i = 1:length(mask_names)

    file_list = cellfun(@(x) fullfile(DIR_BETA, x, sprintf('%s.nii', beta_name)), sbj_list, 'uni', 0);
    beta = cellfun(@niftiread, file_list, 'uni', 0);
    hdr = niftiinfo(file_list{1});

    roi_path = fullfile(mask_path, sprintf('cluster_%s.nii', mask_names{i}));

    roi = spm_read_vols(spm_vol(roi_path), 1);
    roi_idx = find(roi > 0);
    [roi_x, roi_y, roi_z] = ind2sub(size(roi), roi_idx);
    roi_xyz = [roi_x roi_y roi_z]';
    data = cellfun(@(x) nanmean(spm_get_data(x, roi_xyz), 2), file_list);

    %%% - group correlation 
    figure;
    set(gcf, 'Position', [100 100 280 120]);

    subplot(1,11,[1 4.5]);
    [r, p] = corr(reshape(data_behav(group_idx == 0),[],1), reshape(data(group_idx == 0),[],1));
    fig_scatter = scatter(data_behav(group_idx == 0), data(group_idx == 0), 6, 'filled');
    fig_scatter.MarkerFaceColor = plot_color(1, :);
    fig_line = lsline;
    fig_line.LineWidth = 1.2;

    mdl = fitlm(reshape(data_behav(group_idx == 0),[],1), reshape(data(group_idx == 0),[],1));
    x = linspace(min(reshape(data_behav(group_idx == 0),[],1)), max(reshape(data_behav(group_idx == 0),[],1)), 1000);
    [y,ci] = predict(mdl,x');
    hold on
    fill([x,flip(x)]', [[ci(:,1)]', flip(ci(:,2))']', [.5 .5 .5],'FaceAlpha',0.2,'linestyle','none'); 

    if p<0.05; fig_line.Color=[1 0.5 0.5]; end
    xlim([0 1.05]); xticks(0:0.2:1);
    ylabel('beta coefficient'); yticks(-4:2:4); ylim([-4 4.5]);
    set(gca,'XTickLabelRotation',0);
    set(gca,'FontName','Helvetica','FontSize',6, 'FontWeight','bold');
    set(gca,'LineWidth',0.8);
    fprintf('%8s: non - r=%.3f, p=%.3f', mask_names{i}, r, p);

    subplot(1,11,[7.5 11]);
    [r, p] = corr(reshape(data_behav(group_idx == 1),[],1), reshape(data(group_idx == 1),[],1));
    fig_scatter = scatter(data_behav(group_idx == 1), data(group_idx == 1), 6, 'filled');
    fig_scatter.MarkerFaceColor = plot_color(2, :);
    fig_line = lsline;
    fig_line.LineWidth = 1.2;

    mdl = fitlm(reshape(data_behav(group_idx == 1),[],1), reshape(data(group_idx == 1),[],1));
    x = linspace(min(reshape(data_behav(group_idx == 1),[],1)), max(reshape(data_behav(group_idx == 1),[],1)), 1000);
    [y,ci] = predict(mdl,x');
    hold on
    fill([x,flip(x)]', [[ci(:,1)]', flip(ci(:,2))']', [.5 .5 .5],'FaceAlpha',0.2,'linestyle','none'); 

    if p<0.05; fig_line.Color=[1 0.5 0.5]; end
    xlim([0 1.05]); xticks(0:0.2:1);
    ylabel('beta coefficient'); yticks(-4:2:4); ylim([-4 4.5]);
    set(gca,'XTickLabelRotation',0);
    set(gca,'FontName','Helvetica','FontSize',6, 'FontWeight','bold');
    set(gca,'LineWidth',0.8);
    fprintf('\tabused - r=%.3f, p=%.3f\n', r, p);
end
disp(' ');

