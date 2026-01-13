% Episodic memory task for adolescents
% :: GLM first level (SPM12)
% :: code written by Kahyun Choi


%% define parameters

WORKING_DIRECTORY = pwd;

IN_PATH_FUNC = 'C:\fMRI_preprocessed_directory'; % directory of preprocessed fMRI data (nifti)
IN_PATH_REG = 'C:\regressor_directory'; % directory of regressor files

IN_PATH_CONTRAST = 'C:\contrast_weight_directory';
IN_PATH_MOVE = 'C:\movement_info_directory';

d = dir(IN_PATH_FUNC); d(1:2) = [];
sbj_list = {d.name};
sbj_list = sbj_list(contains(sbj_list, 'A'));
sbj_nums = cellfun(@(x) str2double(x(2:4))-100, sbj_list);

DIR_FUNC_LIST = {'M1', 'M2', 'M3', 'M4'};


% load subject list
subject_path = 'C:\subject_info_directory';
load(fullfile(subject_path, 'rejected_subject_list_GLM.mat')); % 'rejected_all_sbj', 'rejected_run_sbj', 'rejected_run_list'

sbj_reject_all = cellfun(@(x) str2double(x(2:end))-100, rejected_all_sbj);
sbj_reject_run = cellfun(@(x) str2double(x(2:end))-100, rejected_run_sbj);
reject_run = rejected_run_list;

sbj_idx = setdiff(sbj_nums, sbj_reject_all);
n_subjects = length(sbj_idx);


% GLM file setting
glm_name = sprintf('1st_level_N%d_CASE%s_txt', n_subjects, CASE_NUM);

OUT_PATH = ['../results/glm_norm_smooth/', glm_name]; % standard space
target_reg_exp = '.*swraf.*';
% OUT_PATH = ['../results/glm_indiv_rough/', glm_name]; % individual space (without normalization)
% target_reg_exp = '\<raf.*';

if ~exist(OUT_PATH, 'dir'); mkdir(OUT_PATH); end


% contrast weight setting
con_names = {'enc1', 'ret1', 'enc-fix1', 'ret-fix1', 'enc-ctrl1', 'ret-ctrl1', 'enc-ctrl-fix1', 'ret-ctrl-fix1', 'reenact1', 'interfere1', ...
             'enc2', 'ret2', 'enc-fix2', 'ret-fix2', 'enc-ctrl2', 'ret-ctrl2', 'enc-ctrl-fix2', 'ret-ctrl-fix2', 'reenact2', 'interfere2', ...
             'ctrl-enc', 'ctrl-ret'};

xls_weights = xlsread(fullfile(IN_PATH_CONTRAST, 'contrast_weight.xlsx'));
xls_weights(isnan(xls_weights)) = 0;

con_weights = arrayfun(@(num) xls_weights(num, :), 1:size(xls_weights, 1), 'uni', 0);
con_weights = cellfun(@(x) [x 0 0 0 0 0 0], con_weights, 'uni',0); % movements + rotation

file_name_qc = fullfile(IN_PATH_MOVE, 'movement_info.mat');
load(file_name_qc); % 'move_all', 'rot_all'


%% run GLM: first-level

error_sbj = []; error_run = [];
for sbj_i = 1:length(sbj_nums)

    % check if rejected sbj
    if ~isempty(find(sbj_reject_all == sbj_nums(sbj_i), 1))
        continue;
    end
    
    for run_i = 1:4   
        
        % check if rejected run
        if ~isempty(find(sbj_reject_run == sbj_nums(sbj_i), 1))
            if ~isempty(find(reject_run{find(sbj_reject_run == sbj_nums(sbj_i), 1)} == run_i, 1))
                continue;
            end
        end
        
        try
            %% setup
            cd(WORKING_DIRECTORY)

            out_path = fullfile(OUT_PATH, sbj_list{sbj_i}, DIR_FUNC_LIST{run_i});
            reg_out_path = fullfile(OUT_PATH, sbj_list{sbj_i}, 'reg');
            if ~exist(out_path, 'dir'); mkdir(out_path); end
            if ~exist(reg_out_path, 'dir'); mkdir(reg_out_path); end

            in_path_func = fullfile(IN_PATH_FUNC, sbj_list{sbj_i}, DIR_FUNC_LIST{run_i});
            file_list = dir(in_path_func);
            file_list = arrayfun(@(x) x.name, file_list, 'uni', 0);
            
            target_idx = cellfun(@(x) regexp(x, target_reg_exp), file_list, 'uni', 0);
            target_idx = cellfun(@isempty, target_idx);
            file_list = file_list(~target_idx);
            
            file_list = cellfun(@(x) fullfile(in_path_func, x), file_list, 'uni', 0);
            
            % movement regressors
            [m_pars] = spm_select('List', in_path_func, '^rp.*\.txt$');
            cont_reg = cellstr([in_path_func,'/', m_pars]);

            % task regressors
            task_reg = fullfile(IN_PATH_REG, sprintf('%s_run%d_trial.mat', sbj_list{sbj_i}, run_i));


            %% batch
            if ~iscell(out_path)
                out_dir = {out_path};
            end
            func = {file_list};
            reg_mat = {task_reg};
            reg_txt = cont_reg;
            p = 0.05;

            file_name_mask = 'C:\mask_directory\mask_ICV.nii';
            
            %
            matlabbatch{1}.spm.stats.fmri_spec.dir = out_dir;
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
            
            %
            if length(func) ~= 1
                for run_i = 1:length(func)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(run_i).scans = func{run_i};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(run_i).multi = {reg_mat{run_i}};
                    if ~isempty(reg_txt)
                        matlabbatch{1}.spm.stats.fmri_spec.sess(run_i).multi_reg = {reg_txt{run_i}};
                    else
                        matlabbatch{1}.spm.stats.fmri_spec.sess(run_i).multi_reg = {''};
                    end
                    matlabbatch{1}.spm.stats.fmri_spec.sess(run_i).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess(run_i).regress = struct('name', {}, 'val', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess(run_i).hpf = 128;
                end
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess.scans = func{1};
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi = reg_mat;
                if ~isempty(reg_txt)
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = reg_txt;
                else
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
                end
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {file_name_mask};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
            
            %
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
            
            %
            matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            for con_i = 1:length(con_names)
                matlabbatch{3}.spm.stats.con.consess{con_i}.tcon.name = con_names{con_i};
                matlabbatch{3}.spm.stats.con.consess{con_i}.tcon.weights = con_weights{con_i};
                matlabbatch{3}.spm.stats.con.consess{con_i}.tcon.sessrep = 'none';
            end
            matlabbatch{3}.spm.stats.con.delete = 0;
            
            %
            matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            for con_i = 1:length(con_names)
                matlabbatch{4}.spm.stats.results.conspec(con_i).titlestr = '';
                matlabbatch{4}.spm.stats.results.conspec(con_i).contrasts = con_i;
                matlabbatch{4}.spm.stats.results.conspec(con_i).threshdesc = 'FWE';
                matlabbatch{4}.spm.stats.results.conspec(con_i).thresh = p;
                matlabbatch{4}.spm.stats.results.conspec(con_i).extent = 0;
                matlabbatch{4}.spm.stats.results.conspec(con_i).conjunction = 1;
                matlabbatch{4}.spm.stats.results.conspec(con_i).mask.none = 1;
            end
            
            matlabbatch{4}.spm.stats.results.units = 1;
            matlabbatch{4}.spm.stats.results.export{1}.png = true;
            
            %% run batch
            batch = matlabbatch;

            spm('defaults','fmri');
            spm_jobman('initcfg');
            spm_jobman('run',batch);

            cd(WORKING_DIRECTORY)

        catch
            cd(WORKING_DIRECTORY)
            error_sbj = [error_sbj sbj_i]; error_run = [error_run run_i];
        end
    end
end

% error check
for i = 1:length(error_sbj)
    fprintf('\n1st level error - sbj%d, run%d\n\n', error_sbj(i), error_run(i));
end

