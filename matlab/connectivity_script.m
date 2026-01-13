% Episodic memory task for adolescents
% :: seed-based functional connectivity script (CONN)
% :: code written by Kahyun Choi


%% define parameters

in_fmri_path = 'C:\fMRI_directory'; % directory of fMRI data (nifti)
in_free_path = 'C:\freesurfer_directory'; % directory of FreeSurfer output data
regressor_path = 'C:\regressor_directory'; % directory of regressor files


% load subject list
subject_path = 'C:\subject_info_directory';
load(fullfile(subject_path, 'subject_list_for_CONN.mat')); % 'bhv_list', 'sbj_list', 'sbj_nums', 'rejected_run_sbj', 'rejected_run_list'

reject_run_sbj = cellfun(@(x) str2double(x(2:end))-100, rejected_run_sbj);
reject_run = rejected_run_list;


% file list
t1_list = cellfun(@(x) [in_free_path x '/mri/T1.mgz'], sbj_list, 'uni', 0);
func_list = cell(1, num_sbj);
for sbj_i = 1:num_sbj
    if ~isempty(find(reject_run_sbj == sbj_nums(sbj_i), 1))
        run_list = setdiff(1:4, reject_run{find(reject_run_sbj == sbj_nums(sbj_i), 1)});
    else
        run_list = 1:4;
    end
    func_list{sbj_i} = cell(1, length(run_list));
    for run_i = 1:length(run_list)
        tmp = [in_fmri_path sbj_list{sbj_i} '/' num2str(run_list(run_i)) '/f*.nii'];
        func_list{sbj_i}{run_i} = cellstr(conn_dir(tmp));
    end
end
STRUCTURAL_FILE = t1_list;
FUNCTIONAL_FILE = func_list;

nsessions = cellfun(@length, func_list);
TR = 2;


% ROI mask files
roi_path = cellfun(@(x) [in_free_path x '/label/lh.aparc.annot'], sbj_list, 'uni', 0);
roi_name = 'FS_atlas';

roi_path2 = cellfun(@(x) [in_free_path x '/mri/aparc+aseg.mgz'], sbj_list, 'uni', 0);
roi_name2 = 'FS_aseg';

roi_name3 = {'LaHPC', 'LpHPC', 'Lbody', 'Ltail', 'LCA1', 'LCA23DG', 'LSub', ...
             'RaHPC', 'RpHPC', 'Rbody', 'Rtail', 'RCA1', 'RCA23DG', 'RSub', ...
             'BaHPC', 'BpHPC', 'Bbody', 'Btail', 'BCA1', 'BCA23DG', 'BSub'};
roi_path3 = cell(1, length(roi_name3));
for roi_i = 1:length(roi_name3)
    mask_path = 'C:\freesurfer_mask_directory';
    roi_path3{roi_i} = cellfun(@(x) [mask_path, x, '/', roi_name3{roi_i}, '.nii'], sbj_list, 'uni', 0);
end


% condition
enc_onsets = cell(1, length(sbj_list)); enc_durations = cell(1, length(sbj_list));
for sbj_i = 1:length(sbj_list)
    curr_sbj = sbj_nums(sbj_i);
    if ~isempty(find(reject_run_sbj == sbj_nums(sbj_i), 1))
        run_num = 4 - length(reject_run{reject_run_sbj == sbj_nums(sbj_i)});
        run_idx = setdiff([1 2 3 4], reject_run{reject_run_sbj == sbj_nums(sbj_i)});
    else
        run_num = 4;
        run_idx = [1 2 3 4];
    end
    count = 0;
%     for run_i = 1:run_num
    for run_i = run_idx
        count = count + 1;
        
        task_reg = fullfile(regressor_path, sprintf('%s_run%d.mat', sbj_list{sbj_i}, run_i));
        reg = load(task_reg);
        
        fixenc_i = 2; enc_i = 3;
        enc_onsets{sbj_i}{count} = reg.onsets{fixenc_i};
        enc_durations{sbj_i}{count} = reg.onsets{enc_i} - reg.onsets{fixenc_i} + reg.durations{enc_i} + 6; % hemodynamic delay
    end
end


%% SETUP

global CONN_gui; CONN_gui.usehighres = true; % to fix "usehighres" error

clear batch;

batch.filename = fullfile('../results/CONN', sprintf('conn_N%d_indiv.mat', num_sbj));
batch.parallel.N = 3;

batch.Setup.nsubjects = num_sbj;
batch.Setup.RT = TR;

batch.Setup.functionals = repmat({{}}, [num_sbj, 1]);
for sbj_i = 1:num_sbj
    for run_i = 1:nsessions(sbj_i)
        batch.Setup.functionals{sbj_i}{run_i} = FUNCTIONAL_FILE{sbj_i}{run_i};
    end
end
batch.Setup.structurals = STRUCTURAL_FILE;

for sbj_i = 1:num_sbj
    for run_i = 1:nsessions(sbj_i)
        batch.Setup.conditions.names={'rest', 'retrieval', 'encoding'};
        batch.Setup.conditions.onsets{1}{sbj_i}{run_i} = 0;
        batch.Setup.conditions.durations{1}{sbj_i}{run_i} = inf;
        batch.Setup.conditions.onsets{2}{sbj_i}{run_i} = enc_onsets{sbj_i}{run_i};
        batch.Setup.conditions.durations{2}{sbj_i}{run_i} = enc_durations{sbj_i}{run_i};
    end
end

batch.Setup.isnew = 1; % new file
batch.Setup.done = 0; % save without running

conn_batch(batch);


%% PREPROCESSING

clear batch;

batch.filename = fullfile('../results/CONN', sprintf('conn_N%d_indiv.mat', num_sbj));
batch.parallel.N = 3; % 3

conn_importaseg;

batch.Setup.isnew = 0; % existing file
batch.Setup.analyses = [1 2 4]; % except voxel-to-voxel
batch.Setup.voxelresolution = 4;

batch.Setup.preprocessing.steps = 'default_ss';
batch.Setup.preprocessing.sliceorder = 'interleaved (bottom-up)'; % 1-3-5-...-2-4-...
batch.Setup.preprocessing.diffusionsteps = 40; % enter number of diffusion steps for smoothing: 40

batch.Setup.done = 1;
batch.Setup.overwrite = 'Yes';

conn_batch(batch);


%% ROI - subject-space connectivity: aparc+aseg.mgz (subject-space)
clear batch;

batch.filename = fullfile('../results/CONN', sprintf('conn_N%d_indiv.mat', num_sbj));

batch.Setup.rois.names = {roi_name2};
batch.Setup.rois.files = {roi_path2};
batch.Setup.rois.dataset = 7; % subject-space (7)
batch.Setup.rois.add = 0; % use 1 to define an additional set of ROIs (to be added to any already-existing ROIs in your project)

batch.Setup.done = 0; % save without running
conn_batch(batch);


%% ROI - subject-space connectivity: freesurfer mask (subject-space)
clear batch;

batch.filename = fullfile('../results/CONN', sprintf('conn_N%d_indiv.mat', num_sbj));
batch.parallel.N = 4;

batch.Setup.rois.names = roi_name3;
batch.Setup.rois.files = roi_path3;
batch.Setup.rois.dataset = ones(1,length(roi_name3))*7; % subject-space
batch.Setup.rois.add = 1; % use 1 to define an additional set of ROIs (to be added to any already-existing ROIs in your project)

batch.Setup.done = 1; % save with running
conn_batch(batch);


%% CONN display: atlas ROI setup
% ROI 1, 2) -> advanced options (Atlas file V / Subject-specific ROI V)

conn
conn('load', fullfile('../results/CONN', sprintf('conn_N%d_indiv.mat', num_sbj)));


%% DENOISING & FIRST-LEVEL ANALYSIS

clear batch;

batch.filename = fullfile('../results/CONN', sprintf('conn_N%d_indiv.mat', num_sbj));
batch.parallel.N = 4;

% DENOISING
batch.Denoising.filter = [0.008 Inf]; % high-pass (Hz)
batch.Denoising.done = 1;
batch.Denoising.overwrite = 'Yes';


% FIRST-LEVEL ANALYSIS
batch.Analysis.name = 'SBC';
batch.Analysis.modulation = 0; % 0: standard weighted GLM
batch.Analysis.conditions = [];

batch.Analysis.done = 1;
batch.Analysis.overwrite = 'Yes';

% RUN
conn_batch(batch);


%% CONN display: for 2nd level analysis (GUI)

conn
conn('load', fullfile('../results/CONN', sprintf('conn_N%d_indiv.mat', num_sbj)));

conn gui_results


