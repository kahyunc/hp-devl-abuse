% Episodic memory task for adolescents
% :: seed-based functional connectivity script (CONN)
% :: code written by Kahyun Choi


%% define parameters

IN_PATH_FMRI = 'C:\fMRI_directory'; % directory of fMRI data (nifti)
IN_PATH_REG = 'C:\regressor_directory'; % directory of regressor files


% load subject list
subject_path = 'C:\subject_info_directory';
load(fullfile(subject_path, 'subject_list_for_GLM.mat')); % 'bhv_list', 'sbj_list', 'sbj_nums', 'rejected_run_sbj', 'rejected_run_list'

reject_run_sbj = cellfun(@(x) str2double(x(2:end))-100, rejected_run_sbj);
reject_run = rejected_run_list;


% file list
t1_list = cellfun(@(x) fullfile(IN_PATH_FMRI, x, 'T1'), sbj_list, 'uni', 0);
func_list = cell(1, num_sbj);
for sbj_i = 1:num_sbj
    if ~isempty(find(reject_run_sbj == sbj_nums(sbj_i), 1))
        run_list = setdiff(1:4, reject_run{find(reject_run_sbj == sbj_nums(sbj_i), 1)});
    else
        run_list = 1:4;
    end
    func_list{sbj_i} = cell(1, length(run_list));
    for run_i = 1:length(run_list)
        tmp = fullfile(IN_PATH_FMRI, sbj_list{sbj_i}, sprintf('M%d', run_list(run_i)), 'f*.nii');
        func_list{sbj_i}{run_i} = cellstr(conn_dir(tmp));
    end
end
STRUCTURAL_FILE = cellstr(cellfun(@(x) conn_dir([x '/s*.nii']), t1_list, 'uni', 0));
FUNCTIONAL_FILE = func_list;

nsessions = cellfun(@length, func_list);
TR = 2;


% ROI mask files
dir_mask = 'C:\mask_directory';
mask_names = {'B_HPC', 'B_aHPC', 'B_pHPC'}; % AAL; AP split at y = -22 mm
mask_paths = cellfun(@(x) fullfile(dir_mask, ['r' x '.nii']), mask_names, 'uni', 0);

% condition
WINDOW_SECONDS = 25; % 25: original; 18: no fixation; 19: no post-encoding extension
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
    for run_i = run_idx
        count = count + 1;
        
        task_reg = fullfile(IN_PATH_REG, sprintf('%s_run%d.mat', sbj_list{sbj_i}, run_i));
        reg = load(task_reg);
        
        fixenc_i = 2; enc_i = 3;
        if WINDOW_SECONDS == 18
            enc_onsets{sbj_i}{count} = reg.onsets{enc_i};
            enc_durations{sbj_i}{count} = reg.durations{enc_i} + 6;
        elseif WINDOW_SECONDS == 19
            enc_onsets{sbj_i}{count} = reg.onsets{fixenc_i};
            enc_durations{sbj_i}{count} = reg.onsets{enc_i} - reg.onsets{fixenc_i} + reg.durations{enc_i};
        else
            enc_onsets{sbj_i}{count} = reg.onsets{fixenc_i};
            enc_durations{sbj_i}{count} = reg.onsets{enc_i} - reg.onsets{fixenc_i} + reg.durations{enc_i} + 6;
        end

    end
end


%% SETUP & PREPROCESSING

global CONN_gui; CONN_gui.usehighres = true; % to fix "usehighres" error

clear batch;

batch.filename = fullfile('../results/CONN', sprintf('conn_N%d_w%d.mat', num_sbj, WINDOW_SECONDS));
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
        batch.Setup.conditions.names={'rest', 'encoding'};
        batch.Setup.conditions.onsets{1}{sbj_i}{run_i} = 0;
        batch.Setup.conditions.durations{1}{sbj_i}{run_i} = inf;
        batch.Setup.conditions.onsets{2}{sbj_i}{run_i} = enc_onsets{sbj_i}{run_i};
        batch.Setup.conditions.durations{2}{sbj_i}{run_i} = enc_durations{sbj_i}{run_i};
    end
end
batch.Setup.rois.names = mask_names;
batch.Setup.rois.files = mask_paths;
batch.Setup.rois.add = 1; %use 1 to define an additional set of ROIs (to be added to any already-existing ROIs in your project)

%
batch.Setup.isnew = 1; % new file
batch.Setup.preprocessing.steps = 'default_mni';
batch.Setup.preprocessing.sliceorder = 'interleaved (bottom-up)'; % 1-3-5-...-2-4-...
batch.Setup.done = 1;
batch.Setup.overwrite = 'Yes';

conn_batch(batch);



%% DENOISING & FIRST-LEVEL ANALYSIS

clear batch;

batch.filename = fullfile('../results/CONN', sprintf('conn_N%d_w%d.mat', num_sbj, WINDOW_SECONDS));
batch.parallel.N = 4;

% DENOISING
% WM/CSF and other confounds use CONN defaults.
batch.Denoising.detrending = 1;
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
conn('load', fullfile('../results/CONN', sprintf('conn_N%d_w%d.mat', num_sbj, WINDOW_SECONDS)));

conn gui_results


