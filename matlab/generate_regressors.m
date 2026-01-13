% Episodic memory task for adolescents
% :: create regressor of each subject
% :: code written by Kahyun Choi


%% define parameters

BHV_RAW_PATH = 'C:\behavior_raw_directory'; % directory of behavioral raw data (mat)
BHV_FILE_PATH = 'C:\behavior_directory'; % directory of behavioral data (mat)
OUT_PATH = '../data/data_regressors/'; % directory to save regressor files (mat)

data_dir = dir(BHV_RAW_PATH); data_dir(1:2) = [];
tmp_idx = arrayfun(@(num) isfolder(fullfile(data_dir(num).folder, data_dir(num).name)), 1:length(data_dir));
data_dir(~tmp_idx) = []; 

sbj_date_list = cellfun(@(x) x(1:6), {data_dir.name}, 'UniformOutput', false); sbj_date_list(1) = [];
sbj_num_list = cellfun(@(x) x(9:11), {data_dir.name}, 'UniformOutput', false); sbj_num_list(1) = [];
sbj_num_array = cellfun(@(x) str2double(x), sbj_num_list);

d = dir(BHV_FILE_PATH);
file_names = {d.name};
file_names = file_names(cellfun(@(x) contains(x, 'bhv'), file_names));
sbj_list = cellfun(@(x) x(9:12), file_names, 'uni', 0); % ####: extracted from file name (bhvData_####.mat)


%% create regressor from behavior file

for sbj_i = 1:length(sbj_list)
    
    sbj_name = sbj_list{sbj_i};
    sbj_idx = find(sbj_num_array == str2double(sbj_name(2:end)));

    sbj_folder_name = [sbj_date_list{sbj_idx} '_A' sbj_num_list{sbj_idx}];
    bhv_file_name = [sbj_date_list{sbj_idx} '_sbj', sbj_num_list{sbj_idx}, '_WCDWWW.mat'];

    if ~exist(fullfile(BHV_RAW_PATH, sbj_folder_name, bhv_file_name), 'file'); continue; end
    load(fullfile(BHV_RAW_PATH, sbj_folder_name, bhv_file_name)); % 'Data'

    if isempty(Data); continue; end
    

    %% extract timing info from raw behavioral data
    
    trial_types = {'inst_enc_1', 'fix_enc_1', 'enc_1', 'inst_int_1', 'interfere_1', ... % 1~5
                   'inst_ret_1', 'fix_ret_1', 'ret_1', 'inst_enact_1', 'enact_1', ...   % 6~10
                   'inst_enc_2', 'fix_enc_2', 'enc_2', 'inst_int_2', 'interfere_2', ...`% 11~15
                   'inst_ret_2', 'fix_ret_2', 'ret_2', 'inst_enact_2', 'enact_2', ...   % 16~20
                   'ctrl_inst_enc', 'ctrl_fix_enc', 'ctrl_enc', 'ctrl_inst_ret', 'ctrl_fix_ret', 'ctrl_ret'}; % 21~26
    
    tmp = {Data.time.label};
    runIdx = find(ismember(tmp, 'run start'));
    instIdx = find(contains(tmp, 'instruction'));
    
    fixEncIdx = find(ismember(tmp, 'enc fixation'));
    fixSimIdx = find(ismember(tmp, 'sim fixation'));
    
    encStartIdx = find(ismember(tmp, 'enc start'));
    encEndIdx = find(ismember(tmp, 'enc end'));
    intStartIdx = find(ismember(tmp, 'int start'));
    intEndIdx = find(ismember(tmp, 'sim instruction'));
    simStartIdx = find(ismember(tmp, 'sim start'));
    simEndIdx = find(ismember(tmp, 'sim end'));
    retStartIdx = find(ismember(tmp, 'ret start'));
    retEndIdx = find(ismember(tmp, 'ret end'));
    
    onsets = {}; durations = {};
    for run_i = 1:length(runIdx)
        if run_i < length(runIdx)
            nextRunIdx = runIdx(run_i + 1);
        else
            nextRunIdx = size(Data.time, 2);
        end
        runTime = Data.time(runIdx(run_i)).time;
        
        tmpOn = []; tmpDur = [];
        for i = 1:length(instIdx) % instruction
            if instIdx(i) > runIdx(run_i) && instIdx(i) < nextRunIdx
                tmpOn = [tmpOn; Data.time(instIdx(i)).time - runTime];
                tmpDur = [tmpDur; Data.time(instIdx(i)+1).time - Data.time(instIdx(i)).time];
            else
                continue;
            end
        end
        onsets{run_i}{1} = tmpOn(1);        durations{run_i}{1} = tmpDur(1);    % inst_enc
        onsets{run_i}{11} = tmpOn(5);       durations{run_i}{11} = tmpDur(5);    
        onsets{run_i}{4} = tmpOn(2);        durations{run_i}{4} = tmpDur(2);    % inst_int
        onsets{run_i}{14} = tmpOn(6);       durations{run_i}{14} = tmpDur(6);
        onsets{run_i}{6} = tmpOn(3);        durations{run_i}{6} = tmpDur(3);    % inst_ret
        onsets{run_i}{16} = tmpOn(7);       durations{run_i}{16} = tmpDur(7);
        onsets{run_i}{9} = tmpOn(4);        durations{run_i}{9} = tmpDur(4);    % inst_enact
        onsets{run_i}{19} = tmpOn(8);       durations{run_i}{19} = tmpDur(8);
        onsets{run_i}{21} = tmpOn(9);       durations{run_i}{21} = tmpDur(9);   % ctrl_inst_enc
        onsets{run_i}{24} = tmpOn(10);      durations{run_i}{24} = tmpDur(10);  % ctrl_inst_ret
        
        tmpOn = []; tmpDur = [];
        for i = 1:length(fixEncIdx) % fixation encoding
            if fixEncIdx(i) > runIdx(run_i) && fixEncIdx(i) < nextRunIdx
                tmpOn = [tmpOn; Data.time(fixEncIdx(i)).time - runTime];
                tmpDur = [tmpDur; Data.time(fixEncIdx(i)+1).time - Data.time(fixEncIdx(i)).time];
            else
                continue;
            end
        end
        onsets{run_i}{2} = tmpOn(1);        durations{run_i}{2} = tmpDur(1);    % fix_enc
        onsets{run_i}{12} = tmpOn(2);       durations{run_i}{12} = tmpDur(2);   % fix_enc
        onsets{run_i}{22} = tmpOn(3);       durations{run_i}{22} = tmpDur(3);   % ctrl_fix_enc
        
        tmpOn = []; tmpDur = [];
        for i = 1:length(fixSimIdx) % fixation retrieval
            if fixSimIdx(i) > runIdx(run_i) && fixSimIdx(i) < nextRunIdx
                tmpOn = [tmpOn; Data.time(fixSimIdx(i)).time - runTime];
                tmpDur = [tmpDur; Data.time(fixSimIdx(i)+1).time - Data.time(fixSimIdx(i)).time];
            else
                continue;
            end
        end
        onsets{run_i}{7} = tmpOn(1);        durations{run_i}{7} = tmpDur(1);    % fix_ret
        onsets{run_i}{17} = tmpOn(2);       durations{run_i}{17} = tmpDur(2);   % fix_ret
        onsets{run_i}{25} = tmpOn(3);       durations{run_i}{25} = tmpDur(3);   % ctrl_fix_ret
        
        tmpOn = []; tmpDur = [];
        for i = 1:length(encStartIdx) % encoding
            if encStartIdx(i) > runIdx(run_i) && encStartIdx(i) < nextRunIdx
                tmpOn = [tmpOn; Data.time(encStartIdx(i)).time - runTime];
                tmpDur = [tmpDur; Data.time(encEndIdx(i)).time - Data.time(encStartIdx(i)).time];
            else
                continue;
            end
        end
        onsets{run_i}{3} = tmpOn(1);        durations{run_i}{3} = tmpDur(1);    % enc
        onsets{run_i}{13} = tmpOn(2);       durations{run_i}{13} = tmpDur(2);   % enc
        onsets{run_i}{23} = tmpOn(3);       durations{run_i}{23} = tmpDur(3);   % ctrl_enc
    
        tmpOn = []; tmpDur = [];
        for i = 1:length(intStartIdx) % interference
            if intStartIdx(i) > runIdx(run_i) && intStartIdx(i) < nextRunIdx
                tmpOn = [tmpOn; Data.time(intStartIdx(i)).time - runTime];
                tmpDur = [tmpDur; Data.time(intStartIdx(i)+1).time - Data.time(intStartIdx(i)).time];
            else
                continue;
            end
        end
        onsets{run_i}{5} = tmpOn(1);        durations{run_i}{5} = tmpDur(1);    % interfere
        onsets{run_i}{15} = tmpOn(2);       durations{run_i}{15} = tmpDur(2);   % interfere
    
        tmpOn = []; tmpDur = [];
        for i = 1:length(simStartIdx) % retrieval
            if simStartIdx(i) > runIdx(run_i) && simStartIdx(i) < nextRunIdx
                tmpOn = [tmpOn; Data.time(simStartIdx(i)).time - runTime];
                tmpDur = [tmpDur; Data.time(simEndIdx(i)).time - Data.time(simStartIdx(i)).time];
            else
                continue;
            end
        end
        onsets{run_i}{8} = tmpOn(1);        durations{run_i}{8} = tmpDur(1);    % ret
        onsets{run_i}{18} = tmpOn(2);       durations{run_i}{18} = tmpDur(2);   % ret
        onsets{run_i}{26} = tmpOn(3);       durations{run_i}{26} = tmpDur(3);   % ctrl_fix_ret
    
        tmpOn = []; tmpDur = [];
        for i = 1:length(retStartIdx) % reenactment
            if retStartIdx(i) > runIdx(run_i) && retStartIdx(i) < nextRunIdx
                tmpOn = [tmpOn; Data.time(retStartIdx(i)).time - runTime];
                tmpDur = [tmpDur; Data.time(retEndIdx(i)).time - Data.time(retStartIdx(i)).time];
            else
                continue;
            end
        end
        onsets{run_i}{10} = tmpOn(1);       durations{run_i}{10} = tmpDur(1);   % enact
        onsets{run_i}{20} = tmpOn(2);       durations{run_i}{20} = tmpDur(2);   % enact
    end
    onsets_all = onsets; durations_all = durations; trial_types_all = trial_types;


    %% save regressors
    
    for run_i = 1:length(onsets_all)
        onsets = onsets_all{run_i};
        durations = durations_all{run_i};
        names = trial_types_all;

        out_file_name = sprintf('%s_run%d_trial.mat', sbj_name, run_i);
        save(fullfile(OUT_PATH, out_file_name), 'names', 'onsets', 'durations');
    end

    clear Data onsets durations names
end
