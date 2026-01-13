% Episodic memory task for adolescents
% :: preprocessing - raw behavioral data
% :: code written by Kahyun Choi


%% define parameters

bhv_raw_path = 'C:\raw_behavior_directory'; % directory of raw behavioral data
bhv_save_path = 'C:\behavior_directory'; % directory to save preprocessed behavioral data

bhv_dir = dir(bhv_raw_path); bhv_dir(1:2) = [];
sbj_date = cellfun(@(x) x(1:6), {bhv_dir.name}, 'uni', 0);
sbj_list = cellfun(@(x) x(9:11), {bhv_dir.name}, 'uni', 0);
sbj_nums = cellfun(@(x) str2double(x), sbj_num_list);


%% preprocess excel data

xlsx_data_path = fullfile(bhv_raw_path, 'sbj_data.xlsx');
sbj_data = readtable(xlsx_data_path);

sbj_age = sbj_data.Birth_age;

svy_list = sbj_data.Properties.VariableNames;
svy_list = svy_list(2:end);

for svy_i = 1:length(svy_list)
    sbj_svy{svy_i} = eval(['sbj_data.', svy_list{svy_i}]);
end

% SES: no response -> NaN
idx = find(ismember(svy_list, 'SES'));
tmp = sbj_svy{idx};
ind = (tmp == 17);
tmp(ind) = NaN;
sbj_svy{idx} = tmp;


%% preprocess behavioral data

for sbj_i = 1:length(sbj_nums)
    sbj_idx = find(sbj_nums == sbj_i + 100);

    sbj_folder_name = [sbj_date{sbj_idx} '_A' sbj_list{sbj_idx}];
    bhv_file_name = [sbj_date{sbj_idx} '_sbj', sbj_list{sbj_idx}, '_WCDWWW.mat'];
    sbj_file_path = fullfile(bhv_raw_path, sbj_folder_name, bhv_file_name);
    
    if ~exist(sbj_file_path, 'file')
        bhv = [];
    else
        load(sbj_file_path); % 'Data'

        %%
        bhv = [];
        
        % sbjInfo
        bhv.sbjInfo.date = Data.sbj.sbjDate;
        bhv.sbjInfo.num = Data.sbj.sbjNum;
        bhv.sbjInfo.age = Data.sbj.sbjAge;
        bhv.sbjInfo.sex = Data.sbj.sbjSex;

        bhv.sbjInfo.prefer = Data.sbj.preferAnimal; % most prefer = first
        bhv.sbjInfo.learned = Data.sbj.learnedAnimal;
        
        % expInfo
        bhv.expInfo.context = Data.trial.trialinfo.context;
        bhv.expInfo.animal = Data.trial.trialinfo.animal;
        bhv.expInfo.location = Data.trial.trialinfo.location;
        if ~interf_type
            bhv.expInfo.inter_back = Data.trial.trialinfo.interference.back;
            bhv.expInfo.inter_stim = Data.trial.trialinfo.interference.stim;
            if str2double(Data.sbj.sbjNum) < 108
                bhv.expInfo.inter_word_one = Data.trial.trialinfo.interference.wordOne;
                bhv.expInfo.inter_word_var = Data.trial.trialinfo.interference.wordVar;
            end
        else
            bhv.expInfo.inter_type = Data.trial.trialinfo.interference.type;
            bhv.expInfo.inter_stim = Data.trial.trialinfo.interference.stim;
        end
        
        % reenactment
        ret_animal = []; ret_location = []; ret_animalAcc = []; ret_locationAcc = [];
        for i = 1:max([Data.trial.retrieval.trial])
            tmp = [Data.trial.retrieval.animal];
            ret_animal = [ret_animal; tmp(4*i-3:4*i)];
            tmp = [Data.trial.retrieval.location];
            ret_location = [ret_location; tmp(4*i-3:4*i)];
            tmp = [Data.trial.retrieval.animalAcc];
            ret_animalAcc = [ret_animalAcc; tmp(4*i-3:4*i)];
            tmp = [Data.trial.retrieval.locationAcc];
            ret_locationAcc = [ret_locationAcc; tmp(4*i-3:4*i)]; 
        end
        bhv.reenact.animal = ret_animal; bhv.reenact.location = ret_location;
        bhv.reenact.animalAcc = ret_animalAcc; bhv.reenact.locationAcc = ret_locationAcc;
        
        %% organized
        % em_acc        [what where what-where what-when where-when www]
        
        enact = bhv.reenact;
        
        what_resp = bhv.reenact.animal; where_resp = bhv.reenact.location;
        ww_resp = what_resp*10 + where_resp;
        what_ans = bhv.expInfo.animal; where_ans = bhv.expInfo.location;
        ww_ans = what_ans*10 + where_ans;
        
        if size(what_resp, 1) ~= size(what_ans, 1)
            what_ans = what_ans(1:size(what_resp, 1), :);
            ww_ans = ww_ans(1:size(ww_resp, 1), :);
        end
        
        what_acc_list = []; where_acc_list = []; ww_acc_list = []; 
        what_err_list = []; where_err_list = []; where_err_animal = [];

        num_block = size(what_resp, 1);
        for block_i = 1:num_block
            if isnan(what_resp(block_i, 1))
                block_what_ans = NaN;
                block_where_ans = NaN;
                block_ww_ans = NaN;
            else
                block_what_ans = length(intersect(what_ans(block_i,:), what_resp(block_i,:)))/4;
                block_where_ans = length(intersect(where_ans(block_i,:), where_resp(block_i,:)))/4;
                block_ww_ans = length(intersect(ww_ans(block_i,:), ww_resp(block_i,:)))/4;
            end
        
            what_acc_list = [what_acc_list;block_what_ans];
            where_acc_list = [where_acc_list;block_where_ans];
            ww_acc_list = [ww_acc_list;block_ww_ans];
        
            what_fail_animal = setdiff(what_ans(block_i, :), what_resp(block_i, :));
            if ~isempty(what_fail_animal)
                what_err_list = [what_err_list what_fail_animal];
            end
            
            where_fail_location = setdiff(where_ans(block_i, :), where_resp(block_i, :));
            fail_event = [];
            if ~isempty(where_fail_location)
                for f_i = 1:length(where_fail_location)
                    fail_event(f_i) = find(where_ans(block_i, :) == where_fail_location(f_i));
                end
                where_err_list = [what_err_list where_fail_location];
                where_err_animal = [where_err_animal what_ans(block_i, fail_event)];
            end
        end
        trial_acc = [what_acc_list'; where_acc_list'; ww_acc_list'; nanmean(enact.animalAcc, 2)'; nanmean(enact.locationAcc, 2)'; ...
                     nanmean((enact.animalAcc + enact.locationAcc) == 2, 2)'];
        
        what_acc = nanmean(what_acc_list);
        where_acc = nanmean(where_acc_list);
        ww_acc = nanmean(ww_acc_list);
        
        what_when_acc = nanmean(enact.animalAcc, 'all'); 
        where_when_acc = nanmean(enact.locationAcc, 'all');
        www_acc = nanmean((enact.animalAcc + enact.locationAcc) == 2, 'all');

        em_acc = [what_acc where_acc ww_acc what_when_acc where_when_acc www_acc];
        
        %
        bhv.organized.em = em_acc;
        bhv.organized.trial_acc = trial_acc;

        %%
        if ~isnan(sbj_age(sbj_i))
            bhv.sbjInfo.age = num2str(sbj_age(sbj_i));
        end
        bhv.sbjInfo.svy_value = cell2mat(cellfun(@(x) x(sbj_i), sbj_svy, 'uni', 0));
        bhv.sbjInfo.svy_name = svy_list;
    end
    
    sbj_file_name = ['bhvData_A' sbj_list{sbj_idx} '.mat'];
    out_file_path = fullfile(bhv_save_path, sbj_file_name);
    save(out_file_path, 'bhv');
    
    fprintf([num2str(sbj_i), '\t']);
end
disp(' ');
disp('preprocessing done');


%% organize behavioral data

tmp_em_acc = {};
for sbj_i = 1:length(bhv_list)
    tmp_em_acc{sbj_i} = bhv_list{sbj_i}.organized.em;
end

em_acc = arrayfun(@(num) cellfun(@(x) x(num), tmp_em_acc), 1:6, 'uni', 0);

sbj_age = cellfun(@(x) x.sbjInfo.age, bhv_list, 'uni', 0);
sbj_age = str2double(sbj_age);

sbj_sex = cellfun(@(x) x.sbjInfo.sex, bhv_list, 'uni', 0);
sbj_sex = ismember(sbj_sex, {'F', 'f'}); % 0: male, 1: female

sbj_svy = {};
for svy_i = 1:length(bhv_list{1}.sbjInfo.svy_name)
    sbj_svy{svy_i} = cell2mat(cellfun(@(x) x.sbjInfo.svy_value(svy_i), bhv_list, 'uni', 0))';
end
svy_label = bhv_list{1}.sbjInfo.svy_name;


% save subject list
subject_path = 'C:\subject_info_directory';
save(fullfile(subject_path, 'behavioral_data.mat'), 'em_acc', 'sbj_age', 'sbj_sex');
save(fullfile(subject_path, 'survey_data.mat'), 'sbj_svy', 'svy_label');
