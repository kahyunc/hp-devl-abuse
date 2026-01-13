% Episodic memory task for adolescents
% :: fMRI data preprocessing (SPM12)
% :: code written by Kahyun Choi


%% define parameters

IN_PATH = 'C:\fMRI_directory'; % directory of fMRI data (nifti)

DIR_T1_NAME = 'T1';
DIR_FUNC_LIST = {'M1', 'M2', 'M3', 'M4'};

sbj_list = dir(IN_PATH);
sbj_list = arrayfun(@(x) x.name, sbj_list(3:end), 'uni', 0);


%% run

fail_list = [];
for sbj_i = 1:length(sbj_list)
    try
        sbj_name = sbj_list{sbj_i};

        dir_t1 = fullfile(IN_PATH, sbj_name, DIR_T1_NAME);
        file_name_t1 = dir(dir_t1);
        file_name_t1 = file_name_t1(3);
        file_name_t1 = fullfile(file_name_t1.folder, file_name_t1.name);

        func_dir_list = cellfun(@(x) fullfile(IN_PATH, sbj_name, x), DIR_FUNC_LIST, 'uni', 0);
        func_list = {};
        for run_i = 1:length(func_dir_list)
            file_list = dir(func_dir_list{run_i});
            file_list = arrayfun(@(x) fullfile(x.folder, x.name), file_list(3:end), 'uni', 0);
            func_list{run_i} = file_list(:);
        end

        if sum(cellfun(@(x) isempty(x), func_list)) ~= 0
            idx = cellfun(@(x) isempty(x), func_list);
            func_list(idx) = [];
        end

        %% preprocessing setup

        func = func_list(:);
        anat = {file_name_t1};
        
        epi_sample = func{1}{1};
        header = niftiinfo(epi_sample);
        nslice = header.ImageSize(3);
        
        time_info = split(header.Description(strfind(header.Description,'TR'):end),'/');
        tr = time_info{contains(time_info,'TR')};
        tr = split(tr,{'=','m'});
        tr = str2double(tr{2}) / 1000;
        
        vox_size = header.PixelDimensions(1); 
        smoothing_factor = 8;
        
        %% 1) slice timing
        
        matlabbatch{1}.spm.temporal.st.scans = func;
        
        % slice timing
        matlabbatch{1}.spm.temporal.st.nslices = nslice;
        matlabbatch{1}.spm.temporal.st.tr = tr;
        matlabbatch{1}.spm.temporal.st.ta = tr-(tr/nslice);
        matlabbatch{1}.spm.temporal.st.so = [1:2:nslice 2:2:nslice]; 
        matlabbatch{1}.spm.temporal.st.refslice = floor(nslice/2);
        matlabbatch{1}.spm.temporal.st.prefix = 'a';
        
        
        %% 2) realignment

        for run_i = 1:length(func)
            matlabbatch{2}.spm.spatial.realign.estwrite.data{run_i}(1) = cfg_dep(sprintf('Slice Timing: Slice Timing Corr. Images (Sess %d)',run_i), ...
                                                                                 substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{run_i}, '.','files'));
        end
        matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
        
        %% 3) coregistration

        matlabbatch{3}.spm.spatial.coreg.estwrite.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
        matlabbatch{3}.spm.spatial.coreg.estwrite.source = anat;
        matlabbatch{3}.spm.spatial.coreg.estwrite.other = {''};
        matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        
        %% 4) segmentation

        matlabbatch{4}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Coregister: Estimate & Reslice: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
        matlabbatch{4}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{4}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{4}.spm.spatial.preproc.channel.write = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(1).tpm = {[spm_path '\tpm\TPM.nii,1']}; % gray
        matlabbatch{4}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{4}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(2).tpm = {[spm_path '\tpm\TPM.nii,2']}; % white
        matlabbatch{4}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{4}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(3).tpm = {[spm_path '\tpm\TPM.nii,3']}; % csf
        matlabbatch{4}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{4}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(4).tpm = {[spm_path '\tpm\TPM.nii,4']}; % soft tissue
        matlabbatch{4}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{4}.spm.spatial.preproc.tissue(4).native = [1 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(5).tpm = {[spm_path '\tpm\TPM.nii,5']}; % bone
        matlabbatch{4}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{4}.spm.spatial.preproc.tissue(5).native = [1 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(6).tpm = {[spm_path '\tpm\TPM.nii,6']}; % other
        matlabbatch{4}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{4}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{4}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{4}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{4}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{4}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{4}.spm.spatial.preproc.warp.affreg = 'eastern';
        matlabbatch{4}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{4}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{4}.spm.spatial.preproc.warp.write = [0 1];
        
        
        %% 5) normalization

        matlabbatch{5}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
        for run_i = 1:length(func)
            matlabbatch{5}.spm.spatial.normalise.write.subj.resample(run_i) = cfg_dep(sprintf('Realign: Estimate & Reslice: Resliced Images (Sess %d)',run_i), ...
                                                                                      substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{run_i}, '.','rfiles'));
        end
        matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                                  78 76 85];
        matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = [vox_size vox_size vox_size];
        matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{5}.spm.spatial.normalise.write.woptions.prefix = 'w';
        
        %% 6) smooth

        matlabbatch{6}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
        matlabbatch{6}.spm.spatial.smooth.fwhm = [smoothing_factor smoothing_factor smoothing_factor];
        matlabbatch{6}.spm.spatial.smooth.dtype = 0;
        matlabbatch{6}.spm.spatial.smooth.im = 0;
        matlabbatch{6}.spm.spatial.smooth.prefix = 's';
        
        %% run

        batch = matlabbatch;

        spm('defaults','fmri');
        spm_jobman('initcfg');
        spm_jobman('run',batch);
        
    catch
        fail_list = [fail_list sbj_i];
    end
end
fprintf('\nfailed sbj_i: %d\n\n', fail_list);



