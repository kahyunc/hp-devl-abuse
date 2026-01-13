% Episodic memory task for adolescents
% :: ROI data extraction (from FreeSurfer output)
% :: code written by Kahyun Choi


%% define parameters

seg_dir = 'C:\freesurfer_directory_path'; % directory of FreeSurfer output files
out_dir = 'C:\directory_to_save'; % directory to save ROI extracted data
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

d = dir(seg_dir);
d_names = {d.name};
sbj_list = d_names(contains(d_names, 'A'));
sbj_nums = cell2mat(arrayfun(@(num) str2double(sbj_list{num}(2:end))-100, 1:length(sbj_list), 'uni', 0));


%% ROI extraction

for sbj_i = 1:length(sbj_list)
    sbj_name = sbj_list{sbj_i};
    
    IN_DIR = fullfile(seg_dir, sbj_name);
    OUT_FILE_NAME = fullfile(out_dir, sprintf('%s.mat', sbj_name));

    try
        %% setup

        sbj = struct();
        
        ctx_seg_file = 'aparc_aseg.nii.gz';
        mask_file = 'brainmask.nii.gz';

        left_hpc_seg_file = 'lh_hippoAmygLabels_T1_v21_FS60_FSvoxelSpace.nii.gz';
        right_hpc_seg_file = 'rh_hippoAmygLabels_T1_v21_FS60_FSvoxelSpace.nii.gz';
        left_hpc_hbt_file = 'lh_hippoAmygLabels_T1_v21_HBT_FSvoxelSpace.nii.gz';
        right_hpc_hbt_file = 'rh_hippoAmygLabels_T1_v21_HBT_FSvoxelSpace.nii.gz';


        %% get segmentation data

        ctx_seg_file = fullfile(IN_DIR, ctx_seg_file);
        sbj.seg.ctx = niftiread(ctx_seg_file);
        sbj.seg.header = niftiinfo(ctx_seg_file);
        
        mask_file = fullfile(IN_DIR, mask_file);
        sbj.seg.mask = niftiread(mask_file);
        
        left_hpc_seg_file = fullfile(IN_DIR, left_hpc_seg_file);
        sbj.seg.hpc_l = niftiread(left_hpc_seg_file);
        
        right_hpc_seg_file = fullfile(IN_DIR, right_hpc_seg_file);
        sbj.seg.hpc_r = niftiread(right_hpc_seg_file);
        
        left_hpc_hbt_file = fullfile(IN_DIR, left_hpc_hbt_file);
        sbj.seg.hpc_hbt_l = niftiread(left_hpc_hbt_file);
        
        right_hpc_hbt_file = fullfile(IN_DIR, right_hpc_hbt_file);
        sbj.seg.hpc_hbt_r = niftiread(right_hpc_hbt_file);
        
        
        %% organize hpc segmentation
        
        % MTL (MTL regions: hippocampus, entorhinal, parahippocampal)
        roi_name_list =  {'Hp','EC','PHC'};
        
        source_mat = sbj.seg.ctx;
        idx_hpc_l = 17; idx_hpc_r = 53;
        idx_ec_l = 1006; idx_ec_r = 2006;
        idx_phc_l = 1016;idx_phc_r = 2016;
        
        ind_list_l = {idx_hpc_l, idx_ec_l, idx_phc_l};
        ind_list_r = {idx_hpc_r, idx_ec_r, idx_phc_r};
        mask_l = cellfun(@(x) source_mat == x, ind_list_l, 'uni', 0);
        mask_r = cellfun(@(x) source_mat == x, ind_list_r, 'uni', 0);
        
        roi_list = [ mask_l, mask_r, cellfun(@(x,y) x|y, mask_l, mask_r, 'uni', 0) ];
        roi_name_list = [ cellfun(@(x) ['Lt.' x], roi_name_list,'UniformOutput',false), ...
                         cellfun(@(x) ['Rt.' x], roi_name_list,'UniformOutput',false), ...
                         cellfun(@(x) ['Bi.' x], roi_name_list,'UniformOutput',false)  ];
        
        sbj.seg.hpc_organized.mtl.roi_list = roi_list;
        sbj.seg.hpc_organized.mtl.roi_name_list = roi_name_list;
        
        
        % HBT (AP axis subregions: head, body, tail, body+tail (posterior))
        roi_name_list = {'head','body','tail','B+T'};
        idx_list = [232, 231, 226]; % head body tail
        
        source_mat = sbj.seg.hpc_hbt_l; 
        mask_l = arrayfun(@(x) source_mat==x, idx_list, 'UniformOutput',false);
        mask_l = [mask_l, mask_l{2}|mask_l{3}];
        
        source_mat = sbj.seg.hpc_hbt_r; 
        mask_r = arrayfun(@(x) source_mat==x, idx_list, 'UniformOutput',false);
        mask_r = [mask_r, mask_r{2}|mask_r{3}];
        
        roi_list = [ mask_l, mask_r, cellfun(@(x,y) x|y, mask_l, mask_r, 'uni', 0) ];
        roi_name_list = [ cellfun(@(x) ['Lt.' x], roi_name_list,'UniformOutput',false), ...
                         cellfun(@(x) ['Rt.' x], roi_name_list,'UniformOutput',false), ...
                         cellfun(@(x) ['Bi.' x], roi_name_list,'UniformOutput',false)  ];
        
        sbj.seg.hpc_organized.hbt.roi_list = roi_list;
        sbj.seg.hpc_organized.hbt.roi_name_list = roi_name_list;
        
        
        % MAIN (subfields: CA23DG, CA1, subiculum)
        roi_name_list = {'CA23DG','CA1','Sub'};
        idx_list = [210, 208, 206, 205];
        
        source_mat = sbj.seg.hpc_l; 
        mask_l = arrayfun(@(x) source_mat==x, idx_list, 'UniformOutput',false);
        mask_l = [mask_l{1}| mask_l{2}, mask_l(3:end)];
        
        source_mat = sbj.seg.hpc_r; 
        mask_r = arrayfun(@(x) source_mat==x, idx_list, 'UniformOutput',false);
        mask_r = [mask_r{1}| mask_r{2}, mask_r(3:end)];
        
        roi_list = [ mask_l, mask_r, cellfun(@(x,y) x|y, mask_l, mask_r, 'uni', 0) ];
        roi_name_list = [ cellfun(@(x) ['Lt.' x], roi_name_list,'UniformOutput',false), ...
                         cellfun(@(x) ['Rt.' x], roi_name_list,'UniformOutput',false), ...
                         cellfun(@(x) ['Bi.' x], roi_name_list,'UniformOutput',false)  ];
        
        sbj.seg.hpc_organized.main.roi_list = roi_list;
        sbj.seg.hpc_organized.main.roi_name_list = roi_name_list;
                
        
        % cortical regions (thalamus, striatum, amygdala)
        source_mat = sbj.seg.ctx;
                
        roi_name_list = {'thalamus','caudate','putamen','pallidum','amygdala','NAcc'};
        ind_list_l = [10,11,12,13,18,26];
        ind_list_r = [49,50,51,52,54,58];
        
        mask_l = arrayfun(@(x) source_mat == x, ind_list_l, 'uni', 0);
        mask_r = arrayfun(@(x) source_mat == x, ind_list_r, 'uni', 0);
        roi_list = [mask_l, mask_r, cellfun(@(x,y) x|y, mask_l, mask_r, 'uni', 0)];
        
        roi_name_list = [ cellfun(@(x) ['Lt.' x], roi_name_list,'UniformOutput',false), ...
                         cellfun(@(x) ['Rt.' x], roi_name_list,'UniformOutput',false), ...
                         cellfun(@(x) ['Bi.' x], roi_name_list,'UniformOutput',false)  ];
        
        sbj.seg.ctx_organized.subc.roi_name_list = roi_name_list;
        sbj.seg.ctx_organized.subc.roi_list = roi_list;
        
        
        %% save
        
        sbj = rmfield(sbj,'seg');
        save(OUT_FILE_NAME, 'sbj');

        fprintf('\n %d done \n', sbj_i);

    catch e
        disp(sbj_i)
        disp(e.message)
    end
end

