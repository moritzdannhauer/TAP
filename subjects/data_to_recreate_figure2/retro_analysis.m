function retro_data = retro_analysis(retro_setup)
%Dependencies:
%subj '/' target_name '_HT' HT '/target.msh' for roi index
%subj '/' subj '_bin.msh' for simnibs fwdsim
%subj '/' subj '_T1fs_conform.nii.gz' (the MR scans)
%mask (e.g.'Subj1PMC.nii.gz')

% computation initialization
global num_of_subj sesh_per_subj blck_per_sesh recover_switch
practice_switch = sum([strcmp(retro_setup.subj,{'Subj1';'Subj2';'Subj3';'Subj4';'Subj5'});...
    strcmp(retro_setup.sesh_names,{'PMC';'PSC';'sham'});...
    strcmp(retro_setup.opt_name,{'sSubj1_PMC_hair0.5mm';'sSubj1_PSC_hair0.5mm';'sSubj1_PMC_hair0.5mm';...
    'sSubj2_PreMotor_hair2.5mm';'sSubj2_PSC_hair2.5mm';'sSubj2_PreMotor_hair2.5mm';...
    'sSubj3_PMC_hair1.5mm';'sSubj3_PSC3_hair2.5mm';'sSubj3_PMC_hair1.5mm';...
    'sSubj4_PMC3_hair2.5mm';'sSubj4_PSC2_hair2.0mm';'sSubj4_PMC3_hair2.5mm';...
    'sSubj5_PMC2_hair2.0mm';'sSubj5_PSC_hair2.0mm';'sSubj5_PMC2_hair3.5mm'})]);
if practice_switch == 23; load('example_data.mat'); end

addpath('matlab');
if isunix;  addpath(char([retro_setup.path2simnibs '/matlab'])); addpath('matlab/additionals/');end
if ispc;    addpath(char([retro_setup.path2simnibs '\matlab'])); addpath('matlab\additionals\');end
addpath(retro_setup.snapped_folder);

num_of_subj = length(retro_setup.subj); sesh_per_subj = length(retro_setup.sesh_names); blck_per_sesh = length(retro_setup.bs_name_s);
sign_check_indx = [1,5,9;2,6,10;3,7,11;4,8,12];
x_axis = [-1,0,0]; y_axis = [0,1,0]; z_axis = [0,0,1];%in BrainSight space

if recover_switch == 1
    load('backup_data.mat')
    last_sesh = data_bckup.last_saved;
    recover_tmp = strsplit(last_sesh,'\t');
    
    subj2start_str = recover_tmp{1};
    sesh2start_str = recover_tmp{2};
    blck2start_str = recover_tmp{4};
    
    disp(['[TAP] Last saved position was ' subj2start_str sesh2start_str ' block' blck2start_str])
    
    % set restart position
    subj2start = find(strcmp(retro_setup.subj,subj2start_str)==1);
    sesh2start = find(strcmp(retro_setup.sesh_names,sesh2start_str)==1);
    blck2start = str2double(blck2start_str);
    
    % recover data
    retro_data = data_bckup;
    retro_data = rmfield(retro_data,'last_saved');% remove the 'last_saved' field
    
    % clean to avoid duplicate data saving
    data_fns = fieldnames(retro_data);
    if retro_setup.E50_switch == 1; data_length = length(data_fns)-1;
    else; data_length = length(data_fns); end
    
    for metric_indx = 1:data_length
        data2clean = retro_data.(data_fns{metric_indx});
        clean_start = find(data2clean{sesh2start,2,subj2start} == blck2start,1);
        data2clean{sesh2start,1,subj2start}(clean_start:end) = [];
        data2clean{sesh2start,2,subj2start}(clean_start:end) = [];
        retro_data.(data_fns{metric_indx}) = data2clean;
    end
    
else
    retro_data.norml_dist = cell(sesh_per_subj,2,num_of_subj);%col1 = data, col2 = stim block
    retro_data.tan_dist = cell(sesh_per_subj,2,num_of_subj);%col1 = data, col2 = stim block
    retro_data.norml_ang = cell(sesh_per_subj,2,num_of_subj);%col1 = data, col2 = stim block
    retro_data.tan_ang = cell(sesh_per_subj,2,num_of_subj);%col1 = data, col2 = stim block
    retro_data.ef_mag = cell(sesh_per_subj,3,num_of_subj);%col1 = data, col2 = stim block , col3 = raw ||E||
    retro_data.ef_dir = cell(sesh_per_subj,2,num_of_subj);%col1 = data, col2 = stim block
    
    retro_data.norml_dist_out_label = cell(sesh_per_subj,1,num_of_subj);
    retro_data.tan_dist_out_label = cell(sesh_per_subj,1,num_of_subj);
    retro_data.norml_ang_out_label = cell(sesh_per_subj,1,num_of_subj);
    retro_data.tan_ang_out_label = cell(sesh_per_subj,1,num_of_subj);
    retro_data.ef_mag_out_label = cell(sesh_per_subj,1,num_of_subj);
    retro_data.ef_dir_out_label = cell(sesh_per_subj,1,num_of_subj);
    
    if retro_setup.E50_switch == 1; retro_data.E50 = zeros(blck_per_sesh,sesh_per_subj,num_of_subj); end
    
    subj2start = 1; sesh2start = 1; blck2start = 1;
end

if exist('retro_tmp','dir')
    if isunix;  system('rm retro_tmp/*');end
    if ispc;  system('rm retro_tmp\*');end
    system('rmdir retro_tmp');
end

% computation
for subj_indx = subj2start:num_of_subj
    
    curt_subj = retro_setup.subj{subj_indx};
    if isunix;  path_to_bin_msh = [curt_subj '/' curt_subj '_bin.msh'];path_to_mr_scans = [curt_subj '/' curt_subj '_T1fs_conform.nii.gz'];end
    if ispc;  path_to_bin_msh = [curt_subj '\' curt_subj '_bin.msh'];path_to_mr_scans = [curt_subj '\' curt_subj '_T1fs_conform.nii.gz'];end
    
    for sesh_indx_wthn_subj = sesh2start:sesh_per_subj
        
        curt_sesh = retro_setup.sesh_names{sesh_indx_wthn_subj};
        sesh_indx_acrs_all = (subj_indx-1)*sesh_per_subj + sesh_indx_wthn_subj;
        curt_HT = retro_setup.HTs(sesh_indx_acrs_all);
        
        sign_opt = retro_setup.sign_opt_all{sesh_indx_acrs_all};
        
        if isunix;  Topt_bs = read_brainsight_file([curt_subj '/' curt_sesh '/' retro_setup.opt_name{sesh_indx_acrs_all} '.txt']);end
        if ispc;  Topt_bs = read_brainsight_file([curt_subj '\' curt_sesh '\' retro_setup.opt_name{sesh_indx_acrs_all} '.txt']);end
        Topt_bs = Topt_bs{1};
        
        opt_center_bs = Topt_bs(1:3,4); opt_z_bs = Topt_bs(1:3,3);opt_sign = sign(Topt_bs);
        %norm(vecA) * cos(ang_vecA_PlaneNormal) * PlaneNormal + vecA_projected = vecA
        %so vecA_projected = vecA - norm(vecA) * cos(ang_vecA_PlaneNormal) * PlaneNormal
        proj_opt = Topt_bs(1:3,2) - norm(Topt_bs(1:3,2)) * cosd( angv1v2(Topt_bs(1:3,2),opt_z_bs) ) * opt_z_bs/norm(opt_z_bs);

        % opt E-field simulation
        curt_mask = [curt_subj curt_sesh '.nii.gz'];
        if isunix;  path_to_target_dot_msh = [curt_subj '/' curt_sesh '/' retro_setup.target_msh_name];
            path2curt_mask = [curt_subj '/' curt_mask];end
        if ispc;  path_to_target_dot_msh = [curt_subj '\' curt_sesh '\' retro_setup.target_msh_name];
            path2curt_mask = [curt_subj '\' curt_mask];end
        
        %check if example data and subject meshes exist
       % data_aval = exist('example_data.mat')+exist(path_to_target_dot_msh,'file')*exist(path_to_bin_msh,'file')*...
        %    exist(path_to_mr_scans,'file')*exist(path2curt_mask,'file');
        
        data_aval = exist(path_to_target_dot_msh,'file')*exist(path_to_bin_msh,'file')*...
            exist(path_to_mr_scans,'file')*exist(path2curt_mask,'file');
        if data_aval >= 16 % if +ve, exist(xyz,'file') = 2, 16 = 2^4, 
            
            Topt_s = Topt_bs;
            incor_signs = [sign(Topt_bs(sign_check_indx)) ~= sign_opt boolean([0;0;0;0])];
            Topt_s(incor_signs) = -Topt_s(incor_signs);%double-check and ensure the signs are correct
            Topt_s_tmp = convert_coord_brainsight_2_simnibs(curt_subj, curt_mask, Topt_bs(1:3,4)');
            Topt_s(1:3,4) = Topt_s_tmp;
            
            Sopt = sim_struct('SESSION');
            Sopt.fnamehead = path_to_bin_msh;%head mesh file, start_pipeling->make_simnibs3_from_2 to generate this
            Sopt.pathfem = 'retro_tmp';
            Sopt.poslist{1} = sim_struct('TMSLIST');
            Sopt.poslist{1}.pos.matsimnibs=Topt_s;
            Sopt.poslist{1}.fnamecoil = [retro_setup.coil_name retro_setup.coil_file_type];
            Sopt.poslist{1}.pos(1).distance = str2double(curt_HT); % distance between coil and scalp, put in HT
            Sopt.open_in_gmsh = retro_setup.gmsh_visualize_switch;
            
            run_simnibs(Sopt);
            
            % extract E-field data in ROI
            target_tmp_s = mesh_load_gmsh4(path_to_target_dot_msh);
            indexroi_s = find(target_tmp_s.element_data{2}.tetdata==1);
            
            if isunix;  optEFresult = mesh_load_gmsh4([Sopt.pathfem '/' curt_subj '_bin_TMS_1-0001_' retro_setup.coil_name '_scalar.msh']);end
            if ispc;  optEFresult = mesh_load_gmsh4([Sopt.pathfem '\' curt_subj '_bin_TMS_1-0001_' retro_setup.coil_name '_scalar.msh']);end
            optE_s = optEFresult.element_data{1}.tetdata;
            optnormE_s = optEFresult.element_data{2}.tetdata;
            optEoi_s = optE_s(indexroi_s,:);
            optnormEoi_s = optnormE_s(indexroi_s);
            if isrow(optnormEoi_s); optnormEoi_s = optnormEoi_s';end
            if retro_setup.E50_switch == 1 %if E50 is needed
                optnormEoi_tmp_s = sort(optnormEoi_s,'descend');
                E50_opt = optnormEoi_tmp_s(50);
            end
            
            % free disk space and end
            if isunix;  system(['rm ' Sopt.pathfem '/*']);end
            if ispc;  system(['rm ' Sopt.pathfem '\*']);end
            system(['rmdir ' Sopt.pathfem]);
            disp(['[TAP] Optimal E-field simulation completed for ' curt_subj '.']);
            
        elseif data_aval < 16 && exist('example_data.mat') && practice_switch == 23 % lack some mesh but have the example_data and is a practice run
            
        elseif data_aval < 2 % no example_data.mat nor meshes
            error('[TAP] Error: Imaging data not found at the specificed directory. For practice runs, ensure example_data.mat exists')
        end
        
        % computation for individual stim block
        for blck_indx = blck2start:blck_per_sesh
            
            % read Brainsight text files
            if isunix;  Tsnapped_bs = read_brainsight_file([retro_setup.snapped_folder '/' curt_subj '/' curt_sesh '_' retro_setup.bs_name_s{blck_indx}]);
                Toriginal_bs = read_brainsight_file([retro_setup.original_folder '/' curt_subj '/' curt_sesh '/' retro_setup.bs_name_o{blck_indx}]);end
            if ispc;  Tsnapped_bs = read_brainsight_file([retro_setup.snapped_folder '\' curt_subj curt_sesh '\' retro_setup.bs_name_s{blck_indx}]);
                Toriginal_bs = read_brainsight_file([retro_setup.original_folder '\' curt_subj curt_sesh '\' retro_setup.bs_name_o{blck_indx}]);end
            
            % delete hotspot samples in block 1 if any
            if blck_indx == 1 && retro_setup.CU_indx_all(sesh_indx_acrs_all) ~= 0
                CU_indx = retro_setup.CU_indx_all(sesh_indx_acrs_all);
                Toriginal_bs(1:CU_indx) = [];
                Tsnapped_bs(1:CU_indx) = [];
            end
            
            % HT extrusion
            for i = 1:length(Tsnapped_bs)
                curt_z_bs = Tsnapped_bs{i}(1:3,3);
                Tsnapped_bs{i}(1:3,4) = Tsnapped_bs{i}(1:3,4) + ...
                    curt_HT*[cosd(angv1v2(curt_z_bs,x_axis)),cosd(angv1v2(curt_z_bs,y_axis)),cosd(angv1v2(curt_z_bs,z_axis))]';
            end
            
            % adjust signs based on cmp between opt and exp bs records
            Toriginal_bs = sign_adj_wthn_bs(opt_sign,Toriginal_bs);
            Tsnapped_bs = sign_adj_wthn_bs(opt_sign,Tsnapped_bs);
            
            % coil position assessments
            norml_dist_dev = [];tan_dist_dev = [];norml_ang_dev = [];tan_ang_dev = [];
            for i = 1:length(Tsnapped_bs)
                norml_dist_dev(i,1) = dot(-Toriginal_bs{i}(1:3,3),Tsnapped_bs{i}(1:3,4)-Toriginal_bs{i}(1:3,4));
                tan_dist_dev(i,1) = norm(Tsnapped_bs{i}(1:3,4)-opt_center_bs);
                norml_ang_dev(i,1) = angv1v2(opt_z_bs,Tsnapped_bs{i}(1:3,3));
                curt_y_bs = Tsnapped_bs{i}(1:3,2);
                %projection of exp coil y on plane normal to opt coil z
                proj_exp = curt_y_bs - norm(curt_y_bs) * cosd( angv1v2(curt_y_bs,opt_z_bs) ) * opt_z_bs/norm(opt_z_bs);
                tan_ang_dev(i,1) = angv1v2(proj_opt,proj_exp);
                if angv1v2(Topt_bs(1:3,1),proj_exp) > 90
                    tan_ang_dev(i,1) = -tan_ang_dev(i,1);
                end
            end
            
            norml_dist_blck = blck_indx*ones(length(norml_dist_dev),1);
            tan_dist_blck = blck_indx*ones(length(tan_dist_dev),1);
            norml_ang_blck = blck_indx*ones(length(norml_ang_dev),1);
            tan_ang_blck = blck_indx*ones(length(tan_ang_dev),1);
            
            retro_data.norml_dist{sesh_indx_wthn_subj,1,subj_indx} =...
                [retro_data.norml_dist{sesh_indx_wthn_subj,1,subj_indx}; norml_dist_dev]; %col1 = data
            retro_data.norml_dist{sesh_indx_wthn_subj,2,subj_indx} =...
                [retro_data.norml_dist{sesh_indx_wthn_subj,2,subj_indx}; norml_dist_blck];%col2 = stim blck
            retro_data.tan_dist{sesh_indx_wthn_subj,1,subj_indx} =...
                [retro_data.tan_dist{sesh_indx_wthn_subj,1,subj_indx}; tan_dist_dev];
            retro_data.tan_dist{sesh_indx_wthn_subj,2,subj_indx} =...
                [retro_data.tan_dist{sesh_indx_wthn_subj,2,subj_indx}; tan_dist_blck];
            retro_data.norml_ang{sesh_indx_wthn_subj,1,subj_indx} =...
                [retro_data.norml_ang{sesh_indx_wthn_subj,1,subj_indx}; norml_ang_dev];
            retro_data.norml_ang{sesh_indx_wthn_subj,2,subj_indx} =...
                [retro_data.norml_ang{sesh_indx_wthn_subj,2,subj_indx}; norml_ang_blck];
            retro_data.tan_ang{sesh_indx_wthn_subj,1,subj_indx} =...
                [retro_data.tan_ang{sesh_indx_wthn_subj,1,subj_indx}; tan_ang_dev];
            retro_data.tan_ang{sesh_indx_wthn_subj,2,subj_indx} =...
                [retro_data.tan_ang{sesh_indx_wthn_subj,2,subj_indx}; tan_ang_blck];
            
            % record outlining samples
            retro_data.norml_dist_out_label{sesh_indx_wthn_subj,1,subj_indx} = ...
                [retro_data.norml_dist_out_label{sesh_indx_wthn_subj,1,subj_indx}; abs(norml_dist_dev)>retro_setup.dist_out_thrsld];
            retro_data.tan_dist_out_label{sesh_indx_wthn_subj,1,subj_indx} = ...
                [retro_data.tan_dist_out_label{sesh_indx_wthn_subj,1,subj_indx}; abs(tan_dist_dev)>retro_setup.dist_out_thrsld ];
            retro_data.norml_ang_out_label{sesh_indx_wthn_subj,1,subj_indx} = ...
                [retro_data.norml_ang_out_label{sesh_indx_wthn_subj,1,subj_indx}; abs(norml_ang_dev)>retro_setup.ang_out_thrsld];
            retro_data.tan_ang_out_label{sesh_indx_wthn_subj,1,subj_indx} = ...
                [retro_data.tan_ang_out_label{sesh_indx_wthn_subj,1,subj_indx}; abs(tan_ang_dev)>retro_setup.ang_out_thrsld];
            
            data_bckup = retro_data;
            data_bckup.last_saved = [curt_subj char(9) curt_sesh char(9) 'block' char(9) num2str(blck_indx)];%char(9)=/tab
            save('backup_data.mat','-V6','data_bckup')
            
            % avg with median
            x_bs = [];y_bs = [];z_bs = [];
            for i = 1:length(Tsnapped_bs)
                x_bs(i) = Tsnapped_bs{i}(1,4);
                y_bs(i) = Tsnapped_bs{i}(2,4);
                z_bs(i) = Tsnapped_bs{i}(3,4);
            end
            coil_median_bs = [median(x_bs) median(y_bs) median(z_bs)];
            dist2median = [];
            for i = 1:length(Tsnapped_bs)
                dist2median(i) = norm(coil_median_bs-Tsnapped_bs{i}(1:3,4));
            end
            [~,min_indx] = min(dist2median);
            if length(min_indx) > 1
                disp('[TAP] Warning: Multiple experimental coil setups have equally minimum distance to a coil center');
                disp('     with medians of experimental X, Y, and Z as its X, Y, and Z (Dannhauer, 2021).')
                disp('     TAP uses the first one of them. ')
            end
            Texp_rep_bs = Tsnapped_bs{min_indx(1)};
            
            if data_aval >= 16 % necessary meshes exist
                % exp E-field simulation
                %transfer matrix prep
                incor_signs = [sign(Texp_rep_bs(sign_check_indx))~=sign_opt boolean([0;0;0;0])];
                Texp_rep_bs(incor_signs) = -Texp_rep_bs(incor_signs);
                coil_center_s = convert_coord_brainsight_2_simnibs(...
                    curt_subj, curt_mask, Texp_rep_bs(1:3,4)');%coil center in simnibs space
                Texp_rep_s = Texp_rep_bs; Texp_rep_s(1:3,4) = coil_center_s;%convert T from bs space to simnibs space
                
                %simulation
                exps = sim_struct('SESSION');
                exps.fnamehead = path_to_bin_msh;
                exps.pathfem = 'retro_tmp';
                exps.poslist{1} = sim_struct('TMSLIST');
                exps.poslist{1}.pos.matsimnibs=Texp_rep_s;
                exps.poslist{1}.fnamecoil = [retro_setup.coil_name retro_setup.coil_file_type];
                exps.poslist{1}.pos(1).distance = str2double(curt_HT);
                exps.open_in_gmsh = retro_setup.gmsh_visualize_switch;
                run_simnibs(exps);
                
                %extract E-field data in ROI and free disk space
                if isunix;  expEFresult = mesh_load_gmsh4([exps.pathfem '/' curt_subj '_bin_TMS_1-0001_' retro_setup.coil_name '_scalar.msh']);end
                if ispc;  expEFresult = mesh_load_gmsh4([exps.pathfem '\' curt_subj '_bin_TMS_1-0001_' retro_setup.coil_name '_scalar.msh']);end
                expE_s = expEFresult.element_data{1}.tetdata;
                expnormE_s = expEFresult.element_data{2}.tetdata;
                expEoi_s = expE_s(indexroi_s,:);
                expnormEoi_s = expnormE_s(indexroi_s);
                if isrow(expnormEoi_s); expnormEoi_s = expnormEoi_s';end
                
                if isunix;  system(['rm ' exps.pathfem '/*']);end
                if ispc;  system(['rm ' exps.pathfem '\*']);end
                system(['rmdir ' Sopt.pathfem]);
                
                % E-field assessments
                %percentage EF mag dev
                ef_mag_dev = 100*(expnormEoi_s-optnormEoi_s)./optnormEoi_s;
                retro_data.ef_mag{sesh_indx_wthn_subj,1,subj_indx} = ...
                    [retro_data.ef_mag{sesh_indx_wthn_subj,1,subj_indx}; ef_mag_dev];
                
                retro_data.ef_mag{sesh_indx_wthn_subj,2,subj_indx} = ...
                    [retro_data.ef_mag{sesh_indx_wthn_subj,2,subj_indx}; blck_indx*ones(length(expnormEoi_s),1)];
                
                retro_data.ef_mag_out_label{sesh_indx_wthn_subj,1,subj_indx} = ...
                [retro_data.ef_mag_out_label{sesh_indx_wthn_subj,1,subj_indx}; abs(ef_mag_dev)>retro_setup.efm_out_thrsld];
                
                %MSO conversion
                if retro_setup.MSO_conversion_swtich == 1
                    curt_MSO = 1.5662*retro_setup.MSOs(sesh_indx_acrs_all)-2.3237;
                    expnormEoi_s = expnormEoi_s*curt_MSO;
                end
                retro_data.ef_mag{sesh_indx_wthn_subj,3,subj_indx} = ...
                    [retro_data.ef_mag{sesh_indx_wthn_subj,3,subj_indx}; expnormEoi_s];
                
                %EF dir dev
                ef_dir_dev = [];
                for i = 1:length(indexroi_s)
                    ef_dir_dev(i,1) = angv1v2(expEoi_s(i,:),optEoi_s(i,:));
                end
                retro_data.ef_dir{sesh_indx_wthn_subj,1,subj_indx} = ...
                    [retro_data.ef_dir{sesh_indx_wthn_subj,1,subj_indx}; ef_dir_dev];
                
                retro_data.ef_dir{sesh_indx_wthn_subj,2,subj_indx} = ...
                    [retro_data.ef_dir{sesh_indx_wthn_subj,2,subj_indx}; blck_indx*ones(length(ef_dir_dev),1)];
                
                retro_data.ef_dir_out_label{sesh_indx_wthn_subj,1,subj_indx} = ...
                    [retro_data.ef_dir_out_label{sesh_indx_wthn_subj,1,subj_indx}; abs(ef_dir_dev)>retro_setup.efd_out_thrsld];
                
                %E50 if it is needed
                if retro_setup.E50_switch == 1
                    expnormEoi_tmp_s = sort(expnormEoi_s,'descend');
                    retro_data.E50(blck_indx,sesh_indx_wthn_subj,subj_indx)...
                        = expnormEoi_tmp_s(50);
                end
                
                data_bckup = retro_data;
                data_bckup.last_saved = [curt_subj char(9) curt_sesh char(9) 'block' char(9) num2str(blck_indx)];
                save('backup_data.mat','-V6','data_bckup')
                
            elseif data_aval < 16 && exist('example_data.mat') && practice_switch == 23
                %assign example data to retro_data
                retro_data.ef_mag = example_data.ef_mag;
                retro_data.ef_dir = example_data.ef_dir;
            elseif data_aval < 2 % no example_data.mat nor meshes
                error('[TAP] Error: Imaging data not found at the specificed directory. For practice runs, ensure example_data.mat exists.')
            end
        end
        disp(['[TAP] Coil and e-field assessments completed for ' curt_subj ' ' curt_sesh '.']);
    end
end

save('retro_data.mat','-V6','retro_data')

end
