%this script will generate all subimages for Figure 2 in the manuscript and
%will perform many SimNIBS TMS simulations and may run for over 80 minutes
%on a regular laptop (e.g., MacBook i5 3.1 GHz, 8 GB RAM)
%however, we do not share the head models for the subjects on github, so
%this SimNIBS computations are omitted
clear all;
abc=tic;
%% analysis setup
retro_setup.path2simnibs='/Users/rothermichk17/Applications/SimNIBS-3.2';
retro_setup.original_folder = 'original';
retro_setup.snapped_folder = 'snapped'; % folder that contains the experiment log files (the BrainSight text files)
retro_setup.coil_name = 'Medtronic_MCF_B65';%for more coil models, https://simnibs.github.io/simnibs/build/html/documentation/coils.html
retro_setup.coil_file_type = '.ccd'; %or '.nii.gz'
retro_setup.gmsh_visualize_switch = 0; % 1 for visualizing and checking interim simulation results, 0 for automated computation
retro_setup.E50_switch = 0; % 1 for calculating E50 per stimulation block, 0 for faster computation

retro_setup.subj = {'Subj1';'Subj2';'Subj3';'Subj4';'Subj5'};% subj{N} should give the name for folder labelled subjN in Figure 1
retro_setup.sesh_names = {'PMC';'sham';'PSC'};% sesh{N} should give the name for folder labelled seshN in Figure 1, and [subj{N} sesh{N}] should form the name for folders labelled subjNseshN in Figure 1
retro_setup.opt_name = {'sSubj1_PMC_hair0.5mm';'sSubj1_PMC_hair0.5mm';'sSubj1_PSC_hair0.5mm';...
    'sSubj2_PreMotor_hair2.5mm';'sSubj2_PreMotor_hair2.5mm';'sSubj2_PSC_hair2.5mm';...
    'sSubj3_PMC_hair1.5mm';'sSubj3_PMC_hair1.5mm';'sSubj3_PSC3_hair2.5mm';...
    'sSubj4_PMC3_hair2.5mm';'sSubj4_PMC3_hair2.5mm';'sSubj4_PSC2_hair2.0mm';...
    'sSubj5_PMC2_hair2.0mm';'sSubj5_PMC2_hair3.5mm';'sSubj5_PSC_hair2.0mm'}; %names of the optimal text files in a subj1sesh1,…, subj1seshN, subj2sesh1, …, subjNseshN order
retro_setup.bs_name_o = {'Hotspot+Block1.txt';'Block2.txt';'Block3.txt';'Block4.txt'};%names of original BrainSight files (requires users to have same names and same blocks of stim in each session)
retro_setup.bs_name_s = {'B1.txt';'B2.txt';'B3.txt';'B4.txt'};%names of snapped experimental BrainSight files
retro_setup.target_msh_name = 'target.msh'; % name of the SimNIBS  generated ROI mesh
retro_setup.dist_out_thrsld = 10;%mm, custom distance outlier threshold 
retro_setup.ang_out_thrsld = 20;%deg, custom angle outlier threshold 
retro_setup.efm_out_thrsld = 50;%percent, custom EF magnitude outlier threshold
retro_setup.efd_out_thrsld = 20;%deg , custom EF direction outlier threshold
retro_setup.MSO_conversion_swtich = 1;%0 for normalized ||E||, 1 for absolute ||E|| in V/m

% The following 3 variables should be in the same subj-sesh order as opt_name
retro_setup.HTs = [0.5;0.5;0.5;...
    2.0;2.5;2.5;...
    1.5;1.5;2.5;...
    2.5;2.0;2.0;...
    2.0;3.5;2.0];%hair thicknesses
retro_setup.CU_indx_all = [141;77;20;...
    108;98;95;...
    112;98;83;...
    144;11;92;...
    87;149;87];%clean up (CU) index, start of the actual stimulation sample, put 0 if no non-stimulation sample prior to actual stimulation
retro_setup.MSOs = [53;52;47;...
    69;76;64;...
    41;39;38;...
    50;46;47;...
    58;55;45];%percentage of maximum stimulator output used in each session
retro_setup.sign_opt_all = {[-1,1,1;1,1,-1;-1,1,-1;0,0,0];[-1,1,1;1,1,-1;-1,1,-1;0,0,0];[1,-1,1;-1,-1,-1;1,-1,-1;0,0,0];...
    [1,-1,1;-1,-1,-1;1,-1,-1;0,0,0];[1,-1,1;-1,-1,-1;1,-1,-1;0,0,0];[1,-1,1;-1,-1,-1;1,-1,-1;0,0,0];...
    [-1,1,1;1,1,-1;-1,1,-1;0,0,0];[-1,1,1;1,1,-1;-1,1,-1;0,0,0];[1,-1,1;-1,-1,-1;1,-1,-1;0,0,0];...
    [-1,1,1;1,1,1;-1,1,-1;0,0,0];[-1,1,1;1,1,1;-1,1,-1;0,0,0];[1,-1,1;-1,-1,1;1,-1,-1;0,0,0];...
    [-1,1,1;1,1,-1;-1,1,-1;0,0,0];[-1,1,1;1,1,-1;-1,1,-1;0,0,0];[1,-1,1;1,-1,-1;1,1,-1;0,0,0]};%sign pattern of the optimal transfer matrices in SimNIBS space (direct output of the prospective dosing)

%% analysis
load retro_data; % add "%" to beginning of line if you got all imaging data requested from authors
%retro_data = retro_analysis(retro_setup); % ... and remove "%" from beginning of this line so retro_analysis can run

 

%% violin prep
vprep.violin_colors = {[1,0,0];[.5,.5,.5];[0,0,1]};%RGB red, gray, blue, check https://www.mathworks.com/help/matlab/ref/plot.html for other colors
vprep.prctile_line_color = [0,1,1]; %RGB cyan
vprep.violin_transparency = .7;
vprep.overall_median_color = [0,0,0];%RGB white, check https://www.mathworks.com/help/matlab/ref/plot.html for more color
vprep.subj_median_markers = ['p';'s';'^';'o';'v'];%check https://www.mathworks.com/help/matlab/ref/plot.html for more marker selection
vprep.scatter_switch = 0;%1 = scatter plot on top of violins, 0 = no scatter
vprep.violin_group_dist = 4;
vprep.individual_violin_dist = 3;

%% violin plot
plot_violins(retro_data,vprep)

%% export data to Excel
export2excel(retro_data);

toc(abc)