clear all;

if (ispc)
    sep='\';
elseif (ismac || isunix)
    sep='/';
end

subj='Ernie_headreco';
subjects_folder='subjects';
target_name='M1';
outputfolder='M1_retrosim'

BS_Trans_mat = read_brainsight_file([ subjects_folder sep subj sep target_name sep 'ernie_headreco_no-conform_M1_hair0.5mm.txt' ]);

S = sim_struct('SESSION');
S.fnamehead = [subjects_folder sep subj sep 'ernie_headreco_no-conform.msh']; % head mesh
S.pathfem = [ subjects_folder sep subj sep outputfolder];
S.poslist{1} = sim_struct('TMSLIST');
S.poslist{1}.fnamecoil = 'Magstim_70mm_Fig8.ccd'; % Choose a coil in the ccd-files folder
S.poslist{1}.anisotropy_type='scalar';
S.poslist{1}.aniso_maxcond=2;
S.poslist{1}.aniso_maxratio=10;
S.poslist{1}.pos(1).distance = 0; % 4 mm distance from coil surface to head surface
S.poslist{1}.pos(1).didt = 1e6; % 4 mm distance from coil surface to head surface
S.poslist{1}.pos(1).matsimnibs=BS_Trans_mat{1};
run_simnibs(S);





