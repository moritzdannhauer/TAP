function coord = get_target(subject, mask_from_noreen)
%addpath('/Users/moritz/Documents/MATLAB/')
%subject='Subj2';
%mask_from_noreen='bia5_Subj2_003_expectedmask.nii.gz';

nii = load_nii([subject '/' mask_from_noreen]);
rascorner = [nii.hdr.hist.qoffset_x nii.hdr.hist.qoffset_y nii.hdr.hist.qoffset_z];
[x y z] = ind2sub(size(nii.img),find(nii.img==1));
p=[257-x-abs(rascorner(1)) 257-y-abs(rascorner(2)) z-abs(rascorner(3))]; 
p=[p(:,1)-1 p(:,2)-1 p(:,3)-1];
p=[-p(:,1) -p(:,2) p(:,3)]; % x and y are switched in seg3d
nii = load_nii([subject '/' subject '_T1fs_conform.nii.gz']);
rascorner = [nii.hdr.hist.qoffset_x nii.hdr.hist.qoffset_y nii.hdr.hist.qoffset_z];
ras2tkreg = [128 -128 -128] - rascorner;
tkregtarget = p + ras2tkreg;

if( any(size(tkregtarget)==[1 3] | size(tkregtarget)==[3 1])==1 )
 coord = tkregtarget;
else
 coord = mean(tkregtarget);
end
