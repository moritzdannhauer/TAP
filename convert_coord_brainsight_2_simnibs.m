function simnibs = convert_coord_brainsight_2_simnibs(subject, mask, brainsight_coord)

p=brainsight_coord;
nii = load_nii([subject '/' mask]);
rascorner = [nii.hdr.hist.qoffset_x nii.hdr.hist.qoffset_y nii.hdr.hist.qoffset_z];
p=[p(:,1)-rascorner(1) p(:,2)-rascorner(2) p(:,3)-abs(rascorner(3))];
p=[-p(:,1) -p(:,2) p(:,3)];
nii = load_nii([subject '/' subject '_T1fs_conform.nii.gz']); %this is a simnibs-generated file
rascorner = [nii.hdr.hist.qoffset_x nii.hdr.hist.qoffset_y nii.hdr.hist.qoffset_z];
ras2tkreg = [128 -128 -128] - rascorner;
simnibs = p + ras2tkreg;


