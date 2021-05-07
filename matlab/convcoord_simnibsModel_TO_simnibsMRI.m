function simnibsModel = simnibsModel2simnibsMRI(subject, mask_name, coordinates)

coordinates(:,1)=-coordinates(:,1);
coordinates(:,2)=-coordinates(:,2);
nii = load_nii([subject '/' subject '_T1fs_conform.nii.gz']);
rascorner = [nii.hdr.hist.qoffset_x nii.hdr.hist.qoffset_y nii.hdr.hist.qoffset_z];
dim=nii.original.hdr.dime.dim(2:4);
offset = [-dim(1)/2 dim(2)/2 dim(3)/2] + rascorner; 
offset(:,3)=-offset(:,3);
simnibsModel = coordinates - offset;

