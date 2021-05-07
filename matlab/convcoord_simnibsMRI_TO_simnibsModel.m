function simnibsModel = simnibsMRI2simnibsModel(subject, mask_name, coordinates)

nii = load_nii([subject '/' subject '_T1fs_conform.nii.gz']);
rascorner = [nii.hdr.hist.qoffset_x nii.hdr.hist.qoffset_y nii.hdr.hist.qoffset_z];
dim=nii.original.hdr.dime.dim(2:4);
offset = [-dim(1)/2 dim(2)/2 dim(3)/2] + rascorner; %simnibs model (0,0,0) is at the volume center, flip all 3 axis
offset(:,3)=-offset(:,3);
simnibsModel = coordinates + offset;
simnibsModel(:,1)=-simnibsModel(:,1);
simnibsModel(:,2)=-simnibsModel(:,2);





