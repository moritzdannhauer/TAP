clear all;
load /Users/rothermichk17/Pipeline/matlab/tmp;
load /Users/rothermichk17/Pipeline/matlab/tmp2;
org=roi;
[tri, nhat, inner_node]=surftri(roi.node',roi.cell');
clear roi;
roi.node=org.node;
roi.face=tri';
roi.field=zeros(1,length(tri));

TargetingNavigator('Brain',brain,'Scalp',scalp,'ROI_size',tms_opt.target_size,'ROI',roi,'ROI_center',...
tms_opt.target,'ROI_normal',-tms_opt.target_direction,'Scalp_normal',scalp_normal);
global abc;
disp(num2str(abc));



