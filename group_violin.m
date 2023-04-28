function vs = group_violin(group_dist,vs,num_in_group,num_of_group)
% vs = group_violin(group_dist,vs,num_in_group,num_of_group)

for i = 2:num_of_group
    for j = (i-1)*num_in_group+1:i*num_in_group
        vs(1,j).ScatterPlot.XData = vs(1,j).ScatterPlot.XData+(i-1)*group_dist;
        vs(1,j).ViolinPlot.XData = vs(1,j).ViolinPlot.XData+(i-1)*group_dist;
        vs(1,j).BoxPlot.XData = vs(1,j).BoxPlot.XData+(i-1)*group_dist;
        vs(1,j).WhiskerPlot.XData = vs(1,j).WhiskerPlot.XData+(i-1)*group_dist;
        vs(1,j).MedianPlot.XData = vs(1,j).MedianPlot.XData+(i-1)*group_dist;
    end
end
end
% for i = num_in_group+1:2*num_in_group
%     vs_in(1,i).ScatterPlot.XData = vs_in(1,i).ScatterPlot.XData+group_dist;
%     vs_in(1,i).ViolinPlot.XData = vs_in(1,i).ViolinPlot.XData+group_dist;
%     vs_in(1,i).BoxPlot.XData = vs_in(1,i).BoxPlot.XData+group_dist;
%     vs_in(1,i).WhiskerPlot.XData = vs_in(1,i).WhiskerPlot.XData+group_dist;
%     vs_in(1,i).MedianPlot.XData = vs_in(1,i).MedianPlot.XData+group_dist;
% end
% for i = 2*num_in_group+1:3*num_in_group
%     vs_in(1,i).ScatterPlot.XData = vs_in(1,i).ScatterPlot.XData+2*group_dist;
%     vs_in(1,i).ViolinPlot.XData = vs_in(1,i).ViolinPlot.XData+2*group_dist;
%     vs_in(1,i).BoxPlot.XData = vs_in(1,i).BoxPlot.XData+2*group_dist;
%     vs_in(1,i).WhiskerPlot.XData = vs_in(1,i).WhiskerPlot.XData+2*group_dist;
%     vs_in(1,i).MedianPlot.XData = vs_in(1,i).MedianPlot.XData+2*group_dist;
% end
% for i = 3*num_in_group+1:4*num_in_group
%     vs_in(1,i).ScatterPlot.XData = vs_in(1,i).ScatterPlot.XData+3*group_dist;
%     vs_in(1,i).ViolinPlot.XData = vs_in(1,i).ViolinPlot.XData+3*group_dist;
%     vs_in(1,i).BoxPlot.XData = vs_in(1,i).BoxPlot.XData+3*group_dist;
%     vs_in(1,i).WhiskerPlot.XData = vs_in(1,i).WhiskerPlot.XData+3*group_dist;
%     vs_in(1,i).MedianPlot.XData = vs_in(1,i).MedianPlot.XData+3*group_dist;
% end