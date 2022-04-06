function vs = group_violin(group_dist,vs,num_in_group,num_of_group)
% vs = group_violin(group_dist,vs,num_in_group,num_of_group)

for i = 2:num_of_group %for violin groups
    for j = (i-1)*num_in_group+1:i*num_in_group %for violins in a group
        vs(1,j).ScatterPlot.XData = vs(1,j).ScatterPlot.XData+(i-1)*group_dist;
        vs(1,j).ViolinPlot.XData = vs(1,j).ViolinPlot.XData+(i-1)*group_dist;
        vs(1,j).BoxPlot.XData = vs(1,j).BoxPlot.XData+(i-1)*group_dist;
        vs(1,j).WhiskerPlot.XData = vs(1,j).WhiskerPlot.XData+(i-1)*group_dist;
        vs(1,j).MedianPlot.XData = vs(1,j).MedianPlot.XData+(i-1)*group_dist;
    end
end
end