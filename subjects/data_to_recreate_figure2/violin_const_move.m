function vs = violin_const_move(move_dist,vs)
% vs = violin_const_move(move_dist,vs)

for i = 1:length(vs)
    vs(1,i).ScatterPlot.XData = vs(1,i).ScatterPlot.XData+move_dist;
    vs(1,i).ViolinPlot.XData = vs(1,i).ViolinPlot.XData+move_dist;
    vs(1,i).BoxPlot.XData = vs(1,i).BoxPlot.XData+move_dist;
    vs(1,i).WhiskerPlot.XData = vs(1,i).WhiskerPlot.XData+move_dist;
    vs(1,i).MedianPlot.XData = vs(1,i).MedianPlot.XData+move_dist;
    
end
end