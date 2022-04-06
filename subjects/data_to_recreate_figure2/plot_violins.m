function figs = plot_violins(retro_data,vprep)

%global num_of_subj sesh_per_subj blck_per_sesh
%% plot preps
% load inputs
data_names = fieldnames(retro_data);
% violin plot parameters
vprep.violin_colors = {[1,0,0];[.5,.5,.5];[0,0,1]};%RGB red, gray, blue, check https://www.mathworks.com/help/matlab/ref/plot.html for other colors
vprep.prctile_line_color = [0,1,1]; %RGB cyan
vprep.violin_transparency = .7;
vprep.overall_median_color = [0,0,0];%RGB white, check https://www.mathworks.com/help/matlab/ref/plot.html for more color
vprep.subj_median_markers = ['p';'s';'^';'o';'v'];%check https://www.mathworks.com/help/matlab/ref/plot.html for more marker selection
vprep.scatter_switch = 0;%1 = scatter plot on top of violins, 0 = no scatter
vprep.violin_group_dist = 4;
vprep.individual_violin_dist = 3;
%size setups
num_of_subj = size(retro_data.(data_names{1}),3);
sesh_per_subj = size(retro_data.(data_names{1}),1);
example_data = retro_data.(data_names{1});
blck_per_sesh = length(unique(example_data{1,2,1}));

violin_width = 1.8*sesh_per_subj/3;
subj_median_size = 1300*num_of_subj/5;
overall_median_size = 1800*subj_median_size/1300;
scatter_size = 20*subj_median_size/1300;
prctile_line_width = 10*violin_width/1.8;
%plot setups
box_dist = 0.05;
num_violin_group = blck_per_sesh;
num_violin_per_group = sesh_per_subj;
violin_count = num_violin_group*num_violin_per_group;
ylabels = {'Normal coil deviation (mm)';'Tangential coil deviation (mm)';...
    ['Normal coil deviation (' char(176) ')'];['Tangential coil deviation (' char(176) ')'];...
    'ROI EF mag. deviation (%)';['ROI EF dir. deviation (' char(176) ')']};

%% violin plotting
num_fig = length(data_names)/2;
%merge outliers across metrics
out_merged=[];
for fn_indx = 1:num_fig-2
    out4plot = retro_data.(data_names{fn_indx+num_fig});
    out4violin = [];
    for sesh_indx = 1:sesh_per_subj
        for subj_indx = 1:num_of_subj
            out4violin = [out4violin;out4plot{sesh_indx,1,subj_indx}];
        end
    end
    if fn_indx == 1
        out_merged = out4violin;
    else
        out_merged(out4violin==1) = 1;
    end
end

out5plot = retro_data.(data_names{1+num_fig});
for fn_indx = 2:num_fig-2
    curt_tmp = retro_data.(data_names{fn_indx+num_fig});
    for sesh_indx = 1:sesh_per_subj
        for subj_indx = 1:num_of_subj
            out5plot{sesh_indx,1,subj_indx}(curt_tmp{sesh_indx,1,subj_indx}==1)=1;
        end
    end
end

%plot violins
for fig_indx = 1:num_fig
    %prepare data
    data2plot = retro_data.(data_names{fig_indx});
    data4violin = []; group4violin = [];
%     if fig_indx >= num_fig-2
        out4plot = retro_data.(data_names{fig_indx+num_fig});out4violin = [];
%     end
    for sesh_indx = 1:sesh_per_subj
        for subj_indx = 1:num_of_subj
            data2plot{sesh_indx,2,subj_indx} = data2plot{sesh_indx,2,subj_indx}+(sesh_indx-1)*box_dist;%space blck number out
            data4violin = [data4violin;data2plot{sesh_indx,1,subj_indx}];%minus sign so outward correction = positive deviation
            group4violin = [group4violin;data2plot{sesh_indx,2,subj_indx}];
            if fig_indx >= num_fig-2; out4violin = [out4violin;out4plot{sesh_indx,1,subj_indx}];end
        end
    end
    % remove outliers
    if fig_indx < num_fig-2
        data4violin(logical(out_merged))=[];
        group4violin(logical(out_merged))=[];
    elseif fig_indx >= num_fig-2
        data4violin(logical(out4violin))=[];
        group4violin(logical(out4violin))=[];
    end
%     if fig_indx == 1
%         out5plot=retro_data.(data_names{fig_indx+num_fig+1});
%     end
    
    %plot violins
    f = figure('WindowState','maximized');
    disp(['[TAP] ' ylabels{fig_indx} ': base violin plotting....'])
    vs = violinplot(data4violin,group4violin,'Width',violin_width);
    ylabel(ylabels{fig_indx})
    for sesh_indx = 1:sesh_per_subj
        for violin_indx = sesh_indx:sesh_per_subj:violin_count-sesh_per_subj+sesh_indx
            vs(1,violin_indx).ViolinColor = vprep.violin_colors{sesh_indx};%assign colors for each session's violins
        end
    end
    for i = 1:violin_count
        vs(1,i).ViolinAlpha = vprep.violin_transparency;
        vs(1,i).MedianPlot.SizeData = overall_median_size; %median circle size
        vs(1,i).BoxPlot.FaceColor = vprep.prctile_line_color;
        vs(1,i).BoxPlot.EdgeColor = vprep.prctile_line_color;
        vs(1,i).WhiskerPlot.Color = vprep.prctile_line_color; %set inter-percentile line color
        vs(i).BoxPlot.LineWidth = prctile_line_width;
        vs(1,i).ShowData = 0;
    end
    
    vs = group_violin(vprep.violin_group_dist,vs,num_violin_per_group,num_violin_group);%sort violins by group (stim block)
    vs = group_violin(vprep.individual_violin_dist,vs,1,violin_count);%add distance between individual violins within a group
    hold on;
    
    subj_median_dist = vs(sesh_per_subj+1).MedianPlot.XData-2;
    f;%normal_center_dev_4plot{session_indx,1,subj_indx}
    disp(['[TAP] ' ylabels{fig_indx} ': subject median plotting....'])
    for subj_indx = 1:num_of_subj
        for sesh_indx = 1:sesh_per_subj
            violin_range = (sesh_indx-1)*blck_per_sesh+1:sesh_indx*blck_per_sesh;
            data2plot{sesh_indx,1,subj_indx}(logical(out4plot{sesh_indx,1,subj_indx}))=[];
            data2plot{sesh_indx,2,subj_indx}(logical(out4plot{sesh_indx,1,subj_indx}))=[];
%             if fig_indx == 1
                data2plot{sesh_indx,1,subj_indx}(logical(out5plot{sesh_indx,1,subj_indx}))=[];
                data2plot{sesh_indx,2,subj_indx}(logical(out5plot{sesh_indx,1,subj_indx}))=[];
%             end
            vs_subj(violin_range) = violinplot(data2plot{sesh_indx,1,subj_indx},data2plot{sesh_indx,2,subj_indx});
            for i = 1:blck_per_sesh
                current_violin = (sesh_indx-1)*blck_per_sesh+i;
                vs_subj(current_violin).ScatterPlot.Marker = vprep.subj_median_markers(subj_indx);
                vs_subj(current_violin).ScatterPlot.SizeData = scatter_size;
                vs_subj(current_violin).ScatterPlot.MarkerFaceAlpha = vprep.scatter_switch;
                
                vs_subj(current_violin).MedianPlot.Marker = vprep.subj_median_markers(subj_indx);
                if vprep.subj_median_markers(subj_indx) == 'p'
                    vs_subj(current_violin).MedianPlot.SizeData = subj_median_size*1.3;
                else
                    vs_subj(current_violin).MedianPlot.SizeData = subj_median_size;
                end
                vs_subj(current_violin).MedianPlot.MarkerFaceColor = vprep.overall_median_color;
                vs_subj(current_violin).MedianPlot.MarkerEdgeAlpha = 1;
                vs_subj(current_violin).MedianPlot.MarkerFaceAlpha = 1;
                vs_subj(current_violin).ViolinPlot.FaceAlpha=0;
                vs_subj(current_violin).ViolinPlot.EdgeAlpha=0;
                vs_subj(current_violin).BoxPlot.EdgeAlpha=0;
                vs_subj(current_violin).BoxPlot.FaceAlpha=0;
                vs_subj(current_violin).WhiskerPlot.Visible=0;
            end
            vs_subj(violin_range) = group_violin(subj_median_dist,vs_subj(violin_range),1,num_violin_group);
            single_violin_dist = (sesh_indx-1)*vprep.individual_violin_dist+sesh_indx-1;
            vs_subj(violin_range) = violin_const_move(single_violin_dist,vs_subj(violin_range));
            hold on
        end
    end
    for sesh_indx = 1:sesh_per_subj
        for violin_indx = sesh_indx : sesh_per_subj : violin_count-sesh_per_subj+sesh_indx
            vs_subj(1,violin_indx).ViolinColor = vprep.violin_colors{sesh_indx};%assign colors for each session's violins
        end
    end
    
    %darken the midline and overall median
    disp(['[TAP] ' ylabels{fig_indx} ': overall medians and interquartile line plotting....'])
    vs2 = violinplot(data4violin,group4violin,'Width',violin_width);
    for i = 1:violin_count
        vs2(i).ShowData = 0; %turn off scatter
        vs2(i).ViolinPlot.Visible = 0; %turn off violin
        vs2(1,i).MedianPlot.SizeData = overall_median_size; % median circle size
        vs2(i).BoxPlot.Visible = 0;
        vs2(i).MedianPlot.MarkerFaceAlpha=0;
    end
    %add distance inbetween blocks
    disp(['[TAP] ' ylabels{fig_indx} ': violin spacing....'])
    vs2 = group_violin(vprep.violin_group_dist,vs2,num_violin_per_group,num_violin_group);% vs = group_violin(group_dist,vs,num_in_group,num_of_group)
    vs2 = group_violin(vprep.individual_violin_dist,vs2,1,violin_count);xticks([]);
end
end