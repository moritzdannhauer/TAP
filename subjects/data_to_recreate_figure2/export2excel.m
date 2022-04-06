function [coil_table,ef_table] = export2excel(retro_data)
%numsesh*2(data+group label)*numsubj
num_subj = size(retro_data.norml_dist,3);
num_sesh = size(retro_data.norml_dist,1);
data_names = fieldnames(retro_data);
num_data_type = length(data_names)/2;

subj_names = {};
for i = 1:num_subj
    subj_names = [subj_names;['subj' num2str(i)]];
end
sesh_names = {};
for i = 1:num_sesh
    sesh_names = [sesh_names;['session' num2str(i)]];
end

title_common = {'subject_number';'session_name';'block_number'};
title_coil = [title_common;data_names(1:4);data_names(7:10)];
title_ef = [title_common;data_names(5:6);data_names(11:12)];

%sort data into single cols
col_data = cell(length(data_names),2);
for i = 1:length(data_names)
    curt_dev_data_col = []; blck_label = [];
    curt_dev_data = retro_data.(data_names{i});
    if i <= num_data_type %dev data
        for subj_indx = 1:num_subj
            for sesh_indx = 1:num_sesh
                curt_dev_data_col = [curt_dev_data_col;curt_dev_data{sesh_indx,1,subj_indx}];
                blck_label = [blck_label;curt_dev_data{sesh_indx,2,subj_indx}];
            end
        end
        col_data{i,1} = curt_dev_data_col;
        col_data{i,2} = blck_label;
        
    elseif i > num_data_type %outlier labels
        for subj_indx = 1:num_subj
            for sesh_indx = 1:num_sesh
                curt_dev_data_col = [curt_dev_data_col;curt_dev_data{sesh_indx,1,subj_indx}];
            end
        end
        col_data{i,1} = curt_dev_data_col;
    end
    
end


cell4coil = cell(length(col_data{1}),11);
%fill in the subj and sesh cols
subj_col_start = 1;
for subj_indx = 1:num_subj
    curt_data = retro_data.norml_dist(:,1,subj_indx);
    curt_subj_length = 0;
    for i = 1:length(curt_data); curt_subj_length = curt_subj_length+length(curt_data{i});end
    subj_col_end = subj_col_start + curt_subj_length - 1;
    for subj_col_indx = subj_col_start:subj_col_end
        cell4coil{subj_col_indx,1} = subj_names{subj_indx};
    end
    
    sesh_col_start = subj_col_start;
    for sesh_indx = 1:num_sesh
        curt_sesh_length = length(retro_data.norml_dist{sesh_indx,1,subj_indx});
        sesh_col_end = sesh_col_start + curt_sesh_length - 1;
        for sesh_col_indx = sesh_col_start:sesh_col_end
            cell4coil{sesh_col_indx,2} = sesh_names{sesh_indx};
        end
        sesh_col_start = sesh_col_end+1;
    end
    
    subj_col_start = subj_col_end+1;
end
%fill in block label, coil data, and outlier label
blck_label = col_data{1,2};
for row_indx = 1:length(cell4coil)
    cell4coil{row_indx,3} = blck_label(row_indx);
end
for dev_type = 1:4
    curt_dev_data = col_data{dev_type,1};
    for row_indx = 1:length(cell4coil)
        cell4coil{row_indx,dev_type+3} = curt_dev_data(row_indx);
    end
end
for out_label = 1:4
    curt_out_label = col_data{out_label+num_data_type,1};
    for row_indx = 1:length(cell4coil)
        cell4coil{row_indx,out_label+7} = curt_out_label(row_indx);
    end
end

cell4ef = cell(length(col_data{1}),7);
%fill in the subj and sesh cols
subj_col_start = 1;
for subj_indx = 1:num_subj
    curt_data = retro_data.ef_mag(:,1,subj_indx);
    curt_subj_length = 0;
    for i = 1:length(curt_data); curt_subj_length = curt_subj_length+length(curt_data{i});end
    subj_col_end = subj_col_start + curt_subj_length - 1;
    for subj_col_indx = subj_col_start:subj_col_end
        cell4ef{subj_col_indx,1} = subj_names{subj_indx};
    end
    
    sesh_col_start = subj_col_start;
    for sesh_indx = 1:num_sesh
        curt_sesh_length = length(retro_data.ef_mag{sesh_indx,1,subj_indx});
        sesh_col_end = sesh_col_start + curt_sesh_length - 1;
        for sesh_col_indx = sesh_col_start:sesh_col_end
            cell4ef{sesh_col_indx,2} = sesh_names{sesh_indx};
        end
        sesh_col_start = sesh_col_end+1;
    end
    
    subj_col_start = subj_col_end+1;
end
%fill in block label, ef data, and outlier label
blck_label = col_data{5,2};
for row_indx = 1:length(cell4ef)
    cell4ef{row_indx,3} = blck_label(row_indx);
end
for dev_type = 1:2
    curt_dev_data = col_data{dev_type+4,1};
    for row_indx = 1:length(cell4ef)
        cell4ef{row_indx,dev_type+3} = curt_dev_data(row_indx);
    end
end
for out_label = 5:6
    curt_out_label = col_data{out_label+num_data_type,1};
    for row_indx = 1:length(cell4ef)
        cell4ef{row_indx,out_label+1} = curt_out_label(row_indx);
    end
end

coil_table=cell2table(cell4coil,'VariableNames',title_coil);
ef_table=cell2table(cell4ef,'VariableNames',title_ef);
writetable(coil_table,'retro_coil.xlsx');
writetable(ef_table,'retro_ef.xlsx');


end