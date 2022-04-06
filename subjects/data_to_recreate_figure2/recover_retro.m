function retro_data = recover_retro(retro_setup)

global recover_switch

% remove simnibs's pathfem if exists
if exist('retro_tmp','dir')
    if isunix;  system('rm retro_tmp/*');end
    if ispc;  system('rm retro_tmp\*');end
    system('rmdir retro_tmp');
end
    
recover_switch = 1;

retro_data = retro_analysis(retro_setup);

end