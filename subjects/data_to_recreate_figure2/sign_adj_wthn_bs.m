function Texp_bs = sign_adj_wthn_bs(opt_sign,Texp_bs)
% Texp_bs = sign_adj_wthn_bs(opt_sign,Texp_bs)

for smp_indx = 1:length(Texp_bs) %sample index
    for el_indx = 1:11 %matrix element index
        sign_tmp(smp_indx,el_indx) = sign(Texp_bs{smp_indx}(el_indx));
    end
end

flip_indx = [];
for el_indx = 1:11
    curt_signs = sign_tmp(:,el_indx);
%     if mode(curt_signs) == 0 % if most samples have 0 for this element
%         disp('[TAP] Warning: Experimental coil locations provided are abnormal, check the text files before continuing.')
%         pause
%     end
    
    % check if this element scatters around 0
    num_of_pos = length(find(curt_signs == 1));
    num_of_neg = length(find(curt_signs == -1));
    curt_el = [];
    for smp_indx = 1:length(Texp_bs)
        curt_el(smp_indx) = Texp_bs{smp_indx}(el_indx);
    end
    
    if abs(num_of_pos - num_of_neg) < .2*length(Texp_bs)...
            && median(abs(curt_el)) < 0.05 %if this element does scatter around 0
        % do nothing, not flipping this element
    elseif sign(opt_sign(el_indx)) ~= mode(curt_signs)
        flip_indx = [flip_indx;el_indx];
    end
    
end

if ~isempty(flip_indx)
    for i = 1:length(Texp_bs)
        Texp_bs{i}(flip_indx) = -Texp_bs{i}(flip_indx);
    end
end

end