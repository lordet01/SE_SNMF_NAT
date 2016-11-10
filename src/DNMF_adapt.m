function B_a = DNMF_adapt(Y, D, B, p)

    %Get H_hat
    p.w_update_ind = false(p.R_x+p.R_d,1);
    p.h_update_ind = true(p.R_x+p.R_d,1);
    p.init_w = B; %given from exemplar basis as initialization
    [~, A_hat] = sparse_nmf(Y, p);    
%     
%     p.w_update_ind = true(p.cluster_buff*p.R_e,1);
%     p.h_update_ind = false(p.cluster_buff*p.R_e,1);
%     p.init_w = B(:,p.R_e+1:p.R_e+p.R_e); %given from exemplar basis as initialization
%     p.init_h = A_hat(p.R_e+1:p.R_e+p.R_e,:); %given from exemplar basis as initialization
%     [B_a, ~] = sparse_nmf(D, p);

    % Adapt noise basis
    p.w_update_ind = true(p.R_d,1);
    p.h_update_ind = false(p.R_d,1);
    p.init_w = B(:,p.R_x+1:p.R_x+p.R_d); %given from exemplar basis as initialization
    p.init_h = A_hat(p.R_x+1:p.R_x+p.R_d,:); %given from exemplar basis as initialization
    [B_a, ~] = sparse_nmf(D, p);
end