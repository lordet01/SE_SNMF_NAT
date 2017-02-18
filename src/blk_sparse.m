function [Q, r_blk_out]=blk_sparse(X, D, r_blk, l, p)

[K,~] = size(X);
% blk_gapN2 = (p.blk_gap-1)/2;
% %Unroll the splice
% e = reshape(e, [p.fftlength/2+1, p.Splice*2+1]);
% d = reshape(d, [p.fftlength/2+1, p.Splice*2+1]);

    % Method1: Q. Hoyer
%     SNR_local = (X.^p.pow) ./ max(X.^p.pow + D.^p.pow, p.nonzerofloor);
    SNR_local = (X.^p.pow) ./ max(D.^p.pow, p.nonzerofloor);
    SNR_local = max(SNR_local, p.regul);
%     disp(SNR_local);

r_blk_out = [r_blk(:,2:p.P_len_l) SNR_local];

Q = [zeros(p.DCbin,1); 0.1*ones(K - p.DCbin,1)];
n = p.P_len_l .* p.P_len_k;
P_len_k2 = floor(p.P_len_k*0.5);
if l> p.P_len_l
    for k = P_len_k2+p.blk_gap + 1 : p.blk_gap : K - P_len_k2
        b = reshape(r_blk_out(k-P_len_k2+1 : k + P_len_k2, :), n, 1);
        %         b = b + abs(min(b));
        b = b ./ max(b);
        SNR_l1 = sum(b, 1);
        SNR_l2 = sqrt(sum(b.^2, 1));
        P_tmp = (sqrt(n) - SNR_l1 ./ SNR_l2) / (sqrt(n)-1); %Method1: Q. Hoyer, 2004
        %         Q(k) = (p.alpha_p) * Q(k-1) + (1-p.alpha_p) * P_tmp .^ p.kappa;
        P_val = (p.alpha_p) * Q(k-1) + (1-p.alpha_p) * P_tmp;
        Q(k - P_len_k2:k + P_len_k2) = (p.beta_p) * Q(k - p.blk_gap - P_len_k2:k - p.blk_gap + P_len_k2) ...
            + (1-p.beta_p) * P_val;
    end
    Q(1:p.P_len_k-1) = Q(p.P_len_k+p.DCbin);
end
%     if mean(Q) < p.Q_sig_c
%         Q = zeros(size(Q));
%     end
Q = sigmf(Q, [p.Q_sig_a p.Q_sig_c]);
Q(1:p.DCbin) = 0;
end
