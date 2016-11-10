function [Q, r_blk_out]=blk_sparse(X, D, r_blk, l, p)

[K,~] = size(X);
blk_gapN2 = (p.blk_gap-1)/2;
% %Unroll the splice 
% e = reshape(e, [p.fftlength/2+1, p.Splice*2+1]);
% d = reshape(d, [p.fftlength/2+1, p.Splice*2+1]);

% Method1: Q. Hoyer
SNR_local = (X) ./ max(D, p.nonzerofloor);
% SNR_local = 10*log10(X ./ max(D, p.nonzerofloor));
SNR_local = SNR_local ./ max(SNR_local);

r_blk_out = [r_blk(:,2:p.P_len_l) SNR_local];

Q = [zeros(p.DCbin,1); 0.1*ones(K - p.DCbin,1)];
n = p.P_len_l .* p.P_len_k;
P_len_k2 = floor(p.P_len_k*0.5);
if l> p.P_len_l
    for k = P_len_k2+p.DCbin : p.blk_gap : K - P_len_k2
        b = reshape(r_blk_out(k-P_len_k2+1 : k + P_len_k2, :), n, 1);
%         b = b + abs(min(b));
%         b = b ./ max(b);
        SNR_l1 = sum(b, 1);
        SNR_l2 = sqrt(sum(b.^2, 1));
        P_tmp = (sqrt(n) - SNR_l1 ./ SNR_l2) / (sqrt(n)-1); %Method1: Q. Hoyer, 2004
%         Q(k) = (p.alpha_p) * Q(k-1) + (1-p.alpha_p) * P_tmp .^ p.kappa;
        P_val = (p.alpha_p) * Q(k-1) + (1-p.alpha_p) * P_tmp;
        Q(k-blk_gapN2:k) = P_val;
        Q(k:k+blk_gapN2) = P_val;
    end
    Q(1:p.P_len_k-1) = Q(p.P_len_k+p.DCbin);
end
% Q = Q ./ max(Q);
% Q = Q .* p.nu;
Q(1:p.DCbin) = 0;
end
