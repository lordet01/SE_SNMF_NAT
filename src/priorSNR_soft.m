function [Sm_out, n_Sm_out, Sm_out_Mel, n_Sm_out_Mel, NPD] = priorSNR_soft(e, d, p, g)

[n,m] = size(e);

F_order = p.F_order;
fftlen2 = round(p.fftlength/2+1);


% Sigmoid parameter
ALPHA = 0.2;
BETA = 4; %dB
%Smoothing Parameter
GAMMA = 2;

% SNR_l1 = (e) ./ (d);
e = e ./ mean(e);
d = d ./ mean(d);
SNR_l2 = 10*log10(max((e.^2) ./ (d.^2), 0.001));
% Sm = SNR_l2;
Sm =  1 ./ (1+exp(-ALPHA * (SNR_l2 - BETA))); %0~1 map using sigmoid function

%Post processing of soft mask (median, averaging)
Sm_med = medfilt2(Sm,[5,p.blk_len_sep]);
h = 1 ./ (GAMMA^2) * ones(GAMMA);
Sm_smooth = filter2(h,Sm_med);
Sm_out = Sm_smooth;
n_Sm_out = 1 - Sm_out;
for k = 1 : p.Splice * 2 + 1
    n_Sm_out(1+(k-1)*fftlen2 : p.DCbin+(k-1)*fftlen2, :) = zeros(p.DCbin,1) + p.nonzerofloor;
end

Sm_out = Sm_out + abs(min(Sm_out));
Sm_out = Sm_out ./ max(Sm_out);

%Generate Mel-domain masks
for k = 1 : p.Splice * 2 + 1
    Sm_out_Mel(1+(k-1)*F_order : k*F_order, :) = ...
        g.melmat*shiftdim(Sm_out(1+(k-1)*fftlen2 : k*fftlen2, :));
    n_Sm_out_Mel(1+(k-1)*F_order : k*F_order, :) = ...
        g.melmat*shiftdim(n_Sm_out(1+(k-1)*fftlen2 : k*fftlen2, :));    
end

Hd_val = mean(mean(Sm_out(p.DCbin+1:end, floor(p.blk_len_sep/4)+1 : p.blk_len_sep - floor(p.blk_len_sep/4)),1));

if Hd_val < p.Hd_THR;
    NPD = 1;
else
    NPD = 0;
end


% %Reroll the splice
% DM = reshape(DM, [(p.fftlength/2+1) * (p.Splice*2+1), 1]);
