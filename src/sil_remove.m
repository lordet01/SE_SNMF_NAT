function s_out = sil_remove(s_in)

%Perform simple VAD to exclude silence part
p.bg_len = 0.01 * p.fs; %50ms
p.min_voiced_len = 0.5 * p.fs; %0.3s
p.min_unvoiced_len = 0.4 * p.fs; %0.3s
p.thr = 0.7;
bg_mean = 0; %initial condition

s_left = s_in;

[~,v_start,v_end, bg_mean] = vadenergy(s_left, bg_mean, p);
s_chop = s_left(v_start:v_end,1);
s_left = s_left(v_end+1:end,1);

end