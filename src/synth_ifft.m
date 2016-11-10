function s_out=synth_ifft(TF_mag, TF_phase, sz, shift, fftlen, win, preemph, pow)


size_crnt = 1;
i = 1;
frame_num = size(TF_mag,2);
s_out = zeros(frame_num*shift + sz-shift,1);
TF_mag = TF_mag ./ pow;
while i < frame_num+1
    TF_mag_sym = [TF_mag(:,i);flipud(TF_mag(2:fftlen/2,i))];
    TF_phase_sym = [TF_phase(:,i);-1*flipud(TF_phase(2:fftlen/2,i))];
    TF = TF_mag_sym .* exp(sqrt(-1).*TF_phase_sym);
    s_proc = real(ifft(TF));
    s_proc = s_proc(1:sz);
    
    s_proc = s_proc .* win;
    
    %de-emphasis
    s_proc = filter(1, [1 -preemph], s_proc);
    s_out(size_crnt : size_crnt-1 + sz) = s_out(size_crnt : size_crnt-1 + sz) + s_proc;
    
    size_crnt = size_crnt + shift;
    i = i+1;
end;