function s_buff=synth_ifft_buff(TF_mag, TF_phase, sz, fftlen, win, preemph, DCbin_back, pow)


i = 1;
[freq_num,frame_num] = size(TF_mag);
%s_out = zeros(frame_num*shift + sz-shift,1);
s_buff = zeros(sz, frame_num);

%Gives zeros to DCbin 
TF_mag(1:DCbin_back,:) = 0;
TF_mag = TF_mag .^ (1/pow);
while i < frame_num+1
    if freq_num == fftlen
        TF = TF_mag(:,i);
    else    
        TF_mag_sym = [TF_mag(:,i);flipud(TF_mag(2:fftlen/2,i))];
        if length(TF_phase) == fftlen
            TF_phase_sym = TF_phase;
        else
            TF_phase_sym = [TF_phase(:,i);-1*flipud(TF_phase(2:fftlen/2,i))];
        end
        TF = TF_mag_sym .* exp(sqrt(-1).*TF_phase_sym);
    end
    s_proc = real(ifft(TF));
    s_proc = s_proc(1:sz);
    
    
    s_proc = s_proc .* win;
    %de-emphasis
    s_proc = filter(1, [1 -preemph], s_proc);
    
    s_buff(:,i) = s_proc;
    %s_out(size_crnt : size_crnt-1 + sz) = s_out(size_crnt : size_crnt-1 + sz) + s_proc;
    
    %size_crnt = size_crnt + shift;
    i = i+1;
end;