function [S_mag, S_phase] = stft_fft( s, sz, shift, fftlen, DCbin, win, preemph)
% Short-time FFT
%
% Inputs:
%  s   input time series (must be row vector), or input complex spectrogram (DC to Nyquist)
%  sz  frame size
%  hp  hop size in samples
%  pd  pad size in samples
%  w   window to use (function name of data vector)
%  ll  highest frequency index to return
%
% Output:
%  f   complex STFT output (only DC to Nyquist components), or time series resynthesis

size_crnt = 1;
i = 1;
frame_num = floor(length(s)/shift);
S_mag = zeros(fftlen/2+1,frame_num);
S_phase = zeros(fftlen/2+1,frame_num);

while size_crnt < length(s) - fftlen
    s_frame = filter([1 -preemph], 1, s(size_crnt:size_crnt-1 + sz));
    s_frame = win.* s_frame;

    s_frame_pad = [s_frame; zeros(fftlen-sz,1)];
    S_frame = fft(s_frame_pad);
    S_frame_mag = abs(S_frame(1:floor(fftlen/2)+1));
    S_phase(:,i) = angle(S_frame(1:floor(fftlen/2)+1));
    
    %Set zero for LPF effect
    S_frame_mag(1:DCbin) = zeros(DCbin,1) + 0.000001;
    
    S_mag(:,i) = S_frame_mag;

    size_crnt = size_crnt + shift;
    i = i+1;
end;