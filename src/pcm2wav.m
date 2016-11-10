function pcm2wav(DIR_wavout, p)

fout_read = fopen(DIR_wavout, 'rb');
out_full = fread(fout_read,inf, 'int16');
% [fname,ext]=strtok(DIR_wavout,'.');
% wavname = [fname,'.wav'];

wavname = DIR_wavout;
out_full = out_full./32767;
wavwrite(out_full, p.fs, 16, wavname);
fclose('all');

end