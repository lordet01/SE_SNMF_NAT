function NTF_sep_event_RT(path_in, path_event, path_noise, path_denoise, B_Mel_x, B_Mel_d, B_DFT_x, B_DFT_d, p)

ch = p.ch;

%multichannel file I/O initialization
for j = 1:ch
    fin(j) = fopen([path_in(j,:)],'rb');
end

for j = 1:ch
    for i = 1:p.EVENT_NUM
        [name, ext] =strtok(path_event(j,:), '.');
        fevent(i,j) = fopen([name,'_',num2str(i),ext],'wb');
    end
    for i =1:p.NOISE_NUM
        [name, ext] =strtok(path_noise(j,:), '.');
        fnoise(i,j) = fopen([name,'_',num2str(i),ext],'wb');
    end
end

for j = 1:ch
    fdenoise(j) = fopen([path_denoise(j,:)],'wb');
end

frame_len = p.framelength;
frame_shift = p.frameshift;

if p.adapt_train_N
    %% Load Recent noise dictionary if it exists
    B_idx=fopen('B_D_u.mat');
    if B_idx > -1
        try
            load('B_D_u.mat');
        catch
            
        end
    end
end
% contour(B_DFT_d);

%% Buffer initialization
if strcmp(p.NMF_algorithm,'NTF') || strcmp(p.NMF_algorithm,'PMWF')  
    g=init_buff_NTF(B_Mel_x, B_Mel_d, B_DFT_x, B_DFT_d, p);
elseif strcmp(p.NMF_algorithm,'SNMF')
    g=init_buff(B_Mel_x, B_Mel_d, B_DFT_x, B_DFT_d, p);
end

y = zeros(ch, frame_len, 1);
d_hat = zeros(p.NOISE_NUM, ch, frame_len, 1);
x_hat = zeros(p.EVENT_NUM, ch, frame_len, 1);
x_tilde = zeros(ch, frame_len,1);

%% Wav processing
header = zeros(ch, 22);
for j = 1:ch
    header(j, :)=fread(fin(j), 22, 'int16'); %Skip wav header parts
end
% end

l = 1;
% g.Fm = 0;
% g.THR_NPD = 0;
% g.THR_cnt = 1;
% BETA_blk_tmp = g.BETA_blk;
cnt_residue = 0;
s_in_sub = zeros(ch, frame_shift);
while (1)
    
    %Check eof
    [~, len] = fread(fin(1), frame_shift, 'int16');
    
    if cnt_residue > p.delay
        break;
    end
    
    if len ~= frame_shift
        cnt_residue = cnt_residue + 1;
        y = zeros(ch, frame_len, 1);
    else
        fseek(fin(1),-2*frame_shift,0); %Rewind file pointer moved by eof check
        for j = 1:ch
            [s_in_sub(j,:), ~] = fread(fin(j), frame_shift, 'int16');
        end
        
        %Frame_wise queing
        y(:, 1:frame_len-frame_shift) = y(:, frame_shift+1:frame_len);
        y(:, frame_len-frame_shift+1:frame_len) = s_in_sub;
    end
    
    %Give initial noise flag for basis update
    if l <= p.init_N_len ...
            % && l >= 2*p.Splice + 1
        g.W = 1;
%         g.BETA_blk = g.BETA_blk_init;
    else
        g.W = 0;
%         g.BETA_blk = BETA_blk_tmp;
    end
    
    %Put frame-wise algorithm here
    if strcmp(p.NMF_algorithm,'NTF')
%         [e_est_frame,n_est_frame, d_frame, g] = bntf_sep_event_RT_ASP(y, g, p);
    elseif strcmp(p.NMF_algorithm,'PMWF')
%         [e_est_frame,n_est_frame, d_frame, g] = bntf_sep_event_RT_PMWF(y, g, p);
    elseif strcmp(p.NMF_algorithm,'SNMF')
%         [e_est_frame,n_est_frame, d_frame, g] = bnmf_sep_event_RT(y, g, p);
        [e_est_frame,n_est_frame, d_frame, g] = bnmf_sep_event_RT_IS16(y, l, g, p);
    end
    
    for j =1:ch
        if l > p.delay
%             for i = 1: p.NOISE_NUM
%                 d_hat(i,j,1:frame_len-frame_shift) = d_hat(i,j,frame_shift+1:frame_len);
%                 d_hat(i,j,frame_len-frame_shift+1:frame_len) = zeros(frame_shift,1);
%                 d_hat(i,j,:) = d_hat(i,j,:) + n_est_frame(i,j,:);
%                 s_noise_sub = d_hat(i,j,1:frame_shift);
%                 fwrite(fnoise(i,j), s_noise_sub, 'int16');
%             end
            
%             for i = 1: p.EVENT_NUM
%                 x_hat(i,j,1:frame_len-frame_shift) = x_hat(i,j,frame_shift+1:frame_len);
%                 x_hat(i,j,frame_len-frame_shift+1:frame_len) = zeros(frame_shift,1);
%                 x_hat(i,j,:) = x_hat(i,j,:) + e_est_frame(i,j,:);
%                 s_event_sub = x_hat(i,j,1:frame_shift);
%                 fwrite(fevent(i,j), s_event_sub, 'int16');
%             end
            x_tilde(j,1:frame_len-frame_shift) = x_tilde(j,frame_shift+1:frame_len);
            x_tilde(j,frame_len-frame_shift+1:frame_len) = zeros(frame_shift,1);
            x_tilde(j,:) = x_tilde(j,:) + d_frame(j,:);
            fwrite(fdenoise(j), x_tilde(j,1:frame_shift), 'int16');
        end
    end
    l = l + 1;
end

%% Save updated noise basis for next use
B_DFT_d=g.B_DFT_d;
B_Mel_d=g.B_Mel_d;
save('B_D_u.mat','B_DFT_d','B_Mel_d');

fclose('all');
% % %%Convert PCM to wav
% % for j = 1:ch
% %     for i = 1:p.EVENT_NUM
% %         [name, ext] =strtok(path_event(j,:), '.');
% %         pcm2wav([name,'_',num2str(i),ext],p);
% %     end
% %     for i =1:p.NOISE_NUM
% %         [name, ext] =strtok(path_noise(j,:), '.');
% %         pcm2wav([name,'_',num2str(i),ext],p);
% %     end
% % end

for j = 1:ch
    pcm2wav(path_denoise(j,:),p);
end

% % %Record noise detection flag
% % if strcmp(p.NMF_algorithm,'NTF')
% %     W_log = W_log(2:end);
% %     lamda_log = lamda_log(2:end);
% %     for j = 1:ch
% %         save([path_denoise(j,:),'.mat'], 'W_log', 'lamda_log'); %Chop off first dummy 0
% %     end
% % end
end