addpath('src');
addpath('settings');

initial_setting_SNMF_NAT;

% fname = 'M04_423C020A_STR.CH6';
fname = 'white_0dB_noisy';
ftype = '.wav';
path_in = ['wav/',fname,ftype];
path_denoise = ['wav/',fname,'_out_v3.9_18_pc',ftype];

% % Load Event basis
% load(['basis/Clean_train_car_alarm/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_20.mat']);
% B_DFT_x = B_DFT_sub; B_Mel_x = B_Mel_sub;
% load(['basis/Clean_train_scream/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_20.mat']);
% B_DFT_x = [B_DFT_x, B_DFT_sub]; B_Mel_x = [B_Mel_x, B_Mel_sub];
% load(['basis/Clean_train_TIMIT_test/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_100.mat']);
% B_DFT_x = [B_DFT_x, B_DFT_sub]; B_Mel_x = [B_Mel_x, B_Mel_sub];

% preemph = 0.92
% load(['basis/Clean_train_TIMIT_test/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_100.mat']);

% preemph = 0.0
load(['basis/Clean_train_TIMIT_test/','TASLP_Splice0-SNMF_p2_DD0','/R_100.mat']);
B_DFT_x = [B_DFT_sub]; B_Mel_x = [B_Mel_sub];

%Load Noise basis
% load(['basis/CHiME3_bgn_ch6/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_20.mat']);
% B_DFT_d = B_DFT_sub; B_Mel_d = B_Mel_sub;

% preemph = 0.92
% load(['basis/CHiME3_bgn_ch6/','IS16_Splice0-SNMF-preemph_p2_DD0','/R_100.mat']);

% preemph = 0.0
load(['basis/CHiME3_bgn_ch6/','TASLP_Splice0-SNMF_p2_DD0','/R_100.mat']);
B_DFT_d = B_DFT_sub; B_Mel_d = B_Mel_sub;

if p.fs == 8000
    B_DFT_x = B_DFT_x(1:257,:);
    B_DFT_d = B_DFT_d(1:257,:);
end

if size(B_DFT_d,2) < p.R_d
    R_tmp = p.R_d - size(B_DFT_d,2);
    B_DFT_d = [B_DFT_d, B_DFT_d(:,1:R_tmp)];
    B_Mel_d = [B_Mel_d, B_Mel_d(:,1:R_tmp)];
end


if strcmp(p.B_sep_mode, 'Mel')
    B1_x = B_Mel_x; B1_d = B_Mel_d;
else
    B1_x = B_DFT_x; B1_d = B_DFT_d;
end
B2_x = B_DFT_x; B2_d = B_DFT_d;


ch = p.ch;

%multichannel file I/O initialization
for j = 1:ch
    fin(j) = fopen([path_in(j,:)],'rb');
end

% for j = 1:ch
%     for i = 1:p.EVENT_NUM
%         [name, ext] =strtok(path_event(j,:), '.');
%         fevent(i,j) = fopen([name,'_',num2str(i),ext],'wb');
%     end
%     for i =1:p.NOISE_NUM
%         [name, ext] =strtok(path_noise(j,:), '.');
%         fnoise(i,j) = fopen([name,'_',num2str(i),ext],'wb');
%     end
% end

for j = 1:ch
    fdenoise(j) = fopen([path_denoise(j,:)],'wb');
end

frame_len = p.framelength;
frame_shift = p.frameshift;

%% Buffer initialization
if strcmp(p.NMF_algorithm,'NTF') || strcmp(p.NMF_algorithm,'PMWF')  
    g=init_buff_NTF(B1_x, B1_d, B2_x, B2_d, p);
elseif strcmp(p.NMF_algorithm,'SNMF')
    g=init_buff(B1_x, B1_d, B2_x, B2_d, p);
end

y = zeros(ch, frame_len, 1);
d_hat = zeros(p.NOISE_NUM, ch, frame_len, 1);
x_hat = zeros(p.EVENT_NUM, ch, frame_len, 1);
x_tilde = zeros(ch, frame_len,1);

%% Wav processing
if strcmp(ftype,'.wav')
    header = zeros(ch, 22);
    for j = 1:ch
        header(j, :)=fread(fin(j), 22, 'int16'); %Skip wav header parts
    end
end

l = 1;
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
        [~,~, d_frame, g] = bnmf_sep_event_RT_IS16(y, l, g, p);
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
%             
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

fclose('all');
% %%Convert PCM to wav
% for j = 1:ch
%     for i = 1:p.EVENT_NUM
%         [name, ext] =strtok(path_event(j,:), '.');
%         pcm2wav([name,'_',num2str(i),ext],p);
%     end
%     for i =1:p.NOISE_NUM
%         [name, ext] =strtok(path_noise(j,:), '.');
%         pcm2wav([name,'_',num2str(i),ext],p);
%     end
% end

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