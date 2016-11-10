function [e_est_frame,n_est_frame, d_frame, g] = PMWF_sep_event_RT_CHiME(s_frame, g, p)


%% Global buffer initialize
TF_blk = g.TF_blk;
TF_phs_blk = g.TF_phs_blk;
TF_E_blk = g.TF_E_blk;
TF_N_blk = g.TF_N_blk;
TF_E_blk_prv = g.TF_E_blk_prv;
TF_N_blk_prv = g.TF_N_blk_prv;
TF_D_blk_prv = g.TF_D_blk_prv;
mean_N_blk_tot_prv = g.mean_N_blk_tot_prv;

Ecov = g.Ecov;
Ncov = g.Ncov;
Ycov = g.Ycov;
Ncov_prv = g.Ncov_prv;
Ycov_prv = g.Ycov_prv;

e_est_blk = g.e_est_blk;
n_est_blk = g.n_est_blk;
d_blk = g.d_blk;

blk_cnt = g.blk_cnt;
B1_e = g.B1_e;
B1_n = g.B1_n;
B2_e = g.B2_e;
B2_n = g.B2_n;
alpha = g.alpha;
W = g.W;

%% set local parameters
[n1,~] = size(B1_n);
[~,r_n] = size(B2_n);
[n2,r_e] = size(B2_e);
n_unit = p.fftlength / 2 + 1;
ch = p.ch;
r = r_e+r_n;
m = p.blk_len_sep;
h = p.blk_hop_sep;
sz = p.framelength;
fftlen = p.fftlength;
fftlen2 = round(fftlen/2+1);
splice = p.Splice;

if blk_cnt > h
    blk_cnt = mod(blk_cnt,h);
end

%% STFT
for j = 1:p.ch
    s_frame(j,:,:) = filter([1 -p.preemph], 1, shiftdim(s_frame(j,:,:)));
end
s_frame = repmat(p.win_STFT, 1,ch)' .* s_frame;
s_frame_pad = [s_frame(:,:), zeros(ch,fftlen-sz)];

S_frame_phs = zeros(ch, fftlen2);
S_frame_mag = zeros(ch, floor(fftlen2));
for j = 1:ch
    S_frame = fft(s_frame_pad(j,:));
    S_frame_phs(j,:) = angle(S_frame(1:floor(fftlen2)));
    S_frame_mag(j,:) = abs(S_frame(1:floor(fftlen2)));
end

%Set zero for LPF effect
S_frame_mag(:, 1:p.DCbin) = zeros(ch, p.DCbin);

%FIt zero to floor value
S_frame_mag(:,:) = S_frame_mag(:,:) + p.nonzerofloor;

if m > 1
    TF_blk(:,:, 1:m-1) = TF_blk(:,:, 2:m); %Block shift
    TF_phs_blk(:,:, 1:m-1) = TF_phs_blk(:,:, 2:m);
end

%Block update
if m > 1
    TF_blk(:,1:n2-n_unit,m) = TF_blk(:,1+n_unit:n2,m-1);
    TF_phs_blk(:,1:n2-n_unit,m) = TF_phs_blk(:,1+n_unit:n2,m-1);
else
    TF_blk(:,1:n2-n_unit,m) = TF_blk(:,1+n_unit:n2,m);
    TF_phs_blk(:,1:n2-n_unit,m) = TF_phs_blk(:,1+n_unit:n2,m);
end
TF_blk(:,n2-n_unit+1:n2,m) = S_frame_mag;
TF_phs_blk(:,n2-n_unit+1:n2,m) = S_frame_phs;

%Splice processing range after the separation
splice_ext = (1+(splice-p.L_PMWF)*fftlen2):((splice+p.L_PMWF)+1)*fftlen2;

if blk_cnt==h
    m_l = 2*p.L_PMWF + 1;
    TF_blk_cur = TF_blk(:,splice_ext,:);
    TF_blk_cur = reshape(TF_blk_cur, [p.ch, n_unit, 2*p.L_PMWF + 1]);
    TF_phs_blk_cur =  TF_phs_blk(:,splice_ext,:);
    TF_phs_blk_cur = reshape(TF_phs_blk_cur, [p.ch, n_unit, m_l]);
    covlen = n_unit;
    
    %% Make Cross Covariance matrix for PMWF
    if p.PMWF        
% %         if W == 1
% %             %Noise cov. matrix
% %             N = zeros(p.ch,covlen,m);
% %             for j = 1:p.ch
% %                 for t = 1:m
% %                 N_mag_sym = TF_blk_cur(j,:,t);
% %                 N_phase_sym = TF_phs_blk_cur(j,:,t);
% %                 N(j,:,t) = N_mag_sym .* exp(sqrt(-1).*N_phase_sym);
% %                 end
% %             end
% %             
% %             Ncov_tmp = PSD_cov_mat(N,covlen,m,p);
% % %             Ncov = p.ALPHA_N_PMWF*Ncov_prv + (1-p.ALPHA_N_PMWF)*Ncov_tmp; 
% %             Ncov = Ncov + Ncov_tmp;
% %             
% %             %Normalize N for a period
% %             if mod(g.cnt,p.init_N_len) == 0
% %                 Ncov = Ncov / (p.init_N_len-1);
% %             end
% %         end
        
        %Noisy cov. matrix
        Y = zeros(p.ch,covlen,m_l);
        for j = 1:p.ch
            for t = 1:m_l
                Y_mag_sym = TF_blk_cur(j,:,t);
                Y_phase_sym = TF_phs_blk_cur(j,:,t);
                Y(j,:,t) = Y_mag_sym .* exp(sqrt(-1).*Y_phase_sym);
            end
        end
        
        
        Ycov_tmp=PSD_cov_mat(Y,covlen,m_l,p);
%         Ycov = p.ALPHA_E_PMWF*Ycov_prv + (1-p.ALPHA_E_PMWF)*Ycov_tmp;
        Ycov = Ycov + Ycov_tmp;
        
        %Normalize Y for a period
        if mod(g.cnt,p.norm_period) == 0
            Ycov = Ycov / (p.norm_period-1);
        end
        
		%Noise cov. matrix
        if W == 1
            Ncov =Ycov;
        end
        
        %Event(speech) cov. matrix
        Ecov = Ycov - Ncov;
        
        %Get trace and lamda of E and N
        NEcov = zeros(p.ch,p.ch, covlen);
        lamda_PMWF = zeros(covlen,1);
        for f = 1:covlen
            NEcov(:,:,f) = (Ncov(:,:,f) + eye(p.ch)*1e-3) \ Ecov(:,:,f);
            lamda_PMWF(f) = trace(NEcov(:,:,f));
        end

        %Construct filter
        H = zeros(p.ch,p.ch, covlen);
        for j = 1:p.ch
            u = zeros(p.ch,1);
            u(j) = 1;
            for f = 1:covlen
                H(j,:,f) = (NEcov(:,:,f) / (p.BETA_PMWF + lamda_PMWF(f) + p.nonzerofloor)) * u;
            end
        end

        %PMWF filtering
        D = zeros(p.ch, covlen, m_l); 
        for j = 1:p.ch
            for f = 1:covlen
                for t = 1:m_l
                    tmp1 = conj(H(j,:,f));
                    tmp2 = Y(:,f,t);
                    D(j,f,t) = tmp1 * tmp2;
                end
            end
        end
    end    
    
    %% smooth enhanced and input magnitudes by configuring alpha
    TF_D_blk = D(:,:,p.L_PMWF+1); 
    
    %% --------block-wise inverse STFT-------------
    %Splice extraction range
%     splice_ext_out = (1+splice*fftlen2):(splice+1)*fftlen2;
    for j = 1:ch
%         for i = 1:p.NOISE_NUM
%             TF_N_blk_ch = shiftdim(TF_N_blk(i,j,:,:));
%             n_est_blk(i,j,:,:)=synth_ifft_buff(TF_N_blk_ch, TF_phs_blk_cur_ch, sz, fftlen, p.win_ISTFT, p.preemph, p.DCbin_back);
%             n_est_blk(i,j,:,:) = n_est_blk(i,j,:,:) * p.overlapscale;
%         end
        
%         for i = 1:p.EVENT_NUM
%             TF_E_blk_ch = shiftdim(TF_E_blk(i,j,:,:));
%             e_est_blk(i,j,:,:)=synth_ifft_buff(TF_E_blk_ch, TF_phs_blk_cur_ch, sz, fftlen, p.win_ISTFT, p.preemph, p.DCbin_back);
%             e_est_blk(i,j,:,:) = e_est_blk(i,j,:,:) * p.overlapscale;
%         end
        TF_phs_blk_cur_ch = shiftdim(TF_phs_blk_cur(j,:,:));
        TF_D_blk_ch = shiftdim(TF_D_blk(j,:,:));
        TF_D_blk_ch_sym = [TF_D_blk_ch; flipud(conj(TF_D_blk_ch(2:covlen-1, :)))];
        d_blk(j,:,:) = synth_ifft_buff(TF_D_blk_ch_sym, TF_phs_blk_cur_ch, sz, fftlen, p.win_ISTFT, p.preemph, p.DCbin_back);
        d_blk(j,:,:) = d_blk(j,:,:) * p.overlapscale;
    end
    
    blk_cnt = 0;
    
    TF_E_blk_prv = TF_E_blk;
    TF_N_blk_prv = TF_N_blk;
    TF_D_blk_prv = TF_D_blk;
    Ncov_prv = Ncov;
    Ycov_prv = Ycov;
end

%% ----------frame signal writing------------
blk_cnt = blk_cnt+1;
e_est_frame = e_est_blk(:,:,:,blk_cnt);
n_est_frame = n_est_blk(:,:,:,blk_cnt);
d_frame = d_blk(:,:,blk_cnt);


%% Global buffer update
g.TF_blk = TF_blk;
g.TF_phs_blk = TF_phs_blk;
g.TF_E_blk_prv = TF_E_blk_prv;
g.TF_N_blk_prv = TF_N_blk_prv;
g.TF_D_blk_prv = TF_D_blk_prv;
g.mean_N_blk_tot_prv = mean_N_blk_tot_prv;

%PMWF buffer
g.Ycov = Ycov;
g.Ncov = Ncov;
g.Ecov = Ecov;
g.Ncov_prv = Ncov_prv;
g.Ycov_prv = Ycov_prv;

g.e_est_blk = e_est_blk;
g.n_est_blk = n_est_blk;
g.d_blk = d_blk;

g.blk_cnt = blk_cnt;
g.B1_e = B1_e;
g.B1_n = B1_n;
g.B2_e = B2_e;
g.B2_n = B2_n;
g.alpha = alpha;
g.W = W;
end