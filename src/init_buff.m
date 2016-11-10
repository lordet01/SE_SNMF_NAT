function [g]=init_buff(B_Mel_x, B_Mel_d, B_DFT_x, B_DFT_d, p)


%set local parameters
[n1,~] = size(B_Mel_x);
[n2,R_x] = size(B_DFT_x);
n2_unit = n2 /(2*p.Splice+1);
[~,R_d] = size(B_DFT_d);
n1_unit = floor(n1 /  (2*p.Splice+1));
n2_unit = floor(n2 /  (2*p.Splice+1));
r = R_x + R_d;
m = p.blk_len_sep;
sz = p.framelength;
EVENT_NUM = p.EVENT_NUM;
NOISE_NUM = p.NOISE_NUM;

%Global buffer update
g.Ym = zeros(n2,m) + p.nonzerofloor;
g.Yp = zeros(n2,m);

g.Xm_Mel = zeros(EVENT_NUM,n1,m);
g.Dm_Mel = zeros(NOISE_NUM,n1,m);

g.Xm_hat = zeros(EVENT_NUM,n2,m);
g.Dm_hat = zeros(NOISE_NUM,n2,m);
g.Xm_tilde = zeros(n2,m);

g.x_tilde = zeros(sz, m);
g.x_hat = zeros(EVENT_NUM,sz, m);
g.d_hat = zeros(NOISE_NUM,sz, m);

g.blk_cnt = 1;
g.B_Mel_x = B_Mel_x;
g.B_Mel_d = B_Mel_d;
g.B_DFT_x = B_DFT_x;
g.B_DFT_d = B_DFT_d;
g.A_d = rand(R_d, m);
g.Ad_blk = rand(p.R_a, p.m_a);
g.lambda_d_blk = zeros(n2, p.m_a);
g.lambda_Gy = zeros(n2,m);
g.update_switch = 1;
g.r_blk = zeros(n2,p.P_len_l);

%Buff from MCRA technique
% g.eta_2term = ones(n2,m);
% g.S_d = zeros(n2,m);
% g.S_y = zeros(n2,m);
% g.lambda_d = zeros(n2,m);
g.lambda_dav = zeros(n2,m);
% g.lambda_dav_long = zeros(n2,m);
% g.St =zeros(n2,m);
% g.Smin = zeros(n2,m);
% g.Smint =zeros(n2,m);
% g.SW = zeros(n2,p.Nwin);
% g.SWt = zeros(n2,p.Nwin);
% g.SMact = zeros(n2,m);
% g.SMactt = zeros(n2,m);
g.l_mod_lswitch = 0;

if strcmp(p.B_sep_mode, 'Mel')
    g.melmat = mel_matrix(p.fs, p.F_order, p.fftlength, 1, p.fs/2)'; %Get Mel Matrix
end

% %PSD Smoothing matrix for block-wise processing with block-size shift
% BETA_blk = zeros(2*m,2*m);
% for j = 1:p.blk_len_sep * 2
%     for i = 1:j-1
%         BETA_blk(i,j) = p.alpha_eta^(j-i);
%         if i > 1
%             BETA_blk(i,j) = p.alpha_eta^(j-i)*(1-p.alpha_eta);
%         end
%     end
%     BETA_blk(j,j) = (1-p.alpha_eta);
% end
% g.BETA_blk = BETA_blk;
% 
% %Beta initial block
% BETA_blk_init = zeros(2*m,2*m);
% for j = 1:p.blk_len_sep * 2
%     for i = 1:j-1
%         BETA_blk_init(i,j) = 0;
%         if i > 1
%             BETA_blk_init(i,j) = 0;
%         end
%     end
%     BETA_blk_init(j,j) = 1;
% end
% g.BETA_blk_init = BETA_blk_init;
% 
% %PSD Smoothing matrix for block-wise processing with block-size shift
% GAMMA_blk = zeros(2*m,2*m);
% for j = 1:p.blk_len_sep * 2
%     for i = 1:j-1
%         GAMMA_blk(i,j) = p.GAMMA^(j-i);
%         if i > 1
%             GAMMA_blk(i,j) = p.GAMMA^(j-i)*(1-p.GAMMA);
%         end
%     end
%     GAMMA_blk(j,j) = (1-p.GAMMA);
% end
% g.GAMMA_blk = GAMMA_blk;

