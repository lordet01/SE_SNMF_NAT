function [g]=init_buff_NTF(B1_e, B1_n, B2_e, B2_n, p)


%set local parameters
[n,r_s] = size(B2_e);
[n,r_n] = size(B2_n);
r = r_s + r_n;
c = p.ch;
m = p.blk_len_sep;
sz = p.framelength;
fftlen = p.fftlength;
fftlen2 = round(fftlen/2+1);
covlen = fftlen2;
EVENT_NUM = p.EVENT_NUM;
NOISE_NUM = p.NOISE_NUM;
splice_len = fftlen2 * (2*p.L_PMWF+1);

%Global buffer update
g.TF_blk = zeros(c,n,m) + p.nonzerofloor;
g.TF_phs_blk = zeros(c,n,m);
g.TF_E_blk = zeros(EVENT_NUM,c,splice_len,m);
g.TF_N_blk = zeros(NOISE_NUM,c,splice_len,m);
g.TF_E_blk_prv = zeros(EVENT_NUM,c,splice_len,m);
g.TF_N_blk_prv = zeros(NOISE_NUM,c,splice_len,m);
g.mean_N_blk_tot_prv = zeros(c,splice_len,m);
g.TF_D_blk_prv = zeros(c,splice_len,m);
g.TF_D_blk = zeros(c,splice_len,m);
g.A_blk_prv = zeros(r, m);

%PMWF buffers
g.Ycov = zeros(p.ch,p.ch,covlen);
g.Ecov = zeros(p.ch,p.ch,covlen);
g.Ncov = zeros(p.ch,p.ch,covlen);
g.Ncov_prv = zeros(p.ch,p.ch,covlen);
g.Ecov_prv = zeros(p.ch,p.ch,covlen);
g.Ycov_prv = zeros(p.ch,p.ch,covlen);

g.d_blk = zeros(c, sz, m);
g.e_est_blk = zeros(EVENT_NUM, c, sz, m);
g.n_est_blk = zeros(NOISE_NUM, c, sz, m);

g.blk_cnt = 0;
g.B1_e = B1_e;
g.B1_n = B1_n;
g.B2_e = B2_e;
g.B2_n = B2_n;
g.alpha = 1;


%PSD Smoothing matrix for block-wise processing with block-size shift
BETA_blk = zeros(2*m,2*m);
for j = 1:p.blk_len_sep * 2
    for i = 1:j-1
        BETA_blk(i,j) = p.BETA^(j-i);
        if i > 1
            BETA_blk(i,j) = p.BETA^(j-i)*(1-p.BETA);
        end
    end
    BETA_blk(j,j) = (1-p.BETA);
end
g.BETA_blk = BETA_blk;

%Beta initial block
BETA_blk_init = zeros(2*m,2*m);
for j = 1:p.blk_len_sep * 2
    for i = 1:j-1
        BETA_blk_init(i,j) = 0;
        if i > 1
            BETA_blk_init(i,j) = 0;
        end
    end
    BETA_blk_init(j,j) = 1;
end
g.BETA_blk_init = BETA_blk_init;

%PSD Smoothing matrix for block-wise processing with block-size shift
GAMMA_blk = zeros(2*m,2*m);
for j = 1:p.blk_len_sep * 2
    for i = 1:j-1
        GAMMA_blk(i,j) = p.GAMMA^(j-i);
        if i > 1
            GAMMA_blk(i,j) = p.GAMMA^(j-i)*(1-p.GAMMA);
        end
    end
    GAMMA_blk(j,j) = (1-p.GAMMA);
end
g.GAMMA_blk = GAMMA_blk;

