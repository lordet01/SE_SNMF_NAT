global p;


% NMF algorithm options
% p.NMF_algorithm = 'NTF';
% p.NMF_algorithm = 'PMWF'; %PMWF with NE based on second order statistics
p.NMF_algorithm = 'SNMF';
% p.NMF_algorithm = 'IMCRA';
% p.NMF_algorithm = 'BNMF_nmoh';
p.useGPU = 0;
p.ForceRewrite = 1; %Forcely rewrite output file
p.ForceRetrain = 0; %Forcely retrain bases
p.ForceRetrain_DNMF = 0; %Forcely retrain DNMF bases

% BNMF block
p.blk_len_sep=1;
p.blk_hop_sep=p.blk_len_sep;
p.Splice = 0;

% Signal Parameters
p.fs = 16000;
p.wintime = 0.040;
p.hoptime = 0.010;
p.ch = 1;
p.framelength =round(p.wintime*p.fs);
p.frameshift  =round(p.hoptime*p.fs);
p.delay = p.Splice + p.blk_len_sep + floor(p.wintime / p.hoptime / 2 + 0.5); %Delay compensated output
% p.delay = 0; %Delayed Output
p.fftlength = 2^ceil(log2(p.framelength));     
p.F_DFT_order = 0.5*p.fftlength+1;
p.F_order = 64;
p.overlapscale = 2*p.frameshift/p.framelength;
p.win_STFT = sqrt(hann(p.framelength,'periodic'));
% p.win_STFT = ones(p.framelength, 1); %Set type of window
p.win_ISTFT = sqrt(hann(p.framelength,'periodic'));
% p.win_ISTFT = hann(p.framelength,'periodic'); %Set type of window
p.pow = 2; %Power coefficient: [1: Mag., 2: Pow]

% NMF parameters
p.EVENT_NUM = 1;
p.EVENT_RANK = [1];
% p.EVENT_RANK = [1, 21, 41];
p.NOISE_NUM = 1;
p.NOISE_RANK = [1];
p.train_Exemplar = 0; %train with SNMF if 0
p.train_DNMF = 0; %Train Basis again using discriminative NMF algorithm
p.cluster_buff = 1; %Maximum rank scale before clustering (1: Turn off clustering)
p.R_x = 100;
p.R_d = 100;
p.clip_subsample = 1;
p.train_file_len_max = p.fs * 60; % 1 min, second unit, inf: turn off this option
p.train_seq_len_max = p.fs * 720; %12 min
p.nonzerofloor = 1e-9;

%Basis update option
p.adapt_train_N = 1;
p.init_N_len = 15; %No. of initial frames used for nosie basis update
p.R_a = 50;
p.m_a = 100; %No. of stacked block for basis adaptation
p.overlap_m_a = 0.1; %No. of stacked block for basis adaptation
p.Ar_up = 2; %Define Ax and Ad ratio for noise dictionary update (Lower: Update frequently, Higher: Update rarely)
 
%Block sparsity options
p.blk_sparse = 1; %block sparsity switch
p.P_len_k = 60; % vertical (frequency bin) size of Block for local sparsity calculation
p.P_len_l = 20; % horizental (time frame index) size of Block for local sparsity calculation
% p.kappa = 1.0;
p.nu = 1.0;
p.alpha_p = 0.4; %DD smoothing factor for P
p.blk_gap = 7; %Blk_gap for complexity issue (1 for ideal), Odd only!

%MDI options
p.MDI_est = 0; %MDI estimation option
p.MDI_est_noise = 0; %MDI estimation option (noise)
p.sparsity_mdi = 5;
p.conv_eps_mdi = 1e-5; 

%PMWF options
p.PMWF = 0;
p.BETA_PMWF = 10; %0: MVDR, >0: PMWF
p.M_PMWF = 2; %spectral neighbor region
p.L_PMWF = 2; %temporal neighbor region
p.ALPHA_E_PMWF = 0.3; %Mixing ratio between current Event PSD-CM vs average one
p.norm_period = p.init_N_len; %Normalization period of PSD cov. matrix
p.Ncov_update = 1;%Update Ncov by cov. mat. from output - input

%Front-end processing
p.preemph = 0.92;
p.DCfreq = 160; %Hz
p.DCbin = floor(p.DCfreq / (p.fs / p.fftlength) + 0.5); %Forcely give 0 to 1~N bin which is not important in speech
%Back-end processing
p.DCbin_back = p.DCbin; %Forcely give 0 to 1~N bin which is not important in speech
%Multi-channel options
p.filegap = p.ch; %No. of file consisting one session

%Run option
p.separation = 1;
p.B_sep_mode = 'DFT'; %['DFT', 'Mel']
p.MelConv = 1; %Use frequency scale conversion. 00: Coupled dictionary

%Sub_options
p.train_VAD = 0;
p.train_ANOT = 0;

%SNMF parameters
p.cf = 'kl';   %  'is', 'kl', 'ed'; takes precedence over setting the beta value
p.sparsity = 5;
p.max_iter = 25; % Stopping criteria
p.conv_eps = 1e-3; 
p.display   = 0; % Display evolution of objective function
p.random_seed = 1; % Random seed: any value over than 0 sets the seed to that value
p.cost_check = 1;
p.basis_update_N = 0;
p.basis_update_E = 0;
p.est_scale = 1.0;
%Single channel enhancement options
p.ENHANCE_METHOD = 'MMSE'; %['Wiener', 'MMSE']
%2.1) Parameters of "Decision-Directed" a Priori SNR Estimate
p.alpha_eta=0.3;	% Recursive averaging parameter
p.eta_min=10^(-18/10);	% Lower limit constraint
%2.2) Smooth over frequency and time
% % % p.w = 2; % Size of frequency smoothing window function=2*w+1
% % % p.b=hanning(2*p.w+1);
% % % p.b=p.b/sum(p.b);     % normalize the window function
% % % p.alpha_s = 0.9; % Recursive averaging parameter for the smoothing operation
% % % %2.3) calculate speech presence and absence probabilites (SPP and SAP)
% % % p.delta_s=1.67;		% Local noise factor
% % % p.Bmin=1.66;
% % % p.delta_y=4.6;		% Local noise factor
% % % p.delta_yt=3;
% % % p.delta_s=1.67;		% Local minimum factor
% % % %2.4) estimate noise magnitues using SAP
p.alpha_d = 0.85; % Recursive averaging parameter for the noise
% % % p.alpha_d_long = 0.99;
% % % p.Vwin = 15;
% % % p.Nwin = 8;
% % % %2.5) Design a noise reduction filter by A-posteriori and A-priori SNRs incorporating SPP and SAP
p.beta = 1.0; %Bias compensation factor in the noise estimator [!!! Give 8 for listening Test!, 2 for objective test]
p.beta_max = 10000.0;

%VAD for speech training
p.speech_train_start = round(p.fs * 0.5);
p.speech_train_end = round(p.fs * 1.5);
p.speech_train_len = p.speech_train_end - p.speech_train_start;

Testname = ['IS16','_Splice',num2str(p.Splice),''];
%Name of the test system
OUTname = [Testname,'_',p.NMF_algorithm,'_A',num2str(p.adapt_train_N),'_M',num2str(p.MDI_est_noise),'_r',num2str(p.R_x),'_p',num2str(p.pow), ...
           '_',p.ENHANCE_METHOD,'_P',num2str(p.blk_sparse),'_Proposed_IS16_20160316'];