% /*********************************************
% 	GIST Source Separation + Enhancement Engine (NMF, NTF) (v3.8)
% 	Human Media Communication & Processing Laboratory (HuCom)
%   Developer: Kwang Myung Jeon
% 	
% 	Change Log.
% 				- [20150715](v1.0) Initial released
%				- [20150715](v1.0)): Support Sparse NMF and NTF algorithm
%				- [20150804](v1.1)): Added Mel-DFT coupled dicrionary
%                                    option [2015_D. Baby, Trans.SALP]
%				- [20150804](v1.1)): Added NMF with GPU option 
%				                     [To choose either SNMF or NMF(GPU)]
%				                     (Very slow due to real-time frameworks
%				                     for now!....)
%				- [20150806](v1.1)): Revised basis storing rule for better
%                                    maintanance per each basis
%				- [20150807](v1.2)): Revised basis generation rule
%				- [20150807](v1.2)): Verified optimal parameters of block size,
%                                    windowing, pre-emphasis options
%				- [20150807](v1.2)): Sync is automatically corrected for
%                                   various window, splice settings
%				- [20150807](v1.2)): Revised basis normalization option for
%                                    both Mel and DFT not to overflow 1.0
%				- [20150807](v1.2)): Revised WIner filter to use smoothed
%                                    PSD and support Block-wise operation correctly
%				- [20150807](v1.2)): Revised MMSE filter to support Block-wise operation
%				- [20150812](v1.3)): Emproved NTF compatibiliy to revised options
%				- [20150812](v1.3)): Added noise detection optino using channel activation of NTF
%				- [20150812](v1.3)): Added Basis update option governed by noise detection
%				- [20150812](v1.3)): Added smoothing between input and enhanced magnitudes by noise detection
%				- [20150814](v1.4)): Revised Mel-magnitude basis generation by applying mel-matrix on unit 
%                                    DFT magnitude instead of applying on super vector (full splice)
%				- [20150825](v1.5)): Added Noise basis update options using initial frames for ASR experiments
%				- [20150825](v1.5)): Solved initial noise reduction probelm by forcing BETA parameter to 0 for initial regions 
%                                    (same as noise update region)
%				- [20150825](v1.6)): Revised NTF to apply sparse constraint
%				- [20151005](v1.7)): Added Parameterized Multichannel non-causal Wiener Filter (PMWF) option
%				- [20151010](v1.8)): Added MDI based spectral restoration option
%				- [20151026](v2.1)): Added NPD and DM estimation module inputting estimated audio and noise spectrums from NMF
%				- [20151208](v2.2)): Revised noise basis update rule to apply it before processing the first TF block of noisy input
%				- [20160104](v2.3)): Fixed DFT, Mel feature compatibility issues
%				- [20160104](v2.3)): Changed noise boost term prior to filtering
%				- [20160209](v2.4)): Added supervised Discriminative NMF option (ref: IS14, F. Weninger)
%               - [20160219](v2.5)): Shuffle training clips for better training performance
%               - [20160226](v3.0)): Simplified Source code to reduce redundancies
%               - [20160310](v3.1)): Redesigned noise basis adapataion by updating basis with fixed activations of previous frames
%               - [20160310](v3.1)): Redesigned Wiener and MMSE filter suits to NMF separation
%               - [20160310](v3.1)): Added p.pow control options to choose mag or pow domain under every steps including training part 
%               - [20160310](v3.2)): Changed NMF operations on smoothed magnitude domain to consider temporal continuity of audio signals
%               - [20160312](v3.3)): Added Block-sparsity estimation to measure degrees of separation in each local TF bins.
%               - [20160314](v3.4)): Added basis update rule that compare power ratios between speech and noise activation befor conduct updates
%               - [20160315](v3.4a)): Prepared to enhance three different events (car alarm, scream, speech) for Techwin demo preparation-1
%               - [20160315](v3.5a)): Optimized OD-NMF options (Length 100, update 10) and added Activation based adaptive nosise floor update
%               - [20160315](v3.6)): Merged GUI, filewise, DB-batch environments for a one core moudules
%               - [20160316](v3.7)): Implemented hybrid noise basis update logic
%               - [20160316](v3.8)): Added Adaptive beta(noise floor) decision rule using the power ratios between speech and noise activation [IS16 submitted!]
% ***********************************************/


clear; clc;
%--------------------------------------------------------------------------
addpath('exp_cond');

%-----Choose one of settings-----
addpath('settings');
% initial_setting_Exemplar;
% initial_setting_SNMF;
% initial_setting_semisupervised;
% initial_setting_Proposed_IS_20160316_Obj_results;
initial_setting_Proposed_IS_20160321_repaired;

dir_DB = 'W:\Interspeech16\DB';
dir_DB_TRAIN = 'W:\CHiME3\DB';
% EVENT_NAME = {'car_alarm', 'scream', 'TIMIT_test'};
EVENT_NAME = {'TIMIT_test'};

%% Train event bases
% DB_train_event1 = [dir_DB,'/','clean_train\',EVENT_NAME{1}];
% DB_train_event2 = [dir_DB,'/','clean_train\',EVENT_NAME{2}];
DB_train_event3 = [dir_DB_TRAIN,'/','isolated\clean\tr05_org'];


%DB_path_list = {DB_train_event1, DB_train_event2, DB_train_event3};
DB_path_list = {DB_train_event3};
dir_Basis = 'basis';

%% Train multiple event bases into one
p.domain_DD = 0; %Run on decision-directed domain (1 only for noise)
if p.train_Exemplar
    B_conf = [Testname, '-Ex-preemph','_p', num2str(p.pow),'_DD',num2str(p.domain_DD)];
else
    B_conf = [Testname, '-SNMF-preemph','_p', num2str(p.pow),'_DD',num2str(p.domain_DD)];
end

dir_Basis_full = {'','',''};
for k = 1 : p.EVENT_NUM
    B_class_event = ['Clean_train_',EVENT_NAME{k}];
    
    dir_Basis_event = [dir_Basis,'/',B_class_event,'/',B_conf];
    mkdir(dir_Basis_event);
    dir_Basis_full{k} = dir_Basis_event;
end

DC_freq_set = [80 80 80];
VAD_set = [1 1 1];
[B_x_DFT, B_x_Mel, A_x_DFT, A_x_Mel] = run_basis_train(DB_path_list, dir_Basis_full, DC_freq_set, ...
    VAD_set, p.EVENT_NUM, B_class_event, p.R_x, p);
    



%% Train noise bases
% set path
DB_train_noise = [dir_DB_TRAIN,'/','backgrounds\CHiME3_background_all\CH6'];
DB_path = {DB_train_noise};
dir_Basis = 'basis';

p.domain_DD = 0; %Run on decision-directed domain (1 only for noise)
if p.train_Exemplar
    B_conf = [Testname, '-Ex-preemph','_p', num2str(p.pow),'_DD',num2str(p.domain_DD)];
else
    B_conf = [Testname, '-SNMF-preemph','_p', num2str(p.pow),'_DD',num2str(p.domain_DD)];
end

B_class_noise = 'CHiME3_bgn_ch6';

dir_Basis_noise = [dir_Basis,'/',B_class_noise,'/',B_conf];
mkdir(dir_Basis_noise);
dir_Basis_full = {dir_Basis_noise};

DC_freq_set = [80];
VAD_set = [0];

[B_d_DFT, B_d_Mel, A_d_DFT, A_d_Mel] = run_basis_train(DB_path, dir_Basis_full, DC_freq_set, ...
                                    VAD_set, p.NOISE_NUM, B_class_noise, p.R_d, p);

%% Perform Discriminative NMF training
if p.train_DNMF
   DNMF_idx = fopen([dir_Basis_event,'/','DNMF_R_',num2str(p.R_x),'.mat']);
   if DNMF_idx == -1 || p.ForceRetrain_DNMF
       %Speech, Noise Read
       x_t = wavread(['basis/',B_class_event,'/','1.wav']);
       d_t = wavread(['basis/',B_class_noise,'/','1.wav']);
       B = [B_x_DFT, B_d_DFT];
       B_hat = run_basis_DNMF(x_t, d_t, B, p);
       B_x_DFT = B_hat(:,1:p.R_x);
       B_d_DFT = B_hat(:,1+p.R_x:p.R_x+p.R_d);
       
       B_Mel = [B_x_Mel, B_d_Mel];
       B_hat = run_basis_DNMF_Mel(x_t, d_t, B_Mel, p);
       B_x_Mel = B_hat(:,1:p.R_x);
       B_d_Mel = B_hat(:,1+p.R_x:p.R_x+p.R_d);
       
       save([dir_Basis_event,'/','DNMF_R_',num2str(p.R_x),'.mat'],'B_x_DFT','B_x_Mel', '-v7.3');
       save([dir_Basis_noise,'/','DNMF_R_',num2str(p.R_x),'.mat'],'B_d_DFT','B_d_Mel', '-v7.3');
   else
       load([dir_Basis_event,'/','DNMF_R_',num2str(p.R_x),'.mat']);
       load([dir_Basis_noise,'/','DNMF_R_',num2str(p.R_x),'.mat']);
   end
end
                                
                                
%% Denoise multiple noise type and recording condition
if p.useGPU
    gpuDevice();
end


B_DFT = [B_x_DFT, B_d_DFT];
B_Mel = [B_x_Mel, B_d_Mel];

target_list = {...
%                  'DKITCHEN', 'DWASHING', 'DLIVING',...
%                  'SCAFE', 'SPSQUARE', 'STRAFFIC', ...
%                  'PSTATION', 'PCAFETER', 'PRESTO', ...
%                  'TBUS', 'TCAR', 'TMETRO', ...
%                  'NFIELD', 'NPARK', 'NRIVER', ...
%                  'OOFFICE', 'OHALLWAY', 'OMEETING', ...=
                  }; 
% target_list = { 'DLIVING', 'STRAFFIC', 'PCAFETER', 'TMETRO', 'NRIVER', 'OOFFICE' };
% SNR_list = {'0', '10','15','5'};
SNR_list = {'10',};
target_list = {'real_CHIME3'};
for k = 1: p.EVENT_NUM
    for target = target_list
        for SNR = SNR_list
            dir_Noisy = ['Noisy_IS16_samples/',EVENT_NAME{k},'/',target{1},'/',SNR{1},'dB'];
            DB_input = [dir_DB,'/',dir_Noisy];
            DB_output = [dir_DB,'/enhanced_IS16/',EVENT_NAME{k},'/',OUTname,'/',target{1},'/',SNR{1},'dB'];
            mkdir(DB_output);


            %NMF family
            if strcmp(p.NMF_algorithm, 'SNMF')
                if strcmp(p.B_sep_mode, 'Mel')
                    run_ntf_sep_RT(DB_input, DB_output, B_x_Mel, B_d_Mel, B_x_DFT, B_d_DFT, p);
                elseif strcmp(p.B_sep_mode, 'DFT')
                    run_ntf_sep_RT(DB_input, DB_output, B_x_DFT, B_d_DFT, B_x_DFT, B_d_DFT, p);
                end
            end

            if strcmp(p.NMF_algorithm, 'IMCRA')
                run_IMCRA(DB_input, DB_output, p);
            end

            if strcmp(p.NMF_algorithm, 'BNMF_nmoh')
                run_BNMF_nmoh(DB_input, DB_output, ...
                    ['basis/',B_class_event,'/','1.wav'], ... 
                    ['basis/',B_class_noise,'/','1.wav'], p);
            end
        end
    end
end


% % target_list = {'dt05_bus_real', 'dt05_bus_simu', 'tr05_bus_real', 'tr05_bus_simu', ...
% %                'dt05_caf_real', 'dt05_caf_simu', 'tr05_caf_real', 'tr05_caf_simu', ...
% %                'dt05_ped_real', 'dt05_ped_simu', 'tr05_ped_real', 'tr05_ped_simu', ...
% %                'dt05_str_real', 'dt05_str_simu', 'tr05_str_real', 'tr05_str_simu'}; 
% target_list = {'ASP_PSTATION_0dB'};
% for target = target_list
%     dir_Noisy = ['isolated/noisy_sample/',target{1}];
%     DB_input = [dir_DB,'/',dir_Noisy];
%     DB_output = [dir_DB,'/enhanced/',OUTname,'/',target{1}];
%     mkdir(DB_output);
%     if strcmp(p.B_sep_mode, 'Mel')
%         run_ntf_sep_RT(DB_input, DB_output, B_x_Mel, B_d_Mel, B_x_DFT, B_d_DFT, p);
%     elseif strcmp(p.B_sep_mode, 'DFT')
%         run_ntf_sep_RT(DB_input, DB_output,  B_x_DFT, B_d_DFT, B_x_DFT, B_d_DFT, p);
%     end 
% end

% B_d_Mel = zeros(size(B_d_Mel,1), size(B_d_Mel,2)) + p.nonzerofloor;
% B_d_DFT = zeros(size(B_d_DFT,1), size(B_d_DFT,2)) + p.nonzerofloor;

% copyfile('initial_setting.m', [DB_output,'/initial_setting-',OUTname,'.m']);
% copyfile('Do_MultiBatch.m', [DB_output,'/Do_MultiBatch-',OUTname,'.m']);
