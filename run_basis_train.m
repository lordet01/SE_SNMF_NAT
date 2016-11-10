function [B_DFT, B_Mel, A_DFT, A_Mel] = run_basis_train(DB_path, dir_Basis_full, DC_freq_set, VAD_set, event_num, B_name, R, p)

addpath('src');

B_DFT = zeros(p.F_DFT_order * (2*p.Splice+1), R * event_num);
B_Mel = zeros(p.F_order * (2*p.Splice+1), R * event_num);
A_DFT = 0;
A_Mel = 0;
DC_bin_set = floor(DC_freq_set ./ (p.fs / p.fftlength) + 0.5);
for l = 1:event_num
    fexist = fopen([dir_Basis_full{l},'/','R_',num2str(R),'.mat']);
    if fexist == -1 || p.ForceRetrain
        disp(['------',num2str(l),'-th event train------']);
        
        %% Make a full training sequence
        EventFileList = dir(DB_path{l});
        ix = randperm(length(EventFileList)-2);
        EventFileList = EventFileList(ix+2);
        s_full = zeros(p.train_seq_len_max + p.train_file_len_max,1);
        s_len_cnt = 0;
        for i=1:p.clip_subsample:length(EventFileList)
            
            [filename, ext] = strtok(EventFileList(i).name, '.');
            disp(['Training:', filename]);
            event = [DB_path{l},'/',filename,ext];
            s = wavread(event);
            s = s .* 32767;
            s_len = length(s);
            
            if VAD_set(l)
                %Perform simple VAD to exclude silence part
                p.bg_len = 0.05 * p.fs; %50ms
                p.min_voiced_len = 0.5 * p.fs; %0.3s
                p.min_unvoiced_len = 0.4 * p.fs; %0.3s
                p.thr = 0.7;
                vad = vadenergy_simple(s, p);
                s = nonzeros(s.*vad);
            elseif p.train_ANOT
                [v_start, v_end] = load_anot(filename,s_len,p);
                s = s(v_start:v_end,1);
            elseif s_len > p.train_file_len_max
                s = s(1:p.train_file_len_max);
            end
            s = s ./ sqrt(var(s));
            s = s ./ max(abs(s)) .* 30000;
            s_len = length(s);
            s_full(1 + s_len_cnt : s_len + s_len_cnt) = s;
            
            s_len_cnt = s_len_cnt + s_len;
            if s_len_cnt > p.train_seq_len_max
                s_full = s_full(1:p.train_seq_len_max);
                break;
            end
        end
        if s_len_cnt <= p.train_seq_len_max
         s_full = s_full(1:s_len_cnt);
        end
        %% Extract feature for NMF
        %Feature1: DFT Magnitude
        [TF_mag, ~] = stft_fft(s_full, p.framelength, p.frameshift, p.fftlength, DC_bin_set(l), p.win_STFT, p.preemph);
        TF_mag = TF_mag(:,any(TF_mag,1)); %Exclude all-zero column
        [TF_mag, ~] = frame_splice(TF_mag,p);
        TF_mag = TF_mag .^ p.pow + p.nonzerofloor;
        if p.domain_DD
            TF_mag_DD = TF_DD(TF_mag, p);
            TF_mag = TF_mag_DD;
        end 
        
        %Feature2: Mel Magnitude
        m = size(TF_mag,2);
        n = p.fftlength/2 + 1;
        melmat = mel_matrix(p.fs, p.F_order, p.fftlength, 1, p.fs/2)'; %Get Mel Matrix
        %     melmat_splice = repmat(melmat,2*p.Splice+1,2*p.Splice+1);
        TF_Mel = zeros(p.F_order*(p.Splice * 2 + 1),m);
        for k = 1 : p.Splice * 2 + 1
            TF_Mel(1+(k-1)*p.F_order : k*p.F_order, :) = ...
                melmat*TF_mag(1+(k-1)*n : k*n, :);
        end
        
        rng('default'); rng(1);
        sample_idx =  randsample(size(TF_mag,2),p.cluster_buff*R);
        B_DFT_init = TF_mag(:,sample_idx);
        B_Mel_init = TF_Mel(:,sample_idx);
        if p.train_Exemplar == 0
            p.w_update_ind = true(p.cluster_buff*R,1);
            p.h_update_ind = true(p.cluster_buff*R,1);
            p.init_w = B_DFT_init; %Given from Exemplar basis as initialization
            [B_DFT_init, A_DFT_init] = sparse_nmf(TF_mag, p);
            
            p.init_w = B_Mel_init; %Given from Exemplar basis as initialization
            [B_Mel_init, A_Mel_init] = sparse_nmf(TF_Mel, p);
        else
            A_DFT_init = 0;
            A_Mel_init = 0;
        end
        wavwrite(s_full./32767, 16000, ['basis/',B_name,'/',num2str(l),'.wav']);
        clear('s_full');
        
%         %Feature2: Mel Magnitude
%         m = size(B_DFT_init,2);
%         n = p.fftlength/2 + 1;
%         melmat = mel_matrix(p.fs, p.F_order, p.fftlength, 1, p.fs/2)'; %Get Mel Matrix
%         %     melmat_splice = repmat(melmat,2*p.Splice+1,2*p.Splice+1);
%         B_Mel_init = zeros(p.F_order*(p.Splice * 2 + 1),m);
%         for k = 1 : p.Splice * 2 + 1
%             B_Mel_init(1+(k-1)*p.F_order : k*p.F_order, :) = ...
%                 melmat*B_DFT_init(1+(k-1)*n : k*n, :);
%         end
%         wn = sqrt(sum(B_Mel_init.^2));
%         B_Mel_init  = bsxfun(@rdivide,B_Mel_init,wn) + 1e-9;
        
        % Basis Normalization
        wn = sqrt(sum(B_DFT_init.^2));
        B_DFT_init  = bsxfun(@rdivide,B_DFT_init,wn) + 1e-9;
        wn = sqrt(sum(B_Mel_init.^2));
        B_Mel_init  = bsxfun(@rdivide,B_Mel_init,wn) + 1e-9;
        
        if p.cluster_buff > 1
            %% Duplicated rank reduction by Clustering based sub-sampling of NMF basis
            [~, ~, ~, D]= kmeans(B_Mel_init', R, 'distance', 'cityblock', ...
                'emptyaction', 'singleto', ...
                'onlinephase', 'off', ...
                'start', 'cluster');
            [~,Dmin_idx]=min(D);
            B_Mel_sub = B_Mel_init(:,Dmin_idx);
            B_DFT_sub = B_DFT_init(:,Dmin_idx);
            A_DFT_sub = A_DFT_init(Dmin_idx,:);
            A_Mel_sub = A_Mel_init(Dmin_idx,:);
        else
            B_Mel_sub = B_Mel_init;
            B_DFT_sub = B_DFT_init;
            A_DFT_sub = A_DFT_init;
            A_Mel_sub = A_Mel_init;
        end
        
        save([dir_Basis_full{l},'/','R_',num2str(R),'.mat'],'B_DFT_sub', 'B_Mel_sub', 'A_DFT_sub', 'A_Mel_sub', '-v7.3');
    else
        load([dir_Basis_full{l},'/','R_',num2str(R),'.mat']);
    end
    
    B_Mel(:,1 + (l-1)*R : l*R) = B_Mel_sub;
    B_DFT(:,1 + (l-1)*R : l*R) = B_DFT_sub;
    if A_DFT == 0
        A_DFT = A_DFT_sub;
        A_Mel = A_Mel_sub;
    else
        A_DFT = [A_DFT, A_DFT_sub];
        A_Mel = [A_Mel, A_Mel_sub];
    end
end


fclose('all');
end