function [x_hat_i,d_hat_i, x_tilde, g] = bnmf_sep_event_RT_IS16(y, l, g, p)


%% Global buffer initialize
% Ym = g.Ym;
% Yp = g.Yp;
Xm_hat = g.Xm_hat;
Dm_hat = g.Dm_hat;
Xm_Mel = g.Xm_Mel;
Dm_Mel = g.Dm_Mel;
% BETA_blk = g.BETA_blk;
% GAMMA_blk = g.GAMMA_blk;
x_hat = g.x_hat;
d_hat = g.d_hat;
x_tilde = g.x_tilde;
blk_cnt = g.blk_cnt;
B_Mel_x = g.B_Mel_x;
B_Mel_d = g.B_Mel_d;
B_DFT_x = g.B_DFT_x;
B_DFT_d = g.B_DFT_d;
Ad_blk = g.Ad_blk;
A_d = g.A_d;
lambda_d_blk = g.lambda_d_blk;
% W = g.W;
% THR_NPD = g.THR_NPD;
% THR_cnt = g.THR_cnt;

% eta_2term = g.eta_2term;
% S_y = g.S_y;
% S_d = g.S_d;
% St = g.St;
% Smin = g.Smin;
% Smint =g.Smint;
% SW = g.SW;
% SWt = g.SWt;
% SMact = g.SMact;
% SMactt = g.SMactt;
% lambda_d = g.lambda_d;
lambda_dav = g.lambda_dav;
% lambda_dav_long = g.lambda_dav_long;
lambda_Gy = g.lambda_Gy;
l_mod_lswitch = g.l_mod_lswitch;
Xm_tilde = g.Xm_tilde;
% lambda_y = g.lambda_y;
r_blk = g.r_blk;

%% set local parameters
[n1,~] = size(B_Mel_d);
[~,R_d] = size(B_DFT_d);
[n2,R_x] = size(B_DFT_x);
n1_unit = floor(n1 /  (2*p.Splice+1));
n2_unit = floor(n2 /  (2*p.Splice+1));
r = R_x+R_d;
m = p.blk_len_sep;
h = p.blk_hop_sep;
sz = p.framelength;
fftlen = p.fftlength;
fftlen2 = round(fftlen/2+1);
splice = p.Splice;
if blk_cnt > h
    blk_cnt = mod(blk_cnt,h);
end
% l = g.cnt; %Time frame index


%% STFT
y = filter([1 -p.preemph], 1, y);
y = p.win_STFT .* y';
y_pad = [y; zeros(fftlen-sz,1)];
Y = fft(y_pad);
Yp = angle(Y(1:floor(fftlen/2)+1));
Ym = abs(Y(1:floor(fftlen/2)+1)) .^ p.pow;

%Set zero for LPF effect
Ym(1:p.DCbin) = zeros(p.DCbin,1);

%FIt zero to floor value
Ym(:,1) = Ym(:,1) + p.nonzerofloor;



% % % % %% -------Estimate Noise magnitudes by IMCRA ------------
% % % % [N_frame_mag]=IMCRA_func(Ym, g, p);

%Block shift
if m > 1
    Ym(:, 1:m-1) = Ym(:, 2:m); %Block shift
    Yp(:, 1:m-1) = Yp(:, 2:m);
end

%Block update
if m > 1
    Ym(1:n2-n2_unit,m) = Ym(1+n2_unit:n2,m-1);
    Yp(1:n2-n2_unit,m) = Yp(1+n2_unit:n2,m-1);
else
    Ym(1:n2-n2_unit,m) = Ym(1+n2_unit:n2,m);
    Yp(1:n2-n2_unit,m) = Yp(1+n2_unit:n2,m);
end
Ym(n2-n2_unit+1:n2,m) = Ym;
Yp(n2-n2_unit+1:n2,m) = Yp;

%Splice processing range after the separation
splice_ext = (1+splice*fftlen2):(splice+1)*fftlen2;
if blk_cnt==h    
    
    %% Feature Frequency Scale Conversion (DFT to Mel)
    if strcmp(p.B_sep_mode, 'Mel')
        Ym_Mel = zeros(n1, m);
        for k = 1 : p.Splice * 2 + 1
            Ym_Mel(1+(k-1)*n1_unit : k*n1_unit, :) = ...
                g.melmat*shiftdim(Ym(1+(k-1)*fftlen2 : k*fftlen2, :));
        end
        
        %Power normaliation of mel_mag_blk matched to Ym
        vn = sqrt(sum(Ym_Mel.^2));
        tn = sqrt(sum(Ym.^2));
        Ym_Mel  = bsxfun(@rdivide,Ym_Mel,vn) + 1e-9;
        Ym_Mel = bsxfun(@times,Ym_Mel,tn);
        Y_sep = Ym_Mel;
    elseif strcmp(p.B_sep_mode, 'DFT')
        Y_sep = Ym;
    end

    %% 1) Perform SNMF separation
    if p.basis_update_N
        p.w_update_ind = [false(R_x,1); true(R_d,1)]; % Semi-supervised:
%         B_Mel = [B_Mel_x, rand(n1,R_d)];
%         B_DFT = [B_DFT_x, rand(n2,R_d)];           
    elseif p.basis_update_E
        p.w_update_ind = [true(R_x,1); false(R_d,1)]; % Semi-supervised:
%         B_Mel = [rand(n1,R_x), B_Mel_d];
%         B_DFT = [rand(n2,R_x), B_DFT_d];         
    elseif p.basis_update_N && p.basis_update_E
        p.w_update_ind = true(r,1);
%         B_Mel = rand(n1,R_x+R_d);
%         B_DFT = rand(n2,R_x+R_d);        
    else
        p.w_update_ind = false(r,1); % Supervised
    end
    B_Mel = [B_Mel_x, B_Mel_d];
    B_DFT = [B_DFT_x, B_DFT_d];
    
    if strcmp(p.B_sep_mode, 'Mel')
        p.init_w = B_Mel; %initialization
    else
        p.init_w = B_DFT; %initialization
    end    
    
    p.h_update_ind = true(r,1);
    if p.useGPU == 0
        [~, A] = sparse_nmf(Y_sep, p);
    else
        [~, A] = sparse_nmf_GPU(Y_sep, p);
    end
    

    
    %Multiclass separation (Event)
    for i = 1:p.EVENT_NUM
        if i == p.EVENT_NUM
            R_x_i = p.EVENT_RANK(i):R_x;
        else
            R_x_i = p.EVENT_RANK(i):p.EVENT_RANK(i+1)-1;
        end
        if strcmp(p.B_sep_mode, 'Mel') && p.MelConv
            tmp_Xm_hat_Mel = B_Mel(:,R_x_i)*A(R_x_i, :);
            %         %Get Last frame from supervector
            Xm_Mel(i,:,:) = tmp_Xm_hat_Mel(:,:);
            for k = 1 : p.Splice * 2 + 1
                Xm_hat(i,1+(k-1)*fftlen2 : k*fftlen2, :) = ...
                    g.melmat'*shiftdim(Xm_Mel(i, 1+(k-1)*n1_unit : k*n1_unit, :));
            end
        else
            tmp_Xm_hat = B_DFT(:,R_x_i)*A(R_x_i, :);
            Xm_hat(i,:,:) = tmp_Xm_hat(:,:);
        end
    end
    
    %Multiclass separation (Noise)
    for i = 1:p.NOISE_NUM
        if i == p.NOISE_NUM
            R_d_i = R_x+p.NOISE_RANK(i):R_x+R_d;
        else
            R_d_i = R_x+p.NOISE_RANK(i):R_x+p.NOISE_RANK(i+1)-1;
        end

        if strcmp(p.B_sep_mode, 'Mel') && p.MelConv
            tmp_Dm_hat_Mel = B_Mel(:,R_d_i)*A(R_d_i, :);
            
            %       %Get Last frame from supervector
            Dm_Mel(i,:,:) = tmp_Dm_hat_Mel;
            for k = 1 : p.Splice * 2 + 1
                Dm_hat(i,1+(k-1)*fftlen2 : k*fftlen2, :) = ...
                    g.melmat'*shiftdim(Dm_Mel(i, 1+(k-1)*n1_unit : k*n1_unit, :));
            end
        else
            tmp_Dm_hat = B_DFT(:,R_d_i)*A(R_d_i, :);
            Dm_hat(i,:,:) = tmp_Dm_hat(:,:);
        end
    end
    Xm_hat_sum = shiftdim(sum(Xm_hat,1));
    Dm_hat_sum = shiftdim(sum(Dm_hat,1));
    
    %Get Ym_Mel_DFT
    Ym_Mel_DFT = Ym;
    if strcmp(p.B_sep_mode, 'Mel') && p.MelConv
        for k = 1 : p.Splice * 2 + 1
            Ym_Mel_DFT(1+(k-1)*fftlen2 : k*fftlen2, :) = ...
                g.melmat'*shiftdim(Ym_Mel(1+(k-1)*n1_unit : k*n1_unit, :));
        end
    end
    
    %% Calculate Block Sparsity
    if p.blk_sparse
        [Q, r_blk]= blk_sparse(Xm_hat_sum, Dm_hat_sum, r_blk, l, p);
    else
        Q = ones(n2_unit, m);
    end
    
    %% 3) Enhancement Filter Construction
    if strcmp(p.ENHANCE_METHOD, 'Wiener') || strcmp(p.ENHANCE_METHOD, 'MMSE')
        %Estimate smoothed noise PSD
        if l == 1
            lambda_dav=Ym_Mel_DFT;
        end
        
        %% A ratio calculation %20160315
        A_d_mag = (sum(A(1+R_x:R_x+R_d,m)) / R_d);
        A_x_mag = (sum(A(1:R_x,m)) / R_x);
        beta = 20*log10(A_d_mag / (A_x_mag));
        beta = beta * p.beta;
        
        if beta < p.beta
            beta = p.beta;
        elseif beta >= p.beta_max
            beta = p.beta_max;
        end
%         disp( beta );
%         beta = p.beta;
      
        lambda_dav=p.alpha_d.*lambda_dav+(1-p.alpha_d).*Dm_hat_sum * beta;
        lambda_d=lambda_dav;  % new version
        
        if strcmp(p.ENHANCE_METHOD, 'Wiener')
            G = Xm_hat_sum ./ (Xm_hat_sum + Dm_hat_sum);
        elseif strcmp(p.ENHANCE_METHOD, 'MMSE')
            eta = (p.alpha_eta * Xm_tilde + ...
                (1-p.alpha_eta) * Xm_hat_sum .* Q) ./ max(lambda_d, p.nonzerofloor);

            
            eta = max(0.0031, eta);
            G =  eta ./ (eta+ones(size(eta)));
        end
        G = min(G, 1);
    end
    if l <= p.init_N_len %Init condition
        G = zeros(n2,m) + p.nonzerofloor;
        A_x_mag = p.nonzerofloor;
    end
    Xm_tilde = G.* Ym;
    
    %% 4) Adaptive basis training using enhanced spectrums
      Q_control = (1-mean(Q)) * p.Ar_up;
%       disp(Q_control);
      if p.adapt_train_N && ( Q_control*A_d_mag > A_x_mag)
          
          if l <= p.init_N_len %Init condition
              D_ref = Ym;
          else
              M_ref =  (1-G);
              M_ref(1:p.DCbin) = zeros(p.DCbin,1) + p.nonzerofloor;
              D_ref = Ym .* M_ref;
          end
          
          %Estimate smoothed noise PSD
          if l == 1
              lambda_Gy=Ym;
          end
          
          lambda_Gy = D_ref;
          lambda_d_blk = [lambda_d_blk(:,2:p.m_a) lambda_Gy];
          %         plot(lambda_d_blk);
          %         contour(lambda_d_blk);
          Ad_blk = [Ad_blk(:,2:p.m_a) A(R_x+1:R_x+p.R_a,:)];
          
          
          r_up = Q_control * mean(Ad_blk,2) > A_x_mag;
          r_up_inv = 1 - r_up;
%           disp(sum(r_up));
          
          Ad_blk_up = Ad_blk .* repmat(r_up, [1 p.m_a]);
          Ad_blk_up = Ad_blk_up(any(Ad_blk_up,2),:); %Exclude all-zero rows
          if g.update_switch == floor(p.overlap_m_a * p.m_a)
              
              if strcmp(p.B_sep_mode, 'Mel')
                  lambda_d_blk_Mel = zeros(n1, p.m_a);
                  for k = 1 : p.Splice * 2 + 1
                      lambda_d_blk_Mel(1+(k-1)*n1_unit : k*n1_unit, :) = ...
                          g.melmat*shiftdim(lambda_d_blk(1+(k-1)*fftlen2 : k*fftlen2, :));
                  end
                  B_Mel_d_up = B_Mel_d(:,1:p.R_a) .* repmat(r_up', [n1 1]);
                  B_Mel_d_up = B_Mel_d_up(:,any(B_Mel_d_up,1)); %Exclude all-zero columns
                  B_Mel_d_rem = B_Mel_d(:,1:p.R_a) .* repmat(r_up_inv', [n1 1]);
                  B_Mel_d_rem = B_Mel_d_rem(:,any(B_Mel_d_rem,1)); %Exclude all-zero columns
                  
                  [~,R_a_up] = size(B_Mel_d_up);
                  B_Mel_d_fix = B_Mel_d(:,p.R_a+1 : end);
                  if R_a_up > 0
                      p.w_update_ind = true(R_a_up,1);
                      p.h_update_ind = false(R_a_up,1);
                      p.init_w = B_Mel_d_up; %given from exemplar basis as initialization
                      p.init_h = Ad_blk_up; %given from exemplar basis as initialization
                      [B_d_tmp, ~] = sparse_nmf(lambda_d_blk_Mel, p);
                      
                      
                      B_Mel_d = [B_Mel_d_rem, B_d_tmp, B_Mel_d_fix];
                  else
                      B_Mel_d = [B_Mel_d_rem, B_Mel_d_fix];
                  end
                      
              else
                  B_DFT_d_up = B_DFT_d(:,1:p.R_a) .* repmat(r_up', [n1 1]);
                  B_DFT_d_up = B_DFT_d_up(:,any(B_DFT_d_up,1)); %Exclude all-zero columns
                  B_DFT_d_rem = B_DFT_d(:,1:p.R_a) .* repmat(r_up_inv', [n1 1]);
                  B_DFT_d_rem = B_DFT_d_rem(:,any(B_DFT_d_rem,1)); %Exclude all-zero columns
                  
                  [~,R_a_up] = size(B_DFT_d_up);
                  B_DFT_d_fix = B_Mel_d(:,p.R_a+1 : end);
%                     disp(R_a_up);
                  if R_a_up > 0
                      p.w_update_ind = true(R_a_up,1);
                      p.h_update_ind = false(R_a_up,1);
                      p.init_w = B_DFT_d_up; %given from exemplar basis as initialization
                      p.init_h = Ad_blk_up; %given from exemplar basis as initialization
                      [B_d_tmp, ~] = sparse_nmf(lambda_d_blk, p);
                      B_DFT_d = [B_DFT_d_rem, B_d_tmp, B_DFT_d_fix];
                  else
                      B_DFT_d = [B_DFT_d_rem, B_DFT_d_fix];
                  end
%                   contour(B_DFT_d);
              end
              
              g.update_switch = 1;
          else
              g.update_switch = g.update_switch + 1;
          end
      end
      
      %% --------block-wise inverse STFT-------------
    for i = 1:p.NOISE_NUM
        tmp_Dm_hat = shiftdim(Dm_hat(i,splice_ext,:));
%         tmp_Dm_hat = Dm_hat_sum;%debug
        d_hat(i,:,:)=synth_ifft_buff(tmp_Dm_hat, Yp(splice_ext,:), sz, fftlen, p.win_ISTFT, p.preemph, p.DCbin_back, p.pow);
        d_hat(i,:,:) = d_hat(i,:,:) * p.overlapscale;
    end
    for i = 1:p.EVENT_NUM
        tmp_Xm_hat = shiftdim(Xm_hat(i,splice_ext,:));
%         tmp_Xm_hat = Xm_hat_sum;%debug
        x_hat(i,:,:)=synth_ifft_buff(tmp_Xm_hat, Yp(splice_ext,:), sz, fftlen, p.win_ISTFT, p.preemph, p.DCbin_back, p.pow);
        x_hat(i,:,:) = x_hat(i,:,:) * p.overlapscale;
    end
    x_tilde = synth_ifft_buff(Xm_tilde(splice_ext,:), Yp(splice_ext,:), sz, fftlen, p.win_ISTFT, p.preemph, p.DCbin_back, p.pow);
    x_tilde = x_tilde * p.overlapscale;
    
    blk_cnt = 0;
end

%% ----------frame signal writing------------
blk_cnt = blk_cnt+1;
x_hat = x_hat(:,:,blk_cnt);
d_hat = d_hat(:,:,blk_cnt);

%Make 2D to 3D matrix for compatibility with NTF functions
tmp = ones(size(x_hat,1),1,size(x_hat,2));
tmp(:,1,:) = x_hat;
x_hat_i = tmp;

tmp_d = ones(1,size(d_hat,1),size(d_hat,2));
tmp_d(:,1,:) = d_hat;
d_hat_i = tmp_d;
x_tilde = x_tilde(:,blk_cnt);
x_tilde = x_tilde';

%% Global buffer update
g.Ym = Ym;
g.Yp = Yp;
g.x_hat = x_hat;
g.d_hat = d_hat;
g.x_tilde = x_tilde;
g.blk_cnt = blk_cnt;
g.B_Mel_x = B_Mel_x;
g.B_Mel_d = B_Mel_d;
g.B_DFT_x = B_DFT_x;
g.B_DFT_d = B_DFT_d;
g.A_d = A_d;
g.Ad_blk = Ad_blk;
g.lambda_d_blk = lambda_d_blk;
% g.W = W;
% g.THR_NPD = THR_NPD;
% g.THR_cnt = THR_cnt;

%MCRA buffers
% g.eta_2term = eta_2term;
% g.S_y = S_y;
% g.S_d = S_d;
% g.St = St;
% g.Smin = Smin;
% g.Smint = Smint;
% g.SW = SW;
% g.SWt = SWt;
% g.SMact = SMact;
% g.SMactt = SMactt;
% g.lambda_d = lambda_d;
g.lambda_dav = lambda_dav;
% g.lambda_dav_long = lambda_dav_long;
g.lambda_Gy = lambda_Gy;
g.l_mod_lswitch = l_mod_lswitch;
g.Xm_tilde = Xm_tilde;
% g.lambda_y = lambda_y;
g.r_blk = r_blk;


end