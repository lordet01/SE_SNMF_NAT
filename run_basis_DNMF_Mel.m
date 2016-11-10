function [B_hat] = run_basis_DNMF_Mel(x, d, B, p)

%Ref: IS14, F. Weninger

if length(x) < length(d)
    d = d(1:length(x));
else
    x = x(1:length(d));
end
y =  x + d;

%Get Spectrograms of signals
[X, ~] = stft_fft(x, p.framelength, p.frameshift, p.fftlength, p.DCbin, p.win_STFT, p.preemph);
X = X(:,any(X,1)); %Exclude all-zero column
[X, ~] = frame_splice(X,p);
X = X .^ p.pow + p.nonzerofloor;
% if p.domain_DD
%     X_DD = TF_DD(X, p);
%     X = X_DD;
% end
%Feature2: Mel Magnitude
m = size(X,2);
n = p.fftlength/2 + 1;
melmat = mel_matrix(p.fs, p.F_order, p.fftlength, 1, p.fs/2)'; %Get Mel Matrix
%     melmat_splice = repmat(melmat,2*p.Splice+1,2*p.Splice+1);
X_Mel = zeros(p.F_order*(p.Splice * 2 + 1),m);
for k = 1 : p.Splice * 2 + 1
    X_Mel(1+(k-1)*p.F_order : k*p.F_order, :) = ...
        melmat*X(1+(k-1)*n : k*n, :);
end

[D, ~] = stft_fft(d, p.framelength, p.frameshift, p.fftlength, p.DCbin, p.win_STFT, p.preemph);
D = D(:,any(D,1)); %Exclude all-zero column
[D, ~] = frame_splice(D,p);
D = D .^ p.pow + p.nonzerofloor;
% % if p.domain_DD
%     D_DD = TF_DD(D, p);
%     D = D_DD;
% % end
%Feature2: Mel Magnitude
% m = size(D,2);
% n = p.fftlength/2 + 1;
% melmat = mel_matrix(p.fs, p.F_order, p.fftlength, 1, p.fs/2)'; %Get Mel Matrix
%     melmat_splice = repmat(melmat,2*p.Splice+1,2*p.Splice+1);
D_Mel = zeros(p.F_order*(p.Splice * 2 + 1),m);
for k = 1 : p.Splice * 2 + 1
    D_Mel(1+(k-1)*p.F_order : k*p.F_order, :) = ...
        melmat*D(1+(k-1)*n : k*n, :);
end


[Y, ~] = stft_fft(y, p.framelength, p.frameshift, p.fftlength, p.DCbin, p.win_STFT, p.preemph);
Y = Y(:,any(Y,1)); %Exclude all-zero column
[Y, ~] = frame_splice(Y,p);
Y = Y .^ p.pow + p.nonzerofloor;
% if p.domain_DD
%     Y_DD = TF_DD(Y, p);
%     Y = Y_DD;
% end
%Feature2: Mel Magnitude
% m = size(X,2);
% n = p.fftlength/2 + 1;
% melmat = mel_matrix(p.fs, p.F_order, p.fftlength, 1, p.fs/2)'; %Get Mel Matrix
%     melmat_splice = repmat(melmat,2*p.Splice+1,2*p.Splice+1);
Y_Mel = zeros(p.F_order*(p.Splice * 2 + 1),m);
for k = 1 : p.Splice * 2 + 1
    Y_Mel(1+(k-1)*p.F_order : k*p.F_order, :) = ...
        melmat*Y(1+(k-1)*n : k*n, :);
end

%Get H_hat, Eq. (6)
p.w_update_ind = false((p.R_x + p.R_d),1);
p.h_update_ind = true((p.R_x + p.R_d),1);
p.init_w = B; %given from exemplar basis as initialization
[~, A_hat] = sparse_nmf(Y_Mel, p);

% GetW_hat, Eq. (7)
p.w_update_ind = true(p.R_x,1);
p.h_update_ind = false(p.R_x,1);
p.init_w = B(:,1:p.R_x); %given from exemplar basis as initialization
p.init_h = A_hat(1:p.R_x,:); %given from exemplar basis as initialization
[B_hat_x, ~] = sparse_nmf(X_Mel, p);

p.w_update_ind = true(p.R_d,1);
p.h_update_ind = false(p.R_d,1);
p.init_w = B(:,p.R_x+1:p.R_x+p.R_d); %given from exemplar basis as initialization
p.init_h = A_hat(p.R_x+1:p.R_x+p.R_d,:); %given from exemplar basis as initialization
[B_hat_d, ~] = sparse_nmf(D_Mel, p);

B_hat = [B_hat_x, B_hat_d];





end