function [y, out]= proc_BNMF_nmoh(fin,fout, fspeech, fnoise, p)

addpath('src/BNMF_nmoh');

% Modified to a function by Kwang Myung Jeon for file-wise comparison with other NR methods 
% This demo provides a framework for supervised and online NMF based noise
% reduction:
% N. Mohammadiha, P. Smaragdis and A. Leijon, “Supervised and Unsupervised Speech 
% Enhancement Using Nonnegative Matrix Factorization,?IEEE Trans. Audio, Speech, and Language Pro-
% cess., vol. 21, no. 10, pp. 2140?151, oct. 2013.
% 
% The implementation is based on a @NMF class. Type 'help NMF' to see the
% help for this class or type 'NMF' to see the properties, then click on the 'methods' in the
% last line to see a list of available methods.
% 
% Nasser Mohammadiha November 2013



method='online';%'online' or 'supervised'

alen=512;ulen=256;%analysis and update length
% alen=p.framelength; ulen=p.frameshift;%analysis and update length
fs=16000;%sampling frequency



%-------------------------read speech-------------------------------
speech=wavread(fspeech);
% speech = speech .* 32767;
speech=speech/sqrt(var(speech));
noise=wavread(fnoise);
% noise = noise .* 32767;
noise=noise/sqrt(var(noise));
mixed = wavread(fin);
mixed = mixed .* 32767;

%learn speech model, use separate training data for test and train. Here we
%use the same data for simplicity.

spec_scale=5; %scale spectrograms to reduce rounding effect. Setting of this parameter is done experimentally.
%For normalized training data, and normalized analysis window (which is
%done here) spec_scale=5 is suitable.

fexist = fopen('src/BNMF_nmoh/B_speech.mat');
% fexist = -1;
if fexist == -1
    num_speech_basis=p.R_x; %number of basis for speech nmf model
    speechTr_Spect=spec_scale*MySpectrogram(speech,alen,ulen); %magnitude spectrogram
    speech_nmf=NMF(speechTr_Spect,num_speech_basis);%construct nmf object for speech model
    setParameters(speech_nmf,'max_it',100,'update_boundFlag',1);
    speech_nmf.train;%train the speech model
    save('src/BNMF_nmoh/B_speech.mat', 'speech_nmf');
else
    load('src/BNMF_nmoh/B_speech.mat');
end
% figure,plot(speech_nmf.LogEvidence),title('Lower bound for the loglikelihood for VB-Training speech basis.')


switch method
    case 'supervised'
        
        fexist = fopen('src/BNMF_nmoh/B_noise.mat');
        fexist= -1;
        if fexist == -1        
            noise_name='my noise';%just give a name
            a_noise=10;%
            %learn separate model for each noise beforehand using some training data
            num_noise_basis=p.R_d;%number of basis for noise nmf model
            noiseTr_Spect=spec_scale*MySpectrogram(noise,alen,ulen);
            noise_nmf=NMF(noiseTr_Spect,num_noise_basis);
            setParameters(noise_nmf,'max_it',100,'update_boundFlag',1);
            noise_nmf.train;
            save('src/BNMF_nmoh/B_noise.mat', 'noise_nmf', 'a_noise', 'noise_name');
        else
           load('src/BNMF_nmoh/B_noise.mat'); 
        end
%         figure,plot(noise_nmf.LogEvidence),title('Lower bound for the loglikelihood for VB-Training speech basis.')
        noise_data=[];
    case 'online'
        %learn noise model online
        noise_name='online';
        num_noise_basis_online=15;
        a_noise=100;% shape parameter for prior of noise activations
        %learn the initial model using the begining of
        %the mixed signal which is supposed to be noise only
        input=mixed(1:15*ulen)/sqrt(var(mixed(1:15*ulen)));%make scale comparable with offline methods
        NoS=spec_scale*MySpectrogram(input, alen, ulen);
        
        noise_data=zeros(ulen+1,50); %buffer n in section III.B in the paper
        if size(NoS,2)>size(noise_data,2)
            noise_data=NoS(:,end-size(noise_data,2)+1:end);
            newNoiseInBuffer=0;
        else
            noise_data(:,end-size(NoS,2)+1:end)=NoS;
            newNoiseInBuffer=size(NoS,2)-size(noise_data,2)+10;
        end
        noise_nmf=NMF(NoS,num_noise_basis_online);
        setParameters(noise_nmf,'max_it',1000,'update_boundFlag',0);
        noise_nmf.train;
        noise_nmf.adjust_ShapeparamBasis(200);%set minimum shape param
        noise_nmf.X=[];
end



%nmf model for the mixed signal
mixed_nmf.OutputDistr(1)=NMF(speech_nmf,noise_nmf,[]);%combine two nmf models
mixed_nmf.UserData{1,1}=[0 a_noise]; %\phi^(s) and \phi^(n) in the paper in section III.C. Should be set emprically for each noise
mixed_nmf.UserData{1,2}=noise_name;

%-----------------------------Enahancement---------------------------------

estimatedSNR=0;%initial value
win=hann(alen,'periodic');
norm_coef=sqrt(sum(win.^2));
win=win/norm_coef; %normalize window
n1 = 1;n2 = alen;
enSpeech=zeros(floor(length(mixed)/ulen)*ulen,1);%memory allocation
for n=1:floor(length(mixed)/ulen)-1
    Y=fft(mixed(n1:n2).*win);Y=Y(1:alen/2+1);
    MagIn=spec_scale*abs(Y);
    
    %run the main estimation function
    if n==1
        [mixed_nmf,Est_MagOut]=BNMF_Factorization_oneFrame(mixed_nmf,MagIn,n,estimatedSNR,spec_scale,noise_data);
    else
        [mixed_nmf,Est_MagOut]=BNMF_Factorization_oneFrame(mixed_nmf,MagIn,n,estimatedSNR,spec_scale);
    end
    
    %reconstruct enhanced signal using overlap-add framework
    X1=1/spec_scale*Est_MagOut.*exp(1i*angle(Y));
    X1(1,:)=real(X1(1,:));    X1(end,:)=real(X1(end,:));%has to be real
    X=[X1; conj(X1(end-1:-1:2,:))]; %coeff 0 and alen/2 are unique and others symmetric
    enSpeech(n1:n2)=enSpeech(n1:n2)+ifft(norm_coef*X);
    
    n1 = n1 + ulen;        n2 = n2 + ulen;
    
    %-----------estimate long term SNR--------------
    % ref:   C. Kim and R. M. Stern, “Robust signal-to-noise ratio estimation based
    % on waveform amplitude distribution analysis,?in Proc. Int. Conf. Spoken
    % Language Process. (Interspeech), 2008, pp. 2598?601.
    %for the sake of demo
    if n>50 %no estimate for first 5o frames. SNR is estimated based on past 50 frames
        noisy_input=mixed((n-50)*ulen+1:n*ulen);
        G=log(mean(abs(noisy_input)))-mean(log(abs(noisy_input)+eps));
        ruts=roots([p(1) p(2) p(3)-G]);
        [vv,ii]=min(abs(ruts));
        estimatedSNR=.998*estimatedSNR+(1-.998)*ruts(ii);
    else
        G_values=[.423 .442 .642 .885];%taken from the above paper
        snrss=[-5 0 10 20];
        p=polyfit(snrss,G_values,2);
        estimatedSNR=0;
    end
    
end
% mixed=mixed(1:length(enSpeech));
% speech=speech(1:length(enSpeech));
% noise=noise(1:length(enSpeech));
% snr_out=10*log10(var(speech)/var(speech-enSpeech))%how well did it go?
% enSpeech = enSpeech / sqrt(var(enSpeech));
wavwrite(enSpeech./32767, 16000, 16, fout);
% wavwrite(mixed, 16000, 16, 'noisy.wav');