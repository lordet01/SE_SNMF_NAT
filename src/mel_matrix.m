function [M,fCen,Mfull]=mel_matrix(fs,NbCh,Nfft,warp,fhigh)
    % [M,fCen,Mfull]=mel_matrix(fs,NbCh,Nfft,warp,fhigh)
    % returns Nfft/2+1 - by - NbCh matrix of triangular weights
    % NbCh filterbanks linearly spaced on a MEL scale till fhigh
    % warp: default 1.0
    % fhigh: default=fs/2
    % Mfull: matrix without truncation to Nfft/2+1 points
    
    if nargin<4,
        warp=1;
    end
    if nargin<5,
        fhigh=fs/2;
    end
    
    LowMel=2595 * log10(1+64/700);
    NyqMel=2595 * log10(1+fhigh/700);
    
    StartMel=LowMel + (0:NbCh-1)/(NbCh+1)*(NyqMel-LowMel);
    fCen=warp*700*(10.^(StartMel/2595)-1);
    StartBin=round(Nfft/fs*fCen)+1;
    
    EndMel=LowMel + (2:NbCh+1)/(NbCh+1)*(NyqMel-LowMel);
    EndBin=round(warp*Nfft/fs*700*(10.^(EndMel/2595)-1))+1;
    
    TotLen=EndBin-StartBin+1;
    
    LowLen=[StartBin(2:NbCh) EndBin(NbCh-1)]-StartBin+1;
    HiLen=TotLen-LowLen+1;
    
    M=sparse([],[],[],ceil(warp*Nfft/2+1),NbCh,sum(TotLen));
    for k=1:NbCh,
        M(StartBin(k):StartBin(k)+LowLen(k)-1,k)=(1:LowLen(k))'/LowLen(k);
        M(EndBin(k)-HiLen(k)+1:EndBin(k),k)=(HiLen(k):-1:1)'/HiLen(k);
    end
    %TotWeight=sum(M,1);
    Mfull=M;
    M=M(1:Nfft/2+1,:);
    %WeightSum=full([sum(M,1);TotWeight]);
    
end