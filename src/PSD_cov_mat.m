function Scov=PSD_cov_mat(S,covlen,m,p)

%Refered from Y.G. Jin, J. W. Shin, N.S. Kim, "Spectro-Temporal FIltering for 
% Multichannel Speech Enhancement in Short-Time Fourier Transform Domain," SPL 2015.

Scov = zeros(p.ch,p.ch,covlen);
for f=1:covlen
    for t=floor(m/2)+1
        Stf = zeros(p.ch, (2*p.M_PMWF+1) * (2*p.L_PMWF+1),1);
        k = 1;
        for m_p = -p.M_PMWF : p.M_PMWF
            for l_p = -p.L_PMWF :p.L_PMWF
                if f-p.M_PMWF < 1 || f+p.M_PMWF > covlen
                    Stf(:,k) = S(:,f,t+l_p);
                else
                    Stf(:,k) = S(:,f+m_p,t+l_p);
                end
                k = k +1;
            end
        end
        Scov(:,:,f)=Scov(:,:,f)+Stf*Stf';
    end
end
end