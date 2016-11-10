function [Feat_splice, n]=frame_splice(Feat, p)


[K,T] =size(Feat);
K_s = (2*p.Splice+1)*K;
Feat_splice = zeros(K_s, T);
s_tmp = zeros(K_s,1);
for t = 1:T
    
    for s = 0:p.Splice
        if (t-s) < 1
            s_tmp(K*(p.Splice+s)+1:K*((p.Splice+s)+1),1) = Feat(:,t+s);
            s_tmp(K*(p.Splice-s)+1:K*((p.Splice-s)+1),1) = zeros(K,1);
        elseif (t+s) > T
            s_tmp(K*(p.Splice+s)+1:K*((p.Splice+s)+1),1) = zeros(K,1);
            s_tmp(K*(p.Splice-s)+1:K*((p.Splice-s)+1),1) = Feat(:,t-s);
        else
            s_tmp(K*(p.Splice+s)+1:K*((p.Splice+s)+1),1) = Feat(:,t+s);
            s_tmp(K*(p.Splice-s)+1:K*((p.Splice-s)+1),1) = Feat(:,t-s);
        end
    end
    Feat_splice(:, t) = s_tmp;
end
n = size(Feat_splice,1);
end