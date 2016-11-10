function [C,A] = GIST_NTF(p, B, S_mag)


C_UPDATE = 1;
% B_UPDATE = 0;
A_UPDATE = 0;

[Channel, ~, M] = size(S_mag);    %size of G
[N, K] = size(B);    %size of B
str =[];

% Initialize C,A
C = rand(Channel, K);
% A = rand(M, K);
A = ones(M,K);

O = ones(Channel, N, M);
Q_size = M;
flr = p.nonzerofloor;

objective = struct;
objective.div = zeros(1,p.max_iter);
objective.cost = zeros(1,p.max_iter);


% Normalize the columns of B and rescale C accordingly
Bn = sqrt(sum(B.^2));
B  = bsxfun(@rdivide,B,Bn);
C  = bsxfun(@times,  C,Bn);

% Initialize X_hat
% for k = 1 : K
%     for h = 1 : Channel
%         Syn = B(:,k)*A(:,k)';
%         Syn_2(h,:,:) = C(h,k).*Syn;
%     end
%     X_hat = X_hat + Syn_2;
% end
X_hat_temp = kr(A,B,C);
X_hat = sum(reshape(X_hat_temp,[Channel N Q_size K]),4);
X_hat = max(X_hat,flr);
P = max(S_mag./X_hat, flr);

% Update C, B, A
for i = 1 : p.max_iter
    %Update for A
    if A_UPDATE == 1
%         for k = 1 : K
%             R(:,:,k) = C(:,k) * B(:,k)';
%         end
        CB = reshape(kr(B,C), [Channel N K]);
%         for r = 1 : Channel
%             RD_temp1 = reshape(CB(r,:,:),N,K);
%             RD_temp2 = reshape(P(r,:,:),N,M);
%             RD = RD + RD_temp2' * RD_temp1;
%             RO_temp1 = reshape(CB(r,:,:),N,K);
%             RO_temp2 = reshape(O(r,:,:),N,M);
%             RO = RO + RO_temp2' * RO_temp1;
%         end
        CB_tmp = reshape(CB, [Channel*N K]);
        P_tmp = reshape(P, [Channel*N Q_size]);
        O_tmp = reshape(O, [Channel*N Q_size]);
        PCB = max(P_tmp' * CB_tmp, flr);
        OCB = max(O_tmp' * CB_tmp, flr);
        
        A = bsxfun(@rdivide, (A.*PCB),OCB + p.sparsity);
        A = max(A, flr);
%         X_hat = zeros(Channel,N,M);
%         for k = 1 : K
%             for h = 1 : Channel
%                 Syn = B(:,k)*A(:,k)';
%                 Syn_2(h,:,:) = C(h,k).*Syn;
%             end
%             X_hat = X_hat + Syn_2;
%         end

% %         %Normalize A
% %         An = sqrt(sum(A.^2));
% %         A = bsxfun(@rdivide,A,An);
% %         C = bsxfun(@times, C, An);
        
        X_hat_temp = kr(A,B,C);
        X_hat = sum(reshape(X_hat_temp,[Channel N Q_size K]),4);
        X_hat = max(X_hat,flr);
        P = max(S_mag./X_hat, flr);
    end
    
    % Update for C
    if C_UPDATE == 1
%         for k = 1 : K
%             test(:,:,k) = B(:,k) * A(:,k)';
%         end
        BA = reshape(kr(A,B),[N Q_size K]);
    
%         for n = 1 : N
%             PD_temp1 = reshape(P(:,n,:),Channel,M);
%             PD_temp2 = reshape(BA(n,:,:),M,K);
%             PD = PD + PD_temp1 * PD_temp2;
%             PO_temp1 = reshape(O(:,n,:),Channel,M);
%             PO_temp2 = reshape(BA(n,:,:),M,K);
%             PO = PO + PO_temp1 * PO_temp2;
%         end
        BA_tmp = reshape(BA, [N*Q_size K]);
        P_tmp = reshape(P, [Channel N * Q_size]);
        O_tmp = reshape(O, [Channel N * Q_size]);
        PBA = max(P_tmp * BA_tmp, flr);
        OBA = max(O_tmp * BA_tmp , flr);
        C = bsxfun(@rdivide, (C.*PBA),OBA + p.sparsity);
        C = max(C,flr);
        %         X_hat = zeros(Channel,N,M);
        %         for k = 1 : K
        %             for h = 1 : Channel
        %                 Syn = B(:,k)*A(:,k)';
        %                 Syn_2(h,:,:) = C(h,k).*Syn;
        %             end
        %             X_hat = X_hat + Syn_2;
        %         end
        
        % %         %Normalize C by channel 1
        %         Cn = sqrt(sum(C(1,:)'.^2));
% %         C = bsxfun(@rdivide,C,C(1,:));
% %         Cn = sqrt(sum(C.^2,2));
% %         C = bsxfun(@rdivide,C,Cn); 
% %         
        X_hat_temp = kr(A,B,C);
        X_hat = sum(reshape(X_hat_temp,[Channel N Q_size K]),4);
        X_hat = max(X_hat,flr);
        P = max(S_mag./X_hat, flr);
    end
    
    div = sum(sum(S_mag(:) .* log(S_mag(:)./X_hat(:)) - S_mag(:) + X_hat(:)));
    cost = div +  sum(sum(p.sparsity .* C));
    objective.div(i)  = div;
    objective.cost(i) = cost;
    
    if p.display ~= 0
        fprintf(repmat('\b',1,length(str)));
        str = sprintf('iteration %d div = %.3e cost = %.3e', i, div, cost);
        fprintf('%s', str);
    end
    
    % Convergence check
    if i > 1 && p.conv_eps > 0
        e = abs(cost - last_cost) / last_cost;
        if (e < p.conv_eps)
            if p.display ~= 0
                disp('\nConvergence reached, aborting iteration\n');
            end
            objective.div = objective.div(1:i);
            objective.cost = objective.cost(1:i);
            break;
        end
    end
    last_cost = cost;
end
if p.display ~= 0
    disp('\nMax Iteration reached, aborting iteration\n');
end


end