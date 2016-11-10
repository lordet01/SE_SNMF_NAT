function [w, h, objective] = sparse_nmf_GPU(v, p)

% SPARSE_NMF Sparse NMF with beta-divergence reconstruction error, 
% L1 sparsity constraint, optimization in normalized basis vector space.
%
% [w, h, objective] = sparse_nmf(v, p)
%
% Inputs:
% v:  matrix to be factorized
% p: optional parameters
%     beta:     beta-divergence parameter (default: 1, i.e., KL-divergence)
%     cf:       cost function type (default: 'kl'; overrides beta setting)
%               'is': Itakura-Saito divergence
%               'kl': Kullback-Leibler divergence
%               'kl': Euclidean distance
%     sparsity: weight for the L1 sparsity penalty (default: 0)
%     max_iter: maximum number of iterations (default: 100)
%     conv_eps: threshold for early stopping (default: 0, 
%                                             i.e., no early stopping)
%     display:  display evolution of objective function (default: 0)
%     random_seed: set the random seed to the given value 
%                   (default: 1; if equal to 0, seed is not set)
%     init_w:   initial setting for W (default: random; 
%                                      either init_w or r have to be set)
%     r:        # basis functions (default: based on init_w's size;
%                                  either init_w or r have to be set)
%     init_h:   initial setting for H (default: random)
%     w_update_ind: set of dimensions to be updated (default: all)
%     h_update_ind: set of dimensions to be updated (default: all)
%
% Outputs:
% w: matrix of basis functions
% h: matrix of activations
% objective: objective function values throughout the iterations
%
%
%
% References: 
% J. Eggert and E. Korner, "Sparse coding and NMF," 2004
% P. D. O'Grady and B. A. Pearlmutter, "Discovering Speech Phones 
%   Using Convolutive Non-negative Matrix Factorisation
%   with a Sparseness Constraint," 2008
% J. Le Roux, J. R. Hershey, F. Weninger, "Sparse NMF ? half-baked or well 
%   done?," 2015
%
% This implementation follows the derivations in:
% J. Le Roux, J. R. Hershey, F. Weninger, 
% "Sparse NMF ? half-baked or well done?," 
% MERL Technical Report, TR2015-023, March 2015
%
% If you use this code, please cite:
% J. Le Roux, J. R. Hershey, F. Weninger, 
% "Sparse NMF ? half-baked or well done?," 
% MERL Technical Report, TR2015-023, March 2015
%   @TechRep{LeRoux2015mar,
%     author = {{Le Roux}, J. and Hershey, J. R. and Weninger, F.},
%     title = {Sparse {NMF} -?half-baked or well done?},
%     institution = {Mitsubishi Electric Research Labs (MERL)},
%     number = {TR2015-023},
%     address = {Cambridge, MA, USA},
%     month = mar,
%     year = 2015
%   }
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2015 Mitsubishi Electric Research Labs (Jonathan Le Roux,
%                                         Felix Weninger, John R. Hershey)
%   Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = size(v, 1);
n = size(v, 2);

if ~exist('p', 'var')
    p = struct;
end

if ~isfield(p, 'max_iter')
    p.max_iter = 100;
end

if ~isfield(p, 'random_seed')
    p.random_seed = 1;
end

if ~isfield(p, 'sparsity')
    p.sparsity = 0;
end

if ~isfield(p, 'conv_eps')
    p.conv_eps = 0;
end

if ~isfield(p, 'cf')
    p.cf = 'kl';
end

switch p.cf
    case 'is'
        p.beta = 0;
    case 'kl'
        p.beta = 1;
    case 'ed'
        p.beta = 2;
    otherwise
        if ~isfield(p, 'beta')
            p.beta = 1;
        end
end

if p.random_seed > 0
    rand('seed', p.random_seed);
end

if ~isfield(p, 'init_w')
    if ~isfield(p, 'r')
        error('Number of components or initialization must be given')
    end
    r = p.r;
    w = rand(m, r);
else
    ri = size(p.init_w, 2);
    w(:, 1:ri) = p.init_w;
    if isfield(p, 'r') && ri < p.r
        w(:, (ri + 1) : p.r) = rand(m, p.r - ri);
        r = p.r;
    else
        r = ri;
    end
end

if ~isfield(p, 'init_h')
    h = rand(r, n);
elseif ischar(p.init_h) && strcmp(p.init_h, 'ones')
    fprintf('sup_nmf: Initalizing H with ones.\n');
    h = ones(r, n);
else
    h = p.init_h;
end

if ~isfield(p, 'w_update_ind')
    p.w_update_ind = true(r, 1);
end

if ~isfield(p, 'h_update_ind')
    p.h_update_ind = true(r, 1);
end

% sparsity per matrix entry
if length(p.sparsity) == 1
    p.sparsity = ones(r, n) * p.sparsity;
elseif size(p.sparsity, 2) == 1
    p.sparsity = repmat(p.sparsity, 1, n);
end

% Normalize the columns of W and rescale H accordingly
wn = sqrt(sum(w.^2));
w  = bsxfun(@rdivide,w,wn);
h  = bsxfun(@times,  h,wn');

% Initialization of GPU Memories
fprintf(1,'Dumping Memories to GPU\n');
w = gpuArray(w);
h = gpuArray(h);
dph = parallel.gpu.GPUArray.zeros(r, n);
dmh = parallel.gpu.GPUArray.zeros(r, n);

if ~isfield(p, 'display') 
    p.display = 0;
end

flr = 1e-9;
lambda = max(w * h, flr);
last_cost = Inf;

objective = struct;
objective.div = zeros(1,p.max_iter);
objective.cost = zeros(1,p.max_iter);

div_beta  = p.beta;
h_ind = p.h_update_ind;
w_ind = p.w_update_ind;
update_h = sum(h_ind);
update_w = sum(w_ind);

fprintf(1,'Performing sparse NMF with beta-divergence, beta=%.1f\n',div_beta);
tic
for it = 1:p.max_iter
    
    % H updates
    if update_h > 0
        switch div_beta
            case 1
                dph = bsxfun(@plus, sum(w(:, h_ind))', p.sparsity);
                dph = max(dph, flr);
                dmh = w(:, h_ind)' * (v ./ lambda);
                h(h_ind, :) = bsxfun(@rdivide, h(h_ind, :) .* dmh, dph);
            case 2
                dph = w(:, h_ind)' * lambda + p.sparsity;
                dph = max(dph, flr);
                dmh = w(:, h_ind)' * v;
                h(h_ind, :) = h(h_ind, :) .* dmh ./ dph;
            otherwise
                dph = w(:, h_ind)' * lambda.^(div_beta - 1) + p.sparsity;
                dph = max(dph, flr);
                dmh = w(:, h_ind)' * (v .* lambda.^(div_beta - 2));
                h(h_ind, :) = h(h_ind, :) .* dmh ./ dph;                
        end
        lambda = max(w * h, flr);
    end

    
    % W updates
    if update_w > 0
        switch div_beta
            case 1
                dpw = bsxfun(@plus,sum(h(w_ind, :), 2)', ...
                    bsxfun(@times, ...
                    sum((v ./ lambda) * h(w_ind, :)' .* w(:, w_ind)), w(:, w_ind)));
                dpw = max(dpw, flr);
                dmw = v ./ lambda * h(w_ind, :)' ...
                    + bsxfun(@times, ...
                    sum(bsxfun(@times, sum(h(w_ind, :),2)', w(:, w_ind))), w(:, w_ind));
                w(:, w_ind) = w(:,w_ind) .* dmw ./ dpw;
            case 2
                dpw = lambda * h(w_ind, :)' ...
                    + bsxfun(@times, sum(v * h(w_ind, :)' .* w(:, w_ind)), w(:, w_ind));
                dpw = max(dpw, flr);
                dmw = v * h(w_ind, :)' + ...
                    bsxfun(@times, sum(lambda * h(w_ind, :)' .* w(:, w_ind)), w(:, w_ind));
                w(:, w_ind) = w(:,w_ind) .* dmw ./ dpw;
            otherwise
                dpw = lambda.^(div_beta - 1) * h(w_ind, :)' ...
                    + bsxfun(@times, ...
                    sum((v .* lambda.^(div_beta - 2)) * h(w_ind, :)' .* w(:, w_ind)), ...
                    w(:, w_ind));
                dpw = max(dpw, flr);
                dmw = (v .* lambda.^(div_beta - 2)) * h(w_ind, :)' ...
                    + bsxfun(@times, ...
                    sum(lambda.^(div_beta - 1) * h(w_ind, :)' .* w(:, w_ind)), w(:, w_ind));
                w(:, w_ind) = w(:,w_ind) .* dmw ./ dpw;
        end
        % Normalize the columns of W
        w = bsxfun(@rdivide,w,sqrt(sum(w.^2)));
        lambda = max(w * h, flr);
    end
    
    
    % Compute the objective function
    switch div_beta
        case 1
            div = sum(sum(v .* log(v ./ lambda) - v + lambda));
        case 2
            div = sum(sum((v - lambda) .^ 2));
        case 0
            div = sum(sum(v ./ lambda - log ( v ./ lambda) - 1)); 
        otherwise
            div = sum(sum(v.^div_beta + (div_beta - 1)*lambda.^div_beta ...
                - div_beta * v .* lambda.^(div_beta - 1))) / (div_beta * (div_beta - 1));
    end
    cost = div + sum(sum(p.sparsity .* h));
    
%     objective.div(it)  = div;
%     objective.cost(it) = cost;
    
    if p.display ~= 0
        fprintf('iteration %d div = %.3e cost = %.3e\n', it, div, cost);
    end
    
    % Convergence check
    if it > 1 && p.conv_eps > 0
        e = abs(cost - last_cost) / last_cost;
        if (e < p.conv_eps)
            disp('Convergence reached, aborting iteration')
%             objective.div = objective.div(1:it);
%             objective.cost = objective.cost(1:it);
            break
        end
    end
    last_cost = cost;
end

fprintf(1,'Take results from GPU memory\n');
w = gather(w);
h = gather(h);
toc
end
