function  [X_DD] = TF_DD(X, p)

[~, L] = size(X);

X_DD = X;
for l = 2:L
    X_DD(:,l) = p.alpha_eta * X_DD(:,l-1) + (1-p.alpha_eta) * X(:,l);
end

end