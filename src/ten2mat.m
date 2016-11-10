function Matrix=ten2mat(Tensor)

[n,r,h] = size(Tensor);

Matrix = zeros(n,r*h);
for i = 1:h
    Matrix(:, (i-1)*r+1:i*r) = Tensor(:,:,i);
end