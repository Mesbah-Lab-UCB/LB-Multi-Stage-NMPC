function Kmat = kernelFn(mat1, mat2, theta)

% Norm of each row of each of the two matrices
mat1SQ = mat1*mat1';
trnorms1 = diag(mat1SQ);
mat2SQ = mat2*mat2';
trnorms2 = diag(mat2SQ);
% Multiply by a vector of ones of appropriate dimension to create two [mat1.shape[0] x mat2.shape[0]] matrices
k1 = trnorms1*ones(size(mat2,1), 1)';
k2 = ones(size(mat1,1),1)*trnorms2';
%Now, having x1^2 and x2^2, we can calculate the kernel function
k = k1 + k2 -2*mat1*mat2';       % k <-- (x1-x2)^2 = x1^2 -2x1x2 +x2^2
%k = k*(-1/(2*l^2));             % k <-- (1/(2l))(x1-x2)^2
%Kmat = (sigma^2)*exp(k);        % k <-- s^2 exp((1/(2l))*(x1-x2)^2)


% Gaussian kernel
h = (1/theta)*sqrt(k);

S = zeros(size(h));
M = size(theta,1);
for j=1:M
     S = S+(0.5*h.^2);
end
R = exp(-S);
% R = (1+sqrt(3)*abs(h)/theta)*exp(-sqrt(3)*abs(h)/theta);
Kmat = R;

end

