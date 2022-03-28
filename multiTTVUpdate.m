function X = multiTTVUpdate(T,U)
%multiTTVUpdate updates a tensor with a singular factor matrix
%
%   X = multiTTVUpdate(T,U) computes an update to tensor T by multiplying
%   the tensor by a singular factor matrix U. This computation should
%   reduce the deminsion of T by exactly one. This reduced form of T will
%   be the return tensor X
%

%% multiTTV update computation

% Tensor Dimensions and Rank Calculation
dims = size(T);
N = length(dims);
R = dims(N);

% Matricization of T and X
Tmat = reshape(double(T), [prod(dims(1:end-1)), R]);
Xmat = zeros(prod(dims(2:N-1)), R);

% iterate over the rank of T
for r = 1:R
    T_block = reshape(Tmat(:,r),[dims(1), prod(dims(2:(N-1)))]);
    X_block = U(:,r)' * T_block;
    Xmat(:,r) = X_block(:);
end

% form matricized Xmat into tensor X
X = tensor(reshape(Xmat, dims(2:N)));

end
