function X = multiTTVUpdate(T,U)
% Inputs: a partialMTTKRP tesnor T, a shortened amount of facotr matrices
% U

% Outputs: X the updated tensor

dims = size(T);
N = length(dims);
R = dims(N);
dims
prod(dims(1:N-1))

Tmat = reshape(T, [prod(dims(1:dims-1)), R]);
Xmat = zeros(prod(dims(2:N-1)), R);


for r = 1:R
    % Reshape as tensor then call ttv in all modes but the first or
    % form kronker product and then multiply on the right
    % get away from all the cases
    T_block = reshape(Tmat(:,R),[dims(1), prod(dims(2:(N-1)))]);
    X_block = U(:,r)' * T_block;
    Xmat(:,r) = X_block(:);
end

X = reshape(Xmat, dims(2:N));

% Old working way
%X = UpdateTensor(T,U,1);

end
