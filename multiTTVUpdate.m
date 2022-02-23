function [X,T_block] = multiTTVUpdate(T,U)
% Inputs: a partialMTTKRP tensor T, a single factor matrix U

% Outputs: X the updated tensor

dims = size(T);
N = length(dims);
R = dims(N);
Tmat = reshape(double(T), [prod(dims(1:end-1)), R]);
Xmat = zeros(prod(dims(2:N-1)), R);


for r = 1:R
    % Reshape as tensor then call ttv in all modes but the first or
    % form kronker product and then multiply on the right
    % get away from all the cases
    
    T_block = reshape(Tmat(:,r),[dims(1), prod(dims(2:(N-1)))]);
    X_block = U(:,r)' * T_block;
    
    Xmat(:,r) = X_block(:);
end

X = tensor(reshape(Xmat, dims(2:N)));

end
