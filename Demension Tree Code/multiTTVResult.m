function V = multiTTVResult(T,U)
%multiTTVResult computes an mttkrp result
%   
%   V = multiTTVResult(T,U) computes a mttkrp result from tesnor T and a
%   shortened amount of factor matrices. T will be multipied by the
%   khatriao product of the factor matrix, which will yield V a singular
%   2-way tensor. V will be the same demensions as the factor matrices 
%   in U.
%

%% MTTKRP result calculation

% Tensor Dimensions and Rank Calculation
N = size(U,2);
R = size(U{1},2)
dims = size(T);


% Matricization of T and preallocation of V
Tmat = reshape(double(T), [prod(dims(1:end-1)), R]);
V = zeros(size(T,1),R);

% iterate efficently over the rank of T
for r = 1:R
    % computes khatriao product of the necesarry column r's of the factor
    % matrices U
    Z = cellfun(@(x) x(:,r), U(1:N), 'UniformOutput', false);
    Kvector = khatrirao(Z(end:-1:1));
    
    % Takes a block of data values from the matricized tensor
    T_block = reshape(Tmat(:,r),[dims(1), prod(dims(2:N+1))]);
    
    % Mutliplication of block from matricized tensor and khatriaoed factor
    % matrices columns
    V(:,r) = T_block * Kvector;
end

end
