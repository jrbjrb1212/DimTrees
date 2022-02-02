function [V,Updated_Tensor, subConst] = multiTTVNEW(T,X,S,n, subConst)
% The multiTTV algortihm should efficiently calculates the matrix product
% of the n-mode matricization of X with the Khatri-Rao product of all
% entries in U, a cell array of matrices, except the nth. The mutliTTV
% takes advtange of partially computed MTTKRP's in the form of reusing them
% in the computation of a specific mttkrp of mode

% T is are temp tensors
% X is cell array of factor matrices
% S is split node
% n is the mode for the mttkrp to compute

N = size(X,2);

% determines size of output MTTKRP product
if (n == 1)
    R = size(X{2},2);
else
    R = size(X{1},2);
end



if n < S-1
    % multittv to find the MTTKRP result for mode 1
    % then use multittv to update tensor
    
    V = zeros(size(T,n),R);
    Updated_Tensor = tenzeros(size(T,1:ndims(T)-1));
    
    for r = 1:R
        % Set up cell array with appropriate vectors for ttv multiplication
        Z = cell(S - n -1,1);
        j = 1;
        for i = n+1:S-1
            list(j) = i - subConst;
            Z{j} = X{i}(:,r);
            j = j + 1;
        end
        % Perform ttv multiplication
        double(ttv(T, Z, list));
        V(:,r) = ans(:,r); 
    end
    Updated_Tensor = UpdateTensor(T,X,n, subConst);
    subConst = subConst + 1;
    
elseif n == S-1 % no need to make another internal node
    V = T.data;
    Updated_Tensor = tenzeros(1);
    subConst = subConst + 1;
    
elseif n < N
    
    V = zeros(size(T,n-subConst),R);
    Updated_Tensor = tenzeros(size(T,1:ndims(T)-1));
    
    for r = 1:R
        % Set up cell array with appropriate vectors for ttv multiplication
        Z = cell(N - n,1);
        j = 1;
        for i = n+1:N
            list(j) = i - subConst;
            Z{j} = X{i}(:,r);
            j = j + 1;
        end
        % Perform ttv multiplication
        double(ttv(T, Z, list));
        V(:,r) = ans(:,r); 
    end
    Updated_Tensor = UpdateTensor(T,X,n,subConst);
    subConst = subConst + 1;
    
    
    
else % n = N case % no need to make another internal node
    V = T.data;
    Updated_Tensor = tenzeros(1);
end

end
