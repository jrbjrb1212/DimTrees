
function T = partialMTTKRPNEW(T,U,S,side)
% Inputs: T is a tensor, U is a cell array, S is the split node of the
% tree, side is variable that describes the side of the tree(1 is left, 2
% is right)

% Output: partialMTTKRP Tensor T

% could also use N = size(U,2);
dims = size(T);
N = length(dims); 
R = dims(N);
% dims = size(T);
% N = ndims(T);

% add tensor to the end of computation

% left partial KRP from N to S
if side == 1
     T = tensor(reshape(double(T),[prod(dims(1:S-1)), prod(dims(S:N))]) * khatrirao(U{N:-1:S}));
     T = tensor(reshape(T, [dims(1:S-1), R]));
     
%      if (mod(N,2) == 0)
%          T = reshape(T, dims(N:-1:S-1));
%      else
%          T = reshape(T, dims(N:-1:S));
%      end

% right partial KRP from S-1 to 1
elseif side == 2
    T = tensor(reshape(double(T),[prod(dims(1:S-1)), prod(dims(S:N))])' * khatrirao(U{S-1:-1:1}));    
    T = tensor(reshape(T, [dims(S:N), R]));
    
%     if (mod(N,2) == 0)
%        T = reshape(T, dims(S:-1:1));
%     else
%        T = reshape(T, dims(S+1:-1:1));
%     end
     
end
end
