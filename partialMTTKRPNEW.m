
function T = partialMTTKRPNEW(Z,U,S,side)
% Inputs: Z is a tensor, U is a cell array, S is the split node of the
% tree, side is variable that describes the side of the tree(1 is left, 2
% is right)

% Output: partialMTTKRP Tensor T

% could also use N = size(U,2);
dims = size(Z);
N = ndims(Z);

% left partial KRP from N to S
if side == 1
    if (mod(N,2) == 0)
        T = tensor(double(reshape(Z,[prod(dims(1:S-1)),prod(dims(S:N))])) * khatrirao(U{N:-1:S}));
        T = reshape(T, dims(N:-1:S-1));
    else
        T = tensor(tenmat(Z,1:S-1) * khatrirao(U{N:-1:S}));
    end
   
% right partial KRP from S-1 to 1
elseif side == 2
    if (mod(N,2) == 0)
       T = tensor(double(reshape(Z,[prod(dims(1:S-1)),prod(dims(S:N))]))' * khatrirao(U{S-1:-1:1}));
       T = reshape(T, dims(S:-1:1));
    else
       T = tensor(tenmat(Z,S:N) * khatrirao(U{S-1:-1:1}));
    end
     
end
end

