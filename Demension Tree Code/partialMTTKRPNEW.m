
function Z = partialMTTKRPNEW(T,U,S,side)
%partialMTTKRPNEW computes a partial tensor of a full tensor 
%
%   Z = partialMTTKRPNEW(T,U,S,side) calculates a partial tensor of a full
%   tennsor by computing all khatriaoproducts on one side of the split node
%   and mixing them into the tensor modes on the opposite side of the split
%   node. T is a full tensor, U is a cell array, S is the split node of 
%   the tree, side is variable that describes the side of the root node 
%   (1 is left, 2 is right)
%
%   Output: partialMTTKRP Tensor Z. Z could be the left or right tensor
%   depending on the passed input of side
%

%% Partial MTTKRP update computation

% Tensor Dimensions and Rank Calculation
dims = size(T);
N = length(dims); 
R = size(U{1}, 2);


% LEFT PARTIAL TENSOR COMPUTATION
% Calculates the khatriao product from N to S and computes with mode
% deminsons from 1 to S-1
if side == 1
     Z = tensor(reshape(double(T),[prod(dims(1:S-1)), prod(dims(S:N))]) * khatrirao(U{N:-1:S}));
     Z = tensor(reshape(Z, [dims(1:S-1), R]));
end

% RIGHT PARTIAL TENSOR COMPUTATION
% Calculates the khatriao product from S-1 to 1 and computes with mode
% deminsons from N to S
if side == 2
    Z = tensor(reshape(double(T),[prod(dims(1:S-1)), prod(dims(S:N))])' * khatrirao(U{S-1:-1:1}));
    Z = tensor(reshape(Z, [dims(S:N), R]));
    
end

end
