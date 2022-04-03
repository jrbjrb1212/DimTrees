% WILL NOT WORK FOR N < 4 %
% Update either LT or RT for output
function [LT, RT] = partialMTTKRP(X,U,S)
% X is a tense, U is a cell array, N is number of nodes
%MTTKRP Matricized tensor times Khatri-Rao product for sparse tensor.
%
%   V = partialMTTKRP(X,U,S,LR) efficiently calculates the matrix product of the
%   n-mode matricization of Z with the Khatri-Rao product of all
%   entries in U, a cell array of matrices, except the nth. The
%   partial MTTKRP is a temporary quantity which can be computed and
%   reused across two distinct MTTKRP computations.

% partialMTTKRP algorithm:
% Z is an N-way tensor and U cell array with N
% factor matrices. S is the split node

% Determines split node as middle node
% S = uint8(1 + (N-1)/2);

N = size(U,2);

% left partial KRP from N to S
% Throw in if statement later to handle out of bounds later
LKRP = khatrirao(U{N}, U{N-1});
for r = N-2:-1:S
    LKRP = khatrirao(LKRP,U{r});
end

% Perform ttv multiplication
LT = tensor(tenmat(X,1:S-1) * LKRP);

% right partial KRP from S-1 to 1
% Throw in if statement later to handle out of bounds later
RKRP = khatrirao(U{S-1}, U{S-2});
for r = S-3:-1:1
   RKRP = khatrirao(RKRP, U{r});
end

% Perform ttv multiplication
RT = tensor(tenmat(X,S:N) * RKRP);



% left partial KRP from N to S
% Throw in if statement later to handle out of bounds
% LP_KRP = khatrirao(U{N}, U{N-1});
% for r = N-2:-1:S
%     LP_KRP = khatrirao(LP_KRP,U{r});
% end
% 
% % right partial KRP from S-1 to 1
% % throw in if statement later to handle out of bounds
% RP_KRP = khatrirao(U{S-1}, U{S-2});
% for r = S-3:-1:1
%    RP_KRP = khatrirao(RP_KRP, U{r});
% end



% left case S to N
% if (S > LR)
%     
%     growSum = khatrirao(U{S}, U{S+1});
%     for r = S+2:N
%        growSum = khatrirao(growSum, U{r});
%     end
%     
% % left case 1 to S-1    
% else
%     growSum = khatrirao(U{1},U{2});
%     for r = 3:S-1
%         growSum = khatrirao(growSum, U{r});
%     end
%  end
% 
% V = growSum;

end

