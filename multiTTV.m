function V = multiTTV(LT,RT,X,S,n,N)
% The multiTTV algortihm should efficiently calculates the matrix product
% of the n-mode matricization of X with the Khatri-Rao product of all
% entries in U, a cell array of matrices, except the nth. The mutliTTV
% takes advtange of partially computed MTTKRP's in the form of reusing them
% in the computation of a specific mttkrp of mode

% LT and RT are temp tensor
% X is cell array of factor matrices
% S is split node
% n is the mode for the mttkrp to compute
% N is the total number of nodes

% determines size of output MTTKRP product
if (n == 1)
    R = size(X{2},2);
else
    R = size(X{1},2);
end


% 1 to S-1 modes
% we use LT as the left partial MTTKRP
if n < S
    % set up output MTTKRP array
    V = zeros(size(LT,n),R);
    for r = 1:R
        % Set up cell array with appropriate vectors for ttv multiplication
        Z = cell(S-2,1);
        j = 1;
        for i = [1:n-1,n+1:S-1]
            list(j) = i;
            Z{j} = X{i}(:,r);
            j = j + 1;
        end
        % Perform ttv multiplication
        double(ttv(LT, Z, list));
        V(:,r) = ans(:,r);
    end
    
% S to N modes
% we use RT as the right partial MTTKRP
else
     % set up output MTTKRP array
     V = zeros(size(RT,n-S+1),R);
    for r = 1:R
        % Set up cell array with appropriate vectors for ttv multiplication
        Z = cell(N-S,1);
        j = 1;
        for i = [S:n-1,n+1:N]
            list(j) = i - S + 1;
            Z{j} = X{i}(:,r);
            j = j + 1;
        end
        % Perform ttv multiplication
        double(ttv(RT, Z, list));
        V(:,r) = ans(:,r);
    end
end
end





%%% old working version
%    currKRP = pKRP;
%    for i = [S-1:-1:n+1,n-1:-1:1]
%          currKRP = khatrirao(currKRP,A{i});  
%     end
%       
% % right case
% % N to S
% else 
%     % conditional statements to 
%     if(n ~= N && n ~= N-1)
%         currKRP = khatrirao(A{N}, A{N-1});
%         for i = [N-2:-1:n+1, n-1:-1:S]
%             currKRP = khatrirao(currKRP, A{i});
%         end
%         
%     elseif (n == N)
%         currKRP = khatrirao(A{N-1}, A{N-2});
%         for i = N-3:-1:S
%             currKRP = khatrirao(currKRP, A{i});
%         end
%         
%     else
%         currKRP = khatrirao(A{N}, A{N-2});
%         for i = N-3:-1:S
%             currKRP = khatrirao(currKRP, A{i});
%         end
%     end
%     
%     currKRP = khatrirao(currKRP, pKRP);
%     %V = tenmat(Z,n) * currKRP;
% end
% V = tenmat(Z,n) * currKRP;
% V = V.data;
% 
% end
