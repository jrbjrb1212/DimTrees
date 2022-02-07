
function V = multiTTVNEW(T,U,S,n,T_OR_R)

% Inputs: a partialMTTKRP tesnor T, a shortened amount of facotr matrices
% U, S the split node for the tree, n the node being calulated on, and
% T_OR_R deciding on either completing a tensor update or mttkrp result

% Outputs: V could be either an updated tensor or MTTKRP result


% mttkrp result return
if T_OR_R == 1
    N = size(U,2);
    
    % determines size of output MTTKRP product
    if (n == 1)
        R = size(U{2},2);
    else
        R = size(U{1},2);
    end
    
    V = zeros(size(T,1),R);
    
    if n < S
        for r = 1:R
            % Set up cell array with appropriate vectors for ttv multiplication
            Z = cell(N-1,1);
            j = 1;
            for i = 2:N
                list(j) = i;
                Z{j} = U{i}(:,r);
                j = j + 1;
            end

            % Perform ttv multiplication
            tempStore = double(ttv(T, Z, list));
            V(:,r) = tempStore(:,r); 
        end
    else
        for r = 1:R
            % Set up cell array with appropriate vectors for ttv multiplication
            Z = cell(N-1,1);
            j = 1;
            for i = 2:N
                list(j) = i;
                Z{j} = U{i}(:,r);
                j = j + 1;
            end
            
            % Perform ttv multiplication
            tempStore = double(ttv(T, Z, list));
            V(:,r) = tempStore(:,r); 
        end
    end
    
% update the tensor 
elseif T_OR_R == 2
    if n < S
        V = UpdateTensor(T,U,n);
    else
        V = UpdateTensor(T,U,n-S+1);
    end
end



% if n < S-1
%     % multittv to find the MTTKRP result for mode 1
%     % then use multittv to update tensor
%     
%     V = zeros(size(T,n),R);
%     Updated_Tensor = tenzeros(size(T,1:ndims(T)-1));
%     
%     for r = 1:R
%         % Set up cell array with appropriate vectors for ttv multiplication
%         Z = cell(S - n -1,1);
%         j = 1;
%         for i = n+1:S-1
%             list(j) = i - subConst;
%             Z{j} = U{i}(:,r);
%             j = j + 1;
%         end
%         % Perform ttv multiplication
%         double(ttv(T, Z, list));
%         V(:,r) = ans(:,r); 
%     end
%     Updated_Tensor = UpdateTensor(T,U,n, subConst);
%     subConst = subConst + 1;
%     
% elseif n == S-1 % no need to make another internal node
%     V = T.data;
%     Updated_Tensor = tenzeros(1);
%     subConst = subConst + 1;
%     
% elseif n < N
%     
%     V = zeros(size(T,n-subConst),R);
%     Updated_Tensor = tenzeros(size(T,1:ndims(T)-1));
%     
%     for r = 1:R
%         % Set up cell array with appropriate vectors for ttv multiplication
%         Z = cell(N - n,1);
%         j = 1;
%         for i = n+1:N
%             list(j) = i - subConst;
%             Z{j} = U{i}(:,r);
%             j = j + 1;
%         end
%         % Perform ttv multiplication
%         double(ttv(T, Z, list));
%         V(:,r) = ans(:,r); 
%     end
%     Updated_Tensor = UpdateTensor(T,U,n,subConst);
%     subConst = subConst + 1;
%     
% else % n = N case % no need to make another internal node
%     V = T.data;
%     Updated_Tensor = tenzeros(1);
% end

end
