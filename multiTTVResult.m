function V = multiTTVResult(T,U)
% Inputs: a partialMTTKRP tesnor T, a shortened amount of facotr matrices
% U

% Outputs: V, the MTTKRP result

N = size(U,2);
    
% determines size of output MTTKRP product
R = size(U{1},2);

V = zeros(size(T,1),R);
dims = size(T);
Tmat = reshape(double(T), [prod(dims(1:end-1)), R]);


for r = 1:R
    Z = cellfun(@(x) x(:,r), U(2:N), 'UniformOutput', false);
    Kvector = khatrirao(Z(end:-1:1));
    
    T_block = reshape(Tmat(:,r),[dims(1), prod(dims(2:N))]);
    
    %tempStore = T_block * Kvector;
    V(:,r) = T_block * Kvector;
end

end

% old info
 % Perform ttv multiplication
    % try to do it without tensortoolbox
    % try khatrioaro
    % 1 R for loop
    % pull out a block from T
    
    %temp_store = reshape(T, prod(dims(1:end-2)), dims(end-1));
    %temp_store
%


% % old working version
% for r = 1:R
%     % Set up cell array with appropriate vectors for ttv multiplication
% %     Z = cell(N-1,1);
% %     j = 1;
% %     
%     % maybe replace with cellfun
%     %Z = cellfun(@(x) x(v(x),:),C)
%    
%     Z = cellfun(@(x) x(:,r), {U{2:N}}, 'UniformOutput', false);
%     
% %     for i = 2:N
% %         list(j) = i;
% %         Z{j} = U{i}(:,r);
% %         j = j + 1;
% %     end
% %     Z{1}
%     
% 
%     % Perform ttv multiplication
%     % try to do it without tensortoolbox
%     % try khatrioaro
%     % 1 R for loop
%     % pull out a block from 
%     % form khatraio of vectors
%     tempStore = double(ttv(T, Z, 2:N));
%     V(:,r) = tempStore(:,r);  
% end
%
