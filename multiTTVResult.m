function V = multiTTVResult(T,U)
% Inputs: a partialMTTKRP tesnor T, a shortened amount of facotr matrices
% U

% Outputs: V, the MTTKRP result

N = size(U,2);
    
% determines size of output MTTKRP product


R = size(U{1},2);


V = zeros(size(T,1),R);

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
