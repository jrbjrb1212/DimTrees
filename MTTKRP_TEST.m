function [W,Y] = MTTKRP_TEST()
% test file

% test case 1
% Z = tensor(rand(2,2,2,2,2,2,2,2,2,2,2,2));
% X = {rand(2,2), rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2)};
% N = 12;

% test case 2
% Z = tensor(rand(50,50,50,50,50));
% X = {rand(50,50), rand(50,50),rand(50,50), rand(50,50), rand(50,50)};
% N = 5;

% test case 3
%Z = tensor(rand(25,25,25,25,25,25));
%X = {rand(25,25),rand(25,25),rand(25,25),rand(25,25),rand(25,25), rand(25,25)};
%N = 6;

% % test case 4
% Z = tensor(rand(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3));
% X = {rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3)};
% N = 17;

%test case 5
% Z = tensor(rand(3,3,3,3,3));
% X = {rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3)};
% N = 5;

% % test case 6
% Z = tensor(rand(3,3,3,3,3,3));
% X = {rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3)};
% N = 6;

% test case 7
% Z = tensor(rand(150,150,150,150));
% X = {rand(150,150),rand(150,150),rand(150,150),rand(150,150)};
% N = 4;

% test case 8
% Z = tensor(rand(30,30,30,30,30));
% X = {rand(30,30), rand(30,30),rand(30,30), rand(30,30), rand(30,30)};
% N = 5;

totalRatio = 0;
numTimes = 1;

for k = 1:numTimes
    tic
    Y = cell(N,1);
    for n = 1:N
        Y{n} = mttkrp(Z,X,n);
    end
    timeOrig = toc;

    % correct strucutre
    tic
    
    % mode dimension of tensors
    dims = size(Z);
    
    % total tensor entries
    total_entries = prod(dims);
    
    % approximate root
    approx_root = sqrt(total_entries);
     
    % do later
    % cumprod from left and cumprod from the right
    % take minium val
    
    % finds split node
    list = find(cumprod(dims) <= approx_root);
    S = list(end) + 1;
    
    
    W = cell(N, 1);
    %MATRIX = cell(N, 1);
    for n = 1:N
        
        if n == 1
            % LEFT KRP TENSOR calculation
            T = partialMTTKRPNEW(Z, X, S, 1);
            % Multittv for MTTKRP result
            W{n} = multiTTVResult(T,X(1:S-1));
            
        elseif n < S-1
            % Multittv for Internal node update
            T = multiTTVUpdate(T,X{n-1});
            % Multittv for MTTKRP result
            W{n} = multiTTVResult(T,X(n:S-1));
          
        elseif n == S-1
            %  Multittv for Internal node update
            T = multiTTVUpdate(T,X{n-1});
            % Multittv for MTTKRP result
            W{n} = double(T);
            
        elseif n == S
            % RIGHT KRP TENSOR calculation
            T = partialMTTKRPNEW(Z, X, S, 2);
            % multittv for MTTKRP result
            W{n} = multiTTVResult(T,X(n:N));
            
        elseif n < N
           % Multittv for Internal node update
           T = multiTTVUpdate(T, X{n-1});
           % Multittv for MTTKRP result
           W{n} = multiTTVResult(T, X(n:N));
           
        else
           %  Multittv for Internal node update
           T = multiTTVUpdate(T,X{n-1});
           % Multittv for MTTKRP result
           W{n} = double(T);
        end
        %MATRIX{n} = T;
        
    end
    timeNew = toc;
    disp("Optimization versus Standard Ratio: " + sprintf("%0.4f",timeOrig/timeNew) + "x speed up");
    totalRatio = totalRatio + (timeOrig/timeNew);
    %disp("New Time: " + timeNew + " is " + (timeOrig - timeNew) + " faster than Original Time: " + timeOrig);
end
disp(" ");
disp("Average speed up " + sprintf("%0.4f",totalRatio/numTimes) + "x");
end
