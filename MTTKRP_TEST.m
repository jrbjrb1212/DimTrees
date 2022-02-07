function [W,Y] = MTTKRP_TEST()
% test file

% test case 1
% Z = tensor(rand(2,2,2,2,2,2,2,2,2,2,2,2));
% X = {rand(2,2), rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2)};
% N = 12;

% test case 2
Z = tensor(rand(50,50,50,50,50));
X = {rand(50,50), rand(50,50),rand(50,50), rand(50,50), rand(50,50)};
N = 5;

% test case 3
% Z = tensor(rand(25,25,25,25,25,25));
% X = {rand(25,25),rand(25,25),rand(25,25),rand(25,25),rand(25,25), rand(25,25)};
% N = 6;

% test case 4
% Z = tensor(rand(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3));
% X = {rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3)};
% N = 17;

% test case 5
% Z = tensor(rand(3,3,3,3,3,3));
% X = {rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3), rand(3,3)};
% N = 6;

for k = 1:10
    tic
    % current implementation
    Y = cell(N,1);
    for n = 1:N
        Y{n} = mttkrp(Z,X,n);
    end
    timeOrig = toc;

    % correct strucutre for dim trees
    tic
    S = uint8(1 + (N-1)/2);
    W = cell(N, 1);
    %MATRIX = cell(N, 1);
    for n = 1:N
        if n == 1
            % LEFT KRP TENSOR calculation
            T = partialMTTKRPNEW(Z, X, S, 1);
            % Multittv for MTTKRP result
            W{n} = multiTTVNEW(T,X(1:S-1),S,n,1);
            
        elseif n < S-1
            % Multittv for Internal node update
            T = multiTTVNEW(T,X(n-1:S-1),S,1,2);
            % Multittv for MTTKRP result
            W{n} = multiTTVNEW(T,X(n:S-1),S,1,1);
          
        elseif n == S-1
            %  Multittv for Internal node update
            T = multiTTVNEW(T,X(n-1:S-1),S,1,2);
            % Multittv for MTTKRP result
            W{n} = T.data;
            
        elseif n == S
            % RIGHT KRP TENSOR calculation
            T = partialMTTKRPNEW(Z, X, S, 2);
            % multittv for MTTKRP result
            W{n} = multiTTVNEW(T,X(n:N),S,S,1);
            
        elseif n < N
           % Multittv for Internal node update
           T = multiTTVNEW(T, X(n-1:N), S, S, 2);
           % Multittv for MTTKRP result
           W{n} = multiTTVNEW(T, X(n:N), S, S, 1);
           
        else
           %  Multittv for Internal node update
           T = multiTTVNEW(T,X(n-1:N),S,S,2);
           % Multittv for MTTKRP result
           W{n} = T.data;
    
        end
        %MATRIX{n} = T;
    end
    timeNew = toc;   
    disp("Optimization versus Standard Ratio: " + sprintf("%0.4f",timeOrig/timeNew) + "x speed up");
    %disp("New Time: " + timeNew + " is " + (timeOrig - timeNew) + " faster than Original Time: " + timeOrig);
end
end
    
    
    
%     tic
%     W = cell(N, 1);
%     S = uint8(1 + (N-1)/2);
%     const = 0;
%     [LT,RT] = partialMTTKRP(Z,X,S);
%     T = LT;
%     for n = 1:S-1
%         [W{n}, T, const] = multiTTVNEW(T,X,S,n,const);
%     end
% 
%     T = RT;
%     for n = S:N
%         [W{n},T,const] = multiTTVNEW(T,X,S,n,const);
%     end
%     timeNew = toc;
    
%     % compare ratio here
%     disp("New Time: " + timeNew + " is " + (timeOrig - timeNew) + " faster than Original Time: " + timeOrig);
%     
% end
% 
% end
