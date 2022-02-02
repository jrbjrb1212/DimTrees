function [W,Y] = MTTKRP_TEST()
% test file
Y = 0;
% test case 1
% Z = tensor(rand(2,2,2,2,2,2,2,2,2,2,2,2));
% X = {rand(2,2), rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2),rand(2,2)};
% N = 12;

% test case 2
% Z = tensor(rand(50,50,50,50,50));
% X = {rand(50,50), rand(50,50),rand(50,50), rand(50,50), rand(50,50)};
% N = 5;

% test case 3
Z = tensor(rand(30,30,30,30,30,30));
X = {rand(30,30),rand(30,30),rand(30,30),rand(30,30),rand(30,30), rand(30,30)};
N = 6;

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
    Y = cell(N,1);
    for n = 1:N
        Y{n} = mttkrp(Z,X,n);
    end
    timeOrig = toc;

    tic
    W = cell(N, 1);
    S = uint8(1 + (N-1)/2);
    const = 0;
    [LT,RT] = partialMTTKRP(Z,X,S);
    T = LT;
    for n = 1:S-1
        [W{n}, T, const] = multiTTVNEW(T,X,S,n,const);
    end

    T = RT;
    for n = S:N
        [W{n},T,const] = multiTTVNEW(T,X,S,n,const);
    end
    timeNew = toc;
    % compare ratio here
    disp("New Time: " + timeNew + " is " + (timeOrig - timeNew) + " faster than Original Time: " + timeOrig);
    
end

end
