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
% Z = tensor(rand(10,10,10,10,10,10));
% X = {rand(10,10),rand(10,10),rand(10,10),rand(10,10),rand(10,10), rand(10,10)};
% N = 6;

% test case 4
% Z = tensor(rand(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3));
% X = {rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3)};
% N = 17;
for k = 1:10
    tic
    Y = cell(N,1);
    for n = 1:N
        Y{n} = mttkrp(Z,X,n);
    end
    timeOrig = toc;

    tic
    W = cell(N, 1);
    [LT,RT,S] = partialMTTKRP(Z,X,N);
    for n = 1:S-1
        W{n} = multiTTV(LT,RT,X,S,n,N);
    end

    for n = S:N
        W{n} = multiTTV(LT,RT,X,S,n,N);
    end
    timeNew = toc;
    disp("New Time: " + timeNew + " is " + (timeOrig - timeNew) + " faster than Original Time: " + timeOrig);
    
end

end
