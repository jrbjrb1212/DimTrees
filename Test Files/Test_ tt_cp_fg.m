%% Test file showing the improvmenet on tt_cp_fg.m function
%   Improvmenet of the function stems from the use of the demension tree
%   strucutre for MTTKRP bottleneck computation.
%
%   Strucuture for efficent MTTKRP computation per work of Eswar, 
%   S., Hayashi, K., Ballard, G., Kannan, R., Matheson, M. A., & Park, H. 
%   (2019). PLANC: Parallel low rank approximation with non-negativity 
%   constraints. arXiv preprint arXiv:1909.01149.


% Input: Array of demension, rank, choice of testing algorithm
% 
% Function computes and returns both DT implementiation time and default
% implementation. Do relative error of two function output
% 
% Output: Time and relative error
%   Either in a print statement or return variables
%% Test Case 1: Dense 6-way tensor, rank 50
Z = tenrand(20,20,20,20,20,20);
X = {rand(20,50),rand(20,50),rand(20,50),rand(20,50),rand(20,50),rand(20,50)};

% Demension Tree
disp("First test case");
disp("****************************************************************");
disp("Testing Demension Tree implementation of tt_cp_fg for " + size(size(Z),2) + " way tensor with rank " + size(X{1},2));
disp(" ");
tic
output = tt_cp_fg(Z,X,'dt', true);
timeNew = toc;
disp("Function output: " + output);

disp(" ");
% defualt implementation 
disp("Testing defualt implementation of tt_cp_fg for with same demensions as above");
disp(" ");
tic
output = tt_cp_fg(Z,X);
timeDefualt = toc;
disp("Function output: " + output);

disp("Improvment versus Defualt Ratio: " + sprintf("%0.4f",timeDefualt/timeNew) + "x speed up");
disp(" ");
disp("****************************************************************");
disp(" ");
disp(" ");

%% Test Case 2: Dense 4-way tensor, rank 5
Z = tenrand(125,125,125,125);
X = {rand(125,5),rand(125,5),rand(125,5),rand(125,5)};

% Demension Tree
disp("Second test case");
disp("****************************************************************");
disp("Testing Demension Tree implementation of tt_cp_fg for " + size(size(Z),2) + " way tensor with rank " + size(X{1},2));
disp(" ");
tic
output = tt_cp_fg(Z,X,'dt', true);
timeNew = toc;
disp("Function output: " + output);

disp(" ");
% defualt implementation 
disp("Testing defualt implementation of tt_cp_fg for with same demensions as above");
disp(" ");
tic
output = tt_cp_fg(Z,X);
timeDefualt = toc;
disp("Function output: " + output);

disp("Improvment versus Defualt Ratio: " + sprintf("%0.4f",timeDefualt/timeNew) + "x speed up");
disp(" ");
disp("****************************************************************");

%% Test Case 3: Dense 5-way tensor, rank 5
Z = tenrand(40,15,25,100,50);
X = {rand(40,5),rand(15,5),rand(25,5),rand(100,5), rand(50,5)};

% Demension Tree
disp("Third test case");
disp("****************************************************************");
disp("Testing Demension Tree implementation of tt_cp_fg for " + size(size(Z),2) + " way tensor with rank " + size(X{1},2));
disp(" ");
tic
output = tt_cp_fg(Z,X,'dt', true);
timeNew = toc;
disp("Function output: " + output);

disp(" ");
% defualt implementation 
disp("Testing defualt implementation of tt_cp_fg for with same demensions as above");
disp(" ");
tic
output = tt_cp_fg(Z,X);
timeDefualt = toc;
disp("Function output: " + output);

disp("Improvment versus Defualt Ratio: " + sprintf("%0.4f",timeDefualt/timeNew) + "x speed up");
disp(" ");
disp("****************************************************************");

%% Test Case 4: Dense 3-way tensor, rank 50
Z = tenrand(150,150,150, 150);
X = {rand(150,50),rand(150,50),rand(150,50), rand(150,50)};
f_1 = norm(Z)^2;

% Demension Tree
disp("Third test case");
disp("****************************************************************");
disp("Testing Demension Tree implementation of tt_cp_fg for " + size(size(Z),2) + " way tensor with rank " + size(X{1},2));
disp(" ");
tic
output = tt_cp_fg(Z,X,f_1,'dt', true);
timeNew = toc
disp("Function output: " + output);

disp(" ");
% defualt implementation 
disp("Testing defualt implementation of tt_cp_fg for with same demensions as above");
disp(" ");
tic
output = tt_cp_fg(Z,X,f_1);
timeDefualt = toc
disp("Function output: " + output);

disp("Improvment versus Defualt Ratio: " + sprintf("%0.4f",timeDefualt/timeNew) + "x speed up");
disp(" ");
disp("****************************************************************");

%% Test Case 5: Dense 10-way Tensor, rank 50
Z = tenrand(5,6,6,8,3,3,3,20,2,5);
X = {rand(5,50),rand(6,50),rand(6,50),rand(8,50), rand(3,50),rand(3,50),rand(3,50),rand(20,50),rand(2,50),rand(5,50)};

% Demension Tree
disp("Fifth test case");
disp("****************************************************************");
disp("Testing Demension Tree implementation of tt_cp_fg for " + size(size(Z),2) + " way tensor with rank " + size(X{1},2));
disp(" ");
tic
output = tt_cp_fg(Z,X,'dt', true);
timeNew = toc;
disp("Function output: " + output);

disp(" ");
% defualt implementation 
disp("Testing defualt implementation of tt_cp_fg for with same demensions as above");
disp(" ");
tic
output = tt_cp_fg(Z,X);
timeDefualt = toc;
disp("Function output: " + output);

disp("Improvment versus Defualt Ratio: " + sprintf("%0.4f",timeDefualt/timeNew) + "x speed up");
disp(" ");
disp("****************************************************************");

% %% Test Case 5: Dense 10-way Tensor, rank 50
% % The demension tree strucutre for calculating the MTTKRP struggles to
% % reach the same performance when the demensions of the tensor are small.
% % However, this may be negliable as the computation will be incredibly fast
% % whether the default or demension tree option is used.
% Z = tenrand(2,3,2,4);
% X = {rand(2,5),rand(3,5),rand(2,5),rand(4,5)};
% 
% % Demension Tree
% disp("Sixth test case");
% disp("****************************************************************");
% disp("Testing Demension Tree implementation of tt_cp_fg for " + size(size(Z),2) + " way tensor with rank " + size(X{1},2));
% disp(" ");
% tic
% output = tt_cp_fg(Z,X,'dt', true);
% timeNew = toc;
% disp("Function output: " + output);
% 
% disp(" ");
% % defualt implementation 
% disp("Testing defualt implementation of tt_cp_fg for with same demensions as above");
% disp(" ");
% tic
% output = tt_cp_fg(Z,X);
% timeDefualt = toc;
% disp("Function output: " + output);
% 
% disp("Improvment versus Defualt Ratio: " + sprintf("%0.4f",timeDefualt/timeNew) + "x speed up");
% disp(" ");
% disp("****************************************************************");
