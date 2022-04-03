function void = tt_cp_fg_TEST(option)

if option == 1
    T = tensor(rand(1000,1000,1000));
    U = {rand(1000,100), rand(1000,100), rand(1000,100)};
elseif option == 2
    T = tensor(rand(75,75,75,75));
    U = {rand(75,20), rand(75,20), rand(75,20), rand(75,20)};
elseif option == 3
    T = tensor(rand(30,30,30,30,30,30));
    U = {rand(30,5), rand(30,5), rand(30,5), rand(30,5), rand(30,5), rand(30,5)};
elseif option == 4
    T = tensor(rand(3,3,3,3,3));
    U = {rand(3,3),rand(3,3),rand(3,3),rand(3,3),rand(3,3)};
else
    disp("Enter a valid option (1-3)");
end


%profile on
% new implementation
tic
X = tt_cp_fg_NEW(T,U);
timeNew = toc;

% old implementation
tic
Y = tt_cp_fg(T,U);
timeOrig = toc;

%profile viewer
%profile off

% X - Y for testing
disp("Optimization versus Standard Ratio: " + sprintf("%0.4f",timeOrig/timeNew) + "x speed up");


end
