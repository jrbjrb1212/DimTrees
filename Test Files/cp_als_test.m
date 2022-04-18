% I use an impossible to reach tolerance to test how much fast the
% optimization is than the default implementation

T = tensor(@rand, [3,3,3,3,3,3]);
R = 50;
profile on
maxiters = 30;
Y = cp_als_new(T,R,'tol',1.0e-10,'maxiters',maxiters,'dt', true);
W = cp_als(T,R,'tol',1.0e-10,'maxiters',maxiters);

profile viewer
profile off
