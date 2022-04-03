% it is working generally
% I use an impossible to reach tolerance to test how much fast the
% optimization is than the default implementation

T = tensor(@rand, [15,15,15,15,15,15]);
R = 50;
profile on
Y = cp_als(T,R,'tol',1.0e-10,'maxiters',30);
W = cp_als_new(T,R,'tol',1.0e-10,'maxiters',30);

profile viewer
profile off
