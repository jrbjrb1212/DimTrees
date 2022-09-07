% test hooi
clear
addpath('tensor_toolbox-master/')

rng(0)
%%
%C = Cauchy_tensor([100,100,100]);

n = 200;
d = [500,500,500];
[ii,ij,ik] = ndgrid(1:d(1),1:d(2),1:d(3));
C = 1./(ii+ij+ik);
C = tensor(C);

%C = tensor(randn(200,200,200));

% m = 100; n = 100; p = 100;
% noise = 1e-4;
% T = create_problem('Type','Tucker','Size',[m n p],'Noise',noise);
% C = T.Data;

% 
% csize = 100;
% tsize = 100;
% seed = 0;
% ratio = .3;
% T = generateRandomTensor(4,csize,tsize,seed,ratio);
% C = tensor(T);


r_init = [100,100,100];
tol = 1e-6;
maxiter = 500;
modes = [3,1,2];


%% STHOSVD
tic; T2 = hosvd(C,tol); t_st = toc;
R_st = size(T2.core);
% 
% relerr2 = norm(tensor(T2)-C)/norm(C)
% 
% 
% %% HOOI with rank from STHOSVD
% 
% T3 = tucker_als(C,R_st);
% 
%  relerr3 = norm(tensor(T3)-C)/norm(C)


%% adaptive HOOI
disp('Adaptive HOOI')
%profile on
[T,iter,times] = hooi_adapt(C,r_init,tol,maxiter,modes);
t_hooi = sum(times);

Rcore = size(T.core);
relerr = norm(tensor(T)-C)/norm(C)

% profile viewer
% profile off

% HOOI with rank from adaptive
% T3 = tucker_als(C,Rcore);
% 
%  relerr3 = norm(tensor(T3)-C)/norm(C)
