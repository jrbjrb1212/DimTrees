clear
maxNumCompThreads(1)
addpath('tensor_toolbox-master/','HOOI')
format short e

files = dir('/home/liz20/HCCI/*.mpi');

T = zeros(672,672,33,626);
for k = 1:length(files)
    str = ['/home/liz20/HCCI/',files(k).name];
    fileID = fopen(str);
    x = fread(fileID,'double');
    
    timestep = reshape(x,672,672,33);
    
    T(:,:,:,k) = timestep;
end

% result is 672 x 672 x 33 x 626
T = tensor(T);


r_init = [10,10,10,10];
tol = 1e-2;
maxiter = 10;

%profile on
tic; [H,iter,times_hooi] = hooi_adapt(T,r_init,tol,maxiter,1:4); t_hooi = toc;
times_hooi
total_hooi = sum(times_hooi)
%profile viewer

tic; T2 = hosvd(T,tol); t_st = toc
r_st = size(T2.core);

% rerun adaptive with new ranks  + 10%
tic; [H2,iter2,times_hnew] = hooi_adapt(T,ceil(1.1*r_st),tol,maxiter,1:4); t_hnew = toc;
times_hnew
total_hnew = sum(times_hnew)


% angles between subspaces
for i = 1:4
    [theta{i},S{i}] = subspace_angle(H.U{i},T2.U{i});
end

save('angles.mat','theta','S')
