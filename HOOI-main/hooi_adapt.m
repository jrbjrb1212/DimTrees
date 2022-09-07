function [T,itercount,times] = hooi_adapt(X,rk,varargin)
% rk: initial guess for ranks
% Optional inputs:
%  - modes: processing order (default 1:d)
%  - err_tol: desired error tolerance (default 1e-4)
%  - maxiter: max number of iterations (default 50)
%  - init: how to initialize factor matrices, choice of 'rand', 'hooi' (one
%  full iter of regular hooi), 'sthosvd', 'rsthosvdkron' (randomized STHOSVD with Kronecker
%  products) (default rand)
%  - print: boolean for displaying iteration info (default true)
%
% times: cell array of iteration times (vector), initialization time, total time 

normXsqr = norm(X)^2;
dims = size(X);
d = length(dims);

% process inputs
params = inputParser;
params.addParameter('err_tol',1e-4,@isscalar);
params.addParameter('maxiter',50,@(x) isscalar(x) & x > 0);
params.addParameter('modes',1:d,@(x) isequal(sort(x),1:d));
params.addParameter('init','rand',@(x) (iscell(x) || ismember(x,{'rand','hooi','rsthosvdkron','sthosvd'})));
params.addParameter('print',true,@islogical);
params.parse(varargin{:});

err_tol = params.Results.err_tol;
maxiter = params.Results.maxiter;
modes = params.Results.modes;
init = params.Results.init;
print = params.Results.print;

% initialize factor matrices
tic;
if strcmp(init,'rand')
    U = cell(1,d);
    for i = modes
        M = randn(dims(i),max(rk));
        [U{i},~] = qr(M,0);
    end
elseif strcmp(init,'hooi')
    T = tucker_als(X,rk,'tol',err_tol,'maxiters',1,dimorder,'modes');
    U = T.U;
elseif strcmp(init,'sthosvd')
    T = hosvd(X,err_tol,'ranks',rk);
    U = T.U;
elseif strcmp(init,'rsthosvdkron')
    T = randSTHOSVDkron(X,rk,'gaussian');
    U = T.U;
end
t_init = toc;
    

itercount = 0;
rk_old = rk;
rkstr = repmat('%i ',1,d);
t_iter = [];

% main iterations
for j = 1:maxiter
    tic;
    itercount = itercount + 1;
    for k = modes
        
        % determine best order for TTM
        ttm_modes = mode_sort(size(X),rk,[1:k-1,k+1:d]);
       
        B = ttm(X,U,ttm_modes,'t');
        disp(B)
        Bk = double(tenmat(B,k));
        
        % compute singular values using Gram
        Gm = Bk*Bk';
        [Uk,S] = eig(Gm);
        [~,p] = sort(diag(S),'descend');
        
        % update factor matrix
        U{k} = Uk(:,p(1:rk(k)));
        
    end 
    
    % compute core after all modes iterated through
    G = ttm(B,U,k,'t'); 
    
    % square entries
    Gsqr = double(G.^2);
    
    % check error
    nrmGsqr = sum(Gsqr(:));
    relerr_sqr = abs((nrmGsqr-normXsqr)/normXsqr);

    if print == true
        outputstr = ['Iter %i: rk =[ ' rkstr '], relerr = %7.6e\n'];
        fprintf(outputstr,j,rk,sqrt(relerr_sqr));
    end 
        
    % determine new 30
    for i = 1:d
        Gsqr = cumsum(Gsqr,i);
    end
    
    if Gsqr(end) < (1-err_tol^2)*normXsqr %if all entries don't satisfy tolerance, take 1.5x the size
        rk = min(dims,ceil(1.5*rk_old));
    else
        temp = Inf;
        % loop through all entries of cumsum to find closest index match
        % for error tolerance
        for q = 1:prod(size(G))
            if Gsqr(q) >= (1-err_tol^2)*normXsqr && Gsqr(q) < temp
                ind = q;
                temp = Gsqr(q);
            end
        end
        % convert to subscripts (better way of doing this?)
        switch d
            case 3
                [rk(1),rk(2),rk(3)] = ind2sub(size(G),ind);
            case 4 
                [rk(1),rk(2),rk(3),rk(4)] = ind2sub(size(G),ind);
            case 5
                [rk(1),rk(2),rk(3),rk(4),rk(5)] = ind2sub(size(G),ind);
        end
    end
    
        
    % stop if rank isn't changing
    if rk_old == rk
        disp('rank stopped changing')
        break
    end

    rk_old = rk;
    t = toc;
    t_iter = [t_iter,t];
end

T = ttensor(G,U);
t_total = sum(t_iter)+t_init;
times = {t_iter,t_init,t_total};
end
