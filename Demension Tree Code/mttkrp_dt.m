function V = mttkrp_dt_all(X,U)
%MTTKRP in all modes of X
%   
%   V = mttkrp_dt_all(X,A) efficently comptues the MTTKRP in each mode of X
%   using the factor matrices of U. The ouput of V will be a cell array
%   containing all the computed MTTKRPs. The size of V will be equal to the
%   amount of factor matrices in U.
%
%

%% Variable Declarations
N = size(U,2);
V = cell(1,N);
dims = size(X);

% finds split node
approx_root = sqrt(prod(dims));
list = find(cumprod(dims) <= approx_root);
S = list(end) + 1;


%% Main Demension Tree Computations
for n = 1:N
    if n == 1
        % LEFT PARTIAL TENSOR 
        T = partialMTTKRPNEW(X, U, S, 1);
        % Multittv for first MTTKRP result
        V{n} = multiTTVResult(T,U(1:S-1));

    elseif n < S-1
        % Internal node update
        T = multiTTVUpdate(T,U{n-1});
        V{n} = multiTTVResult(T,U(n:S-1));

    elseif n == S-1
        T = multiTTVUpdate(T,U{n-1});
        V{n} = double(T);

    elseif n == S
        % RIGHT PARTIAL TENSOR 
        T = partialMTTKRPNEW(X, U, S, 2);
        % multiTTV for MTTKRP result
        V{n} = multiTTVResult(T,U(n:N));

    elseif n < N
       % Internal node update
       T = multiTTVUpdate(T, U{n-1});
       V{n} = multiTTVResult(T, U(n:N));

    else
       T = multiTTVUpdate(T,U{n-1});
       V{n} = double(T);
    end
end

end
