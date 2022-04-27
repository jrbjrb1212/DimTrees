function V = mttkrp_dt_all(X,U)
%MTTKRP in all modes of X
%   
%   V = mttkrp_dt_all(X,A) efficently comptues the MTTKRP in each mode of X
%   using the factor matrices of U. The ouput of V will be a cell array
%   containing all the computed MTTKRPs. The size of V will be equal to the
%   amount of factor matrices in U.
%
%   Strucuture for efficent MTTKRP computation per work of Eswar, 
%   S., Hayashi, K., Ballard, G., Kannan, R., Matheson, M. A., & Park, H. 
%   (2019). PLANC: Parallel low rank approximation with non-negativity 
%   constraints. arXiv preprint arXiv:1909.01149.

%   Also add more comments

%% Variable Declarations
N = size(U,2);
V = cell(1,N);
dims = size(X);

% Efficently finds split node
% A split node is a mode of the N demension that approximately splits the
% number of tensor entries into two by modes
approx_root = sqrt(prod(dims));
list = find(cumprod(dims) <= approx_root);
S = list(end) + 1;


%% Main Demension Tree Computations
for n = 1:N
    if n == 1
        % LEFT PARTIAL TENSOR 
        T = partialMTTKRPNEW(X, U, S, 1);
        % Multittv for first MTTKRP result
        V{n} = multiTTVResult(T,U(n+1:S-1));
        
    elseif n < S-1
        % Internal node update for n-1 node
        T = multiTTVUpdate(T,U{n-1});
        % Multittv for MTTKRP result
        V{n} = multiTTVResult(T,U(n+1:S-1));

    elseif n == S-1
        % Last mode of left half of tree
        T = multiTTVUpdate(T,U{n-1});
        % Convert the 2-D tensor T into a matrix as it is a computed MTTKRP
        V{n} = double(T);

    elseif n == S
        % RIGHT PARTIAL TENSOR 
        T = partialMTTKRPNEW(X, U, S, 2);
        % multiTTV for MTTKRP result
        V{n} = multiTTVResult(T,U(n+1:N));

    elseif n < N
       % Internal node update for n-1 node
       T = multiTTVUpdate(T, U{n-1});
       % multiTTV for MTTKRP result
       V{n} = multiTTVResult(T, U(n+1:N));

    else
       % Last mode of right half of tree
       T = multiTTVUpdate(T,U{n-1});
       % Convert the 2-D tensor T into a matrix as it is a computed MTTKRP
       V{n} = double(T);
    end
end

end
