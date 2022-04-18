% does not work for unbalanced trees

function [V, temp_vals]= mttkrp_dt(X,U,split_node,n,temp_vals)
%MTTKRP in mode n of full tensor X
%   
%   [V, temp_vals]= mttkrp_dt(X,U,S,n,temp_vals) efficently comptues a MTTKRP 
%   in mode n of x using the factor matrices of U and previously comptued
%   values temp_vals in the demension tree. The ouput of V will be matrix
%   of size demension of X by rank of U and an updated cell array of the
%   temp_vals. At most, each call to mttkrp_dt will add two new entries
%   into the temp_vals cell array.
%
%   Strucuture for efficent MTTKRP computation per work of Eswar, 
%   S., Hayashi, K., Ballard, G., Kannan, R., Matheson, M. A., & Park, H. 
%   (2019). PLANC: Parallel low rank approximation with non-negativity 
%   constraints. arXiv preprint arXiv:1909.01149.
%

%% Variable Declarations

N = size(U,2);

%% Main Demension Tree Computations

if n == 1
    % LEFT PARTIAL TENSOR 
    T = partialMTTKRPNEW(X, U, split_node, 1);
    temp_vals{n} = T;
    
    % Multittv for first MTTKRP result
    V = multiTTVResult(T,U(2:split_node-1));

elseif n < split_node-1
    % Internal node update
    T = multiTTVUpdate(temp_vals{n-1},U{n-1});
    temp_vals{n} = T;
    
    % multiTTV for MTTKRP result
    V = multiTTVResult(T,U(n+1:split_node-1));

elseif n == split_node-1
    % multiTTV for MTTKRP result
    T = multiTTVUpdate(temp_vals{n-1},U{n-1});
    % Convert the 2-D tensor T into a matrix as it is a computed MTTKRP
    V = double(T);

elseif n == split_node
    % RIGHT PARTIAL TENSOR 
    T = partialMTTKRPNEW(X, U, split_node, 2);
    temp_vals{n} = T;
    
    % multiTTV for MTTKRP result
    V = multiTTVResult(T,U(n+1:N));

elseif n < N
   % Internal node update
   T = multiTTVUpdate(temp_vals{n-1}, U{n-1});
   temp_vals{n} = T;
   
   % multiTTV for MTTKRP result
   V = multiTTVResult(T, U(n+1:N));

else
   % Internal node update
   T = multiTTVUpdate(temp_vals{n-1},U{n-1});
   % Convert the 2-D tensor T into a matrix as it is a computed MTTKRP
   V = double(T);
end


end
