function [V, comp_vals]= mttkrp_dt(X,U,S,n,comp_vals)
%MTTKRP in mode n of full tensor X
%   
%   [V, comp_vals]= mttkrp_dt(X,U,S,n,comp_vals) efficently comptues a MTTKRP 
%   in mode n of x using the factor matrices of U and previously comptued
%   values comp_vals in the demension tree. The ouput of V will be matrix
%   of size demension of X by rank of U and an updated cell array of the
%   comp_vals. At most, each call to mttkrp_dt will add two new entries
%   into the comp_vals cell array.
%
%

%% Variable Declarations

% move outside of funciton ASAP
N = size(U,2);
% dims = size(X);
% 
% % finds split node
% approx_root = sqrt(prod(dims));
% list = find(cumprod(dims) <= approx_root);
% S = list(end) + 1;

%% Main Demension Tree Computations

if n == 1
    % LEFT PARTIAL TENSOR 
    T = partialMTTKRPNEW(X, U, S, 1);
    comp_vals{1} = T;
    
    % Multittv for first MTTKRP result
    V = multiTTVResult(T,U(1:S-1));
    comp_vals{2} = V;

elseif n < S-1
    % Internal node update
    T = multiTTVUpdate(comp_vals{end-1},U{n-1});
    comp_vals{end+1} = T;
    
    % multiTTV for MTTKRP result
    V = multiTTVResult(T,U(n:S-1));
    comp_vals{end+1} = V;

elseif n == S-1
    % multiTTV for MTTKRP result
    T = multiTTVUpdate(comp_vals{end-1},U{n-1});
    V = double(T);
    comp_vals{end+1} = V;

elseif n == S
    % RIGHT PARTIAL TENSOR 
    T = partialMTTKRPNEW(X, U, S, 2);
    comp_vals{end+1} = T;
    
    % multiTTV for MTTKRP result
    V = multiTTVResult(T,U(n:N));
    comp_vals{end+1} = V;

elseif n < N
   % Internal node update
   T = multiTTVUpdate(comp_vals{end-1}, U{n-1});
   comp_vals{end+1} = T;
   
   % multiTTV for MTTKRP result
   V = multiTTVResult(T, U(n:N));
   comp_vals{end+1} = V;

else
   T = multiTTVUpdate(comp_vals{end-1},U{n-1});
   V = double(T);
   comp_vals{end+1} = V;
end


end
