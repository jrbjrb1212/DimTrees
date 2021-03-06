function [f,G] = tt_cp_fg(Z,A,Znormsqr,varargin)
%TT_CP_FG Computes function and gradient of the CP function.
%
%   [F,G] = TT_CP_FG(Z,A) calculates F = (1/2) ||Z - ktensor(A)||^2 where
%   Z is an N-way tensor and A is a ktensor or a cell array with N
%   factor matrices. It also calculates the gradient of the CP fit
%   function where Z is an N-way tensor and A is a ktensor or a
%   cell array with N factor matrices. The result is also a cell
%   array with N factor matrices corresponding to the gradients; in
%   other words, G{n}(:,r) is the partial derivative of the fit
%   function with respect to A{n}(:,r). 
%
%   [F,G] = TT_CP_FG(Z,A,NORMZSQR) also passes in the pre-computed
%   norm of Z, which makes the computations faster. 
%
%   See also CP_OPT, TT_CP_FUN.
%
%MATLAB Tensor Toolbox. Copyright 2018, Sandia Corporation.


%% Set algorithm parameters from input or by using defaults
params = inputParser;
params.addParameter('dt',false,@islogical);
params.parse(varargin{:});


%% Set-up
% if ~isa(Z,'tensor') && ~isa(Z,'sptensor')
%     error('Z must be a tensor or a sptensor');
% end
N = ndims(Z);

if ~iscell(A) && ~isa(A,'ktensor');
    error('A must be a cell array or ktensor');
end

if isa(A,'ktensor')
    A = tocell(A);
end
R = size(A{1},2);

%% Upsilon and Gamma
Upsilon = cell(N,1);
for n = 1:N
    Upsilon{n} = A{n}'*A{n};
end

Gamma = cell(N,1);
Gamma{1} = 1;
for n = 2:n
    Gamma{n} = Gamma{n-1} .* Upsilon{n-1};
end

right = 1;
for n=N:-1:1
    Gamma{n} = Gamma{n} .* right;
    right = right .* Upsilon{n};
end


%% Calculation

%F1
if exist('Znormsqr','var')
    f_1 = Znormsqr;
else
    f_1 = norm(Z)^2;
end

%% Calculate gradient and F2
G = cell(N,1);

% Default Implementation
if ~params.Results.dt
    U = mttkrp(Z,A,1, 2);
    V = A{1} .* U;
    f_2 = sum(V(:));
    G{1} = -U + A{1}*Gamma{1};
    for n = 2:N
        U = mttkrp(Z,A,n, 2);
        G{n} = -U + A{n}*Gamma{n};
    end
    
% Demension Tree Option for better performance
else
    % Computes MTTKRP in all modes and stores in cell array
    % Index within cell is same as mode in tensor
    U = mttkrp_dt_all(Z,A);

    V = A{1} .* U{1};
    f_2 = sum(V(:));
    G{1} = -U{1} + A{1}*Gamma{1};

    % Compute gradient using precomputed MTTKRPs
    for n = 2:N
        G{n} = -U{n} + A{n}*Gamma{n};
    end
end

%F3
W = Gamma{1} .* Upsilon{1};
f_3 = sum(W(:));

%SUM
f = 0.5 * f_1 - f_2 + 0.5 * f_3;


