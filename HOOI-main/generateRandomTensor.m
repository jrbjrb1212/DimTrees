function [T, sv] = generateRandomTensor(d, coreDim, tensorDim, seed, ratio)

% Generate a cubic random tensor with a cubic core 
% ratio: If specified the singular values of the unfolding on each mode 
%   will be the same and decrease geometrically at a rate specified by ratio.
%   If ratio is not specified the core is random.
% d: The number of modes of the output tensor.
% coreDim: The size of the core tensor as a scalar. 
% tensorDim: The size of the output tensor as a scalar.

    rng(seed);
    assert(coreDim<=tensorDim, "coreDim must be less than or equal to tensorDim.");
    coreSize = ones(1,d).*coreDim;

    if nargin < 4
        % random core
        core = randn(coreSize);
    else
        % diagonal core with specified s-values
        core = zeros(coreSize);
        sv = (1/ratio)*ratio.^(1:coreDim);
        for i=1:coreDim
          c = num2cell(ones(1,d).*i);
          core(c{:}) = sv(i); 
        end
    end
    
    % random orthogonal factor matrices
    U = cell(1,d);
    for i = 1:d
        [Utemp,~] = qr(randn(tensorDim,tensorDim));
        U{i} = Utemp(:, 1:coreDim);
    end 
    T = ttensor(tensor(core), U);
end