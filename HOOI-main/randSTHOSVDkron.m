function T = randSTHOSVDkron(X,ranks,randmat)
% runs randomized range finder in each mode in sequence (like STHOSVD)
% using Kronecker product of random sketch matrices (chosen 
% independently in each mode)

    % store dimensions and properties
    dims = size(X);
    ndims = length(dims);
    
    % compute subranks
%     subranks = ceil(prod(ranks)^(1/(ndims-1)) ./ ranks);
    subranks = getSubRanks(ranks);
    % make sure subranks are powers of 2 for SRHT matrix
    if strcmp(randmat,'srht')
        ready = false;
        while ~ready
            if mod(subranks,2) == zeros(1,ndims)
                ready = true;
            end
                idx = find(mod(subranks,2));
                subranks(idx) = subranks(idx)+1;
        end
    end
    
    % initialize factor matrices
    U = cell(ndims,1);
    
    % initialize random matrices
    Phi = cell(ndims,1);
    
    G = X;
    for n = 1:ndims
        % generate random matrices (new ones for each mode)
        for j = 1:ndims
            if j < n
                %numrows = ranks(j);
                numrows = prod(subranks([1:j-1,j+1:ndims]));
                if numrows > dims(j)
                    numrows = dims(j);
                end
            else
                numrows = dims(j);
            end
            switch randmat
                case 'gaussian'
                    Phi{j} = randn(numrows, subranks(j));
                case 'srht'
                    Phi{j} = srht(numrows, subranks(j));
            end

        end
        % compute sketch
        Y = ttm(G,Phi,-n,'t'); 
        % compute factor matrix
        [U{n},~] = qr(double(tenmat(Y,n)),0);
        % truncate in mode n
        G = ttm(G,U{n},n,'t');
    end
    
    % return Tucker tensor
    T = ttensor(G,U);

end