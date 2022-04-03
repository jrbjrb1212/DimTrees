function T = UpdateTensor(Z,X,n)

if (n == 1)
    R = size(X{2},2);
else
    R = size(X{1},2);
end


N = ndims(Z);
Updated_Tensor = tenzeros(size(Z,1:ndims(Z)-1));

if N == 3
    for r = 1:R
        Update = ttv(Z, X{n}(:,r),1);
        Updated_Tensor(:,r) = Update(:,r);
       
    end
elseif N == 4
    for r = 1:R
        Update = ttv(Z, X{n}(:,r),n);
        Updated_Tensor(:,:,r) = Update(:,:,r);
    end
elseif N == 5
    for r = 1:R
        Update = ttv(Z, X{n}(:,r),n);
        Updated_Tensor(:,:,:,r) = Update(:,:,:,r);
    end
elseif N == 6
    for r = 1:R
        Update = ttv(Z, X{n}(:,r),n);
        Updated_Tensor(:,:,:,:,r) = Update(:,:,:,:,r);
    end
elseif N == 7
    for r = 1:R
        Update = ttv(Z, X{n}(:,r),n);
        Updated_Tensor(:,:,:,:,:,r) = Update(:,:,:,:,:,r);
    end     
elseif N == 8
    for r = 1:R
        Update = ttv(Z, X{n}(:,r),n);
        Updated_Tensor(:,:,:,:,:,:,r) = Update(:,:,:,:,:,:,r);
    end
elseif N == 9
    for r = 1:R
        Update = ttv(Z, X{n}(:,r),n);
        Updated_Tensor(:,:,:,:,:,:,:,r) = Update(:,:,:,:,:,:,:,r);
    end
elseif N == 10
    for r = 1:R
        Update = ttv(Z, X{n}(:,r),n);
        Updated_Tensor(:,:,:,:,:,:,:,:,r) = Update(:,:,:,:,:,:,:,:,r);
    end
elseif N == 11
    for r = 1:R
        Update = ttv(Z, X{n}(:,r),n);
        Updated_Tensor(:,:,:,:,:,:,:,:,:,r) = Update(:,:,:,:,:,:,:,:,:,r);
    end
elseif N == 12
    for r = 1:R
        Update = ttv(Z, X{n}(:,r),n);
        Updated_Tensor(:,:,:,:,:,:,:,:,:,:,r) = Update(:,:,:,:,:,:,:,:,:,:,r);
    end
elseif N == 13
    for r = 1:R
        Update = ttv(Z, X{n}(:,r),n);
        Updated_Tensor(:,:,:,:,:,:,:,:,:,:,:,r) = Update(:,:,:,:,:,:,:,:,:,:,:,r);
    end

elseif N == 14
    for r = 1:R
        Update = ttv(Z, X{n}(:,r),n);
        Updated_Tensor(:,:,:,:,:,:,:,:,:,:,:,:,r) = Update(:,:,:,:,:,:,:,:,:,:,:,:,r);
    end

elseif N == 15
    for r = 1:R
        Update = ttv(Z, X{n}(:,r),n);
        Updated_Tensor(:,:,:,:,:,:,:,:,:,:,:,:,:,r) = Update(:,:,:,:,:,:,:,:,:,:,:,:,:,r);
    end
else
     T = null; 
     return
end

T = Updated_Tensor;
end
