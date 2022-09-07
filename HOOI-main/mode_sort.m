function order = mode_sort(r,k,modes)
% determines order for multi TTM: Y = G x U
% r: size of G
% k: size of Y

d = length(modes);
if d <= 1
    order = modes;
    return
else

    j = modes(randi(d));
    inds = modes;
    lhs = k(inds).*r(inds).*r(j)+k(inds).*r(j).*k(j);
    rhs = r(inds).*r(j).*k(j)+r(inds).*k(inds).*k(j);
    
    L = mode_sort(r,k,modes(lhs < rhs)); % modes to go before
    M = modes(lhs == rhs);
    R = mode_sort(r,k,modes(lhs > rhs)); % modes to go after
    
    order = [L, M, R];
end


