function M = findExact(x,A)
n = length(x); M = zeros(n,1);
for i = 1 : n
    a = A(i); 
    if x(i) <= 2.1, a_s = 1; else, a_s = 1.4531; end
    g = @(m) ((1/1.2*(1+0.2*m.^2)).^6) ./ (m.^2) - (a/a_s)^2;
    
    if x(i)<=1.5, Me=0.5; elseif x(i)<2.1, Me=2; else Me=0.5; end
    M(i) = fzero(g,Me);
end
end

