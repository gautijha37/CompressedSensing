function [PRBS_mat, rnd] = construct_PRBS_mat_11b(N, CF)

% Constructs M by N measurement matrix with +-1 bernoulli entries

M = N/CF;

n = 11;  %% LFSR bit length
ps = zeros(1,n);
ns = zeros(1,n); 
end_len = min(2^n - 1, N);
rnd = zeros(1,end_len);

for i = 1:end_len
    ns(2:n) = ps(1:n-1) ;
    ns(1) = ~xor(ps(11), ps(9));
    rnd(i) = 2*ps(1) - 1; % This is correct! verified with LTspice !!
    ps = ns;
end

PRBS_mat = zeros(M, N); j = 1;
for i = 1 : M
    if (j + CF - 1 > end_len)
        PRBS_mat(i,(i - 1)*CF + 1 : (i - 1)*CF + CF) = [rnd(j : end) , rnd(1 : j + CF - 1 - end_len)];
        j = j + CF - end_len;
    else
        PRBS_mat(i,(i - 1)*CF + 1 : (i - 1)*CF + CF) = rnd(j : j + CF - 1);
        j = j + CF;
    end
end

end

