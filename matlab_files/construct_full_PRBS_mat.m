function [mat] = construct_full_PRBS_mat(M,N, n)
% n = PRBS bit length
% Constructs fully filled M by N measurement matrix with +-1 bernoulli entries

% n = 7;  %% LFSR bit length
% ps = zeros(1,n);
% ns = zeros(1,n); 
% end_len = 2^n - 1;
% rnd = zeros(1,end_len);
% 
% for i = 1:end_len
%     ns(2:n) = ps(1:n-1) ;
%     ns(1) = ~xor(ps(n),ps(1));
%     rnd(i) = 2*ps(1) - 1; % This is correct! verified with LTspice !!
%     ps = ns;
% end

range = [-1, 1];
band = [0 1];
rnd = idinput(2^n - 1, 'prbs', band, range);

mat = zeros(M, N); 
for i = 1 : M
    for j = 1 : N
        mat(i, j) = rnd( mod((i - 1)*N + j, length(rnd)) + 1);
    end
end

end

