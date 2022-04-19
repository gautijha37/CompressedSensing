function [PRBS_mat] = construct_PRBS_mat(N, M, product_length, rnd)

% Constructs M by N measurement matrix with +-1 bernoulli entries
% product length = no. of consecutive nonzero terms in a row of PRBS
% matrix, preferably a power of 2

% M = N/CF;

% % n = 7;  %% LFSR bit length
% ps = zeros(1,n);
% ns = zeros(1,n); 
end_len = min(length(rnd) , N);
% rnd = zeros(1,end_len);
% 
% for i = 1:end_len
%     ns(2:n) = ps(1:n-1) ;
%     ns(1) = ~xor(ps(n),ps(1));
%     rnd(i) = 2*ps(1) - 1; % This is correct! verified with LTspice !!
%     ps = ns;
% end
% rnd = repmat(rnd, [1, M * product_length/length(rnd)]);


PRBS_mat = zeros(M, N); j = 1; k = 0;
for i = 1 : M
%     end_coord = min((i - 1)*product_length + product_length, N);
    if (j + product_length - 1 > end_len)
        PRBS_mat(i,k + 1 : k + product_length) = [rnd(j : end_len) , rnd(1 : j + product_length - 1 - end_len)];
        j = j + product_length - end_len;
    else
        PRBS_mat(i,k + 1 : k + product_length) = rnd(j : j + product_length - 1);
        j = j + product_length;
    end
    
    k = mod(k + product_length, N);
end

end

