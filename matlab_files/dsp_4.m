function [sparse_indices,coeff] = dsp_4(phi_id, phi,y,sparsity)
    
    r = y; sparse_indices = zeros(1,sparsity); 
    [m, n] = size(phi);
    A_dash = zeros(m,n); coeff = zeros(1,n);
    
    for j = 1:sparsity
        [~,mx_id] = max(abs(phi_id' * r));
        sparse_indices(j) = mx_id ;
        A_dash(:,j) = phi(:,mx_id) ;
        A = A_dash(:,[1:j]);
        P = (inv(A' * A))*(A') ;
        estimate = P*y;
        r = y - A*estimate ;
      
    end
    
    for i = 1:sparsity
        coeff(sparse_indices(i)) = estimate(i);
    end
    
end