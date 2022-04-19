function [theta] = matrix_construct6(rest_sig, M, rnd, fs)
    N = length(rest_sig);
    theta = zeros(M, N);
    t = linspace(0, (N - 1)/fs, N);
    for l = 1:N
        f = ((l-1)/N)*fs;
        x = sin(2*pi*f*t);
%         ideal_x = x(1:fs/freq_PNS:end-1);
        
        prod = zeros(1, N);  %%Discrete time Modulated Signal
        i = 1;
     
        while(i <= N)
            prod(i) = x(i)*rnd(1+mod(i,length(rnd)));
            i = i+1;
        end
        
        obs_samples = zeros(1,M);CF = N/M;
        for i = 1:M
            obs_samples(i) = sum(prod((i-1)*CF+1:CF*i)); 
        end
        
%         obs_samples = (max(abs(obs_samples))/(2^(Resolution-1)))*round((2^(Resolution-1))*obs_samples/max(abs(obs_samples)));
        theta(:,l) = obs_samples';
    end
end

