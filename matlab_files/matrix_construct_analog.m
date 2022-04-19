function [theta_mat] = matrix_construct_analog(fs,T,freq_PNS,rnd,N,M,Resolution)
    t = 0:1/fs:T ; theta_mat = zeros(M,N);
    for l = 1:N
        f = ((l-1)/N)*freq_PNS;
        x = sin(2*pi*f*t);
        
        prod = zeros(1,length(t));  %% Analog Modulated Signal
        i = 1;
        while(i <= length(t))
            j = 1;
         while(j <= length(rnd) && i <= length(t))
            for k = 1 : fs/freq_PNS
                if(i > length(t))
                    break;
                end
                prod(i) = x(i)*rnd(j);
                i = i+1;
            end
            j = j+1;
         end
        end
        
        obs_samples = zeros(1,M);
        i = 1; j = 1;
        while (i < length(t)) % length(t) = T*fs + 1
            obs_samples(j) = sum(prod(i : i + T*fs/M - 1));
            j = j + 1;
            i = i + T*fs/M ;
        end
        obs_samples = obs_samples/fs;
        obs_samples = (max(abs(obs_samples))/(2^(Resolution-1)))*round((2^(Resolution-1))*obs_samples/max(abs(obs_samples)));
        theta_mat(:,l) = obs_samples';
    end
    
end