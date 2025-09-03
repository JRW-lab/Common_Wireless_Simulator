function Gamma_MN = gen_Gamma_MN(M,N)

I = eye(M*N);
Gamma_m = zeros(M*N,N);
Gamma_MN = zeros(M*N);
for m = 0:M-1
    for n = 0:N-1
        Gamma_m(:,n+1) = I(:,m+(n*M)+1);
    end
    Gamma_MN(:,m*N+1:(m+1)*N) = Gamma_m(:,:);
end