function result = sinc_trunc(t,Ts,q)
% Generates a truncated sinc pulse from -qTs to qTs of amplitude sqrt(1/Ts)

tol = 0.0001 * Ts;

cond1 = abs(t) - q*Ts > tol;
cond2 = not(cond1);

result = 0 * t;
result(cond1) = 0;
result(cond2) = sqrt(1/Ts) * sinc(t(cond2)./Ts);