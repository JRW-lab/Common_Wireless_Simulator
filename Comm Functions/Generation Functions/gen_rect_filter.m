function result = gen_rect_filter(t,Ts,~)
% Generates a rectangular pulse of amplitude sqrt(1/Ts) from 0 to Ts
% (the third input can be taken off, I'm just too lazy to change all my
% code)

tol = 0.0001 * Ts;

cond1 = t < tol | (t - Ts) >= tol;
cond2 = not(cond1);

result(cond1) = 0;
result(cond2) = sqrt(1/Ts);
