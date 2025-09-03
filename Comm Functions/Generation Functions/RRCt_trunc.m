function result = RRCt_trunc(t,alpha,Ts,q)
% Generates a VERSION OF the RRC pulse (there are other definitions but
% this is valid)

tol = 0.0001 * Ts;

result = 0 * t;
R = 0.5/Ts;

cond1 = abs(t) < tol;
cond2 = abs(abs(t) - (64 * alpha^2 * R^2)^(-1/2)) < tol;
cond3 = (abs(t) - q*Ts) > tol;
cond4 = not(cond1 | cond2 | cond3);

result(cond1) = sqrt(2.*R) .* (1 - alpha + 4.*alpha./pi); % check this is correct
result(cond2) = (sqrt(2.*R) ./ (2.*pi.*R.*(1 - 192.*(alpha.^2).*(R.^2).*(t(cond2).^2)))).*(2.*pi.*R.*(1-alpha).*cos(2.*pi.*R.*(1-alpha).*t(cond2)) + 8.*R.*alpha.*(cos(2.*pi.*R.*(1+alpha).*t(cond2)) - 2.*pi.*R.*(1+alpha).*t(cond2).*sin(2.*pi.*R.*(1+alpha).*t(cond2))));
result(cond3) = 0;
result(cond4) = (sqrt(2.*R) ./ (1 - 64.*alpha.^2.*R.^2.*t(cond4).^2)) .* (sin(2.*pi.*R.*(1-alpha).*t(cond4))./(2.*pi.*R.*t(cond4)) + 4.*alpha./pi .* cos(2.*pi.*R.*(1+alpha).*t(cond4)));
