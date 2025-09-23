function a_t = gen_pulse(t_range, shape, Ts, Q, alpha)
    switch shape
        case "rect"
            a_t = gen_filter_rect(t_range, Ts);
        case "sinc"
            a_t = gen_filter_sinc(t_range, Ts, Q, true);
        case "rrc"
            a_t = gen_filter_RRCt(t_range, Ts, Q, true, alpha);
        otherwise
            error('Unknown pulse shape');
    end
end

function result = gen_filter_rect(t,Ts)
% Generates a rectangular pulse of amplitude sqrt(1/Ts) from 0 to Ts
%
% JRW 3/7/25

tol = 0.0001 * Ts;

cond1 = t < tol | (t - Ts) >= tol;
cond2 = not(cond1);

result = 0 * t;
result(cond1) = 0;
result(cond2) = sqrt(1/Ts);

end

function result = gen_filter_sinc(t,Ts,Q,causal)
% Generates a truncated sinc pulse from -qTs to qTs of amplitude sqrt(1/Ts)
% The time support is shifted to only include t >=0 when causal is TRUE
%
% JRW 3/7/25

if causal
    t = t - (Q*Ts/2);
end

tol = 0.0001 * Ts;

cond1 = abs(t) - Q*Ts/2 > tol;
cond2 = not(cond1);

result = 0 * t;
result(cond1) = 0;
result(cond2) = sqrt(1/Ts) * sinc(t(cond2)./Ts);

end

function result = gen_filter_RRCt(t,Ts,Q,causal,alpha)
% Generates a VERSION OF the RRC pulse
%     (there are other definitions but this is valid)
% The time support is shifted to only include t >=0 when causal is TRUE
%
% JRW 3/7/25

if causal
    t = t - (Q*Ts/2);
end

tol = 0.0001 * Ts;

result = 0 * t;
R = 0.5/Ts;

cond1 = abs(t) < tol;
cond2 = abs(abs(t) - (64 * alpha^2 * R^2)^(-1/2)) < tol;
cond3 = (abs(t) - Q*Ts/2) > tol;
cond4 = not(cond1 | cond2 | cond3);

result(cond1) = sqrt(2.*R) .* (1 - alpha + 4.*alpha./pi); % check this is correct
result(cond2) = (sqrt(2.*R) ./ (2.*pi.*R.*(1 - 192.*(alpha.^2).*(R.^2).*(t(cond2).^2)))).*(2.*pi.*R.*(1-alpha).*cos(2.*pi.*R.*(1-alpha).*t(cond2)) + 8.*R.*alpha.*(cos(2.*pi.*R.*(1+alpha).*t(cond2)) - 2.*pi.*R.*(1+alpha).*t(cond2).*sin(2.*pi.*R.*(1+alpha).*t(cond2))));
result(cond3) = 0;
result(cond4) = (sqrt(2.*R) ./ (1 - 64.*alpha.^2.*R.^2.*t(cond4).^2)) .* (sin(2.*pi.*R.*(1-alpha).*t(cond4))./(2.*pi.*R.*t(cond4)) + 4.*alpha./pi .* cos(2.*pi.*R.*(1+alpha).*t(cond4)));

end