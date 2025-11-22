function result = ambig(t,f,Ts,shape,alpha)
% This function results the result of the ambiguity function of two
% of the same shaped pulse filters. Supported shapes are "rect", "sinc" and
% "rrc". (alpha is only used for "rrc" and can usually be left out)
%
% Instructions:
% 1. t:         enter an m-by-n matrix of time values (in seconds)
% 2. f:         enter an m-by-n matrix of frequency values (in Hertz)
% 3. Ts:        enter symbol period (in seconds)
% 4. shape:     enter shape of pulse filters ("rect"/"sinc"/"rrc")
% 5. alpha:     enter value of roll-off factor
%                   (unused if "rrc" not selected)
%
% Output: an m-by-n sized matrix of ambiguity values, with each value
%   calculated element-wise from the given t and f matrices
%
% Coded by Jeremiah Rhys Wimer, 4/23/2024

% Tolerance for exceptions
tol = .001 * Ts;

% Premake result, reset t and f, set tolerance
result = t * 0;

% Make result for each shape
switch shape
    case "rect"
        % Define f=0 condition
        f_zero = abs(f) < tol;

        % Typical condition definitions
        cond1n = not(f_zero) & (t>0 & t<=Ts);
        cond2 = not(f_zero) & (t<=0 & t>=-Ts);

        % Limit condition definitions
        cond3n = f_zero & (t>0 & t<=Ts);
        cond4n = f_zero & (t<=0 & t>=-Ts);

        % Positive and negative definitions
        result(cond1n) = (exp(1j.*2.*pi.*f(cond1n).*Ts) - exp(1j.*2.*pi.*f(cond1n).*t(cond1n))) ./ (1j.*2.*pi.*f(cond1n).*Ts);
        result(cond2) = (exp(1j.*2.*pi.*f(cond2).*(Ts+t(cond2))) - 1) ./ (1j.*2.*pi.*f(cond2).*Ts);

        % Positive and negative limits as f goes to 0
        result(cond3n) = 1 - t(cond3n) / Ts;
        result(cond4n) = 1 + t(cond4n) / Ts;

    case "sinc"
        % Define f=0 condition
        t_zero = abs(t) < tol;

        % Typical condition definitions
        cond1n = not(t_zero) & (f>0 & f<=1/Ts);
        cond2 = not(t_zero) & (f<=0 & f>=-1/Ts);

        % Limit condition definitions
        cond3n = t_zero & (f>0 & f<=1/Ts);
        cond4n = t_zero & (f<=0 & f>=-1/Ts);
        
        % Positive and negative definitions
        result(cond1n) = (Ts ./ (pi.*t(cond1n))) .* sin(((pi.*t(cond1n)) ./ Ts) - pi.*f(cond1n).*t(cond1n)) .* exp(1j.*pi.*f(cond1n).*t(cond1n));
        result(cond2) = (Ts ./ (pi.*t(cond2))) .* sin(((pi.*t(cond2)) ./ Ts) + pi.*f(cond2).*t(cond2)) .* exp(1j.*pi.*f(cond2).*t(cond2));
        
        % Positive and negative limits as t goes to 0
        result(cond3n) = (1 - f(cond3n).*Ts);
        result(cond4n) = (1 + f(cond4n).*Ts);

    case "rrc"
        % Swap f for definition
        f = -f;
        t = t / Ts;
        f = f * Ts;
        Ts = 1;

        % RRC definitions
        W0 = 1/(2*Ts);
        W = (1+alpha)*W0;

        % Internal definitions for RRC function
        x1 = 2*W0-W;
        x2 = W;
        a = (pi/4) / (W-W0);
        phi1 = pi*(W-2*W0)/(4*(W-W0));
        phi2 = pi*(-f+W-2*W0)/(4*(W-W0));
        phi3 = pi*(f+W-2*W0)/(4*(W-W0));
        r = zeros(4,1);
        if alpha <= 0.5
            r(1) = x2-x1;
            r(2) = 2*x1;
        else
            r(1) = 2*x1;
            r(2) = x2-x1;
        end
        r(3) = x2+x1;
        r(4) = 2*x2;

        % Typical condition definitions
        cond1n = -r(1) <= f & f < 0;
        cond21n = -r(2) <= f & f < -r(1) & alpha >=0.5;
        cond22n = -r(2) <= f & f < -r(1) & alpha <0.5;
        cond3n = -r(3) <= f & f < -r(2);
        cond4n = -r(4) <= f & f <= -r(3);
        cond1p = 0 <= f & f < r(1);
        cond21p = r(1) <= f & f < r(2) & alpha >=0.5;
        cond22p = r(1) <= f & f < r(2) & alpha <0.5;
        cond3p = r(2) <= f & f < r(3);
        cond4p = r(3) <= f & f <= r(4);

        % Negative range definitions
        result(cond1n) = h2(t(cond1n),-x1+f(cond1n),-x2,-a,-a,phi1,phi3(cond1n))...
            + h1(t(cond1n),-x1,-x1+f(cond1n),-a,phi1)...
            + h0(t(cond1n),x1+f(cond1n),-x1)...
            + h1(t(cond1n),x1,x1+f(cond1n),a,phi2(cond1n))...
            + h2(t(cond1n),x2+f(cond1n),x1,a,a,phi1,phi2(cond1n));
        result(cond21n) = h2(t(cond21n),-x1+f(cond21n),-x2,-a,-a,phi1,phi3(cond21n))...
            + h1(t(cond21n),x1+f(cond21n),-x1+f(cond21n),-a,phi1)...
            + h2(t(cond21n),-x1,x1+f(cond21n),-a,a,phi1,phi2(cond21n))...
            + h1(t(cond21n),x1,-x1,a,phi2(cond21n))...
            + h2(t(cond21n),x2+f(cond21n),x1,a,a,phi1,phi2(cond21n));
        result(cond22n) = h1(t(cond22n),-x1,-x2,-a,phi1)...
            + h0(t(cond22n),x1+f(cond22n),-x1)...
            + h1(t(cond22n),x2+f(cond22n),x1+f(cond22n),a,phi2(cond22n));
        result(cond3n) = h1(t(cond3n),x1+f(cond3n),-x2,-a,phi1)...
            + h2(t(cond3n),-x1,x1+f(cond3n),-a,a,phi1,phi2(cond3n))...
            + h1(t(cond3n),x2+f(cond3n),-x1,a,phi2(cond3n));
        result(cond4n) = h2(t(cond4n),x2+f(cond4n),-x2,-a,a,phi1,phi2(cond4n));

        % Positive range definitions
        result(cond1p) = h2(t(cond1p),-x1,-x2+f(cond1p),-a,a,phi1,-phi3(cond1p))...
            + h1(t(cond1p),-x1+f(cond1p),-x1,a,-phi3(cond1p))...
            + h0(t(cond1p),x1,-x1+f(cond1p))...
            + h1(t(cond1p),x1+f(cond1p),x1,a,phi1)...
            + h2(t(cond1p),x2,x1+f(cond1p),a,-a,phi1,-phi2(cond1p));
        result(cond21p) = h2(t(cond21p),-x1,-x2+f(cond21p),-a,a,phi1,-phi3(cond21p))...
            + h1(t(cond21p),x1,-x1,a,-phi3(cond21p))...
            + h2(t(cond21p),-x1+f(cond21p),x1,a,a,phi1,-phi3(cond21p))...
            + h1(t(cond21p),x1+f(cond21p),-x1+f(cond21p),a,phi1)...
            + h2(t(cond21p),x2,x1+f(cond21p),a,-a,phi1,-phi2(cond21p));
        result(cond22p) = h1(t(cond22p),-x1+f(cond22p),-x2+f(cond22p),a,-phi3(cond22p))...
            + h0(t(cond22p),x1,-x1+f(cond22p))...
            + h1(t(cond22p),x2,x1,a,phi1);
        result(cond3p) = h1(t(cond3p),x1,-x2+f(cond3p),a,-phi3(cond3p))...
            + h2(t(cond3p),-x1+f(cond3p),x1,a,a,phi1,-phi3(cond3p))...
            + h1(t(cond3p),x2,-x1+f(cond3p),a,phi1);
        result(cond4p) = h2(t(cond4p),x2,-x2+f(cond4p),a,a,phi1,-phi3(cond4p));

    otherwise
        error("Unsupported shape selected! Please choose: 'rect' / 'sinc' / 'rc'")
end

end

function result = h0(t,f_pos,f_neg)
% Integrates exp(-1j*2*pi*f*t), from f1 to f2

% Set tolerance
tol = 1e-6;

% Define conditions
cond1 = abs(t) > tol;
cond2 = not(cond1);

% Define expression for usual cases
result_pos_vec = (1 ./ (-1j.*2.*pi.*t)) .* exp(-1j.*2.*pi.*f_pos.*t);
result_neg_vec = (1 ./ (-1j.*2.*pi.*t)) .* exp(-1j.*2.*pi.*f_neg.*t);
result_vec = result_pos_vec - result_neg_vec;

% Set results based on conditions
result(cond1) = result_vec(cond1);

% Define expression for abnormal cases
result_vec = f_pos - f_neg;
result(cond2) = result_vec(cond2);

end

function result = h1(t,f_pos,f_neg,w,phi)
% Integrates exp(-1j*2*pi*f*t)*cos(w*f+phi), from f1 to f2

% Set tolerance
tol = 1e-6;

% Define conditions
cond1 = abs(w.^2 - (2.*pi.*t).^2) > tol;
cond2 = not(cond1);

% Define expression for usual cases
result_pos_vec = exp(-1j.*2.*pi.*f_pos.*t) ./ (w.^2 - (2.*pi.*t).^2) .* (w.*sin(w.*f_pos+phi) - 1j.*2.*pi.*t.*cos(w.*f_pos+phi));
result_neg_vec = exp(-1j.*2.*pi.*f_neg.*t) ./ (w.^2 - (2.*pi.*t).^2) .* (w.*sin(w.*f_neg+phi) - 1j.*2.*pi.*t.*cos(w.*f_neg+phi));
result_vec = result_pos_vec - result_neg_vec;

% Set results based on usual conditions
result(cond1) = result_vec(cond1);

% Define expression for abnormal cases
result_neg_vec = (4.*1i.*pi.^2.*f_pos.*(2.*1i.*pi.*f_pos.*1j.*cos(f_pos.*w+phi).*t + 1i.*f_pos.*w.*sin(f_pos.*w+phi) - 2.*1j.*cos(f_pos.*w+phi)).*exp(-2.*1i.*pi.*f_pos.*t)) ./ (-8.*pi.^2);
result_pos_vec = (4.*1i.*pi.^2.*f_neg.*(2.*1i.*pi.*f_neg.*1j.*cos(f_neg.*w+phi).*t + 1i.*f_neg.*w.*sin(f_neg.*w+phi) - 2.*1j.*cos(f_neg.*w+phi)).*exp(-2.*1i.*pi.*f_neg.*t)) ./ (-8.*pi.^2);
result_vec = result_pos_vec - result_neg_vec;

% Set results based on abnormal conditions
result(cond2) = result_vec(cond2);

end

function result = h2(t,f_pos,f_neg,w1,w2,phi1,phi2)
% Integrates exp(-1j*2*pi*f*t)*cos(w1*f+phi1)*cos(w2*f+phi2), from f1 to f2

% Use other h1 function for this (cos(a)cos(b) splits into 2 parts)
part1 = h1(t,f_pos,f_neg,w1+w2,phi1+phi2);
part2 = h1(t,f_pos,f_neg,w1-w2,phi1-phi2);
result = 0.5 * (part1 + part2);

end