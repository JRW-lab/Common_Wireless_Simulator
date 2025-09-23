function result = xambig(t,f,T,select_pulse_shape)
% This function results the result of the cross-ambiguity function of two
% pulses of similar shape. Supported pulse shapes are rectangular ("rect")
% and sinc ("sinc").
%
% Coded by Jeremiah Rhys Wimer, 2/21/2024
tol = .001 * T;
result = zeros(size(t));
switch select_pulse_shape
    case "rect"
        % Set conditions: 1 and 2 are upper and lower bounds
        %                 3 and 4 set the limits as f->0
        %                 5 sets to the otherwise statement
        cond1 = abs(f) >= tol & (-T <= t & t < 0);
        cond2 = abs(f) >= tol & (0 <= t & t <= T);
        cond3 = abs(f) < tol & (-T <= t & t < 0);
        cond4 = abs(f) < tol & (0 <= t & t <= T);

        % Select value based on condition
        result(cond1) = (exp(1j.*2.*pi.*f(cond1).*t(cond1)) - exp(-1j.*2.*pi.*f(cond1).*T)) ./ (1j.*2.*pi.*f(cond1).*T);
        result(cond2) = (1 - exp(1j.*2.*pi.*f(cond2).*(t(cond2)-T))) ./ (1j.*2.*pi.*f(cond2).*T);
        result(cond3) = (1 + t(cond3)./T);
        result(cond4) = (1 - t(cond4)./T);
    case "sinc"
        % Set conditions: 1 and 2 are upper and lower bounds
        %                 3 and 4 set the limits as f->0
        %                 5 sets to the otherwise statement
        cond1 = abs(t) >= tol & (-1/T <= f & f < 0);
        cond2 = abs(t) >= tol & (0 <= f & f <= 1/T);
        cond3 = abs(t) < tol & (-1/T <= f & f < 0);
        cond4 = abs(t) < tol & (0 <= f & f <= 1/T);

        % Select value based on condition
        result(cond1) = (T ./ (pi.*t(cond1))) .* sin(((pi.*t(cond1)) ./ T) + pi.*f(cond1).*t(cond1)) .* exp(1j.*pi.*f(cond1).*t(cond1));
        result(cond2) = (T ./ (pi.*t(cond2))) .* sin(((pi.*t(cond2)) ./ T) - pi.*f(cond2).*t(cond2)) .* exp(1j.*pi.*f(cond2).*t(cond2));
        result(cond3) = (1 + f(cond3).*T);
        result(cond4) = (1 - f(cond4).*T);
    case "ideal"
        % Convert to discrete numerics
        n = t / T;
        m = f * T;

        % Set conditions: 1 is for 
        cond1 = (n == 0 & m == 0);
        cond2 = ~cond1;

        % Select value based on condition
        result(cond1) = 1;
        result(cond2) = 0;
    otherwise
        % Spit error if selecting different pulse shape
        error('Select a supported pulse shape: rect or sinc')
end
end