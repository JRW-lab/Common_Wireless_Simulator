function [chan_coef, delays, Dopplers] = channel_generation(car_fre, max_speed)
% function [chan_coef, delays,Dopplers] = channel_generation(car_fre, max_speed)
%
% physical channel generation with an EVA delay profile
%
% car_fre: carrier frequency
% max_speed: max UE speed in km/hr
%
% Yanjun Pan
% University of Arkansas
% 10/18/2022
%

% rng('default');

% EVA delay profile
% excess tap delay
delays = [0 30 150 310 370 710 1090 1730 2510]*10^(-9);     % in sec
% relative power
pdp = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9];        % in dB 

pow_prof = 10.^(pdp/10);                                    % dB to watt
pow_prof = pow_prof/sum(pow_prof);                          % normalization                        
chan_coef = sqrt(pow_prof) .* (sqrt(1/2) * ...
    (randn(1,length(delays)) + 1i*randn(1,length(delays)) ) );          % chan coefficient

max_UE_speed = max_speed * (1000/3600);                                 % max UE speed in m/s
max_Doppler = (max_UE_speed*car_fre) / (physconst('LightSpeed'));
Dopplers = (max_Doppler*cos(2*pi*rand(1,length(delays))));              % Dopplers
end