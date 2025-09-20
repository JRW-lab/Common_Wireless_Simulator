function [fading] = generate_fading(Ns, L, Fd, Ts, t_offset)
% FUNCTION FADING = RAYLEIGH1I_quick(Ns, Fd, Ts)
%    Fd : Maxmimum Doppler Frequency
%    Ts : Sampling Period
%    Ns : the number of samples in the fading
%
%    the same as Rayleigh1I.m, only a little faster
%    not proper for generating long Rayleigh fading sequences
%    Generate Rayleigh fading with the Zheng's new simulation model
%
%    Author: Jingxian Wu
%    Organization: University of Missouri - Columbia
%    Construction Date: Feb. 28, 2002
%
%    Rev 1.1  (Apr. 10, 2002)
%

% INPUTS - FOR DEBUGGING
% clear; clc;
% Ns = 2;
% Fd = 1;
% Ts = 1;
% M = 16;

% WU'S METHOD-------------------------------------------------------------
M = 16;
n = 1:M;

fading = zeros(L,Ns);
for l = 1:L
    theta = (2*rand-1)*pi;
    phi = (2*rand(1,M)-1)*pi;
    psi = (2*rand(1,M)-1)*pi;
    alpha = (2*pi*n-pi+theta)/(4*M);
    omega_d = 2*pi*Fd;
    t = (0:Ns-1)*Ts + t_offset;

    omega_t = omega_d*cos(alpha)'*t;
    phase_c = omega_t + phi.'*ones(1, Ns);
    phase_s = omega_t + psi.'*ones(1, Ns);

    Uc = sum(cos(phase_c))/sqrt(M);
    Us = sum(sin(phase_s))/sqrt(M);

    fading(l,:) = Uc + 1j*Us;
end

% % RHYS' METHOD------------------------------------------------------------
% R = zeros(Ns);
% for m = 1:Ns
%     for n = 1:Ns
%         R(m,n) = besselj(0, 2*pi*Fd*(m-n)*Ts);
%     end
% end
% 
% [v,u,~] = svd(R);
% R_half = v * sqrt(u);
% 
% fading = zeros(L,Ns);
% for l = 1:L
%     uncorr_vec = sqrt(1/2) * (randn(Ns,1) + 1j*randn(Ns,1));
%     fading(l,:) = (R_half * uncorr_vec).';
% end

% % MOBIN'S METHOD----------------------------------------------------------
% M = 1000;
% DopplerShifts = Fd * cos(rand([L M]) * 2*pi);
% RandomPhase = rand([L M]);
% IndexDelayTaps = 1:L;
% fading3 = zeros(L,Ns);
% for i = 1:Ns
%     t = (i - 1) * Ts;
%     fading3(IndexDelayTaps,i) = sum(exp(1j * 2*pi * (RandomPhase + DopplerShifts * t)), 2) ./ sqrt(M);
% end

end