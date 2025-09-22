clc; clear;
% Settings
addpath(fullfile(pwd, 'Meta Functions'));
addpath(fullfile(pwd, 'Comm Functions'));
    addpath(fullfile(pwd, 'Comm Functions/Custom Functions'));
    addpath(fullfile(pwd, 'Comm Functions/Generation Functions'));
    addpath(fullfile(pwd, 'Comm Functions/OFDM Functions'));
    addpath(fullfile(pwd, 'Comm Functions/OTFS Functions'));
    addpath(fullfile(pwd, 'Comm Functions/OTFS-DD Functions'));
    addpath(fullfile(pwd, 'Comm Functions/ODDM Functions'));
    addpath(fullfile(pwd, 'Comm Functions/TODDM Functions'));
    addpath(fullfile(pwd, 'Comm Functions/TX RX Functions'));
T = 1/15000;
res_length = 101;
res_int = 10000;
% shape = "rect";
shape = "sinc";
q = 40;
alpha = 1;

% Set up ranges
t_range = linspace(-T,T,res_length);
f_range = linspace(-1/T,1/T,res_length);

xambig_vals = zeros(res_length);
ambig_vals = zeros(res_length);
for t_idx = 1:length(t_range)
    for f_idx = 1:length(f_range)

        fprintf("Testing element (%d,%d)\n",t_idx,f_idx)

        % Select instance inputs
        t = t_range(t_idx);
        f = f_range(f_idx);

        % Generate both versions of ambiguity values
        xambig_vals(t_idx,f_idx) = xambig(t, f, T, shape);
        ambig_vals(t_idx,f_idx) = exp(1j*2*pi*f*t) * ambig_direct(t,-f,T,shape,alpha,q,res_int);
        % ambig_vals(t_idx,f_idx) = ambig_direct(t,f,T,shape,alpha,q,res);


    end
end

figure(1)
subplot(3,1,1)
mesh(abs(xambig_vals))
subplot(3,1,2)
mesh(real(xambig_vals))
subplot(3,1,3)
mesh(imag(xambig_vals))
figure(2)
subplot(3,1,1)
mesh(abs(ambig_vals))
subplot(3,1,2)
mesh(real(ambig_vals))
subplot(3,1,3)
mesh(imag(ambig_vals))
norm(xambig_vals - ambig_vals)