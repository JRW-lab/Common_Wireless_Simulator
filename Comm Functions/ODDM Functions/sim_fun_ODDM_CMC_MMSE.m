function [metrics,frame_data] = sim_fun_ODDM_CMC_MMSE(new_frames,parameters)
% Native DELAY-DOMAIN CMC-MMSE receiver for CP-Free ODDM (References/
% CP-Free ODDM.pdf, Eq. 19-27), via the FRESH, from-the-paper
% equalizer_CMC_MMSE.m (see that file's own header for why this is NOT a
% reuse of the questionable equalizer_CMC_MMSE_AWGN.m). Mirrors
% sim_fun_ODDM_SIC_MMSE.m's exact frame/guard/noise conventions
% (same data generation, same zero-padding-at-the-end guard, same
% (N+2)/N noise scaling) so this is a directly paired, apples-to-apples
% comparison against that receiver -- built 2026-09-14 specifically to
% test whether the paper's own native CMC-MMSE receiver closes any of the
% still-open v=500/high-EbN0 gap that equalizer_SIC_MMSE.m (adapted from
% a DIFFERENT paper, SIC-MMSE.pdf/OTFS) leaves unexplained -- see
% notes/AMBIGUITY_TABLE_AUDIT.md.
%
% KEY SIMPLIFICATION vs sim_fun_ODDM_SIC_MMSE.m: equalizer_CMC_MMSE.m
% operates DIRECTLY on gen_HDD_direct.m's native (delay-major) output --
% no Kronecker-DFT similarity transform, no permutation to Doppler-major,
% no time-domain conversion at all. The paper's own Eq. 17-19 system model
% IS this native MxM-grid-of-NxN-blocks layout.

% Make parameters
fields = fieldnames(parameters);
for i = 1:numel(fields)
    eval([fields{i} ' = parameters.(fields{i});']);
end

if isfield(parameters,'csi_settings') && parameters.csi_settings.method ~= "none"
    error("sim_fun_ODDM_CMC_MMSE:unsupportedCsi", ...
        "sim_fun_ODDM_CMC_MMSE only supports perfect CSI (csi_settings.method=='none') currently.");
end

if CP
    error("CMC-MMSE for ODDM requires CP-Free mode - zero-padding (not a cyclic prefix) is what makes the channel causal here.")
end
if ~exist('cmc_iters','var')
    if exist('N_iters','var')
        cmc_iters = N_iters;
    else
        cmc_iters = 8;
    end
end

% Define parameters (identical to sim_fun_ODDM_SIC_MMSE.m)
res = 10;
Es = 1;
syms_per_f = M*N;
Ts = T / M;
L1 = Q + 1;
L2 = Q + 1 + floor(2510*10^(-9) / Ts);
L = L1 + L2;
if L >= M
    error("Settings can not satisfy the ambiguity assumption (zero-padding length L >= M).")
end
Eb = Es / log2(M_ary);
N0 = Eb / (10^(EbN0 / 10)) * ((N+2)/N); % same CP-free noise convention as sim_fun_ODDM_v3/sim_fun_ODDM_SIC_MMSE

if shape == "rect"
    Q = 1;
    alpha = 1;
elseif shape == "sinc"
    alpha = 1;
end

if M_ary == 2
    bit_order = [0;1];
    alphabet_set = linspace(1,M_ary,M_ary)';
    S = sqrt(Es) .* exp(-1j * 2*pi .* (alphabet_set) ./ M_ary);
elseif M_ary == 4
    bit_order = [0,0;0,1;1,0;1,1];
    S = zeros(4,1);
    S(1) = (sqrt(2)/2) + (1j*sqrt(2)/2);
    S(2) = (sqrt(2)/2) - (1j*sqrt(2)/2);
    S(3) = -(sqrt(2)/2) + (1j*sqrt(2)/2);
    S(4) = -(sqrt(2)/2) - (1j*sqrt(2)/2);
    S = sqrt(Es) .* S;
end

[Ambig_Table.vals,Ambig_Table.t_range,Ambig_Table.f_range] = gen_DD_cross_ambig_table(N,M,T,Fc,vel,shape,alpha,Q,res);

% The last L delay bins are the zero-padding guard, exactly matching
% sim_fun_ODDM_SIC_MMSE.m's own convention (verified to already be the
% empirically-better "guard at the end" placement -- see CLAUDE.md-
% adjacent notes/AMBIGUITY_TABLE_AUDIT.md "zero-padding placement" entry).
zero_syms = N*L;
change_map_vert = false(M*N,1);
change_map_vert((M-L)*N+1:end) = true;
known_mask = false(M,1); known_mask((M-L+1):M) = true;   % delay-layer mask (Mx1), for equalizer_CMC_MMSE's own known_mask arg

bit_errors = zeros(new_frames,1);
sym_errors = zeros(new_frames,1);
frm_errors = zeros(new_frames,1);
iters_vec = zeros(new_frames,1);
t_RXfull_vec = zeros(new_frames,1);

for frame = 1:new_frames
    tStartRXfull = tic;

    % Generate data (native delay-major), force the zero-padding guard
    % bins to zero -- NO Doppler-major reindex, NO time-domain transform:
    % equalizer_CMC_MMSE.m consumes xDD/HDD in exactly this native layout.
    [TX_bit,TX_sym,xDD] = gen_data(bit_order,S,syms_per_f);
    TX_bit(change_map_vert,:) = -1;
    TX_sym(change_map_vert) = -1;
    xDD(change_map_vert) = 0;

    % Generate the ODDM DD-domain channel, native delay-major -- used
    % directly, no permutation or Kronecker transform.
    t_offset = max_timing_offset * Ts;
    HDD_native = gen_HDD_direct(T,N,M,Fc,vel,Q,Ambig_Table,t_offset,false);

    % Generate noise and received signal, directly in the native DD domain.
    w = sqrt(N0/2) * (randn(syms_per_f,1) + 1j*randn(syms_per_f,1));
    y = HDD_native*xDD + w;

    % 2026-09-14: switched to equalizer_CMC_MMSE_AWGN.m (the pre-existing
    % implementation, now fixed to handle known/guard layers) - verified
    % bit-for-bit identical to equalizer_CMC_MMSE_native.m across 300
    % paired frames (v=500/EbN0=16dB) after that fix, so this is a pure
    % consolidation onto one canonical implementation, not a behavior
    % change. Argument order note: this function takes (Lp,Ln), the
    % REVERSE of L1,L2's own (backward,forward) sense - Lp=L2 (forward
    % bound), Ln=-L1 (negative backward bound) - see its own header.
    [x_hat,iters_vec(frame)] = equalizer_CMC_MMSE_AWGN(y,HDD_native,N,M,L2,-L1,Es,N0,S,cmc_iters,known_mask,xDD);

    % Hard detection for final x_hat (x_hat's data-layer entries are
    % already hard-decided inside equalizer_CMC_MMSE.m; this re-slice is
    % a no-op for those and correctly classifies known-layer entries too,
    % matching sim_fun_ODDM_SIC_MMSE.m's own final dist/min pattern)
    dist = abs(x_hat.' - S).^2;
    [~,min_index] = min(dist);
    RX_sym = min_index.' - 1;
    RX_bit = bit_order(RX_sym+1,:);

    bit_error_vec = TX_bit(~change_map_vert,:) ~= RX_bit(~change_map_vert,:);
    sym_error_vec = TX_sym(~change_map_vert) ~= RX_sym(~change_map_vert);
    bit_errors(frame) = sum(bit_error_vec(:));
    sym_errors(frame) = sum(sym_error_vec);
    if any(bit_error_vec(:))
        frm_errors(frame) = 1;
    end

    t_RXfull_vec(frame) = toc(tStartRXfull);
end

frame_duration = N * T;
bandwidth_hz = M / T;

metrics.BER = sum(bit_errors,"all") / (new_frames*(syms_per_f-zero_syms)*log2(M_ary));
metrics.SER = sum(sym_errors,"all") / (new_frames*(syms_per_f-zero_syms));
metrics.FER = sum(frm_errors,"all") / (new_frames);
metrics.Thr = (log2(M_ary) * (syms_per_f-zero_syms) * (1 - metrics.FER)) / (frame_duration * bandwidth_hz);
metrics.RX_iters = mean(iters_vec);
metrics.t_RXfull = mean(t_RXfull_vec);

frame_data.bit_errors = bit_errors;
frame_data.sym_errors = sym_errors;
frame_data.frm_errors = frm_errors;
frame_data.t_RXfull = t_RXfull_vec;
frame_data.bits_per_frame = (syms_per_f-zero_syms) * log2(M_ary);
frame_data.syms_per_frame = syms_per_f-zero_syms;
end
