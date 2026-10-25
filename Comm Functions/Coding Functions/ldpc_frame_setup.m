function C = ldpc_frame_setup(parameters, nDataBits, default_N)
%LDPC_FRAME_SETUP  The LDPC layer of one frame for the MMSE sim functions.
%
%   C = ldpc_frame_setup(parameters, nDataBits, default_N)
%
%   PURPOSE. Builds (or reuses) the LDPC code a config's ldpc_settings
%   describes and lays whole codewords over a frame's data bits, so that
%   sim_fun_OFDM_v2, sim_fun_OTFS and sim_fun_ODDM_v3 share one coding
%   layer: ldpc_frame_encode fills the frame, symbol_llr and
%   ldpc_frame_decode recover it, ldpc_frame_metrics reports it. The code
%   and interleaver are those of the ODDM PT-MMSE path (ldpc_construct with
%   seed 12345, ldpc_interleaver), so a coded curve means the same code on
%   every system (LDPC_MMSE notes S1).
%
%   INPUTS
%     parameters  merged point struct; ldpc_settings is read through
%                 resolve_ldpc_settings, so its defaults and sweep aliases
%                 (ldpc_rate, ldpc_N, ...) apply unchanged.
%     nDataBits   data bits the frame carries (CP positions excluded).
%     default_N   codeword length when ldpc_settings sets no N.
%
%   OUTPUT
%     C  struct. Always: enable, R_eff (K/N, or 1 when uncoded). When
%        enabled also: Lc (the code), numCW, perm (TX interleaving,
%        identity when interleave is off), bp_iters.
%
%   LIMITS. The codeword length must divide nDataBits exactly; anything
%   else is refused, as on the PT-MMSE path.

%% ===== 1. SETTINGS =====================================================
% Resolves ldpc_settings; an uncoded config returns at once and touches no
% random stream.

% Resolved settings; uncoded configs stop here.
persistent cache
ls = resolve_ldpc_settings(parameters, default_N);
C = struct('enable', ls.enable, 'R_eff', 1);
if ~ls.enable, return; end

%% ===== 2. CODE AND LAYOUT ==============================================
% Builds the code once per (N, method, rate, column weight), checks that
% whole codewords fill the frame, and draws the interleaver.

% The code, cached across calls.
if isempty(cache), cache = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
key = sprintf('%d_%s_%g_%d', ls.N, ls.method, ls.rate, ls.colWeight);
if ~isKey(cache, key), cache(key) = ldpc_construct(ls.N, ls.rate, ls.colWeight, 12345, ls.method); end
Lc = cache(key);

% Whole codewords must fill the data bits exactly.
numCW = nDataBits / Lc.N;
if numCW < 1 || numCW ~= fix(numCW)
    error("ldpc_frame_setup:sizeMismatch", ...
        "LDPC codeword length %d does not divide the frame's %d data bits evenly -- choose ldpc_settings.N to divide it.", ...
        Lc.N, nDataBits);
end

% Interleaver over the whole frame (private stream), or the identity.
if ls.interleave
    perm = ldpc_interleaver(nDataBits);
else
    perm = (1:nDataBits).';
end

% The layer's description.
C.Lc = Lc;
C.numCW = numCW;
C.perm = perm;
C.bp_iters = ls.bp_iters;
C.R_eff = Lc.K / Lc.N;
end
