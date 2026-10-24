function ctx = ldpc_encode_sp(Lc, u, sp)
%LDPC_ENCODE_SP  Encode with shortening and puncturing applied.
%
%   ctx = ldpc_encode_sp(Lc, u, sp) takes (sp.Kx x numCW) information bits,
%   pads them with sp.s leading ZEROS (the shortened bits, known to both
%   ends and never transmitted), encodes with the full mother code, then
%   drops the shortened information bits and the sp.p trailing parity bits.
%   Returns the (sp.Nx x numCW) bits actually sent.
%
%   See ldpc_shorten_puncture.m for why s and p are equal and what that
%   buys: one mother code, one decoder, one rate across every frame
%   geometry, so a coded comparison between schemes is not also a
%   comparison between codes.

arguments
    Lc struct
    u  double
    sp struct
end
if size(u,1) ~= sp.Kx
    error("ldpc_encode_sp:badInfoLength", ...
        "u has %d rows, expected sp.Kx = %d.", size(u,1), sp.Kx);
end
numCW = size(u,2);
uFull = [zeros(sp.s, numCW); u];           % shortened bits are known zeros
cFull = ldpc_encode_manual(Lc, uFull);     % N x numCW, systematic [u; parity]
ctx = cFull(sp.s+1 : Lc.N - sp.p, :);      % drop shortened info + punctured parity
end
