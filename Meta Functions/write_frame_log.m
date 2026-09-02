function frame_idx_start = write_frame_log(log_dir, param_hash, frame_idx_start, frame_data)
% WRITE_FRAME_LOG Append per-frame simulation results to a binary .mat log.
%
%   frame_idx_start = write_frame_log(log_dir, param_hash, frame_idx_start, frame_data)
%
%   frame_data is a struct with fields (all column vectors of length N):
%     bit_errors    — per-frame bit error counts
%     sym_errors    — per-frame symbol error counts
%     frm_errors    — per-frame frame error indicator (0/1)
%     t_RXfull      — per-frame receiver runtime (seconds)
%
%   The log is stored at Logs/<table>/<param_hash>.mat and holds frame_idx,
%   bit_errors, sym_errors, frm_errors, t_RXfull. Frames are appended using
%   the provided frame_idx_start as the offset (global frame counter).
%
%   Returns the new total frame count (last frame index written).

log_path = fullfile(log_dir, param_hash + ".mat");

N = length(frame_data.frm_errors);
new_idx = (frame_idx_start + 1 : frame_idx_start + N).';

if isfile(log_path)
    existing = load(log_path);
    frame_idx = [existing.frame_idx; new_idx];
    bit_errors = [existing.bit_errors; frame_data.bit_errors];
    sym_errors = [existing.sym_errors; frame_data.sym_errors];
    frm_errors = [existing.frm_errors; frame_data.frm_errors];
    t_RXfull = [existing.t_RXfull; frame_data.t_RXfull];
else
    frame_idx = new_idx;
    bit_errors = frame_data.bit_errors;
    sym_errors = frame_data.sym_errors;
    frm_errors = frame_data.frm_errors;
    t_RXfull = frame_data.t_RXfull;
end

% Use compact numeric types to minimize disk footprint
if ~isempty(bit_errors) && max(bit_errors) <= intmax('uint16')
    bit_errors = uint16(bit_errors);
end
if ~isempty(sym_errors) && max(sym_errors) <= intmax('uint16')
    sym_errors = uint16(sym_errors);
end
if ~isempty(frm_errors)
    frm_errors = uint8(frm_errors);
end
if ~isempty(t_RXfull)
    t_RXfull = single(t_RXfull);
end

save(log_path, 'frame_idx', 'bit_errors', 'sym_errors', 'frm_errors', 't_RXfull', '-v7.3');

frame_idx_start = frame_idx_start + N;
end
