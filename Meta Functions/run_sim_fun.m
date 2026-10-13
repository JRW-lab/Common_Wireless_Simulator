function [metrics, frame_data] = run_sim_fun(new_frames, parameters)
%RUN_SIM_FUN  Run new_frames frames of the simulation function that a
%point's parameters select. The single definition of that choice.
%
%   [metrics, frame_data] = run_sim_fun(new_frames, parameters)
%
%   Factored out of sim_save.m on 2026-10-13, unchanged, so that every
%   caller that needs "simulate this point" -- sim_save for collection,
%   collectors/cost_head.m for the back-to-back runtime benchmark -- reaches
%   the same sim_fun by the same rules. A second copy of this switch is how
%   a benchmark ends up timing different code from what was collected.
%
%   Writes nothing: sim_save does the database side.
switch parameters.system_name
    case "TODDM"
        [metrics, frame_data] = sim_fun_TODDM_v3(new_frames,parameters);
    case "ODDM"
        if isfield(parameters,'receiver_name') && parameters.receiver_name == "SIC-MMSE"
            % Time-domain SIC-MMSE via Kronecker-DFT channel recovery,
            % CP-Free only. Perfect-CSI vs. channel-estimation (e.g.
            % DD-RELAX) is chosen INSIDE this function via the
            % parameters.csi_settings field, not via receiver_name -
            % see sim_fun_ODDM_SIC_MMSE.m's own dispatch.
            [metrics, frame_data] = sim_fun_ODDM_SIC_MMSE(new_frames,parameters);
        elseif isfield(parameters,'receiver_name') && parameters.receiver_name == "PT-MMSE"
            % Partitioned Time-Domain MMSE -- the DD-RELAX paper's OWN
            % receiver (its Sec. V / Algorithm 2, published in the
            % 2026-09-16 revision). Added so the sibling `ODDM
            % estimation paper` project's Figure 3 can source its
            % dashed perfect-CSI reference from a detector MATCHING
            % the one producing its solid curve - see
            % sim_fun_ODDM_PTMMSE.m's own header for why a SIC-MMSE or
            % CMC-MMSE reference is not interchangeable here.
            % Perfect-CSI only currently.
            [metrics, frame_data] = sim_fun_ODDM_PTMMSE(new_frames,parameters);
        elseif isfield(parameters,'receiver_name') && parameters.receiver_name == "CMC-MMSE"
            % Native delay-domain CMC-MMSE (References/CP-Free ODDM.pdf
            % Eq. 19-27), fresh implementation (2026-09-14, see
            % equalizer_CMC_MMSE.m's own header for why this is NOT a
            % reuse of the questionable equalizer_CMC_MMSE_AWGN.m).
            % Perfect-CSI only currently - see sim_fun_ODDM_CMC_MMSE.m.
            [metrics, frame_data] = sim_fun_ODDM_CMC_MMSE(new_frames,parameters);
        else
            [metrics, frame_data] = sim_fun_ODDM_v3(new_frames,parameters);
        end
    case "OTFS"
        [metrics, frame_data] = sim_fun_OTFS(new_frames,parameters); % Common method in literature
    case "OTFS-DD"
        if isfield(parameters,'channel_estimation_method') && parameters.channel_estimation_method ~= "none"
            [metrics, frame_data] = sim_fun_OTFS_MUSIC(new_frames,parameters); % MUSIC channel estimation
        else
            [metrics, frame_data] = sim_fun_OTFS_DD_v3(new_frames,parameters); % Dr. Jingxian Wu's design, perfect CSI
        end
    case "OFDM"
        [metrics, frame_data] = sim_fun_OFDM_v2(new_frames,parameters);
    otherwise
        error("Invalid system selected.")
end
end
