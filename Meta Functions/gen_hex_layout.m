function gen_hex_layout(save_data,conn,table_name,default_parameters,configs,figure_data)

clf
% Figure settings
figures_folder = 'Figures';
line_val = 2;
font_val = 16;
hex_radius = 1;
color_sel = "parula";
tile_scale = 1.6;
m_offset = 0.1;
u_linelength = 1.5;
u_offset = 0.1;

% Load data from DB and set new frame count
switch save_data.priority
    case "mysql"
        if save_data.save_mysql
            T = mysql_load(conn,table_name,"*");
        elseif save_data.save_excel
            try
                T = readtable(save_data.excel_path, 'TextType', 'string');
            catch
                T = table;
            end
        end
    case "local"
        if save_data.save_excel
            try
                T = readtable(save_data.excel_path, 'TextType', 'string');
            catch
                T = table;
            end
        elseif save_data.save_mysql
            T = mysql_load(conn,table_name,"*");
        end
end

% Import settings
data_type = figure_data.data_type;
primary_var = figure_data.primary_var;
primary_vals = figure_data.primary_vals;
save_sel = figure_data.save_sel;

% Load loop
var_names = fieldnames(default_parameters);
results_vec = zeros(length(primary_vals),length(configs));
for primvar_sel = 1:length(primary_vals)

    % Set primary variable
    primval_sel = primary_vals(primvar_sel);

    % Go through each settings profile
    for sel = 1:length(configs)

        % Populate parameters from configs or default_parameters
        params_inst = struct();
        for var_sel = 1:length(var_names)
            var_name = var_names{var_sel};
            if isfield(configs{sel}, var_name)
                value = configs{sel}.(var_name);
            elseif ~strcmp(var_name, primary_var)
                value = default_parameters.(var_name);
            else
                value = primval_sel;
            end
            params_inst.(var_name) = value;
        end

        % Make number of symbols per frame
        syms_per_f = params_inst.M*params_inst.N*params_inst.U;
        U = params_inst.U;

        % Load data from DB
        [~,paramHash] = jsonencode_sorted(params_inst);
        sim_result = T(string(T.param_hash) == paramHash, :);

        % Select data to extract
        results_inst = jsondecode(sim_result.metrics{1});
        switch data_type
            case "BER"
                results_val = results_inst.BER;
            case "SER"
                results_val = results_inst.SER;
            case "FER"
                results_val = results_inst.FER;
            case "Thr"
                frame_duration = params_inst.N * params_inst.T;
                bandwidth_hz = params_inst.M * U / params_inst.T;
                results_val = (log2(params_inst.M_ary) * syms_per_f * (1 - results_inst.FER)) / (frame_duration * bandwidth_hz);
            case "Cap"
                results_val = results_inst.C / 1000;

        end

        % Add result to stack
        if isempty(results_val)
            results_vec(primvar_sel,sel) = NaN;
        else
            results_vec(primvar_sel,sel) = results_val;
        end

    end
end

val_list = zeros(length(configs),3);
for i = 1:length(configs)
    val_list(i,1) = configs{i}.M;
    val_list(i,2) = configs{i}.N;
    val_list(i,3) = configs{i}.U;
end

% Get starting values and coordinates
starting_vals = min(val_list,[],1);
coord_mat = log2(val_list) - log2(starting_vals) + 1;
starting_exp = log2(starting_vals);

% Extract normalized coordinates
m_idx = coord_mat(:,1);
n_idx = coord_mat(:,2);
u_idx = coord_mat(:,3);

% Determine changes in movement
dx = sqrt(3) * hex_radius * cos(pi/6);
dy = sqrt(3) * hex_radius;

% Get hex coordinates
x_coords = dx * u_idx;
y_coords = dy * m_idx + dy*(u_idx-1)/2;

% Get u-label parameters
height = max(y_coords);
u_level = height + u_linelength;

% Startup figure
figure(1);
hold on;
axis equal off;
cmap = eval(color_sel + "(256)");

% Avoid log(0)
epsilon = min(results_vec(results_vec > 0)) / 10;
% epsilon = 1e-15;
log_results = log10(results_vec + epsilon);

% Define the log tick range (orders of magnitude to display)
log_tick_min = floor(min(log_results));
log_tick_max = ceil(max(log_results));
log_ticks = log_tick_min:log_tick_max;

% Normalize results based on log range
norm_results = (log_results - log_tick_min) / (log_tick_max - log_tick_min);
color_idx = round(norm_results * 255) + 1;

% Clamp to valid colormap index range
color_idx = max(min(color_idx, 256), 1);

% Generate hex for each point
m_count = 0;
u_count = 0;
m_coords = zeros(length(unique(m_idx)),2);
n_coords = zeros(length(unique(n_idx)),2);
u_coords = zeros(length(unique(u_idx)),2);
for i = 1:length(results_vec)

    if i ~= 1
        last_u = u_idx(i-1);
    else
        last_u = 1;
    end

    % Get coordinates
    x = x_coords(i);
    y = y_coords(i);

    % Get needed m coordinate
    if u_idx(i) == 1
        m_count = m_count + 1;
        angle = 7*pi/6;
        radius = sqrt(3)/2 + m_offset;
        m_coords(m_count,1) = x + radius * cos(angle);
        m_coords(m_count,2) = y + radius * sin(angle);
    end

    % Get needed n coordinate
    for j = 1:length(n_coords)
        if n_idx(i) == j && m_idx(i) == 1
            angle = -pi/6;
            radius = sqrt(3)/2 + m_offset;
            n_coords(j,1) = x + radius * cos(angle);
            n_coords(j,2) = y + radius * sin(angle);
        end
    end

    % Get needed u coordinate
    if n_idx(i) == 1
        u_count = u_count + 1;
        u_coords(u_count,1) = x;
        u_coords(u_count,2) = y + hex_radius * sqrt(3) / 2;
    end

    % Get color for this hex
    c = cmap(color_idx(i), :);

    % Draw hexagon
    draw_hexagon(x, y, hex_radius, c);

    % Compute luminance (using standard perceptual weights)
    luminance = 0.299*c(1) + 0.587*c(2) + 0.114*c(3);

    % Choose black or white based on brightness
    if luminance > 0.5
        contrast_color = [0, 0, 0];  % dark text on light background
    else
        contrast_color = [1, 1, 1];  % light text on dark background
    end

    % Write value of hex
    text_string = sprintf("%1.1d",results_vec(i));
    text(x, y, text_string, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        "Color",contrast_color,"FontSize",font_val/tile_scale);

end

% Flip coordinates
% m_coords(:,1) = flip(m_coords(:,1));
% m_coords(:,2) = flip(m_coords(:,2));
n_coords(:,1) = flip(n_coords(:,1));
n_coords(:,2) = flip(n_coords(:,2));
u_coords(:,1) = flip(u_coords(:,1));
u_coords(:,2) = flip(u_coords(:,2));

% Get unique indices
unique_m = unique(m_idx);
unique_n = unique(n_idx);
unique_u = unique(u_idx);

% Plot M indices
for i = 1:length(unique_m)

    % Add the label "M=..."
    text(m_coords(i,1), m_coords(i,2), sprintf('M=%d', 2^(unique_m(i)-1+starting_exp(1))), ...
        'HorizontalAlignment', 'right', ...
        'FontSize', font_val, ...
        'Rotation',30);

end

% Plot N indices
for i = 1:length(unique_n)

    % Add the label "N=..."
    text(n_coords(i,1), n_coords(i,2), sprintf('N=%d', 2^(unique_n(length(unique_n)-i+1)-1+starting_exp(2))), ...
        'HorizontalAlignment', 'left', ...
        'FontSize', font_val, ...
        'Rotation',-30);
end

% Plot U indices
for i = 1:length(unique_u)

    % Add the label "U=..."
    text(u_coords(i,1), u_level, sprintf('U=%d', 2^(unique_u(i)-1+starting_exp(3))), ...
        'HorizontalAlignment', 'right', ...
        'FontSize', font_val, ...
        'Rotation', -45);

    % Plot U index line
    plot([u_coords(i,1), u_coords(i,1)], [u_level-u_offset, u_coords(i,2)], 'k--', 'LineWidth', line_val);

end

% Create color bar
cb = colorbar(fontsize=font_val);
% cb.Position = colorbar_loc;
colormap(cb.Parent, eval(color_sel)); %#ok<EVLCS>
cb.Ticks = (log_ticks - log_tick_min) / (log_tick_max - log_tick_min);
cb.TickLabels = arrayfun(@(x) sprintf('10^{%d}', x), log_ticks, 'UniformOutput', false);
xlabel(cb, "BER","FontSize",font_val);

% Create folders if they don't exist
subfolder = fullfile(figures_folder, ['/' char(data_type)]);
subsubfolder = fullfile(subfolder, ['/' char(primary_var)]);
if ~exist(figures_folder, 'dir')
    mkdir(figures_folder);
end
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end
if ~exist(subsubfolder, 'dir')
    mkdir(subsubfolder);
end

% Save figure
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
timestamp_str = char(timestamp);
figure_filename = fullfile(subsubfolder, "Figure_" + timestamp_str + ".png");
if save_sel
    saveas(figure(1), figure_filename);
end

end

function draw_hexagon(x, y, radius, face_color)
angles = linspace(0, 2*pi, 7);
x_hex = x + radius * cos(angles);
y_hex = y + radius * sin(angles);
fill(x_hex, y_hex, face_color, 'EdgeColor', 'k');
end