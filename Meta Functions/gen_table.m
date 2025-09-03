function gen_table(save_data,conn,table_name,hash_cell,configs,figure_data)

% Import settings
data_type = figure_data.data_type;
primary_var = figure_data.primary_var;
primary_vals = figure_data.primary_vals;
% legend_vec = figure_data.legend_vec;
% line_styles = figure_data.line_styles;
% line_colors = figure_data.line_colors;
% save_sel = figure_data.save_sel;

% % Figure settings
% figures_folder = 'Figures';
% line_val = 2;
% mark_val = 10;
% font_val = 16;

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


% Load loop
results_mat = cell(length(primary_vals),length(configs));
col_names = cell(length(configs),1);
for primvar_sel = 1:size(hash_cell,1)

    % Go through each settings profile
    for sel = 1:size(hash_cell,2)

        % Load data from DB
        paramHash = hash_cell{primvar_sel,sel};
        sim_result = T(string(T.param_hash) == paramHash, :);

        % Select data to extract
        metrics_loaded = jsondecode(sim_result.metrics{1});
        results_val = metrics_loaded.(data_type);

        % Add result to stack
        if isempty(results_val)
            results_mat{primvar_sel,sel} = NaN;
        else
            results_mat{primvar_sel,sel} = results_val;
        end

        if primvar_sel == 1
            % Add row name
            params = configs{sel};
            col_names{sel} = jsonencode_sorted(params);
        end

    end
end

% Create results table
row_names = cellstr(string(primary_vals));
data_table = cell2table(results_mat);
data_table.Properties.DimensionNames{1} = char(primary_var);
data_table.Properties.RowNames = row_names;
data_table.Properties.VariableNames = col_names;

% Display table
disp(primary_var + " sweep table:");
disp(data_table)