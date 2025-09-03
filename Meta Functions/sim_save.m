function sim_save(save_data,conn,table_name,current_frames,parameters,paramHash)

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
try
    sim_result = T(string(T.param_hash) == paramHash, :);
catch
    sim_result = [];
end

if ~isempty(sim_result)
    % Find new frame count to simulate
    if sim_result.frames_simulated < current_frames
        new_frames = current_frames - sim_result.frames_simulated;
        run_flag = true;
    else
        run_flag = false;
    end
else
    % Simulate given frame count
    new_frames = current_frames;
    run_flag = true;
end

% Run if needed
if run_flag

    % Simulate needed system
    switch parameters.system_name
        case "TODDM"
            metrics_add = sim_fun_TODDM_v3(new_frames,parameters);
        case "ODDM"
            metrics_add = sim_fun_ODDM_v3(new_frames,parameters);
        case "OTFS"
            metrics_add = sim_fun_OTFS_v3(new_frames,parameters);
        case "OFDM"
            metrics_add = sim_fun_OFDM_v2(new_frames,parameters);
    end

    % Write to database
    switch save_data.priority
        case "mysql"
            if save_data.save_mysql
                mysql_write(conn,table_name,parameters,new_frames,metrics_add);
            end
            if save_data.save_excel
                T = mysql_load(conn,table_name,"*");
                excel_path = save_data.excel_path;
                writetable(T, excel_path);
            end
        case "local"
            if save_data.save_excel
                excel_path = save_data.excel_path;
                local_write(excel_path,parameters,new_frames,metrics_add);
            end
            if save_data.save_mysql
                mysql_write(conn,table_name,parameters,new_frames,metrics_add);
            end
    end

end