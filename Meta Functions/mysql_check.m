function flag_val = mysql_check(conn,flag_sel)

% Run commands
sqlquery = sprintf("SELECT * FROM system_flags WHERE id = '%d'", ...
    flag_sel);
flag_row = fetch(conn, sqlquery);
flag_val = flag_row.flag_value;