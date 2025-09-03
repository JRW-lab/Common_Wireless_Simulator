function mysql_set(conn,flag_sel)

% Run commands
sqlquery = sprintf("UPDATE system_flags SET flag_value=%d WHERE id=%d", ...
    true, flag_sel);
exec(conn, sqlquery);