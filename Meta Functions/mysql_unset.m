function mysql_unset(conn,flag_sel)

% Run commands
sqlquery = sprintf("UPDATE system_flags SET flag_value=%d WHERE id=%d", ...
    false, flag_sel);
exec(conn, sqlquery);