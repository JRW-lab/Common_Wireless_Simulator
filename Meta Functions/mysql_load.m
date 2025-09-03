function sim_result = mysql_load(conn,table_name,paramHash)

if paramHash ~= "*"
    sqlquery_results = sprintf("SELECT * FROM %s WHERE param_hash = '%s'", ...
        table_name,paramHash);
else
    sqlquery_results = sprintf("SELECT * FROM %s", ...
        table_name);
end
sim_result = fetch(conn, sqlquery_results);