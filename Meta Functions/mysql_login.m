function conn = mysql_login(dbname)

driver = 'com.mysql.cj.jdbc.Driver';
port = 3306;

% Resolve this machine's MySQL credentials (prompts once on first use;
% cached locally afterward). Returns [] if MySQL is disabled.
creds = get_mysql_credentials();
if isempty(creds)
    error("MySQL is not configured for this project - delete mysql_local.json to re-run setup, or disable 'Enable MySQL' to use local Excel storage only.");
end

% Auto-provision the schema the first time this account connects to a
% fresh MySQL server (a plain `database(dbname,...)` call below fails
% outright if dbname doesn't exist yet).
try
    admin_conn = database('information_schema', creds.user, creds.password, driver, ...
        sprintf('jdbc:mysql://%s:%d/information_schema', creds.host, port));
    if isopen(admin_conn)
        execute(admin_conn, "CREATE DATABASE IF NOT EXISTS " + dbname);
        close(admin_conn);
    end
catch
    % Schema may already exist, or this account may lack CREATE privileges -
    % fall through and let the real connection attempt below report the
    % actual problem.
end

dburl = sprintf('jdbc:mysql://%s:%d/%s', creds.host, port, dbname);
conn = connectWithRetry(dbname, creds.user, creds.password, driver, dburl);

% Connection Check
if ~isopen(conn)
    error("Failure to form connection to MySQL database...");
end
