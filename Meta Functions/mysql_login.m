function conn = mysql_login(dbname)

% Parameters
users = ["root","remote_user"];
servers = ["localhost","192.168.42.2"];
password = "Jumpydoll3";
port     = 3306;
driver = 'com.mysql.cj.jdbc.Driver';

% Connect through all options
for i = 1:length(users)
    dburl = ['jdbc:mysql://' servers(i) ':' char(string(port)) '/' dbname];
    if users(i) == "root"
        conn = database(dbname, users(i), password, driver, strjoin(dburl));
    else
        conn = connectWithRetry(dbname, users(i), password, driver, strjoin(dburl));
    end
    if isopen(conn)
        break
    end
end

% Connection Check
if ~isopen(conn)
    error("Failure to form connection to MySQL database...")
end