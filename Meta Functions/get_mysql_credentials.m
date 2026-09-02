function creds = get_mysql_credentials()
% GET_MYSQL_CREDENTIALS Resolve the {host, user, password} this project
% uses to reach MySQL, walking a first-time user through setup.
%
% The result is cached in mysql_local.json at the project root (an
% untracked, git-ignored file - never committed). On first call, this
% detects whether a local MySQL server is reachable and prompts once for
% credentials; every call after that reads the cached file silently.
%
% Returns [] if the user opts out of MySQL entirely, meaning callers
% should fall back to local (Excel) storage instead of connecting.

cred_file = fullfile(pwd, 'mysql_local.json');

if isfile(cred_file)
    creds = jsondecode(fileread(cred_file));
    return
end

fprintf('\nNo local MySQL configuration found for this project yet.\n');

driver = 'com.mysql.cj.jdbc.Driver';
port = 3306;

if is_port_open('localhost', port, 1.5)
    fprintf('Local MySQL server detected on localhost:%d.\n', port);
    host = 'localhost';
    user_in = strtrim(input('MySQL username [root]: ', 's'));
else
    fprintf('No local MySQL server detected on localhost:%d.\n', port);
    host = strtrim(input('Enter a MySQL host to use, or press Enter to use local Excel storage only: ', 's'));
    if isempty(host)
        creds = [];
        return
    end
    user_in = strtrim(input('MySQL username [root]: ', 's'));
end
if isempty(user_in)
    user_in = 'root';
end

password = input(sprintf('MySQL password for ''%s''@''%s'': ', user_in, host), 's');

% Verify the credentials actually work before saving them
try
    test_conn = database('information_schema', user_in, password, driver, ...
        sprintf('jdbc:mysql://%s:%d/information_schema', host, port));
    if ~isopen(test_conn)
        error('Connection object did not open.');
    end
    close(test_conn);
catch ME
    error(['Could not connect to MySQL with those credentials (%s).\n' ...
        'Re-run to try again, or leave the host blank next time to use local Excel storage only.'], ME.message);
end

creds = struct('host', host, 'user', user_in, 'password', password);
fid = fopen(cred_file, 'w');
fprintf(fid, '%s', jsonencode(creds));
fclose(fid);
fprintf('Saved local MySQL configuration to %s\n(this file is git-ignored and never leaves your machine).\n\n', cred_file);

end

function tf = is_port_open(host, port, timeout_s)
tf = false;
try
    socket = java.net.Socket();
    socket.connect(java.net.InetSocketAddress(host, port), round(timeout_s * 1000));
    tf = true;
    socket.close();
catch
end
end
