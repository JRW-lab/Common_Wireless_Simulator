function conn = connectWithRetry(dbname, username, password, driver, url)
    % CONNECTWITHRETRY: Keeps retrying until database connection succeeds
    %
    % Example:
    %   conn = connectWithRetry('mydb', 'user', 'pass', ...
    %       'org.postgresql.Driver', 'jdbc:postgresql://192.168.42.1:5432/mydb');

    maxRetries = 50;      % Max number of attempts
    waitTime   = 30;      % Seconds to wait between retries

    for attempt = 1:maxRetries
        try
            conn = database(dbname, username, password, driver, url);

            % If connected successfully
            if isopen(conn)
                fprintf('✅ Connected on attempt %d\n', attempt);
                return;
            else
                fprintf('⚠️ Attempt %d: connection object created but not open.\n', attempt);
            end

        catch ME
            fprintf('❌ Attempt %d failed: %s\n', attempt, ME.message);
        end

        % Wait before retry
        pause(waitTime);
    end

    error('[%s] Failed to connect to database after %d attempts.', datetime('now'), maxRetries);
end