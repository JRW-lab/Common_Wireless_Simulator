function addMysqlJarOnce()
    persistent jarAdded
    if isempty(jarAdded)
        javaaddpath('mysql-connector-j-8.4.0.jar');
        jarAdded = true;
    end
end
