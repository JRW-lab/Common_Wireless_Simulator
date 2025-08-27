function [profile_sel,num_frames,delete_sel] = profile_select(profile_names,frame_sel)

% Parameters
flag = true;

% Set up interface to MySQL server
clc
while flag
    clc;
    fprintf("___________________________________________________\n")
    fprintf("|__________________________________________________|\n")
    fprintf("||                                                ||\n")
    fprintf("|| Saved Profiles:                                ||\n")
    fprintf("||________________________________________________||\n")
    for i = 1:length(profile_names)
        fprintf("||                                                ||\n")
        row = sprintf("|| %d:",i);
        fprintf("%s",row)
        for j = 1:(50-length(char(row)))
            fprintf(" ")
        end
        fprintf("||\n")
        row = sprintf("|| %s",profile_names{i});
        fprintf("%s",row)
        for j = 1:(50-length(char(row)))
            fprintf(" ")
        end
        fprintf("||\n")
    end
    fprintf("||________________________________________________||\n")
    fprintf("||________________________________________________||\n")

    % Decision logic
    fprintf("\n")
    profile_sel = str2double(input(' > Select profile: ', 's'));
    if profile_sel <= length(profile_names) && profile_sel > 0
        flag = false;
    else
        fprintf("\n   ❌ Enter valid profile number!\n");
        input('', 's')
    end
    if frame_sel
        num_frames = str2double(input(' > Enter number of frames: ', 's'));
    else
        num_frames = 0;
    end
    if num_frames == 0
        delete_sel = false;
    else
        delete_sel = str2double(input(' > Delete selected data? ', 's'));
    end

end
clc