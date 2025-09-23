function [render_figure,save_sel,render_delay] = figure_settings()



clc;
                         fprintf("_______________________\n")
                         fprintf("|_____________________|\n")
                         fprintf("||                   ||\n")
                         fprintf("||  Figure Settings  ||\n")
                         fprintf("||___________________||\n")


% Enter values
fprintf("||\n")
render_figure =    str2num(input('|| ❓ Render figure? > ', 's'));

% Settings
if ~render_figure
    render_delay = Inf;
    save_sel = false;
else
    render_delay = 5;
    save_sel = true;
end
clc