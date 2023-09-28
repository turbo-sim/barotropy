function set_plot_options
% Printing options
FontName = 'monospace';
FontSize = 11;
AxLineWidth = 0.75;
LineWidth = 1.00;
MarkerSize = 5;

% Setting LaTeX as text interpreter
set(groot,'DefaultTextInterpreter','latex', ...
    'DefaultAxesTickLabelInterpreter','latex', ...
    'DefaultLegendInterpreter','latex', ...
    'DefaultColorbarTickLabelInterpreter','latex', ...
    'DefaultConstantlineInterpreter', 'latex', ...
    'DefaultGraphplotInterpreter', 'latex', ...
    'DefaultPolaraxesTickLabelInterpreter', 'latex', ...
    'DefaultTextarrowshapeInterpreter', 'latex', ...
    'DefaultTextboxshapeInterpreter', 'latex');%, ...
%     'DefaultColorbarLabelInterpreter','latex');

% Setting several default options for the graphics
set(groot, 'DefaultFigureColor','White' , ...
'DefaultFigurePaperType', 'a4letter', ...
'DefaultAxesColor', 'white', ...
'DefaultAxesFontUnits', 'points',...
'DefaultAxesFontSize', FontSize, ...
'DefaultAxesFontAngle', 'normal', ...
'DefaultAxesGridLineStyle', '-', ...
'DefaultAxesGridAlpha', 0.25, ...
'DefaultAxesGridColor', [0, 0, 0], ...
'DefaultAxesInterruptible', 'on', ...
'DefaultAxesLayer', 'Bottom', ...
'DefaultAxesNextPlot', 'replace', ...
'DefaultAxesUnits', 'normalized', ...
'DefaultAxesXcolor', [0, 0, 0], ...
'DefaultAxesYcolor', [0, 0, 0], ...
'DefaultAxesZcolor', [0, 0, 0], ...
'DefaultAxesVisible', 'on', ...
'DefaultAxesLineWidth', AxLineWidth, ...
'DefaultLineLineWidth', LineWidth, ...
'DefaultLineMarkerSize', MarkerSize, ...
'DefaultTextColor', [0, 0, 0], ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontName', FontName, ...    % Not applied when LaTeX is used
'DefaultAxesFontName', FontName, ...    % Not applied when LaTeX is used
'DefaultTextFontSize', FontSize, ...
'DefaultLegendFontSize', FontSize, ...
'DefaultLegendFontSizeMode','manual', ...
'DefaultTextVerticalAlignment', 'middle', ...
'DefaultTextHorizontalAlignment', 'left', ...
'DefaultLegendItemTokenSize', [15,18], ...
'DefaultColorbarFontsize', FontSize, ...
'DefaultColorbarTickLength', 0.015, ...
'DefaultAxesTickLength', [0.02 0.025]);

% Setting a position for the figure
set(groot,'DefaultFigurePosition', [360   198   560   420]);

end


% Small script to print all the available fields to a file
% s = get(groot, 'factory');
% file_handle = fopen('factory_settings.txt', 'wt');
% field_names = fieldnames(s);
% for i = 1:numel(field_names)
%     fprintf(file_handle, '%s\n', field_names{i});
% end
% fclose(file_handle);
