
% Define the base directory
baseDir = '/scratch.global/kuma0458/';

% Define the layout dimensions
rows = 3;
cols = 2;

% Initialize the figure
%figure('Name', 'Flowrate Analysis', 'Position', [100, 100, 1000, 900]);
f=figure;
% Create a 3x2 tiled layout with compact spacing
tl = tiledlayout(rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');

% Cell arrays containing directory names and tile titles
% Formatted as {row1_col1, row1_col2; row2_col1, row2_col2; row3_col1, row3_col2}
directories = {
    'c-2ak1_re180', 'c-2ak2_re180';
    'c0ak1_re180',  'c0ak2_re180';
    'c2ak1_re180',  'c2ak2_re180'
};

titles = {
    'ak=0.1, c=-2', 'ak=0.2, c=-2';
    'ak=0.1, c=0',  'ak=0.2, c=0';
    'ak=0.1, c=2',  'ak=0.2, c=2'
};

% Define y-limits for each row (corresponding to c=-2, c=0, c=2)
yLimits = {
    [15, 17]; % Row 1 (c = -2)
    [13, 15]; % Row 2 (c = 0)
    [11, 13]  % Row 3 (c = 2)
};

% Loop through rows and columns to generate plots
for r = 1:rows
    for col = 1:cols
        % Advance to the next tile (fills across rows left-to-right by default)
        nexttile;
        
        % Build the full file path including baseDir and the 'run' folder
        filepath = fullfile(baseDir, directories{r, col}, 'run', 'flowrate.mat');
        
        % Load and plot the data if the file exists
        if exist(filepath, 'file')
            % Loads variables 't' and 'flowrate' into a struct 'data'
            data = load(filepath, 't', 'flowrate'); 
            
            % Plot the variables
            plot(data.t, data.flowrate, 'b-', 'LineWidth', 1.5);
        else
            % Display warning and a placeholder text on the plot if missing
            warning('File not found: %s', filepath);
            text(0.5, 0.5, 'Data File Missing', 'HorizontalAlignment', 'center');
        end
        
        % Add the specific title for this tile
        title(titles{r, col}, 'FontWeight', 'bold');
        
        % Apply the specific y-axis limits for this row
        ylim(yLimits{r});
        
        grid on;
    end
end

% Add shared overarching labels and title to the entire tiled layout
xlabel(tl, 't', 'FontWeight', 'bold', 'FontSize', 12);
ylabel(tl, 'flowrate', 'FontWeight', 'bold', 'FontSize', 12);
title(tl, 'Flowrate vs. t across different (ak, c) configurations', 'FontSize', 14);


saveas(f,'flowrate_plots.fig')