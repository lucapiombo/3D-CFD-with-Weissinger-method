function plotFigures(plotData, geom, savefig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot aerodynamic distributions for vortex lattice analysis
% Inputs:
%   plotData: Structure containing aerodynamic coefficients
%   geom:     Geometry structure
%   savefig:  Boolean flag to save figures (default: false)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set default for savefig if not provided
    if nargin < 3
        savefig = false;
    end
    
    % Extract commonly used parameters
    [N, M] = size(plotData.Cp);
    quarter_chord_idx = round(N/4);
    mid_span_idx = min(6, M);  % Use 6th span station or max available
    
    % Create figure handles for better control
    fig_handles = gobjects(6, 1);
    
    % 1. Cp distribution at quarter chord
    fig_handles(1) = figure('Name', 'Cp Quarter Chord');
    plotCpQuarterChord(plotData, geom, quarter_chord_idx);
    saveFigureIfRequested(fig_handles(1), 'Cp single wing', savefig);
    
    % 2. Cp distribution chordwise
    fig_handles(2) = figure('Name', 'Cp Chordwise');
    plotCpChordwise(plotData, geom);
    saveFigureIfRequested(fig_handles(2), 'Cp chord', savefig);
    
    % 3. CL distribution spanwise
    fig_handles(3) = figure('Name', 'CL Spanwise');
    plotCLSpanwise(plotData, geom, mid_span_idx);
    saveFigureIfRequested(fig_handles(3), 'CL single wing', savefig);
    
    % 4. CD distribution spanwise
    fig_handles(4) = figure('Name', 'CD Spanwise');
    plotCDSpanwise(plotData, geom);
    saveFigureIfRequested(fig_handles(4), 'CD single wing', savefig);
    
    % 5. Gamma distribution 3D surface
    fig_handles(5) = figure('Name', 'Gamma Distribution 3D');
    plotGamma3D(plotData, geom);
    saveFigureIfRequested(fig_handles(5), 'GAMMA wing', savefig);
    
    % 6. Gamma distribution spanwise
    fig_handles(6) = figure('Name', 'Gamma Spanwise');
    plotGammaSpanwise(plotData, geom);
    saveFigureIfRequested(fig_handles(6), 'GAMMA spanwise', savefig);
end


%% Helper Functions for Individual Plots
function plotCpQuarterChord(plotData, geom, quarter_chord_idx)
    plot(geom.centroids.y(quarter_chord_idx, :), ...
         plotData.Cp(quarter_chord_idx, :), ...
         '--ob', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    
    formatPlot('C_p distribution at 1/4 chord', ...
               'Span', 'C_p', true);
end


function plotCpChordwise(plotData, geom)
    % Use color cycling for better distinction
    colors = lines(size(plotData.Cp, 2));
    
    for i = 1:size(plotData.Cp, 2)
        plot(geom.centroids.x(:, i), plotData.Cp(:, i), ...
             'Color', colors(i, :), 'LineWidth', 1.5);
        hold on;
    end
    hold off;
    
    formatPlot('C_p distribution along chord', ...
               'Chord', 'C_p', true);
    
    % Add legend for first few span stations
    if size(plotData.Cp, 2) <= 10
        legend(arrayfun(@(x) sprintf('Station %d', x), 1:size(plotData.Cp, 2), ...
               'UniformOutput', false), 'Location', 'best');
    end
end


function plotCLSpanwise(plotData, geom, mid_span_idx)
    plot(geom.centroids.y(mid_span_idx, :), ...
         plotData.C_L(mid_span_idx, :), ...
         '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
         'MarkerFaceColor', [0.2, 0.6, 1]);
    
    formatPlot('C_L distribution along the span', ...
               'Span', 'C_L', true);
end


function plotCDSpanwise(plotData, geom)
    
    for i = 1:size(plotData.C_D, 1)
        plot(geom.centroids.y(i, :), plotData.C_D(i, :), 'LineWidth', 1.2);
        hold on;
    end
    hold off;
    
    formatPlot('C_D distribution along the span', ...
               'Span', 'C_D', true);
end


function plotGamma3D(plotData, geom)
    % Create 3D surface plot
    surf(geom.panels.x, geom.panels.y, geom.panels.z, plotData.gamma, ...
         'EdgeColor', 'k', 'FaceAlpha', 0.9);
    
    colorbar;
    title('Circulation Distribution \Gamma', 'FontSize', 15);
    xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
    zlabel('z', 'FontSize', 12, 'FontWeight', 'bold');
    axis equal;
    view(3);
    grid on;
    box on;
    
    % Set z-limits dynamically based on data
    z_limits = [min(geom.panels.z(:)) - 0.1, max(geom.panels.x(:)) + 0.1];
    zlim(z_limits);
end


function plotGammaSpanwise(plotData, geom)
    % Use distinctive markers for each chord station
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};
    
    hold on;
    for i = 1:min(size(plotData.gamma, 1), length(markers))
        marker_idx = mod(i-1, length(markers)) + 1;
        plot(geom.centroids.y(i, :), plotData.gamma(i, :), ...
             ['-', markers{marker_idx}], ...
             'LineWidth', 1.2, 'MarkerSize', 6);
    end
    hold off;
    
    formatPlot('\Gamma distribution along the span', ...
               'Span', '\Gamma', true);
    
    % Add legend for first few chord stations
    if size(plotData.gamma, 1) <= 10
        legend(arrayfun(@(x) sprintf('Chord %d', x), 1:size(plotData.GAMMA, 1), ...
               'UniformOutput', false), 'Location', 'best');
    end
end


%% Utility Functions
function formatPlot(plotTitle, xLabel, yLabel, addGrid)
    title(plotTitle, 'FontSize', 15);
    xlabel(xLabel, 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(yLabel, 'FontSize', 12, 'FontWeight', 'bold');
    
    if addGrid
        grid on;
        grid minor;
    end
    
    set(gca, 'FontSize', 11);
    set(gcf, 'Color', 'white');
end


function saveFigureIfRequested(figHandle, filename, saveflag)
    if saveflag
        % Create figures directory if it doesn't exist
        if ~exist('figures', 'dir')
            mkdir('figures');
        end
        
        % Save in high resolution
        saveas(figHandle, fullfile('figures', filename), 'png');
        
        % Also save as MATLAB figure for editing
        savefig(figHandle, fullfile('figures', [filename, '.fig']));
        
        fprintf('Saved figure: %s\n', filename);
    end
end