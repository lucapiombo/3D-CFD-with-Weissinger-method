function plotFigures(plotData, geom, savefig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot aerodynamic distributions for vortex lattice analysis
% Inputs:
%   plotData: Structure or cell array containing aerodynamic coefficients
%             For single geometry: struct
%             For two geometries: {struct1, struct2}
%   geom:     Geometry structure or cell array of structures
%   savefig:  Boolean flag to save figures (default: false)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set default for savefig if not provided
    if nargin < 3
        savefig = false;
    end
    
    % Determine if we have single or multiple geometries
    if iscell(plotData)
        nGeom = length(plotData);
        isMultipleGeom = (nGeom > 1);
    else
        nGeom = 1;
        isMultipleGeom = false;
        % Convert to cell for consistent handling
        plotData = {plotData};
        geom = {geom};
    end
    
    % Create figure handles for better control
    fig_handles = gobjects(6, 1);
    
    % 1. Cp distribution at quarter chord
    fig_handles(1) = figure('Name', 'Cp Quarter Chord');
    if isMultipleGeom
        plotCpQuarterChordMultiple(plotData, geom);
        saveFigureIfRequested(fig_handles(1), 'Cp multiple wings', savefig);
    else
        plotCpQuarterChord(plotData{1}, geom{1});
        saveFigureIfRequested(fig_handles(1), 'Cp single wing', savefig);
    end
    
    % 2. Cp distribution chordwise
    fig_handles(2) = figure('Name', 'Cp Chordwise');
    if isMultipleGeom
        plotCpChordwiseMultiple(plotData, geom);
        saveFigureIfRequested(fig_handles(2), 'Cp chord multiple', savefig);
    else
        plotCpChordwise(plotData{1}, geom{1});
        saveFigureIfRequested(fig_handles(2), 'Cp chord', savefig);
    end
    
    % 3. CL distribution spanwise
    fig_handles(3) = figure('Name', 'CL Spanwise');
    if isMultipleGeom
        plotCLSpanwiseMultiple(plotData, geom);
        saveFigureIfRequested(fig_handles(3), 'CL multiple wings', savefig);
    else
        plotCLSpanwise(plotData{1}, geom{1});
        saveFigureIfRequested(fig_handles(3), 'CL single wing', savefig);
    end
    
    % 4. CD distribution spanwise
    fig_handles(4) = figure('Name', 'CD Spanwise');
    if isMultipleGeom
        plotCDSpanwiseMultiple(plotData, geom);
        saveFigureIfRequested(fig_handles(4), 'CD multiple wings', savefig);
    else
        plotCDSpanwise(plotData{1}, geom{1});
        saveFigureIfRequested(fig_handles(4), 'CD single wing', savefig);
    end
    
    % 5. Gamma distribution 3D surface
    fig_handles(5) = figure('Name', 'Gamma Distribution 3D');
    if isMultipleGeom
        plotGamma3DMultiple(plotData, geom);
        saveFigureIfRequested(fig_handles(5), 'GAMMA multiple wings', savefig);
    else
        plotGamma3D(plotData{1}, geom{1});
        saveFigureIfRequested(fig_handles(5), 'GAMMA wing', savefig);
    end
    
    % 6. Gamma distribution spanwise
    fig_handles(6) = figure('Name', 'Gamma Spanwise');
    if isMultipleGeom
        plotGammaSpanwiseMultiple(plotData, geom);
        saveFigureIfRequested(fig_handles(6), 'GAMMA spanwise multiple', savefig);
    else
        plotGammaSpanwise(plotData{1}, geom{1});
        saveFigureIfRequested(fig_handles(6), 'GAMMA spanwise', savefig);
    end
end


%% Helper Functions for Single Geometry (original functions, slightly modified)
function plotCpQuarterChord(plotData, geom)
    [N, ~] = size(plotData.Cp);
    quarter_chord_idx = round(N/4);
    
    plot(geom.centroids.y(quarter_chord_idx, :), ...
         plotData.Cp(quarter_chord_idx, :), ...
         '--ob', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    
    formatPlot('C_p distribution at 1/4 chord', ...
               'Span', 'C_p', true);
end

function plotCpQuarterChordMultiple(plotData, geom)
    colors = lines(length(plotData));
    hold on;
    
    for i = 1:length(plotData)
        [N, ~] = size(plotData{i}.Cp);
        quarter_chord_idx = round(N/4);
        
        plot(geom{i}.centroids.y(quarter_chord_idx, :), ...
             plotData{i}.Cp(quarter_chord_idx, :), ...
             '--o', 'LineWidth', 1.5, 'MarkerFaceColor', colors(i,:), ...
             'Color', colors(i,:), 'DisplayName', sprintf('Geometry %d', i));
    end
    hold off;
    
    formatPlot('C_p distribution at 1/4 chord', ...
               'Span', 'C_p', true);
    legend('show', 'Location', 'best');
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

function plotCpChordwiseMultiple(plotData, geom)
    nGeom = length(plotData);
    subplotRows = ceil(sqrt(nGeom));
    subplotCols = ceil(nGeom / subplotRows);
    
    for g = 1:nGeom
        subplot(subplotRows, subplotCols, g);
        
        colors = lines(size(plotData{g}.Cp, 2));
        for i = 1:size(plotData{g}.Cp, 2)
            plot(geom{g}.centroids.x(:, i), plotData{g}.Cp(:, i), ...
                 'Color', colors(i, :), 'LineWidth', 1.5);
            hold on;
        end
        hold off;
        
        formatPlot(sprintf('C_p distribution (Geometry %d)', g), ...
                   'Chord', 'C_p', true);
    end
end

function plotCLSpanwise(plotData, geom)
    [~, M] = size(plotData.Cp);
    mid_span_idx = min(6, M);
    
    plot(geom.centroids.y(mid_span_idx, :), ...
         plotData.C_L(mid_span_idx, :), ...
         '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
         'MarkerFaceColor', [0.2, 0.6, 1]);
    
    formatPlot('C_L distribution along the span', ...
               'Span', 'C_L', true);
end

function plotCLSpanwiseMultiple(plotData, geom)
    nGeom = length(plotData);
    subplotRows = ceil(sqrt(nGeom));
    subplotCols = ceil(nGeom / subplotRows);
    
    for g = 1:nGeom
        subplot(subplotRows, subplotCols, g);
        
        [~, M] = size(plotData{g}.Cp);
        mid_span_idx = min(6, M);
        
        plot(geom{g}.centroids.y(mid_span_idx, :), ...
             plotData{g}.C_L(mid_span_idx, :), ...
             '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
             'MarkerFaceColor', [0.2, 0.6, 1]);
        
        formatPlot(sprintf('C_L distribution (Geometry %d)', g), ...
                   'Span', 'C_L', true);
    end
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

function plotCDSpanwiseMultiple(plotData, geom)
    nGeom = length(plotData);
    subplotRows = ceil(sqrt(nGeom));
    subplotCols = ceil(nGeom / subplotRows);
    
    for g = 1:nGeom
        subplot(subplotRows, subplotCols, g);
        
        for i = 1:size(plotData{g}.C_D, 1)
            plot(geom{g}.centroids.y(i, :), plotData{g}.C_D(i, :), 'LineWidth', 1.2);
            hold on;
        end
        hold off;
        
        formatPlot(sprintf('C_D distribution (Geometry %d)', g), ...
                   'Span', 'C_D', true);
    end
end

function plotGamma3D(plotData, geom)
    % Create 3D surface plot for single geometry
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
    z_limits = [min(geom.panels.z(:)) - 0.1, max(geom.panels.z(:)) + 0.1];
    zlim(z_limits);
end


function plotGamma3DMultiple(plotData, geom)
    % Create 3D surface plot with multiple geometries shown together
    
    hold on;
    for i = 1:length(plotData)
        surf(geom{i}.panels.x, geom{i}.panels.y, geom{i}.panels.z, ...
             plotData{i}.gamma, 'EdgeColor', 'k');
    end
    hold off;
    
    colorbar;
    title('Circulation Distribution \Gamma (Multiple Geometries)', 'FontSize', 15);
    xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
    zlabel('z', 'FontSize', 12, 'FontWeight', 'bold');
    axis equal;
    view(3);
    grid on;
    box on;
    
    % Set appropriate limits for combined geometries
    all_x = [];
    all_y = [];
    all_z = [];
    for i = 1:length(geom)
        all_x = [all_x; geom{i}.panels.x(:)];
        all_y = [all_y; geom{i}.panels.y(:)];
        all_z = [all_z; geom{i}.panels.z(:)];
    end
    xlim([min(all_x)-0.1, max(all_x)+0.1]);
    ylim([min(all_y)-0.1, max(all_y)+0.1]);
    zlim([min(all_z)-0.1, max(all_z)+0.1]);
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
        legend(arrayfun(@(x) sprintf('Chord %d', x), 1:size(plotData.gamma, 1), ...
               'UniformOutput', false), 'Location', 'best');
    end
end


function plotGammaSpanwiseMultiple(plotData, geom)
    % Plot gamma for multiple geometries together in one plot
    colors = lines(length(plotData));
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};
    
    hold on;
    for g = 1:length(plotData)
        % For each geometry, plot only the first few chord stations
        nStations = min(size(plotData{g}.gamma, 1), 3); % Limit to 3 stations for clarity
        
        for i = 1:nStations
            marker_idx = mod(i-1, length(markers)) + 1;
            plot(geom{g}.centroids.y(i, :), plotData{g}.gamma(i, :), ...
                 ['-', markers{marker_idx}], ...
                 'LineWidth', 1.2, 'MarkerSize', 6, ...
                 'Color', colors(g,:), ...
                 'DisplayName', sprintf('Geom%d, Chord%d', g, i));
        end
    end
    hold off;
    
    formatPlot('\Gamma distribution along the span (Multiple Geometries)', ...
               'Span', '\Gamma', true);
    legend('show', 'Location', 'best', 'NumColumns', 2);
end



%% Utility Functions (unchanged)
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