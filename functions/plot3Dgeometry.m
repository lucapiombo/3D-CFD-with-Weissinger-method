function plot3Dgeometry(geom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the geometry of the wing(s)
% Input:
%   geom: Geometry structure or cell array of structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Determine if we have single or multiple geometries
    if iscell(geom)
        nGeom = length(geom);
        isMultipleGeom = (nGeom > 1);
    else
        nGeom = 1;
        isMultipleGeom = false;
        % Convert to cell for consistent handling
        geom = {geom};
    end
    
    % Create figure
    figure();
    
    % Calculate overall limits for consistent axes
    all_x = [];
    all_y = [];
    all_z = [];
    
    hold on;
    
    for g = 1:nGeom
        % Extract data for current geometry
        px = geom{g}.panels.x;
        py = geom{g}.panels.y;
        pz = geom{g}.panels.z;
        
        cx = geom{g}.centroids.x;
        cy = geom{g}.centroids.y;
        cz = geom{g}.centroids.z;
        
        vx = geom{g}.vortices.x;
        vy = geom{g}.vortices.y;
        vz = geom{g}.vortices.z;
        
        % Plot panels
        surf(px, py, pz, 'FaceColor', [0.8 0.8 0.8], ...
             'EdgeColor', 'k', 'FaceAlpha', 0.5);
        
        % Plot centroids
        plot3(cx(:), cy(:), cz(:), 'ob', ...
              'MarkerFaceColor', 'b', ...
              'MarkerSize', 2);
        
        % Plot vortices
        plot3(vx(:), vy(:), vz(:), 'xr', ...
              'MarkerFaceColor', 'r', ...
              'MarkerSize', 3);
        
        % Collect all points for axis limits
        all_x = [all_x; px(:); cx(:); vx(:)];
        all_y = [all_y; py(:); cy(:); vy(:)];
        all_z = [all_z; pz(:); cz(:); vz(:)];
    end
    
    hold off
    
    % Set labels and properties
    xlabel('x', 'FontSize', 12, 'fontweight', 'bold');
    ylabel('y', 'FontSize', 12, 'fontweight', 'bold');
    zlabel('z', 'FontSize', 12, 'fontweight', 'bold');
    
    % Set title
    if isMultipleGeom
        title(sprintf('Wing Geometry (%d wings)', nGeom), 'FontSize', 14, 'FontWeight', 'bold');
    else
        title('Wing Geometry', 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    % Set axis properties
    axis equal;
    view(3);
    grid on;
    box on;
    
    % Set consistent axis limits with some padding
    if ~isempty(all_x)
        x_padding = (max(all_x) - min(all_x)) * 0.1;
        y_padding = (max(all_y) - min(all_y)) * 0.1;
        z_padding = (max(all_z) - min(all_z)) * 0.1;
        
        if x_padding > 0
            xlim([min(all_x)-x_padding, max(all_x)+x_padding]);
        end
        if y_padding > 0
            ylim([min(all_y)-y_padding, max(all_y)+y_padding]);
        end
        if z_padding > 0
            zlim([min(all_z)-z_padding, max(all_z)+z_padding]);
        end
    end
    
    % Add lighting for better 3D visualization
    camlight headlight;
    lighting gouraud;
    material dull;
end