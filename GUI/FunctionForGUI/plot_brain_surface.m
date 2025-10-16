function plot_brain_surface(ax1, ax2, AtlasState)

    % Validate required fields
    if ~isfield(AtlasState, 'pialsurf') || ...
       ~isfield(AtlasState.pialsurf, 'mesh') || ...
       ~isfield(AtlasState.pialsurf.mesh, 'vertices') || ...
       ~isfield(AtlasState.pialsurf.mesh, 'faces')
        error('AtlasState.pialsurf.mesh must contain vertices and faces fields');
    end
    brain_vertices = AtlasState.pialsurf.mesh.vertices;
    faces = AtlasState.pialsurf.mesh.faces;
    int_at_pos = zeros(size(brain_vertices, 1), 1); % Initialize intensity values
    plot_intensity(ax1, [-1 1],faces, brain_vertices, int_at_pos, 'L')
    plot_intensity(ax2, [-1 1],faces, brain_vertices, int_at_pos, 'R')
end