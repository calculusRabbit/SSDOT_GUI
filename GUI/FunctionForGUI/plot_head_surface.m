function plot_head_surface(ax1, ax2, AtlasState)

    % Validate required fields
    if ~isfield(AtlasState, 'headsurf') || ...
       ~isfield(AtlasState.headsurf, 'mesh') || ...
       ~isfield(AtlasState.headsurf.mesh, 'vertices') || ...
       ~isfield(AtlasState.headsurf.mesh, 'faces')
        error('AtlasState.headsurf.mesh must contain vertices and faces fields');
    end
    brain_vertices = AtlasState.headsurf.mesh.vertices;
    faces = AtlasState.headsurf.mesh.faces;
    int_at_pos = zeros(size(brain_vertices, 1), 1); % Initialize intensity values
    plot_intensity(ax1, [-1 1],faces, brain_vertices, int_at_pos, 'L')
    plot_intensity(ax2, [-1 1],faces, brain_vertices, int_at_pos, 'R')
end

