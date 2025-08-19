function G = Make_G_matrix(At_file, M, threshold_brain, threshold_scalp, sigma_brain, sigma_scalp, Adot_orig, Adot_scalp_orig)

atlasViewer = load(At_file);

brain_vertices = atlasViewer.fwmodel.mesh.vertices; %20004*3
scalp_vertices = atlasViewer.fwmodel.mesh_scalp.vertices; %9563*3

mask_brain = M.mask_brain;
mask_scalp = M.mask_scalp;

brain_vertices_masked = brain_vertices(mask_brain,:);
scalp_vertices_masked = scalp_vertices(mask_scalp,:);

[brain_vertices_new] = down_sample_vertices(brain_vertices_masked, threshold_brain);
[scalp_vertices_new] = down_sample_vertices(scalp_vertices_masked, threshold_scalp);


G_brain = make_kernel_matrix_gpu(brain_vertices_new, brain_vertices_masked,sigma_brain);
G_scalp = make_kernel_matrix_gpu(scalp_vertices_new, scalp_vertices_masked,sigma_scalp);

G.G_brain = G_brain;
G.G_scalp = G_scalp;