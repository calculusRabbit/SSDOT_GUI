function A = Make_A_matrix(Adot, Adot_scalp, E, M)

if ~isempty(M)
    mask = M.mask_brain;
    mask_scalp = M.mask_scalp;

    Adot = Adot(:, mask, :);
    Adot_scalp = Adot_scalp(:, mask_scalp, :);
end

A.Amatrix_brain.w1_HbO = squeeze(Adot(:,:,1))*E(1,1);
A.Amatrix_brain.w1_HbR = squeeze(Adot(:,:,1))*E(1,2);
A.Amatrix_brain.w2_HbO = squeeze(Adot(:,:,2))*E(2,1);
A.Amatrix_brain.w2_HbR = squeeze(Adot(:,:,2))*E(2,2);

A.Amatrix_scalp.w1_HbO = squeeze(Adot_scalp(:,:,1))*E(1,1);
A.Amatrix_scalp.w1_HbR = squeeze(Adot_scalp(:,:,1))*E(1,2);
A.Amatrix_scalp.w2_HbO = squeeze(Adot_scalp(:,:,2))*E(2,1);
A.Amatrix_scalp.w2_HbR = squeeze(Adot_scalp(:,:,2))*E(2,2);