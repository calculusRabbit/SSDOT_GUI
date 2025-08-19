function M = Make_mask(threshold, Adot, Adot_scalp)
intensity = log10(sum(Adot(:,:,1),1));
mask = find(intensity>threshold);
intensity_scalp = log10(sum(Adot_scalp(:,:,1),1));
mask_scalp = find(intensity_scalp>threshold);
M.mask_brain = mask;
M.mask_scalp = mask_scalp;




