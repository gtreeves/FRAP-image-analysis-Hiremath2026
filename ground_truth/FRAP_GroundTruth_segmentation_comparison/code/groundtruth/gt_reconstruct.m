function gt_reconstruct(statefile, idx_bleach, nframes, outdir, id)
%GT_RECONSTRUCT Reconstruct the bleached nucleus with the pipeline's own
%   reconstruct_missing_nucleus (copy that also returns the pre-erosion
%   labels), using the saved state from gt_run_pipeline.
%   Writes <outdir>\<id>_pipeline_masks_reconstructed.mat
S = load(statefile);
[L_nuc, ~, ~, ~, ~, label] = reconstruct_missing_nucleus_track(S.mask_data, S.label, S.nuc_data, ...
    S.I_dl, S.I_nuc, S.t, idx_bleach, nframes, S.datum, S.filenameshort);
close all
T = size(S.I_nuc, 3); L_track = zeros(size(S.I_nuc));
for i = 1:T, L_track(:,:,i) = label{i,1}; end
I_nuc = S.I_nuc; I_dl = S.I_dl; t = S.t; scalings = S.scalings; filenameshort = S.filenameshort;
side = S.side; filenum = S.filenum; reconstructed = true; recon_frames = nframes;
save(fullfile(outdir, [id '_pipeline_masks_reconstructed.mat']), 'L_track', 'L_nuc', 'I_nuc', 'I_dl', 't', ...
    'scalings', 'filenameshort', 'side', 'filenum', 'reconstructed', 'recon_frames', '-v7.3');
fprintf('%s: reconstructed bleached nucleus %d over frames 2-%d\n', id, idx_bleach, nframes);
end
